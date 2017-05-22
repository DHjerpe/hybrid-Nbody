#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "main.h"
#include <assert.h>

#define EPSILON 1e-3

int NPROCS;

/* Five input arguments:
 1: N = number of stars                                         (argv[1])
 2: nsteps = number of time steps                               (argv[2])
 3: delta_t = size of time step                                 (argv[3])
 4: VERIFY, 0 = NO, 1 = YES (Default)                           (argv[4])
 */
int main(int argc, char** argv) {
    
    if (argc > 4) { // run only if number of arguments is correct
        
        double wtime_end;
        double wtime_start = MPI_Wtime();
        int my_rank;
        int N = atoi(argv[1])*5;
        int nsteps = atoi(argv[2]);
	int verify = atoi(argv[4]);
        double delta_t = atof(argv[3]);
        
        double * p = NULL;
        p = (double *)malloc(sizeof(double)*N);
        
        
        
        generateStars(p, N/5);
        
        
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
    
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
        
        
        forceCalc(p,N,delta_t,nsteps);
        
        
        MPI_Finalize();

        
        if (verify) {
            
            double * p_ref = NULL;
            p_ref = (double *)malloc(sizeof(double)*N);
            generateStars(p_ref, N/5);
            
            forceCalcVerify(p_ref,N,delta_t,nsteps);
            
            for (int i = 0; i<N/5; i++) {
                assert(p[i] == p_ref[i]);
                assert(p[i+1] == p_ref[i+1]);
                assert(p[i+2] == p_ref[i+2]);
                assert(p[i+3] == p_ref[i+3]);
                assert(p[i+4] == p_ref[i+4]);
            }
            
            printf("Verify OK! \n");
            
        }
        
        
        free(p);
        p = NULL;
        
        wtime_end = MPI_Wtime();
        
        if (my_rank == 0)

            printf("Total time: %.4f \n", wtime_end - wtime_start);
        
        
    } else {printErrorMsg();}
    
    
    return 0;
}

void forceCalc(double * p, int N,double delta_t, int nsteps) {
  
    N = N/5;
    double G = 100.0/N;
    
    int upperLimit = N;
    
    int my_rank;
    int interval = upperLimit/NPROCS;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
    int lowerLimit = my_rank * interval;
    
        if (my_rank <  NPROCS - 1) {
            
            upperLimit = my_rank * interval + interval;
            
        }
 
    
    int numberToSend = (upperLimit - lowerLimit)*5;


    double * p_buffer = NULL;
    p_buffer = (double *)malloc(sizeof(double)*N*5);
    memcpy(p_buffer,p,sizeof(double)*N*5);
    
    
    for (int n = 0; n<nsteps; n++) {
        
#pragma omp parallel
        {
            double r; double mass;
            double x_diff, y_diff;
            double a_x, a_y;
            double u_x, u_y;
            
            
        #pragma omp for 
        for (int i = lowerLimit; i<upperLimit; i++) {
            double sumX = 0;
            double sumY = 0;
            
            for (int j = 0; j<N; j++) {
                
                if (i != j) {
                    x_diff = (p[i*5]-p[j*5]);     // rij (vector), x-component
                    y_diff = (p[i*5+1]-p[j*5+1]); // rij (vector), y-component
                    
                    r = sqrt((x_diff*x_diff) + (y_diff*y_diff));
                    
                    mass = p[j*5+2];
                    mass /= ((r+EPSILON)*(r+EPSILON)*(r+EPSILON)); // not mass. total constant
                    
                    x_diff *= mass;
                    y_diff *= mass;
                    
                    sumX += x_diff;
                    sumY += y_diff;
                }
            }
            
            
            a_x = -G*sumX;
            a_y = -G*sumY;
            
            // veloceties
            u_x = p[i*5+3] + delta_t*a_x;
            u_y = p[i*5+4] + delta_t*a_y;
            
            // update veloceties
            p_buffer[i*5+3] = u_x;
            p_buffer[i*5+4] = u_y;
            
            // update position
            p_buffer[i*5] = p[i*5] + delta_t*u_x;
            p_buffer[i*5+1] = p[i*5+1] + delta_t*u_y;
        }
        
        }
        
        MPI_Allgather(MPI_IN_PLACE, numberToSend, MPI_DOUBLE, p_buffer, numberToSend, MPI_DOUBLE, MPI_COMM_WORLD);
        
        
        // copy the new values to p
        memcpy(p,p_buffer,sizeof(double)*N*5);
        
    }
    
    free(p_buffer);
    p_buffer = NULL;
}


void forceCalcVerify(double * p, int N,double delta_t, int nsteps) {
    
    N = N/5;
    double G = 100.0/N;
    
    double * p_buffer = NULL;
    p_buffer = (double *)malloc(sizeof(double)*N*5);
    memcpy(p_buffer,p,sizeof(double)*N*5);
    
    
    for (int n = 0; n<nsteps; n++) {

            double r; double mass;
            double x_diff, y_diff;
            double a_x, a_y;
            double u_x, u_y;
    

            for (int i = 0; i<N; i++) {
                double sumX = 0;
                double sumY = 0;
                
                for (int j = 0; j<N; j++) {
                    
                    if (i != j) {
                        x_diff = (p[i*5]-p[j*5]);     // rij (vector), x-component
                        y_diff = (p[i*5+1]-p[j*5+1]); // rij (vector), y-component
                        
                        r = sqrt((x_diff*x_diff) + (y_diff*y_diff));
                        
                        mass = p[j*5+2];
                        mass /= ((r+EPSILON)*(r+EPSILON)*(r+EPSILON)); // not mass. total constant
                        
                        x_diff *= mass;
                        y_diff *= mass;
                        
                        sumX += x_diff;
                        sumY += y_diff;
                    }
                }
            
                a_x = -G*sumX;
                a_y = -G*sumY;
                
                // veloceties
                u_x = p[i*5+3] + delta_t*a_x;
                u_y = p[i*5+4] + delta_t*a_y;
                
                // update veloceties
                p_buffer[i*5+3] = u_x;
                p_buffer[i*5+4] = u_y;
                
                // update position
                p_buffer[i*5] = p[i*5] + delta_t*u_x;
                p_buffer[i*5+1] = p[i*5+1] + delta_t*u_y;
            }
    
        // copy the new values to p
        memcpy(p,p_buffer,sizeof(double)*N*5);
    }
    
    free(p_buffer);
    p_buffer = NULL;
}



/**
 p is a vector containing information about all stars/particles in the system
 p contains: 
 - x position
 - y position
 - mass
 - x velocity
 - y velocity
 
 N is the total number of stars/particles
 */
void generateStars(double * p, int N) {
    
    int randLoc = 256;
    srand(randLoc); // choose a locatoin in the "random numbers" array
    
    for (int i = 0; i<N; i++) {
        
        p[i*5] = (double)rand() / (double)((unsigned)RAND_MAX + 1); // pos x
        ++randLoc; srand(randLoc);
        p[i*5+1] = (double)rand() / (double)((unsigned)RAND_MAX + 1); // pos y
        ++randLoc; srand(randLoc);
        p[i*5+2] = (double)rand() / (double)((unsigned)RAND_MAX + 1); // mass
        ++randLoc; srand(randLoc);
        p[i*5+3] = (double)rand() / (double)((unsigned)RAND_MAX + 1); // velocity x
        ++randLoc; srand(randLoc);
        p[i*5+4] = (double)rand() / (double)((unsigned)RAND_MAX + 1); // velocity y
        ++randLoc; srand(randLoc);
    }
}

void printErrorMsg() {
    
    printf("Error, wrong number of input arguments.\nThe input arguments should have the following form:  \n");
    printf("1: The number of stars \n");
    printf("2: The number of time steps \n");
    printf("3: The size of time step  \n");
    printf("4: Verify ON/OFF (1/0)  \n");
}

