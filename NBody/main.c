#include <stdio.h>
#include <stdlib.h>
#include "file_operations.h"
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>


#define EPSILON 1e-3

void printErrorMsg();
void forceCalc(double * p, int N,double delta_t, int nsteps);

int NPROCS;
int NUMTHREADS;


/* Five input arguments:
 1: N = number of stars                                         (argv[1])
 2: filename = file to read initial configuration of the system (argv[2])
 3: nsteps = number of time steps                               (argv[3])
 4: delta_t = size of time step                                 (argv[4])
 5: graphics is either 1/0 (on/off)                             (argv[5])
 6: number of threads in each node                              (argv[6])
 */
int main(int argc, char** argv) {
    
    if (argc > 5) { // run only if number of arguments is correct
        
        int N = atoi(argv[1])*5;
        int nsteps = atoi(argv[3]);
        double delta_t = atof(argv[4]);
        
        NUMTHREADS = atoi(argv[6]);
        
        
        double * p = NULL;
        p = (double *)malloc(sizeof(double)*N);
        
        read_doubles_from_file(N, p, argv[2]);
        
        char *outp = "result.gal";
        
        
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
        
        
        forceCalc(p,N,delta_t,nsteps);
        
        
        MPI_Finalize();
        
        write_doubles_to_file(N,p,outp);
        
        free(p);
        p = NULL;
        
    } else {printErrorMsg();}
    
    
    return 0;
}
// TODO: Optimise
void forceCalc(double * p, int N,double delta_t, int nsteps) {
  
    N = N/5;
    double G = 100.0/N;
    
    double r; double mass;
    double x_diff, y_diff;
    double a_x, a_y;
    double u_x, u_y;
    
    
    int upperLimit = N;
    
    //args_p args[NUMTHREADS];
   // pthread_t threads[NUMTHREADS];
    
    
    int my_rank;
    int interval = upperLimit/NPROCS;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
   // for (int i = 0; i< NPROCS; i++) {
        
       // args[i].lowerLimit = i*interval;
    
    int lowerLimit = my_rank * interval;
    
        if (my_rank <  NPROCS - 1) {
          //  args[i].upperLimit = i*interval + interval;
            
            upperLimit = my_rank * interval + interval;
            
        }
    //    else {
           // args[i].upperLimit = upperLimit;
      //  }
        
       // args[i].p_buffer = p_buffer;
       // args[i].p = p;
      // args[i].quad_tree = quad_tree;
    //}
    
    int numberToSend = (upperLimit - lowerLimit)*5;
    //int numberToSend = (upperLimit - lowerLimit);
    

    double * p_buffer = NULL;
    p_buffer = (double *)malloc(sizeof(double)*N*5);
    memcpy(p_buffer,p,sizeof(double)*N*5);
    
    
    omp_set_num_threads(NUMTHREADS);
    
    
    for (int n = 0; n<nsteps; n++) {
        
#pragma omp parallel
        {
        
        #pragma omp parallel for
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
            
            //printf("forceX: %lf \t forceY: %lf \n",sumX,sumY);
            
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
        
        /* int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
         void *recvbuf, int recvcount, MPI_Datatype recvtype,
         MPI_Comm comm) */
        
        }
        
        MPI_Allgather(MPI_IN_PLACE, numberToSend, MPI_DOUBLE, p_buffer, numberToSend, MPI_DOUBLE, MPI_COMM_WORLD);
        
        
        // copy the new values to p

        memcpy(p,p_buffer,sizeof(double)*N*5);
        
        
        

        
    }
    
    free(p_buffer);
    p_buffer = NULL;
}



void printErrorMsg() {
    
    printf("Error, wrong number of input arguments.\nThe input arguments shoulc hav the form:  \n");
    printf("1: number of stars \n");
    printf("2: file to read initial configuration of the system from \n");
    printf("3: number of time steps \n");
    printf("4: size of time step  \n");
    printf("5: graphics is either 1/0 (on/off) \n");
}
