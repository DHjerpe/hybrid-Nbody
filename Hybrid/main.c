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

/**
 struct to store information about each particle
 */
typedef struct particle
{
    double x_pos; // x-component (position)
    double y_pos; // y-component (position)
    double mass;  // mass of the particle
    double x_vel; // x-component (velocity)
    double y_vel; // y-component (velocity)
} particle;

void forceCalc(particle * particles, int N,double delta_t, int nsteps) {
void forceCalcVerify(particle * particles, int N,double delta_t, int nsteps);
void generateStars(particle * particles, int size);
void printErrorMsg();
    
/**
 Main function of the program. That simulates the N-body problem.
 Initialises and finalises MPI.
 Params:
    - N: number of stars                                         (argv[1])
    - nsteps: number of time steps                               (argv[2])
    - delta_t: = size of time step                               (argv[3])
    - verify: 0 = NO, 1 = YES                                    (argv[4])
 */
int main(int argc, char** argv) {
    
    if (argc > 4) { // run only if number of arguments is correct
        
        double wtime_end;
        double wtime_start = MPI_Wtime();
        int my_rank;
        //int N = atoi(argv[1])*5;
        int N = atoi(argv[1]);
        int nsteps = atoi(argv[2]);
        int verify = atoi(argv[4]);
        double delta_t = atof(argv[3]);
        
      //  double * p = NULL;
     //   p = (double *)malloc(sizeof(double)*N);
       
        
        particle * particle = NULL;
        particle = (particle *)malloc(sizeof(particle)*N);
        
        generateStars(particles, N);
        
        //generateStars(p, N/5);
        
        
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
    
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
        
        
       // forceCalc(p,N,delta_t,nsteps);
        
        forceCalc(particles, N, delta_t, nsteps);
        
        
        MPI_Finalize();

        
        if (verify && my_rank == 0) {
            
//            double * p_ref = NULL;
//            p_ref = (double *)malloc(sizeof(double)*N);
//            generateStars(p_ref, N/5);
            
            double * particles_ref = NULL;
            particles_ref = (particle *)malloc(sizeof(particle));
            generateStars(particles_ref, N);
            
           // forceCalcVerify(p_ref,N,delta_t,nsteps);
            forceCalcVerify(particles_ref, N, delta_t, nsteps);
            
            for (int i = 0; i<N/5; i++) {
//                assert(p[i] == p_ref[i]);
//                assert(p[i+1] == p_ref[i+1]);
//                assert(p[i+2] == p_ref[i+2]);
//                assert(p[i+3] == p_ref[i+3]);
//                assert(p[i+4] == p_ref[i+4]);
                
                assert(particles.x_pos == particles_ref.x_pos);
                assert(particles.y_pos == particles_ref.y_pos);
                assert(particles.mass == particles_ref.mass);
                assert(particles.x_vel == particles_ref.x_vel);
                assert(particles.y_vel == particles_ref.y_vel);
                
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

//void forceCalc(double * p, int N,double delta_t, int nsteps) {
    void forceCalc(particle * particles, int N,double delta_t, int nsteps) {

   // N = N/5;
    double G = 100.0/N;
    
    int upperLimit = N;
    
    int my_rank;
    int interval = upperLimit/NPROCS;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
    int lowerLimit = my_rank * interval;
    
        if (my_rank <  NPROCS - 1) {
            
            upperLimit = my_rank * interval + interval;
            
        }
 
    
   // int numberToSend = (upperLimit - lowerLimit)*5;
    int numberToSend = (upperLimit - lowerLimit);

    double * particles_buffer = NULL;
    particles_buffer = (double *)malloc(sizeof(double)*N*5);
    memcpy(particles_buffer,particles,sizeof(particle)*N);
        
    
    
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
                    x_diff = (particles[i].x_pos - particles[j].x_pos);
                    y_diff = (particles[i].y_pos - particles[j].y_pos);
                    
                    r = sqrt((x_diff*x_diff) + (y_diff*y_diff));
                    
                    mass = p[j].mass;
                    mass /= ((r+EPSILON)*(r+EPSILON)*(r+EPSILON)); // apply plummer constant
                    
                    x_diff *= mass;
                    y_diff *= mass;
                    
                    sumX += x_diff;
                    sumY += y_diff;
                }
            }
            
            
            a_x = -G*sumX;
            a_y = -G*sumY;
            
            // calculate new veloceties
            u_x = particles[i].x_vel + delta_t*a_x;
            u_y = particles[i].y_vel + delta_t*a_y;
            
            // update new veloceties
            particles_buffer[i].x_vel = u_x;
            particles_buffer[i].y_vel = u_y;
            
            // update new position
            particles_buffer[i].x_pos = particles[i].x_pos + delta_t*u_x;
            particles_buffer[i].y_pos = particles[i].y_pos + delta_t*u_y;
        }
        
        }
        
        MPI_Allgather(MPI_IN_PLACE, numberToSend, MPI_DOUBLE, particles_buffer, numberToSend, MPI_DOUBLE, MPI_COMM_WORLD);
        
        
        // copy the new values to p
        memcpy(particles,particles_buffer,sizeof(particle)*N);
        
    }
    
    //free(p_buffer);
  //  p_buffer = NULL;
        
        free(particles_buffer);
        particles_buffer = NULL;
}

/**
 A vanilla function (serial), that computes a reference 
 to the above parallel implementation
 */
//void forceCalcVerify(double * p, int N,double delta_t, int nsteps) {
void forceCalc(particle * particles, int N,double delta_t, int nsteps) {

   // N = N/5;
    double G = 100.0/N;
    
//    double * p_buffer = NULL;
//    p_buffer = (double *)malloc(sizeof(double)*N*5);
//    memcpy(p_buffer,p,sizeof(double)*N*5);
//    
    
    double * particles_buffer = NULL;
    particles_buffer = (double *)malloc(sizeof(double)*N*5);
    memcpy(particles_buffer,particles,sizeof(particle)*N);
    
    
    
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
                        x_diff = (particles[i].x_pos - particles[j].x_pos);
                        y_diff = (particles[i].y_pos - particles[j].y_pos);
                        
                        r = sqrt((x_diff*x_diff) + (y_diff*y_diff));
                        
                        mass = p[j].mass;
                        mass /= ((r+EPSILON)*(r+EPSILON)*(r+EPSILON)); // apply plummer constant
                        
                        x_diff *= mass;
                        y_diff *= mass;
                        
                        sumX += x_diff;
                        sumY += y_diff;
                    }
                }
            
                a_x = -G*sumX;
                a_y = -G*sumY;
                
                // calculate new veloceties
                u_x = particles[i].x_vel + delta_t*a_x;
                u_y = particles[i].y_vel + delta_t*a_y;
                
                // update new veloceties
                particles_buffer[i].x_vel = u_x;
                particles_buffer[i].y_vel = u_y;
                
                // update new position
                particles_buffer[i].x_pos = particles[i].x_pos + delta_t*u_x;
                particles_buffer[i].y_pos = particles[i].y_pos + delta_t*u_y;
            }
    
        // copy the new values to p
       // memcpy(p,p_buffer,sizeof(double)*N*5);
         memcpy(p,p_buffer,sizeof(particle)*N);
    }
    
    //free(p_buffer);
    //p_buffer = NULL;
    
    free(particles_buffer);
    particles_buffer = NULL;
    
}



/**
 Generate particles (stars) at random position, with random masses,
 random positions, and random velocities
 
 Params:
    - p: vector containing information about all stars/particles in the system
 p contains: 
 - x position
 - y position
 - mass
 - x velocity
 - y velocity
 
    - N: the total number of stars/particles
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


/**
 Function that prints error message if 
 the number of input arguments (from the command line)
 to the program is wrong
 */
void printErrorMsg() {
    
    printf("Error, wrong number of input arguments.\nThe input arguments should have the following form:  \n");
    printf("1: The number of stars \n");
    printf("2: The number of time steps \n");
    printf("3: The size of time step  \n");
    printf("4: Verify ON/OFF (1/0)  \n");
}

