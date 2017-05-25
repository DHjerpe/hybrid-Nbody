#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
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

void forceCalc(particle * particles, int N,double delta_t, int nsteps);
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
        int N = atoi(argv[1]);
        int nsteps = atoi(argv[2]);
        int verify = atoi(argv[4]);
        double delta_t = atof(argv[3]);
        
        
        particle * particles = NULL;
        particles = (particle *)malloc(sizeof(particle)*N);
        
        generateStars(particles, N);
        
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
    
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
        
        forceCalc(particles, N, delta_t, nsteps);
        
    
        MPI_Finalize();

        
        if (verify && my_rank == 0) {
            
            particle * particles_ref = NULL;
            particles_ref = (particle *)malloc(sizeof(particle)*N);
            generateStars(particles_ref, N);
            
            
            forceCalcVerify(particles_ref, N, delta_t, nsteps);
            
            for (int i = 0; i<N/5; i++) {
                
                assert(particles[i].x_pos == particles_ref[i].x_pos);
                assert(particles[i].y_pos == particles_ref[i].y_pos);
                assert(particles[i].mass == particles_ref[i].mass);
                assert(particles[i].x_vel == particles_ref[i].x_vel);
                assert(particles[i].y_vel == particles_ref[i].y_vel);
            }
            printf("Verify OK! \n");
            free(particles_ref);
            particles_ref = NULL;
        }
        
        
        free(particles);
        particles = NULL;
        
        wtime_end = MPI_Wtime();
        
        if (my_rank == 0)
            printf("Total time: %.4f \n", wtime_end - wtime_start);
        
    } else {printErrorMsg();}
    
    
    return 0;
}

/**
 Hybrid function utilising both MPI and OpenMP.
 Computes the new position and velocities of the system
 
 Params:
 - particles: array of structs, each struct containing
 information about a particle
 - N: the total number of stars/particles
 - delta_t: size of the time step
 - nsteps: number of time steps
 */
    void forceCalc(particle * particles, int N,double delta_t, int nsteps) {

    double G = 100.0/N;
    
    int upperLimit = N;
    
    int my_rank;
    int interval = upperLimit/NPROCS;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
    int lowerLimit = my_rank * interval;
    
        if (my_rank <  NPROCS - 1) {
            
            upperLimit = my_rank * interval + interval;
            
        }
 
    int numberToSend = (upperLimit - lowerLimit);

    particle * particles_buffer = NULL;
    particles_buffer = (particle *)malloc(sizeof(particle)*N);
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
                    
                    mass = particles[j].mass;
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
        free(particles_buffer);
        particles_buffer = NULL;
}

/**
 A vanilla function (serial), that computes a reference 
 to the above parallel implementation
 
 Params:
    - particles: array of structs, each struct containing
    information about a particle
    - N: the total number of stars/particles
    - delta_t: size of the time step
    - nsteps: number of time steps
 */
void forceCalcVerify(particle * particles, int N,double delta_t, int nsteps) {

    double G = 100.0/N;
    
    particle * particles_buffer = NULL;
    particles_buffer = (particle *)malloc(sizeof(particle)*N);
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
                        
                        mass = particles[j].mass;
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
         memcpy(particles,particles_buffer,sizeof(particle)*N);
    }

    free(particles_buffer);
    particles_buffer = NULL;
}



/**
 Generate particles (stars) at random position, with random masses,
 random positions, and random velocities
 
 Params:
    - particles: array of structs, each struct containing
    information about a particle
    - N: the total number of stars/particles
 */
void generateStars(particle * particles, int N) {
    
    int randLoc = 256;
    srand(randLoc); // choose a locatoin in the "random numbers" array
    
    for (int i = 0; i<N; i++) {
        
        particles[i].x_pos = (double)rand() / (double)((unsigned)RAND_MAX + 1); // pos x
        ++randLoc; srand(randLoc);
        particles[i].y_pos = (double)rand() / (double)((unsigned)RAND_MAX + 1); // pos y
        ++randLoc; srand(randLoc);
        particles[i].mass = (double)rand() / (double)((unsigned)RAND_MAX + 1); // mass
        ++randLoc; srand(randLoc);
        particles[i].x_vel = (double)rand() / (double)((unsigned)RAND_MAX + 1); // velocity x
        ++randLoc; srand(randLoc);
        particles[i].y_vel = (double)rand() / (double)((unsigned)RAND_MAX + 1); // velocity y
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

