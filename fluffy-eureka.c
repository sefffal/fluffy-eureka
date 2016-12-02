#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "vectors.h"
//#include "fftw-3.3.5/fftw3.h"
#include <execinfo.h>
#include <signal.h>

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

#define GRID_SIZE 100
#define PARTICLE_COUNT 3000
#define PARTICLE_MASS 1.0
#define RADIUS 1.0

/* Define an alias, Grid3, for a 3D array of fixed size. */
typedef double Grid3[GRID_SIZE][GRID_SIZE][GRID_SIZE];

/* Create a struct to represent a particle, with position and velocity. */
typedef struct {
    Vector3 position;
    Vector3 velocity;
} Particle;

void initialize(Particle[], double*, double*, double*, int, char**);
void integrate(Particle[], double, double, double);
void output_positions(Particle[], double);
void calculate_bariness(Particle[], double[3][3]);
void shift_all_by(Particle[], Vector3);
Vector3 center_of_mass(Particle[]);

int main(int argc, char** argv) {
    signal(SIGSEGV, handler);   // install our handler

    Particle particles[PARTICLE_COUNT]; 
    double dt;
    double maxt;
    double grid_space;

    initialize(particles, &dt, &maxt, &grid_space, argc, argv);
    integrate(particles, dt, maxt, grid_space);

    shift_all_by(particles, smulv(-1, center_of_mass(particles)));
    output_positions(particles, maxt);

    return 0;
}

void initialize(Particle particles[], double* dt, double* tmax, double* grid_space, int argc, char** argv) {

    // Check if we got all required arguments.
    // We need at least four, but the 0th element is the name of the program.
    if (argc < 1) {
        // If we didn't receive enough, print error message to STDERR and exit with error.
        fprintf(stderr, "Error: Please provide arguments for fraction_to_flip\n");
        exit(1);
    }

    
    double G = 2*M_PI;
    // Parse the command line arguments and use them as our initial conditions.
    double fraction_to_flip = atof(argv[1]);

    /* Calculate the characteristic time scale of the system -- use to find dt and tmax */
    double t_charac = RADIUS/sqrt(G*PARTICLE_MASS*PARTICLE_COUNT/RADIUS);
    *dt = 0.007 * t_charac;
    *tmax = t_charac*60;
    *grid_space = RADIUS*8/GRID_SIZE;

    fprintf(stderr, "dt: %lf\n", *dt);    
    fprintf(stderr, "tmax: %lf\n", *tmax);
    fprintf(stderr, "fraction_to_flip: %lf\n", fraction_to_flip);    

    srand( (unsigned) time(NULL) * getpid());
    // srand((unsigned)time(NULL));
    // srand(100);

    /* Distribute particles uniformly in a unit sphere */ 
    int i;
    i=0;
    while (i<=PARTICLE_COUNT) {
        
        particles[i].position.x = (2.0*RADIUS*(double)rand()/(double)RAND_MAX)-RADIUS;
        particles[i].position.y = (2.0*RADIUS*(double)rand()/(double)RAND_MAX)-RADIUS;
        particles[i].position.z = (2.0*RADIUS*(double)rand()/(double)RAND_MAX)-RADIUS;

        /* Reject particles outside of a unit sphere */
        if (vmag(particles[i].position) < RADIUS && vmag(particles[i].position) > RADIUS*0.01) {
            i++;
        }
    }

    /* Give each particle an oribital velocity */
    i=0;
    while (i<PARTICLE_COUNT) {
        
        double rmag = vmag(particles[i].position);
        // double Mr = (4.0/3.0) * M_PI * pow(rmag, 3) * density;

        /* Calculate the mass within this sphere of this radius */
        double Mr = 0.0;
        for (int j=0; j<PARTICLE_COUNT; j++) {
            if (vmag(particles[j].position) <= rmag) {
                Mr += PARTICLE_MASS;
            }
        }

        Vector3 r = particles[i].position;        

        /* Rotate into a primed coordinate system where
         * the r vector is along the z' axis.
         * `angle` is then the angle about the z'-axis,
         * which we project into y and x components using
         * the sines and cosines of the angle.
         * Finally, we rotate back into the unprimed
         * coordinate system.
         */
        double speed = sqrt(G*Mr/rmag);

        /* Pick a random angle between 0 and 2*pi */
        double angle = 2.0*M_PI*(double)rand()/(double)RAND_MAX;

        double phi = -atan2(r.y, r.x);
        double theta = -acos(r.z/rmag);

        Vector3 alpha = {
            .x = cos(phi)*cos(theta),
            .y = -cos(theta)*sin(phi),
            .z = sin(theta)
        };

        Vector3 beta = {
            .x = sin(phi),
            .y = cos(phi),
            .z = 0
        };

        // Vector3 gamma = {
        //     .x = -sin(theta)*cos(phi),
        //     .y = sin(theta)*sin(phi),
        //     .z = cos(theta)
        // };

        Vector3 vprime = {
            .x = cos(angle)*speed,
            .y = sin(angle)*speed,
            .z = 0.0
        };

        
        particles[i].velocity = vadd(smulv(vprime.x, alpha), smulv(vprime.y, beta));

        // Check if angular momentum is negative
        if (vcross(r, particles[i].velocity).z <0) {
            // Flip some fraction of orbits from retrograde to prograde
            if ((double)rand()/(double)RAND_MAX < fraction_to_flip) {
                particles[i].velocity.x *= -1;
                particles[i].velocity.y *= -1;
            }
        }

        i+=1;
    }


    /* Set the total translational momentum to zero so that it stays put. */
    Vector3 net_velocity = {.x=0.0, .y=0.0, .z=0.0};
    for (i=0; i<PARTICLE_COUNT; i++) {
        net_velocity = vadd(net_velocity, particles[i].velocity);
    }
    Vector3 ave_velocity = smulv(1/PARTICLE_COUNT, net_velocity);
    fprintf(stderr, "Initial ave. velocity: %lf %lf %lf\n", ave_velocity.x, ave_velocity.y, ave_velocity.z);
    for (i=0; i<PARTICLE_COUNT; i++) {
         particles[i].velocity = vsub(particles[i].velocity, ave_velocity);
    }

    Vector3 cent = center_of_mass(particles);
    fprintf(stderr, "initial center of mass: %lf %lf %lf\n", cent.x, cent.y, cent.z);
    
    /* Shift the origin so that everything is positive.
     * This will make it easier to match up the grid.
     * A shift of 8 times the radius will put everything firmly positive.
     * The grid can then range from 0 to 16 times the radius.
     */
    Vector3 shift = {.x=8*RADIUS,.y=8*RADIUS,.z=8*RADIUS};
    shift_all_by(particles, shift);

    return;
}

void output_positions(Particle particles[], double t) {
    double inertia_tensor[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    calculate_bariness(particles, inertia_tensor);
    for (int i=0; i<PARTICLE_COUNT; i++) {
        printf("%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n",
             t,
             particles[i].position.x,
             particles[i].position.y,
             particles[i].position.z,
             inertia_tensor[0][0],
             inertia_tensor[0][1],
             inertia_tensor[0][2],
             inertia_tensor[1][0],
             inertia_tensor[1][1],
             inertia_tensor[1][2],
             inertia_tensor[2][0],
             inertia_tensor[2][1],
             inertia_tensor[2][2]
            );
    }
    printf("\n\n");
}

void integrate(Particle particles[], double dt, double tmax, double grid_space) {

    Particle old[PARTICLE_COUNT];
    Particle new[PARTICLE_COUNT];
    
    /* Copy initial conditions into particles array: */
    memcpy(new, particles, sizeof(new));

    Vector3 tempPosition[PARTICLE_COUNT];

    double G = 2*M_PI;
    double softening_constant = 0.02*RADIUS * pow(PARTICLE_COUNT, -1.0/3.0); /* Ave particle spacing, prop. to Radius * N^1/3 */
    fprintf(stderr, "Softening: %lf\n", softening_constant);

    int i=0;
    for (double t=0.0; t<tmax; t+=dt, i++) {

        /* Flip old and new temp grids */
        memcpy(particles, old, sizeof(Particle)*PARTICLE_COUNT);
        memcpy(old, new, sizeof(Particle)*PARTICLE_COUNT);
        memcpy(new, particles, sizeof(Particle)*PARTICLE_COUNT);
    
        if (i%25==0) {
            output_positions(old, t);
            fprintf(stderr, "%2.2f%% complete\n", t/tmax*100);
        }
        /* For the half step forward, just calculate once */
        double dtover2 = dt/2;
        
        /* Calculate the position of each particle one half step forward */
        for (int i=0; i<PARTICLE_COUNT; i++) {
            tempPosition[i].x = old[i].position.x + old[i].velocity.x * dtover2;
            tempPosition[i].y = old[i].position.y + old[i].velocity.y * dtover2;
            tempPosition[i].z = old[i].position.z + old[i].velocity.z * dtover2;
        }

        /* Calculate the acceleration of each particle for each other particle,
         * using the temporary positions calculated earlier.
         */
        
        
        Vector3 accelerations[PARTICLE_COUNT];
        /* Zero out initial accelerations */
        for (int i=0; i<PARTICLE_COUNT; i++) {
            accelerations[i].x = 0.0;
            accelerations[i].y = 0.0;
            accelerations[i].z = 0.0;
        }

        /* Go through each combination of particles to add the forces */
        for (int i=0; i<PARTICLE_COUNT; i++) {
            for (int j=i+1; j<PARTICLE_COUNT; j++) {
                Vector3 difference = vsub(tempPosition[i], tempPosition[j]);
                double r = vmag(difference) + softening_constant;
                double field = -G * pow(PARTICLE_MASS, 2) /pow(r, 3);

                accelerations[i].x += field*difference.x;
                accelerations[i].y += field*difference.y;
                accelerations[i].z += field*difference.z;
                accelerations[j] = vadd(accelerations[j], smulv(-1, accelerations[i]));
            } /* End of inner (j) particle loop */
        }

        /* Now that we have all the accelerations, advance one timestep */
        for (int i=0; i<PARTICLE_COUNT; i++) {
            /* Update the velocity in the new grid */
            new[i].velocity.x = old[i].velocity.x + accelerations[i].x * dt;          
            new[i].velocity.y = old[i].velocity.y + accelerations[i].y * dt;          
            new[i].velocity.z = old[i].velocity.z + accelerations[i].z * dt;

            /* Update position in the new grid */
            new[i].position.x = tempPosition[i].x + new[i].velocity.x * dtover2;
            new[i].position.y = tempPosition[i].y + new[i].velocity.y * dtover2;
            new[i].position.z = tempPosition[i].z + new[i].velocity.z * dtover2;

        } /* End outer (i) particle loop */
    } /* End t+=dt time loop */
}

// void calculate_potential(Particle particles[], Grid3 potential, double grid_space, double G) {
//     double c = 4*M_PI*G;
//     Grid3 density;
//     int i;
//     int j;
//     int k;


//     // Zero out density grid
//     for (i=0; i<GRID_SIZE; i++) {
//         for (j=0; j<GRID_SIZE; j++) {
//             for (k=0; k<GRID_SIZE; k++) {
//                 density[i][j][k] = 0;
//             }
//         }
//     }

//     // Distribute mass to grid points using nearest grid point method

//     for (i=0; i<PARTICLE_COUNT; i++) {
//         Vector3 scaled = smulv(1/grid_space, particles[i].position);
//         // there are eight grid points it could be close to
//         double dist1 = vmag(vsub(scaled, ))
//         particles[i].position
//     }

//     for (i=0; i<GRID_SIZE; i++) {
//         for (j=0; j<GRID_SIZE; j++) {
//             for (k=0; k<GRID_SIZE; k++) {
                
//                 if 
                
//                 density[i][j][k] = 0;
//             }
//         }
//     }


//     // Calculate potential at each grid point using direct summation
    
//     //
// }


void calculate_bariness(Particle particles[], double inertia_tensor[3][3]) {

    /* Shift the galaxy center back to the origin for this calculation */
    Vector3 shift = smulv(-1, center_of_mass(particles));
    shift_all_by(particles, shift);

    for (int i=0; i<PARTICLE_COUNT; i++) {
        inertia_tensor[0][0] += particles[i].position.x*particles[i].position.x;
        inertia_tensor[0][1] += particles[i].position.x*particles[i].position.y;
        inertia_tensor[0][2] += particles[i].position.x*particles[i].position.z;

        inertia_tensor[1][0] += particles[i].position.y*particles[i].position.x;
        inertia_tensor[1][1] += particles[i].position.y*particles[i].position.y;
        inertia_tensor[1][2] += particles[i].position.y*particles[i].position.z;

        inertia_tensor[2][0] += particles[i].position.z*particles[i].position.x;
        inertia_tensor[2][1] += particles[i].position.z*particles[i].position.y;
        inertia_tensor[2][2] += particles[i].position.z*particles[i].position.z;
    }

    // fprintf(stderr, "{{%lf,%lf,%lf},\n", inertia_tensor[0][0], inertia_tensor[0][1], inertia_tensor[0][2]);
    // fprintf(stderr, " {%lf,%lf,%lf},\n", inertia_tensor[1][0], inertia_tensor[1][1], inertia_tensor[1][2]);
    // fprintf(stderr, " {%lf,%lf,%lf}}\n", inertia_tensor[2][0], inertia_tensor[2][1], inertia_tensor[2][2]);

    /* Shift the galaxy center back to where it was */
    shift_all_by(particles, smulv(-1, shift));

}




void shift_all_by(Particle particles[], Vector3 shift) {
    for (int i=0; i<PARTICLE_COUNT; i++) {
        particles[i].position = vadd(particles[i].position, shift);
    }
    fprintf(stderr, "Shift all by: %lf %lf %lf\n", shift.x, shift.y, shift.z);    
}


Vector3 center_of_mass(Particle particles[]) {
    Vector3 position;
    for (int i=0; i<PARTICLE_COUNT; i++) {
        position = vadd(particles[i].position, position);
    }
    return smulv(1.0/PARTICLE_COUNT, position);
}
