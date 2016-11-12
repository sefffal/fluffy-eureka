#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "vectors.h"

#define GRID_SIZE 100
#define PARTICLE_COUNT 10

/* Define an alias, Grid3, for a 3D array of fixed size. */
typedef double Grid3[GRID_SIZE][GRID_SIZE][GRID_SIZE];

/* Create a struct to represent a particle, with position and velocity. */
typedef struct {
    Vector3 position;
    Vector3 velocity;
} Particle;

/* Define an alias, ParticleArray, for a fixed array of particles. */
/* typedef Particle ParticleArray[PARTICLE_COUNT]; */

void initialize(Particle[], double*, double*);
void integrate(Particle[], double, double);
void output_positions(Particle[], double);

int main(/*int arc, char** argv*/) {

    Particle particles[PARTICLE_COUNT]; 
    double dt;
    double maxt;
    
    initialize(particles, &dt, &maxt);
    output_positions(particles, 1.0);

    integrate(particles, dt, maxt);

    output_positions(particles, maxt);

    return 0;
}

void initialize(Particle particles[], double* dt, double* tmax) {

    *dt = 0.1;
    *tmax = 1.0;

    /* Allocate particles, that's the hard bit. */
    for (int i=0; i<PARTICLE_COUNT; i++) {
        
        particles[i].position.x = i*50.0;
        particles[i].position.y = i*50.0;
        particles[i].position.z = i*50.0;

        particles[i].velocity.x = i*50.0;
        particles[i].velocity.y = i*50.0;
        particles[i].velocity.z = i*50.0;  
    }

    return;
}

void output_positions(Particle particles[], double t) {
    for (int i=0; i<PARTICLE_COUNT; i++) {
        printf("%.18g\t%.18g\t%.18g\t%.18g\n", t, particles[i].position.x, particles[i].position.y, particles[i].position.z);
    }
    printf("\n\n");
}

void integrate(Particle particles[], double dt, double tmax) {

    Particle old[PARTICLE_COUNT];
    Particle new[PARTICLE_COUNT];
    
    /* Copy initial conditions into particles array: */
    memcpy(new, particles, sizeof(new));

    Vector3 tempPosition[PARTICLE_COUNT];

    double GM = 2*M_PI;
    double softening_constant = 1.0;


    for (double t=0.0; t<tmax; t+=dt) {

        /* Flip old and new temp grids */
        memcpy(particles, old, sizeof(Particle)*PARTICLE_COUNT);
        memcpy(old, new, sizeof(Particle)*PARTICLE_COUNT);
        memcpy(new, particles, sizeof(Particle)*PARTICLE_COUNT);

        // output_positions(old, t);

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
        for (int i=0; i<PARTICLE_COUNT; i++) {

            Vector3 acceleration;
            acceleration.x = 0.0;
            acceleration.y = 0.0;
            acceleration.z = 0.0;

            /* Go through each other particle to add the forces */
            for (int j=0; j<PARTICLE_COUNT; j++) {
                /* Don't add self-forces */
                if (i==j)
                    continue;
                Vector3 difference = vsub(tempPosition[i], tempPosition[j]);
                double r = vmag(difference) + softening_constant;
                double field = -GM/pow(r, 3);

                acceleration.x += field*tempPosition[j].x;
                acceleration.y += field*tempPosition[j].y;
                acceleration.z += field*tempPosition[j].z;
            }

            /* Update the velocity in the new grid */
            new[i].velocity.x = old[i].velocity.x + acceleration.x * dt;          
            new[i].velocity.y = old[i].velocity.y + acceleration.y * dt;          
            new[i].velocity.z = old[i].velocity.z + acceleration.z * dt;

            /* Update position in the new grid */
            new[i].position.x = tempPosition[i].x + new[i].velocity.x * dtover2;
            new[i].position.y = tempPosition[i].y + new[i].velocity.y * dtover2;
            new[i].position.z = tempPosition[i].z + new[i].velocity.z * dtover2;
        }
    }
}


