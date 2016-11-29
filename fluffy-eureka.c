#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "vectors.h"

#define GRID_SIZE 100
#define PARTICLE_COUNT 2000
#define PARTICLE_MASS 1.0
#define RADIUS 1.0

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
    // output_positions(particles, 1.0);
    // return 0;

    integrate(particles, dt, maxt);

    // output_positions(particles, maxt);

    return 0;
}

void initialize(Particle particles[], double* dt, double* tmax) {
    
    double G = 2*M_PI;
    double fraction_to_flip = 1.0;

    // Calculate the characteristic time scale of the system -- use to find dt and tmax
    double t_charac = RADIUS/sqrt(G*PARTICLE_MASS*PARTICLE_COUNT/RADIUS);
    *dt = 0.01 * t_charac;
    *tmax = t_charac*30;//100;

    fprintf(stderr, "dt: %lf\n", *dt);    
    fprintf(stderr, "tmax: %lf\n", *tmax);
    fprintf(stderr, "fraction_to_flip: %lf\n", fraction_to_flip);    


    // Naive density calculation, assuming smooth distribution
    // double density = PARTICLE_COUNT*PARTICLE_MASS/((4.0/3.0)*M_PI*pow(RADIUS, 3));

    // srand((unsigned)time(NULL));
    srand(100);

    /* Distribute particles uniformly in a unit sphere */ 
    int i;
    i=0;
    while (i<=PARTICLE_COUNT) {
        
        particles[i].position.x = (2.0*(double)rand()/(double)RAND_MAX)-1.0;
        particles[i].position.y = (2.0*(double)rand()/(double)RAND_MAX)-1.0;
        particles[i].position.z = (2.0*(double)rand()/(double)RAND_MAX)-1.0;

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
            // Flip some percentage of orbits from retrograde to prograde
            if ((double)rand()/(double)RAND_MAX < fraction_to_flip) {
                particles[i].velocity.x *= -1;
                particles[i].velocity.y *= -1;
            }
        }

        i+=1;
    }

    return;
}

void output_positions(Particle particles[], double t) {
    for (int i=0; i<PARTICLE_COUNT; i++) {
        printf("%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n",
             t,
             particles[i].position.x,
             particles[i].position.y,
             particles[i].position.z,
             particles[i].velocity.x,
             particles[i].velocity.y,
             particles[i].velocity.z
            );
    }
    // printf("\n\n");
}

void integrate(Particle particles[], double dt, double tmax) {

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
                double field = -G * pow(PARTICLE_MASS, 2) /pow(r, 3);

                acceleration.x += field*difference.x;
                acceleration.y += field*difference.y;
                acceleration.z += field*difference.z;
            } /* End of inner (j) particle loop */

            /* Update the velocity in the new grid */
            new[i].velocity.x = old[i].velocity.x + acceleration.x * dt;          
            new[i].velocity.y = old[i].velocity.y + acceleration.y * dt;          
            new[i].velocity.z = old[i].velocity.z + acceleration.z * dt;

            /* Update position in the new grid */
            new[i].position.x = tempPosition[i].x + new[i].velocity.x * dtover2;
            new[i].position.y = tempPosition[i].y + new[i].velocity.y * dtover2;
            new[i].position.z = tempPosition[i].z + new[i].velocity.z * dtover2;

        } /* End outer (i) particle loop */
    } /* End t+=dt time loop */
}


