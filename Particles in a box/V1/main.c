#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<float.h>
#include "functions.h"
#include <time.h>


int main()
{
    srand(time(NULL)); // Seed init for random numbers 

    int nr_particles = 2500;
    int nr_steps     = 500; 
    double box_width = 20.;
    double dt        = 0.075;

    gsl_matrix* positions    = gsl_matrix_alloc(nr_particles,2);
    gsl_matrix* velocities   = gsl_matrix_alloc(nr_particles,2);
    gsl_matrix* attributes   = gsl_matrix_alloc(nr_particles,2);

    Particles my_particles;
    my_particles.attributes = attributes;
    my_particles.positions  = positions;
    my_particles.velocities = velocities;

    gsl_vector* times_vector       = gsl_vector_alloc(2*nr_particles); // Enough entries for all collisions during 1 dt 
    gsl_vector* indices_vector     = gsl_vector_alloc(2*nr_particles); // Enough entries for all collisions during 1 dt
    gsl_vector* sorted_coordinates = gsl_vector_alloc(nr_particles);
    
    Simulation my_simulation;
    my_simulation.times_vector       = times_vector;
    my_simulation.indices_vector     = indices_vector;
    my_simulation.sorted_coordinates = sorted_coordinates;

    initialize_particles(my_particles, box_width);

    write_out(my_particles.positions,"coordinates.txt");
    write_out(my_particles.velocities, "velocities.txt");
    write_out(my_particles.attributes, "attributes.txt");
    for(int i = 0; i < nr_steps; i++)
    {
        check_n_step(my_particles,my_simulation,box_width,dt);
        append_out(my_particles.positions,"coordinates.txt");
        append_out(my_particles.positions,"velocities.txt");
        printf("i = %d\n",i+1);

    }


    gsl_matrix_free(positions); 
    gsl_matrix_free(velocities);
    gsl_matrix_free(attributes);

    gsl_vector_free(times_vector);
    gsl_vector_free(indices_vector);
    gsl_vector_free(sorted_coordinates);

    system("python Animate.py -particles 2500 -steps 500 -bw 20.0");
}
