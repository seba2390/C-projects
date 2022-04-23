#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>
#include<float.h>
#include<time.h>
#include<string.h>
#include "functions.h"

typedef struct Particles  {gsl_matrix* positions; gsl_matrix* velocities; gsl_matrix* attributes;} Particles;
typedef struct Simulation {gsl_vector* collision_times; gsl_vector* collision_indices;} Simulation;

//////// FUNCTIONS FOR CCD /////////
double time_to_particle_collision(Particles my_particles, int idx1, int idx2, double dt);
double time_to_wall_collision(double x_i, double y_i, double vx_i, double vy_i, double x_i_1, double y_i_1, double R, double box_width, double dt);
void check_n_step(Particles my_particles, Simulation my_simulation, double box_width, double dt);
void update_velocity_wall(Particles my_particles, int idx, double box_width);
void update_velocity_particle(Particles my_particles, int idx1, int idx2);


/////// SIMULATION FUNCTIONS /////////
void initialize_particles(Particles my_particles, double upper_bound);
void step_particle(Particles my_particles, int idx, double dt);
void step_all_particles(Particles my_particles, double dt);


/////// WRITE OUT /////////
void append_out(gsl_matrix* matrix, char filename[]);
void write_out(gsl_matrix* matrix,char filename[]);


/////// HELP FUNCTIONS /////////
void arr_subtract(double* arr1, double* arr2, double* arr_res);
double inner_prod(double* arr1, double* arr2);
double fRand(double fMin, double fMax);


/////// PRINTING /////////
void print_permut(char s[], gsl_permutation* V);
void print_matrix(char s[], gsl_matrix* M);
void print_vector(char s[], gsl_vector* V);

#endif