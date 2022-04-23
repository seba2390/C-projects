#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdlib.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<float.h>
#include<time.h>
#include<string.h>
#include "functions.h"



/////// SIMULATION FUNCTIONS FOR CCD /////////

/* Function for computing time until particle-particle collision */
double time_to_particle_collision(Particles my_particles, int idx1, int idx2, double dt)
{
	double r1[2] = {gsl_matrix_get(my_particles.positions,idx1,0), gsl_matrix_get(my_particles.positions,idx1,1)};
	double r2[2] = {gsl_matrix_get(my_particles.positions,idx2,0), gsl_matrix_get(my_particles.positions,idx2,1)};
	double v1[2] = {gsl_matrix_get(my_particles.velocities,idx1,0), gsl_matrix_get(my_particles.velocities,idx1,1)};
	double v2[2] = {gsl_matrix_get(my_particles.velocities,idx2,0), gsl_matrix_get(my_particles.velocities,idx2,1)};
	double c1[2]; arr_subtract(v2,v1,c1); // v2 - v1
	double c2[2]; arr_subtract(r2,r1,c2); // r2 - r1
	double R1 = gsl_matrix_get(my_particles.attributes,idx1,1);
	double R2 = gsl_matrix_get(my_particles.attributes,idx2,1);

	double A = inner_prod(c1,c1);
	double B = 2.*inner_prod(c1,c2);
	double C = inner_prod(c2,c2)-pow((R1+R2),2.);
	double D = pow(B,2.)-4.*A*C;
	//printf("D = %g\n",D);
	double time1 = (-B+sqrt(D))/(2.*A);
	double time2 = (-B-sqrt(D))/(2.*A);
	double final_time = time1;
	
	if(final_time > time2)
	{
		final_time = time2;
	}
	//printf("time1, time2, dt: %g, %g, %g\n",time1, time2, dt);
	//printf("final_time = %g \n",final_time);
	return final_time;

}

/* Function for computing time until particle-wall collision */
double time_to_wall_collision(double x_i, double y_i, double vx_i, double vy_i, double x_i_1, double y_i_1, double R, double box_width, double dt)
{
	double time;
	if(x_i_1 - R <= 0){time = (R-x_i)/(vx_i);} 				        // Left
	if(y_i_1 - R <= 0){time = (R-y_i)/(vy_i);} 			            // Bottom
	if(x_i_1 + R >= box_width){time = -(R+x_i-box_width)/(vx_i);}   // Right
	if(y_i_1 + R >= box_width){time = -(R+y_i-box_width)/(vy_i);}   // Top
	//printf("Calculated wall time: %g\n",time);
	if(time <= 0){printf("Time is %g\n",time);printf("Radius, xnext, ynext = %g, %g, %g\n",R,x_i_1,y_i_1);printf("Radius, x, y = %g, %g, %g\n",R,x_i,y_i);}
	assert(time >= 0);
	return time;
}

/* Function for advancing a single particle dt */
void step_particle(Particles my_particles, int idx, double dt)
{
	// Getting
	double x  = gsl_matrix_get(my_particles.positions,idx,0),   y = gsl_matrix_get(my_particles.positions,idx,1);
	double vx = gsl_matrix_get(my_particles.velocities,idx,0), vy = gsl_matrix_get(my_particles.velocities,idx,1);

	// Stepping
	x = x + vx * dt;
	y = y + vy * dt;

	// Setting
	gsl_matrix_set(my_particles.positions,idx,0,x);
	gsl_matrix_set(my_particles.positions,idx,1,y);
}

/* Function for advancing all particles dt */
void step_all_particles(Particles my_particles, double dt)
{
	int nr_particls = (*my_particles.positions).size1;

	for(int idx = 0; idx < nr_particls; idx++)
	{
		// Getting
		double x  = gsl_matrix_get(my_particles.positions,idx,0),   y = gsl_matrix_get(my_particles.positions,idx,1);
		double vx = gsl_matrix_get(my_particles.velocities,idx,0), vy = gsl_matrix_get(my_particles.velocities,idx,1);

		// Stepping
		x = x + vx * dt;
		y = y + vy * dt;

		// Setting
		gsl_matrix_set(my_particles.positions,idx,0,x);
		gsl_matrix_set(my_particles.positions,idx,1,y);
	}
}

/* Function for updating the velocity vector of a particle after wall collision */
void update_velocity_wall(Particles my_particles, int idx, double box_width)
{
	// Getting
	double vx = gsl_matrix_get(my_particles.velocities,idx,0);
	double vy = gsl_matrix_get(my_particles.velocities,idx,1);
	double x = gsl_matrix_get(my_particles.positions,idx,0);
	double y = gsl_matrix_get(my_particles.positions,idx,1);
	double R = gsl_matrix_get(my_particles.attributes,idx,1);

	// Setting
	if(x + R >= box_width || x - R <= 0){gsl_matrix_set(my_particles.velocities,idx,0,-vx);}
	if(y + R >= box_width || y - R <= 0){gsl_matrix_set(my_particles.velocities,idx,1,-vy);}
}

/* Function for updating velocity after particle collision */
void update_velocity_particle(Particles my_particles, int idx1, int idx2)
{
	double m1 = gsl_matrix_get(my_particles.attributes,idx1,0), m2 = gsl_matrix_get(my_particles.attributes,idx2,0);
	double r1[2] = {gsl_matrix_get(my_particles.positions,idx1,0), gsl_matrix_get(my_particles.positions,idx1,1)};
	double r2[2] = {gsl_matrix_get(my_particles.positions,idx2,0), gsl_matrix_get(my_particles.positions,idx2,1)};
	double v1[2] = {gsl_matrix_get(my_particles.velocities,idx1,0), gsl_matrix_get(my_particles.velocities,idx1,1)};
	double v2[2] = {gsl_matrix_get(my_particles.velocities,idx2,0), gsl_matrix_get(my_particles.velocities,idx2,1)};
	double r1_r2[2]; arr_subtract(r1,r2,r1_r2);
	double r2_r1[2]; arr_subtract(r2,r1,r2_r1);
	double v1_v2[2]; arr_subtract(v1,v2,v1_v2);
	double v2_v1[2]; arr_subtract(v2,v1,v2_v1);

	double mfac1  = 2.*m2/(m1+m2), mfac2 = 2.*m1/(m1+m2);

	double v1fx = v1[0]-2.*m2/(m1+m2)*inner_prod(v1_v2,r1_r2)/inner_prod(r1_r2,r1_r2)*r1_r2[0];
	double v1fy = v1[1]-2.*m2/(m1+m2)*inner_prod(v1_v2,r1_r2)/inner_prod(r1_r2,r1_r2)*r1_r2[1];

	double v2fx = v2[0]-2.*m1/(m1+m2)*inner_prod(v2_v1,r2_r1)/inner_prod(r2_r1,r2_r1)*r2_r1[0];
	double v2fy = v2[1]-2.*m1/(m1+m2)*inner_prod(v2_v1,r2_r1)/inner_prod(r2_r1,r2_r1)*r2_r1[1];

	gsl_matrix_set(my_particles.velocities,idx1,0,v1fx);
	gsl_matrix_set(my_particles.velocities,idx1,1,v1fy);

	gsl_matrix_set(my_particles.velocities,idx2,0,v2fx);
	gsl_matrix_set(my_particles.velocities,idx2,1,v2fy);
}

/* Function for doing all */
void check_n_step(Particles my_particles, Simulation my_simulation, double box_width, double dt)
{	 
	if(dt > 0)
	{
		gsl_vector_set_all(my_simulation.collision_times, 2.*dt);
		int nr_collisions    = 0;
		int nr_particles     = (*my_particles.positions).size1;

		/* CHECKING FOR WALL COLLISIONS AND STORING TIME 2 COLLISION */
		int wall_index;
		double x,y,vx,vy,x_next,y_next,time_2_wall, R;
		for(int particle_i = 0; particle_i < nr_particles; particle_i++)
		{
			R = gsl_matrix_get(my_particles.attributes,particle_i,1);
			x = gsl_matrix_get(my_particles.positions,particle_i,0);   y = gsl_matrix_get(my_particles.positions,particle_i,1);
			vx = gsl_matrix_get(my_particles.velocities,particle_i,0); vy = gsl_matrix_get(my_particles.velocities,particle_i,1);
			x_next = x + vx * dt; y_next = y + vy * dt;
			if(x_next - R <= 0 || x_next + R >= box_width || y_next - R <= 0 || y_next + R >= box_width)
			{
				time_2_wall = time_to_wall_collision(x,y,vx,vy,x_next,y_next,R,box_width,dt);
				if(time_2_wall > 0)
				{
					gsl_vector_set(my_simulation.collision_times,nr_collisions,time_2_wall);
					gsl_vector_set(my_simulation.collision_indices,nr_collisions,particle_i);
					nr_collisions++;
				}
				if(time_2_wall == 0)
				{
					update_velocity_wall(my_particles,particle_i,box_width); // Special case of particle is exactly at wall
				}
			}
		}
		
		/* BRUTE FORCE METHOD; CHECKING DISTANCE BETWEEN ALL PARTICLES O(N^2) */
		int particle_index_1, particle_index_2;
		double x1, x2, y1, y2, vx1, vy1, vx2, vy2, R1, R2, time_2_particle;
		double x1_next, y1_next, x2_next, y2_next, remaining_time;
		for(int particle1 = 0; particle1 < nr_particles-1; particle1++)
		{
			R1 = gsl_matrix_get(my_particles.attributes,particle1,1);
			x1  = gsl_matrix_get(my_particles.positions,particle1,0);  y1 = gsl_matrix_get(my_particles.positions,particle1,1);
			vx1 = gsl_matrix_get(my_particles.velocities,particle1,0);vy1 = gsl_matrix_get(my_particles.velocities,particle1,1);
			x1_next = x1 + vx1 * dt; y1_next = y1 + vy1 * dt;

			for(int particle2 = particle1+1; particle2 < nr_particles; particle2++)
			{
				R2 = gsl_matrix_get(my_particles.attributes,particle2,1);
				x2  = gsl_matrix_get(my_particles.positions,particle2,0);  y2 = gsl_matrix_get(my_particles.positions,particle2,1);
				vx2 = gsl_matrix_get(my_particles.velocities,particle2,0);vy2 = gsl_matrix_get(my_particles.velocities,particle2,1);
				x2_next = x2 + vx2 * dt; y2_next = y2 + vy2 * dt;

				double dist = sqrt((x2_next-x1_next)*(x2_next-x1_next)+(y2_next-y1_next)*(y2_next-y1_next));
				if(dist <= R1+R2)
				{
					time_2_particle = time_to_particle_collision(my_particles, particle1, particle2, dt);
					if(time_2_particle > 0)
					{
							gsl_vector_set(my_simulation.collision_times,nr_collisions,time_2_particle);
							gsl_vector_set(my_simulation.collision_indices,nr_collisions,particle1);
							nr_collisions++;

							gsl_vector_set(my_simulation.collision_times,nr_collisions,time_2_particle);
							gsl_vector_set(my_simulation.collision_indices,nr_collisions,particle2);
							nr_collisions++;
					}
					if (time_2_particle == 0)
					{
						update_velocity_particle(my_particles,particle1,particle1);
					}
							
				}
			}
		}


		if(nr_collisions > 0)
		{ 
			gsl_sort_vector2(my_simulation.collision_times, my_simulation.collision_indices);
			double smallest_time = gsl_vector_get(my_simulation.collision_times,0);
			double remaining_time = dt - smallest_time;
			if(smallest_time == gsl_vector_get(my_simulation.collision_times,1))               // Particle -> Particle collision
			{ 	
				double particle_idx1 = gsl_vector_get(my_simulation.collision_indices,0);
				double particle_idx2 = gsl_vector_get(my_simulation.collision_indices,1);
				step_all_particles(my_particles,smallest_time);
				update_velocity_particle(my_particles, particle_idx1, particle_idx2);	       // Updating velocity of colliding particles       
				check_n_step(my_particles, my_simulation, box_width, remaining_time); 		   // Recursive call for remainder of time step
			}
			else																			   // Particle -> wall collision
			{
				double particle_idx1 = gsl_vector_get(my_simulation.collision_indices,0);
				step_all_particles(my_particles,smallest_time);
				update_velocity_wall(my_particles,particle_idx1,box_width);				       // Updating velocity of the two colliding particles 
				check_n_step(my_particles, my_simulation, box_width, remaining_time);          // Recursive call for remainder of time step
			}
		}
		else
		{
			step_all_particles(my_particles,dt);
		}
	}
}

/* Initializing given nr of particles in square array */
void initialize_particles(Particles my_particles, double upper_bound)
{
	int nr_particles      = (*my_particles.positions).size1;
	int particles_pr_axis = (int)ceil(sqrt((double)nr_particles));
	double particle_spacing  = upper_bound/((double)particles_pr_axis);

	// Setting initial position (grid)
	int done = 0;
	for(int row = 0; row < particles_pr_axis; row++)
	{
		if(done == 1){break;}
		for(int col = 0; col < particles_pr_axis; col++)
		{
			int n = row*particles_pr_axis+col;
			if(n < nr_particles)
			{
				double x = (1+col)*particle_spacing-0.5*particle_spacing;
				double y = (1+row)*particle_spacing-0.5*particle_spacing;
				gsl_matrix_set(my_particles.positions,n,0,x);
				gsl_matrix_set(my_particles.positions,n,1,y);
			}
			else{done = 1;break;}
		}
	}

	for(int row = 0; row < nr_particles; row++)				// Setting initial velocities (gauss distribution) 
	{
		double vx = fRand(-0.75,0.75); 						// Random double in (-0.5;0.5)
		double vy = fRand(-0.75,0.75);						// Random double in (-0.5;0.5)
		gsl_matrix_set(my_particles.velocities,row,0,vx);
		gsl_matrix_set(my_particles.velocities,row,1,vy);
	}

	
	for(int row = 0; row < nr_particles; row++)				// Setting pseudorandom particle radius and particle mass (constant density)
	{
		double radius  = fRand(0.0755,0.155); 				// Random double in (0.075;0.325)
		double density = 1.;
		double mass    = density*4./3.*M_PI*pow(radius,3.); // Volume of sphere = 4/3*pi*radius^3 and density = mass/volume
		gsl_matrix_set(my_particles.attributes,row,0,mass);
		gsl_matrix_set(my_particles.attributes,row,1,radius);
	}
 
	write_out(my_particles.attributes,"attributes.txt"); 				// writing out attributes of particles (mass and radia)
}


/////// HELP FUNCTIONS /////////
/* Function for generatinrg pseudorandom gaussian distributed double in [fmin;fmix] */
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

/* Function for computing inner product between two 2d vectors */
double inner_prod(double* arr1, double* arr2)
{
	double inner_prod = 0.0;
	for(int i = 0; i < 2; i++)
	{
		inner_prod = inner_prod + arr1[i]*arr2[i];
	}
	return inner_prod;
}

/* Function for subtracting arr_res = arr1 - arr2 */
void arr_subtract(double* arr1, double* arr2, double* arr_res)
{
	for(int i = 0; i < 2; i++)
	{
		arr_res[i] = arr1[i] - arr2[i];
	}
}



/////// WRITING OUT /////////
/* Function for writing out entries of (nx2) gsl_matrix */
void write_out(gsl_matrix* matrix,char filename[])
{
	FILE* my_out_stream = fopen(filename,"w");
	for(int i = 0; i < (*matrix).size1; i++)
	{
		fprintf(my_out_stream, "%g	%g\n",gsl_matrix_get(matrix,i,0),gsl_matrix_get(matrix,i,1));
	}
	fclose(my_out_stream);
}

/* Function for appending entries of (nx2) gsl_matrix to end of already existing file */
void append_out(gsl_matrix* matrix, char filename[])
{
	FILE* my_out_stream = fopen(filename,"a");
	/* fopen() return NULL if unable to open file in given mode. */
    if (my_out_stream == NULL)
    {
        /* Unable to open file hence exit */
        printf("\nUnable to open '%s' file.\n", filename);
        printf("Please check whether file exists and you have write privilege.\n");
        exit(EXIT_FAILURE);
    }
	else
	{
		for(int row = 0; row < (*matrix).size1; row++)
		{
			double x = gsl_matrix_get(matrix,row,0);
			double y = gsl_matrix_get(matrix,row,1);
			fprintf(my_out_stream,"%g	%g\n",x,y);
		}
		fclose(my_out_stream);
	}

}



/////// PRINTING /////////
/* defining function for printing matrix */
void print_matrix(char s[], gsl_matrix* M)
{
    printf("%s\n",s);
    printf("\n");
	int nr_columns = (*M).size2, nr_rows = (*M).size1;
	for(int i = 0; i < (*M).size1; i++)
	{
		for(int j = 0; j < (*M).size2; j++)
		{
			printf("|%13g", gsl_matrix_get(M,i,j));
		}
		printf("|");
		printf("\n");
	}
    printf("\n");
}

/* Defining function for printing vector */
void print_vector(char s[], gsl_vector* V)
{
	printf("%s\n",s);
	for(int i = 0; i < (*V).size; i++)
	{
		printf(".----------.\n");
		printf("|%10g|\n",gsl_vector_get(V,i));
	}
	printf(".----------.\n");
}

/* Defining function for entries in permutation object as vector */
void print_permut(char s[], gsl_permutation* V)
{
	
	printf("%s\n",s);
	for(int i = 0; i < (*V).size; i++)
	{
		printf(".----------.\n");
		printf("|%10lu|\n",gsl_permutation_get(V,i));
	}
	printf(".----------.\n");
}