#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <omp.h>
             
#include "../celular_model/celular_model.h"
#include "../config/config.h"
#include "../stimuli/stimuli.h"
#include "../plot/plot.h"
#include "../utils/timer.h"
#include "../utils/utils.h"

//#define OUTPUT_VTK			// Flag used to control if we will be writing the output to a VTK file (Paraview visualization)

#define UM2_TO_CM2 0.00000001f
#define CELL_LENGTH 0.01		// {cm}

struct monodomain_solver;
struct tridiagonal_matrix_thomas;

struct monodomain_solver
{
	// Type of cable
	bool use_homogenous_conductivity;

	// Initial condition configuration
	bool use_steady_state;
	std::string sst_filename;

	// Finite difference parameters configuration
	double dx;
	double dt;
	double tmax;
	double lmax;

	int Ncell;
	int Niter;

	// State vector for each cell
	double *sv;

	// Stimulus current for each cell
	double *stim_current;

	// Cell link mask
	bool *link_mask; 

	// Transmembrane potential
	double *vm;

	struct tridiagonal_matrix_thomas *the_matrix;
};

struct tridiagonal_matrix_thomas
{
	uint32_t n;

	double *a;
	double *b;
	double *c;
	double *d;

	double *c_star;
	double *d_star;
};

struct monodomain_solver* new_monodomain_solver ();
struct tridiagonal_matrix_thomas* new_tridiagonal_matrix_thomas ();
void free_monodomain_solver (struct monodomain_solver *solver);
void free_tridiagonal_matrix_thomas (struct tridiagonal_matrix_thomas *matrix);


void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter);

void solve_diffusion (const double *sv, double *vm, bool *cell_mask, const double beta, const double cm, const double sigma, const double G, const double dx, const double d, const double dt, const int ncell, const int nodes);
void solve_diffusion_v2 (struct monodomain_solver *solver);

void solve_linear_system_using_thomas (const double *a, const double *b, const double *c, const double *d,\
			  double *c_star, double *d_star, double *vms, const uint32_t ncell);

void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes);
void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt);
void assembly_discretization_matrix (struct monodomain_solver *solver,\
			 const double sigma_c, const double G_gap, const double sigma_bar, const double alpha, const double beta, const double cm,\
			 const double d, const double dx, const double dt,\
			 const bool use_homogenous_conductivity);
void assembly_load_vector (struct monodomain_solver *solver, const double alpha);

double calc_homogenous_conductivity (const double sigma, const double G);

void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options);
void configure_cell_mask (struct monodomain_solver *solver);
void print_monodomain_solver (const struct monodomain_solver *solver,\
		const double d, const double sigma, const double G, const double sigma_bar,const double beta,const double cm, const bool use_homogenous_conductivity);


#endif
