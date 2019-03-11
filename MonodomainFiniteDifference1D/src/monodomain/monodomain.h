#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "../celular_model/celular_model.h"
#include "../config/config.h"
#include "../stimuli/stimuli.h"
#include "../plot/plot.h"
#include "../utils/timer.h"
#include "../utils/utils.h"

#define OUTPUT_VTK			// Flag used to control if we will be writing the output to a VTK file (Paraview visualization)

struct monodomain_solver
{
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

	// Transmembrane potential
	double *vm;
};

struct monodomain_solver* new_monodomain_solver ();
void free_monodomain_solver (struct monodomain_solver *solver);


void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter);

void solve_diffusion (double *sv, const double dx, const double dt,\
 			const int np, const int nodes,\
			double *vm);
void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes);
void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt);

void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options);
void print_monodomain_solver (const struct monodomain_solver *solver);


#endif
