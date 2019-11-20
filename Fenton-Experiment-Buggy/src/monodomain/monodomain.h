#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <omp.h>
             
//#include <Eigen/Sparse>
//#include <Eigen/SparseLU>

#include "../celular_model/celular_model.h"
#include "../config/config.h"
#include "../stimuli/stimuli.h"
#include "../plot/plot.h"
#include "../utils/timer.h"
#include "../utils/utils.h"

#define OUTPUT_VTK			// Flag used to control if we will be writing the output to a VTK file (Paraview visualization)

#define UM2_TO_CM2 0.00000001f
#define UM_PER_MS_TO_M_PER_S 0.001f

struct monodomain_solver
{
	// OpenMP configuration
	int num_threads;

	// Initial condition configuration
	bool use_steady_state;
	std::string sst_filename;

	// Finite difference parameters configuration
	double dx;
	double dt;
	double tmax;
	double lmax;
	double diameter;

	double sigma_c;
	double G_gap;

	int Ncell;
	int Niter;
	int num_cell_div;

	// State vector for each cell
	double *sv;

	// Stimulus current for each cell
	double *stim_current;

	// Transmembrane potential
	double *vm;

	// Cells configuration mask
	bool *cell_mask;

	// Variables used for the propagation velocity calculus
	bool calc_activation_time;
	double *activation_time;
	double *max_dvdt;
};

struct monodomain_solver* new_monodomain_solver ();
void free_monodomain_solver (struct monodomain_solver *solver);


void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter);

void solve_diffusion (const double *sv, double *vm, const bool *mask, const double alpha, const double gamma, const double G_gap, const int ncell, const int nodes);

void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes);
void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt);
//void assembly_matrix (Eigen::SparseMatrix<double> &A, const double h, const double dt,const int ncell);
//void assembly_load_vector (Eigen::VectorXd &b, const double *sv,const double h, const double dt, const int ncell, const int nodes);

void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options);
void configure_cell_mask (bool cell_mask[], const int ncell, const int num_cell_div);
void configure_activation_time(double *activation_time, double *max_dvdt, const int ncell);

double calculate_flux (const bool type, const double gamma, const double G_gap);
void calculate_derivative (double *at, double *max_dvdt, double *v_old, double *v_new, const int ncell, const int nodes, const double dt, const double cur_time);
void calculate_propagation_velocity (const double *at, const int *plot_ids, const double dx, const int ncell);

void print_monodomain_solver (const struct monodomain_solver *solver);


#endif
