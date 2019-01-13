// -------------------------------------------------------------------------------------------------------------------------
// Author: Lucas Berg
//
// This program solves the 1D Monodomain equation using a Finite Difference approach and the Noble celular model from 1962
// 
// The objective of this project is to reproduce the results from the papers
// 	- Model-based control of cardiac alternans in Purkinje fibers from Alejandro Garzón and Roman O. Grigoriev (2011)
//	- Spatiotemporal control of cardiac alternans from Blas Echebarria, and Alain Karma (2002)
//
// The system of ODE's from the celular model has been solved using the Standart Euler Scheme and by using OpenMP.
//
// *************************************************************************************************************************
// HOW TO USE:
//
// To run a simulation using this program the user needs to provide a configuration file with the parameters related to the
// solution of the monodomain equation (discretization sizes, stimulus, output rates, domain size, simulation time).
//
// In order to have a reliable calculation of the output data, we perform 2 simulations:
//	- The first one has the objective of reaching the steady-state of the system (200 stimuli). In the last timestep we
// dump the state vectors of each cell of the fiber in a output file.
//	- The second simulation will actually calculate the output data (conduction velocity and APD's). This is done by
// loading the steady-state file from the previous simulation as the initial condition of the system and we consider only
// 3 periods for the stimulus. 
// 	- The output files will be stored in the 'output' folder.
//	- To process the data and analyze its results the user can run several Python scripts that are stored in the 'scripts'
// folder. (Action potential plot, APD plot, BCL x APD diagram, ...)
// *************************************************************************************************************************  
// EXAMPLE:
//	1) First simulation:
//		./bin/FDMMonodomain1D examples/sst_sample.ini
//	2) Second simulation:
//		./bin/FDMMonodomain1D examples/simple_sample.ini
//	3) Plotting the Action Potential from the cell id 0:
//		cd scripts; python plot_potential.py ../output/sv-0.dat
// -------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <omp.h>

#include "../include/timer.h"
#include "../include/config.h"
#include "../include/utils.h"
#include "../include/plot.h"
#include "../include/stimuli.h"
#include "../include/monodomain.h"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc-1 != 1)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	// Initialize all the structures 
	struct user_options *options;
	options = new_user_options();
	
	struct monodomain_solver *solver;
	solver = new_monodomain_solver();

	struct plot_config *plotter;
	plotter = new_plot_config();

	struct stim_config *stim;
	stim = new_stim_config();

	// Reading input configuration file
	read_input_file(options,argv[1]);

	// Parse the user input into the structures
	configure_solver_from_options(solver,options);
	configure_plot_from_options(plotter,solver,options);
	configure_stimulus_from_options(stim,options);

	// OpenMP configuration
	omp_set_dynamic(0);
    	omp_set_num_threads(2);

	print_configuration_parameters(solver->dx,solver->dt,solver->tmax,solver->lmax,\
					solver->Ncell,solver->Niter,Nodes,\
					plotter->plot_cell_ids,plotter->print_rate,plotter->sst_rate);

	// Call the solver function for the monodomain equation
	solve_monodomain(solver,stim,plotter);
	
	// Free memory
	free_user_options(options);
	free_monodomain_solver(solver);
	free_plot_config(plotter);
	free_stim_config(stim);

	return 0;  
}
