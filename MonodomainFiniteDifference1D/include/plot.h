#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "utils.h"
#include "config.h"
#include "monodomain.h"

struct plot_config
{
	std::ofstream *plot_files;
	int *plot_cell_ids;	
	int print_rate;
	int sst_rate;
};

struct plot_config* new_plot_config ();
void free_plot_config (struct plot_config *plotter);

void configure_plot_from_options (struct plot_config *plotter, struct monodomain_solver *solver, struct user_options *options);

//void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options);
//void print_monodomain_solver (const struct monodomain_solver *solver);



#endif
