#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../utils/utils.h"
#include "../config/config.h"
#include "../monodomain/monodomain.h"

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

#endif
