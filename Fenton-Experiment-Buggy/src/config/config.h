#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

struct user_options
{
	// Main configuration
	int num_threads;

	// Finite difference configuration section
	double dx;	
	double cell_length;
	int num_cell_div;
	double dt;
	double tmax;
	double lmax;
	double sigma_c;
	double G_gap;
	double diameter;

	// Plot configuration section
	std::ofstream plot_files[5];
	int plot_cell_ids[5];	
	int print_rate;
	int sst_rate;
	bool use_steady_state;
	bool calc_activation_time;

	// Stimulus configuration section
	double stim_start;
	double stim_duration;
	double start_period;
	double end_period;
	double period_step;
	int n_cycles;

	// Steady-State configuration section
	std::string sst_filename;	

};

struct user_options* new_user_options ();
void free_user_options (struct user_options *options);

void read_input_file(struct user_options *options, const char filename[]);

void print_user_options (const struct user_options *options);

#endif
