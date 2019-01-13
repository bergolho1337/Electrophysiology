#ifndef STIMULI_H
#define STIMULI_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.h"
#include "utils.h"

struct stim_config
{
	double stim_start;
	double stim_duration;
	double start_period;
	double end_period;
	double period_step;
	int n_cycles;
};

struct stim_config* new_stim_config ();
void free_stim_config (struct stim_config *stim);

void configure_stimulus_from_options (struct stim_config *stim, struct user_options *options);
void print_stim_config (struct stim_config *stim);

double get_spatial_stim_currents (const double x);
void compute_stimulus (struct stim_config *stim, double *stims, const double cur_time, const int np, const double dx);

double gaussian (const double x);

#endif
