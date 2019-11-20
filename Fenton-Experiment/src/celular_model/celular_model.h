#ifndef CELULAR_MODEL_H
#define CELULAR_MODEL_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

static const int Nodes = 4;

void compute_initial_conditions (double *sv, const int Ncell, const int Nodes);
void read_initial_conditions_from_file (double *sv, const int Ncell, const int Nodes, const std::string filename);
double dvdt (const double V, const double m, const double h, const double n, const double stim_current);
double dmdt (const double V, const double m);
double dhdt (const double V, const double h);
double dndt (const double V, const double n);

double dmdt_RL (const double V, const double m, const double dt);
double dhdt_RL (const double V, const double h, const double dt);
double dndt_RL (const double V, const double n, const double dt);

#endif
