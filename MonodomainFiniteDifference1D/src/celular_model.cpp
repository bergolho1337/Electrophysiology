// Here we are considering the Noble model for Purkinje fibers from 1962.
// TODO: Implement a generic way to include other celular models (see New-Fenton-Experiment code ...)

#include "../include/celular_model.h"

void compute_initial_conditions (double *sv, const int Ncell, const int Nodes)
{
	printf("[Celular Model] Using Noble 1962 celular model\n");
	printf("[Celular Model] Initializing cells with default initial conditions\n");
	//#pragma omp parallel
	for (int i = 0; i < Ncell; i++)
	{
		sv[i*Nodes+0] = -87.0;   	// V
		sv[i*Nodes+1] = 0.01;     	// m
		sv[i*Nodes+2] = 0.8;     	// h       
		sv[i*Nodes+3] = 0.01;     	// n
	}
}

void read_initial_conditions_from_file (double *sv, const int Ncell, const int Nodes, const std::string filename)
{
	printf("[Celular Model] Using Noble 1962 celular model\n");
	printf("[Celular Model] Initializing cells with steady-state file '%s'\n",filename.c_str());
	
	FILE *file = fopen(filename.c_str(),"r");
	for (int i = 0; i < Ncell; i++)
	{
		for (int j = 0; j < Nodes; j++)
		{
			fscanf(file,"%lf",&sv[i*Nodes+j]);
		}
	}
	fclose(file);
}


double dvdt (const double V, const double m, const double h, const double n, const double stim_current)
{
	// Parameters
	static const double g_Na_max = 4.0e+02;
	static const double E_Na = 4.0e+01;
	static const double g_L = 7.5e-02;
	static const double E_L = -6.0e+01;
	static const double CM = 1.2e+01;

	double gna = ((pow(m,3.0e+00)*h*g_Na_max));
	double ina = (gna+1.4e-01)*(V-E_Na);		// Original
	//double ina = (gna+1.2e-01)*(V-E_Na);		// Turn to excitable ...
	double gk1 = (((1.2*exp((((-V)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V+9.0e+01)/6.0e+01)))));
	double gk2 = ((1.2*pow(n,4.0e+00)));
	double ik = (((gk1+gk2)*(V+1.0e+02)));
	double ileak = ((g_L*(V-E_L)));
	double istim = stim_current;

	return ((-(ina+ik+ileak)+istim)/CM);
}

double dmdt (const double V, const double m)
{
	double beta_m = (((1.2e-01*(V+8.0e+00))/(exp(((V+8.0e+00)/5.0e+00))-1.0e+00)));
    	double alpha_m = (((1.0e-01*((-V)-4.8e+01))/(exp((((-V)-4.8e+01)/1.5e+01))-1.0e+00)));

	return ((alpha_m*(1.0e+00-m))-(beta_m*m));
}

double dhdt (const double V, const double h)
{
	double beta_h = ((1.0/(1.0e+00+exp((((-V)-4.2e+01)/1.0e+01)))));
    	double alpha_h = ((1.7e-01*exp((((-V)-9.0e+01)/2.0e+01))));

	return ((alpha_h*(1.0e+00-h))-(beta_h*h));
}

double dndt (const double V, const double n)
{
	double beta_n = ((2.0e-03*exp((((-V)-9.0e+01)/8.0e+01))));
    	double alpha_n = (((1.0e-04*((-V)-5.0e+01))/(exp((((-V)-5.0e+01)/1.0e+01))-1.0e+00)));

	return ((alpha_n*(1.0e+00-n))-(beta_n*n));
}


