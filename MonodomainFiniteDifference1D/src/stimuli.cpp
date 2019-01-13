#include "../include/stimuli.h"

struct stim_config* new_stim_config ()
{
	struct stim_config *stim = (struct stim_config*)malloc(sizeof(struct stim_config));

	return stim;
}

void free_stim_config (struct stim_config *stim)
{
	free(stim);
}

void configure_stimulus_from_options (struct stim_config *stim, struct user_options *options)
{
	stim->stim_start = options->stim_start;
	stim->stim_duration = options->stim_duration;
	stim->start_period = options->start_period;
	stim->end_period = options->end_period;
	stim->period_step = options->period_step;
	stim->n_cycles = options->n_cycles;

	print_stim_config(stim);
}

double get_spatial_stim_currents (const double x)
{
	static const double xmin = 0.0;
	static const double xmax = 0.2;

	static const double stim_current = 200.0;

	// Old stimulus code ...
	/*
	if (x >= xmin && x <= xmax)
		return stim_current;
	else
		return 0.0;
	*/

	// New stimulus code using a Gaussian as from the paper
	static const double xp = 0.25;		// Position of the electrode

	return stim_current*gaussian(x-xp);
}

void compute_stimulus (struct stim_config *stim, double *stims, const double cur_time, const int np, const double dx)
{
	double stim_start = stim->stim_start;
	double stim_duration = stim->stim_duration;
	double start_period = stim->start_period;
	double end_period = stim->end_period;
	double period_step = stim->period_step;
	int n_cycles = stim->n_cycles;

	double time = cur_time;
	double new_time, stim_period;
	
	new_time = 0.0;
        // New Jhonny stimulus protocol for alternans simulations ...
        for (double new_period = start_period; new_period >= end_period; new_period -= period_step)
        {
            if ( time >= new_time && (time < new_time + n_cycles*new_period || new_period == end_period) )
            {
                stim_period = new_period;
                time -= new_time;
                break;
            }
            new_time += n_cycles*new_period;

        }
        if( (time-floor(time/stim_period)*stim_period>=stim_start) && ( time - floor(time/stim_period)*stim_period <= stim_start + stim_duration ) )
        {
            //#pragma omp parallel for
            for (int i = 0; i < np; i++) 
            {
		double x = i*dx;
                stims[i] = get_spatial_stim_currents(x);
            }
        }
	else
	{
	    //#pragma omp parallel for
	    for (int i = 0; i < np; i++)
		stims[i] = 0.0;
	}
        time = cur_time;
}

double gaussian (const double x)
{
	static const double sigma = 0.1;	// Gaussian width

	return (1.0/(sigma*sqrt(2.0*M_PI)))*exp(-0.5*pow(x/sigma,2.0));
}

void print_stim_config (struct stim_config *stim)
{
	printf("%s\n",PRINT_LINE);
	printf("[Stimulus] stim_start = %.10lf\n",stim->stim_start);
	printf("[Stimulus] stim_duration = %.10lf\n",stim->stim_duration);
	printf("[Stimulus] start_period = %.10lf\n",stim->start_period);
	printf("[Stimulus] end_period = %.10lf\n",stim->end_period);
	printf("[Stimulus] period_step = %.10lf\n",stim->period_step);
	printf("[Stimulus] n_cycles = %d\n",stim->n_cycles);
	printf("%s\n",PRINT_LINE);
}
