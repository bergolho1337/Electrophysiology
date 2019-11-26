#include "plot.h"

struct plot_config* new_plot_config ()
{
	struct plot_config *plotter = (struct plot_config*)malloc(sizeof(struct plot_config));

	plotter->plot_files = NULL;
	plotter->plot_cell_ids = NULL;

	return plotter;
}

void free_plot_config (struct plot_config *plotter)
{
	close_plot_files(plotter->plot_files);

	delete [] plotter->plot_files;
	delete [] plotter->plot_cell_ids;

	free(plotter);
}

void configure_plot_from_options (struct plot_config *plotter, struct monodomain_solver *solver, struct user_options *options)
{
	// Allocate memory 
	plotter->plot_files = new std::ofstream[5];
	plotter->plot_cell_ids = new int[5];

	// Get input information to build the output folder
	double lmax = options->lmax;
	double bcl = options->start_period;

	int offset = nearbyint(solver->lmax / solver->dx) / 5;
	for (int i = 0; i < 5; i++)
	{
		std::stringstream filename;
		plotter->plot_cell_ids[i] = i*offset;	
	
		//filename << "output/sv-" << plotter->plot_cell_ids[i] << ".dat";
		filename << std::setprecision(2) << "results/" << lmax << "cm/" << int(bcl) << "ms/sv-" << plotter->plot_cell_ids[i] << ".dat";
		//printf("%s\n",filename.str().c_str());
		
		plotter->plot_files[i].open(filename.str().c_str());		
	}
	
	plotter->print_rate = options->print_rate;
	plotter->sst_rate = options->sst_rate;

	printf("[Plot] Output state-vector files will be saved inside the 'results' folder !\n");
}
