#include "config.h"

struct user_options* new_user_options ()
{
	struct user_options *options = (struct user_options*)malloc(sizeof(struct user_options));
	
	// Default values
	options->num_threads = 1;
	options->print_rate = 100;
	options->sst_rate = 600000;
	options->use_steady_state = false;
	options->dx = 100.0;
	options->dt = 0.01;
	options->tmax = 600.0;
	options->lmax = 50000.0;
	options->diameter = 100.0;
	options->sigma_c = 0.0004;
	options->G_gap = 0.628;
	options->num_cell_div = 1;

	return options;
}

void free_user_options (struct user_options *options)
{
	free(options);
}

void read_input_file(struct user_options *options, const char filename[])
{
	std::string str, str2, str3;
	std::ifstream file(filename);
		
	// Reading main section
	while (getline(file,str) && str != "[main]");
	file >> str >> str2 >> options->num_threads;
	file >> str >> str2 >> options->dx;
	file >> str >> str2 >> options->dt;
	file >> str >> str2 >> options->tmax;
	file >> str >> str2 >> options->lmax;
	file >> str >> str2 >> options->diameter;
	file >> str >> str2 >> options->num_cell_div;
	file >> str >> str2 >> options->sigma_c;
	file >> str >> str2 >> options->G_gap;
	file >> str >> str2 >> options->print_rate;
	file >> str >> str2 >> options->sst_rate;
	file >> str >> str2 >> str3;
	if (str3 == "true" || str3 == "yes")
		options->use_steady_state = true;
	else
		options->use_steady_state = false;
	file >> str >> str2 >> str3;
	options->sst_filename = str3;
	file >> str >> str2 >> str3;
	if (str3 == "true" || str3 == "yes")
		options->calc_activation_time = true;
	else
		options->calc_activation_time = false;		
	
	// Reading stimulus section
	while (getline(file,str) && str != "[stimulus]");
	file >> str >> str2 >> options->stim_start;
	file >> str >> str2 >> options->stim_duration;
	file >> str >> str2 >> options->start_period;
	file >> str >> str2 >> options->end_period;
	file >> str >> str2 >> options->period_step;
	file >> str >> str2 >> options->n_cycles;;

	// DEBUG
	//print_user_options(options);

	file.close();

}

void print_user_options (const struct user_options *options)
{
	printf("[MAIN]\n");
	printf("[Config] num_threads = %d\n",options->num_threads);
	printf("[Config] cell_length = %.10lf\n",options->cell_length);
	printf("[Config] num_cell_div = %d\n",options->num_cell_div);
	printf("[Config] dx = %.10lf\n",options->dx);
	printf("[Config] dt = %.10lf\n",options->dt);	
	printf("[Config] tmax = %.10lf\n",options->tmax);
	printf("[Config] lmax = %.10lf\n",options->lmax);
	printf("[Config] diameter = %.10lf\n",options->diameter);
	printf("[Config] sigma_c = %.10lf\n",options->sigma_c);
	printf("[Config] G_gap = %.10lf\n",options->G_gap);
	printf("[Config] print_rate = %d\n",options->print_rate);
	printf("[Config] sst_rate = %d\n",options->sst_rate);
	printf("[Config] use_steady_state = %d\n",options->use_steady_state);
	printf("[Config] steady_state_filename = %s\n",options->sst_filename.c_str());
	
	printf("\n[STIMULUS]\n");
	printf("[Config] stim_start = %.10lf\n",options->stim_start);
	printf("[Config] stim_duration = %.10lf\n",options->stim_duration);
	printf("[Config] start_period = %.10lf\n",options->start_period);
	printf("[Config] end_period = %.10lf\n",options->end_period);
	printf("[Config] period_step = %.10lf\n",options->period_step);
	printf("[Config] n_cycles = %d\n",options->n_cycles);
	
}
