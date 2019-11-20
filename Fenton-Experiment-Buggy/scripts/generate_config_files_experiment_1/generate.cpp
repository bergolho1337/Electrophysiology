// ========================================================================================================
// This simulation will be the FIRST experiment of the project
// We consider a cable with homogenous conductivity
// ========================================================================================================

#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;

const int MAX_FILENAME_SIZE = 500;

const int NUM_THREADS = 6;
const double DX = 100.0;
const double DT = 0.001;
const double DIAMETER = 100.0;
const int NUM_CELL_DIVISIONS = 0;
const double G_GAP = 0.628;

const int N_CYCLES = 20;
const double PERIOD = 280.0;

const double SST_TMAX = N_CYCLES*PERIOD;
const unsigned int SST_PRINT_RATE = SST_TMAX/DT + 1;
const double TMAX = PERIOD;
const unsigned int PRINT_RATE = 1000;
const unsigned int SST_RATE = SST_TMAX/DT;

const double STIM_START = 0.0;
const double STIM_DURATION = 2.0;
const double PERIOD_STEP = 100.0;

void write_sst_config_file (const double cable_length, const double sigma)
{
    char filename[MAX_FILENAME_SIZE];
    sprintf(filename,"files/sst-cable-%.1lfum-sigma-%lf.ini",cable_length,sigma);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads = %d\n",NUM_THREADS);
    fprintf(file,"dx = %lf\n",DX);
    fprintf(file,"dt = %lf\n",DT);
    fprintf(file,"tmax = %lf\n",SST_TMAX);
    fprintf(file,"lmax = %lf\n",cable_length);
    fprintf(file,"diameter = %lf\n",DIAMETER);
    fprintf(file,"num_cell_div = %d\n",NUM_CELL_DIVISIONS);
    fprintf(file,"sigma_c = %lf\n",sigma);
    fprintf(file,"G_gap = %lf\n",G_GAP);
    fprintf(file,"print_rate = %u\n",SST_PRINT_RATE);
    fprintf(file,"sst_rate = %u\n",SST_RATE);
    fprintf(file,"use_steady_state = no\n");
    fprintf(file,"sst_filename = steady_state/sst-cable-%.1lfum-sigma-%lf.sst\n",cable_length,sigma);
    fprintf(file,"calc_activation_time = no\n");
    
    fprintf(file,"\n\n");
    
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = %lf\n",STIM_START);
    fprintf(file,"stim_duration = %lf\n",STIM_DURATION);
    fprintf(file,"start_period = %lf\n",PERIOD);
    fprintf(file,"end_period = %lf\n",PERIOD);
    fprintf(file,"period_step = %lf\n",PERIOD_STEP);
    fprintf(file,"n_cycles = %d\n",N_CYCLES);
    
    fclose(file);
}

void write_simulation_config_file (const double cable_length, const double sigma)
{

    char filename[MAX_FILENAME_SIZE];
    sprintf(filename,"files/cable-%.1lfum-sigma-%lf.ini",cable_length,sigma);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads = %d\n",NUM_THREADS);
    fprintf(file,"dx = %lf\n",DX);
    fprintf(file,"dt = %lf\n",DT);
    fprintf(file,"tmax = %lf\n",TMAX);
    fprintf(file,"lmax = %lf\n",cable_length);
    fprintf(file,"diameter = %lf\n",DIAMETER);
    fprintf(file,"num_cell_div = %d\n",NUM_CELL_DIVISIONS);
    fprintf(file,"sigma_c = %lf\n",sigma);
    fprintf(file,"G_gap = %lf\n",G_GAP);
    fprintf(file,"print_rate = %u\n",PRINT_RATE);
    fprintf(file,"sst_rate = %u\n",SST_RATE);
    fprintf(file,"use_steady_state = yes\n");
    fprintf(file,"sst_filename = steady_state/sst-cable-%.1lfum-sigma-%lf.sst\n",cable_length,sigma);
    fprintf(file,"calc_activation_time = yes\n");
    
    fprintf(file,"\n\n");
    
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = %lf\n",STIM_START);
    fprintf(file,"stim_duration = %lf\n",STIM_DURATION);
    fprintf(file,"start_period = %lf\n",PERIOD);
    fprintf(file,"end_period = %lf\n",PERIOD);
    fprintf(file,"period_step = %lf\n",PERIOD_STEP);
    fprintf(file,"n_cycles = %d\n",N_CYCLES);
    
    fclose(file);
}

void write_config_files (const double cable_length, const double sigma)
{
    write_sst_config_file(cable_length,sigma);
    write_simulation_config_file(cable_length,sigma);
}

int main ()
{
    double start_cable_length = 10000.0;
    double end_cable_length = 50000.0;
    double cable_length_step = 5000.0;
    
    double start_sigma = 0.0001;
    double end_sigma = 0.001;
    double sigma_step = 0.0001;

    double cable_length = start_cable_length;
    double sigma = start_sigma;

    while (cable_length <= end_cable_length)
    {
        sigma = start_sigma;

        while (sigma <= end_sigma)
        {
            printf("Writing 'cable_length = %lf' and 'sigma = %lf' files ...\n",cable_length,sigma);
            write_config_files(cable_length,sigma);

            sigma += sigma_step;
        }

        cable_length += cable_length_step;
    }
    
    return 0;
}
