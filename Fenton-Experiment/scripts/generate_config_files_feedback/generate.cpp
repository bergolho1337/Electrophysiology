#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;

const int MAX_SIZE = 500;
const double DT = 0.01;                 // {ms}
const double N_CYCLES = 200.0;           
const double CABLE_LENGTH = 1.0;        // {cm}

int SST_RATE = nearbyint(N_CYCLES/DT);

void write_sst_config_file (const int period, const int start_period, const int step_period)
{
    double tmax = N_CYCLES*period;
    int sst_rate = nearbyint(tmax/DT);
    int print_rate = sst_rate + 1;

    char filename[MAX_SIZE];
    sprintf(filename,"files/prepacing:cable-%g_bcl-%d.ini",CABLE_LENGTH,period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"dx = 0.01\n");
    fprintf(file,"dt = %g\n",DT);
    fprintf(file,"tmax = %g\n",tmax);
    fprintf(file,"lmax = %g\n",CABLE_LENGTH);
    fprintf(file,"print_rate = %d\n",print_rate);
    fprintf(file,"sst_rate = %d\n",sst_rate);
    fprintf(file,"use_steady_state = no\n");
    fprintf(file,"sst_filename = teste.sst\n");
    /*
    if (period == start_period)
    {
	fprintf(file,"use_steady_state = no\n");
        fprintf(file,"sst_filename = teste.sst\n");
    }
    else
    {
	fprintf(file,"use_steady_state = yes\n");
        fprintf(file,"sst_filename = steady_state/cable-5cm-%dms.sst\n",period+step_period);
    }
    */
    fprintf(file,"\n\n");
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 2.0\n");
    fprintf(file,"start_period = %d\n",period);
    fprintf(file,"end_period = %d\n",period);
    fprintf(file,"period_step = 100\n");
    fprintf(file,"n_cycles = 300\n");
    
    fclose(file);
}

void write_simulation_config_file (const int period, const int start_period)
{
    double tmax = 3*period;
    int print_rate = 100;
    int sst_rate = __INT32_MAX__;

    char filename[MAX_SIZE];
    sprintf(filename,"files/cable-%g_bcl-%d.ini",CABLE_LENGTH,period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"dx = 0.01\n");
    fprintf(file,"dt = %g\n",DT);
    fprintf(file,"tmax = %g\n",tmax);
    fprintf(file,"lmax = %g\n",CABLE_LENGTH);
    fprintf(file,"print_rate = %d\n",print_rate);
    fprintf(file,"sst_rate = %d\n",sst_rate);
    fprintf(file,"use_steady_state = yes\n");
    fprintf(file,"sst_filename = steady_state/experiment-1/cable-%g_bcl-%d.sst\n",CABLE_LENGTH,period);
    fprintf(file,"\n\n");
    fprintf(file,"[stimulus]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 2.0\n");
    fprintf(file,"start_period = %d\n",period);
    fprintf(file,"end_period = %d\n",period);
    fprintf(file,"period_step = 100\n");
    fprintf(file,"n_cycles = 20\n");
    
    fclose(file);
}

void write_config_files (const int period, const int start_period, const int step_period)
{
    write_sst_config_file(period,start_period,step_period);
    write_simulation_config_file(period,start_period);
}

int main ()
{
    int start_period = 280;
    int end_period = 150;
    int period_step = 5;
    for (int period = start_period; period >= end_period; period -= period_step)
    {
        printf("Writing %dms files ...\n",period);
        write_config_files(period,start_period,period_step);
    }
    return 0;
}
