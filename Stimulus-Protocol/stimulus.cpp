#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

// Mitchell-Shaeffer variables
double v0 = 0.00000820413566106744f;
double h0 = 0.8789655121804799f;
double IstimStart = 0.0f;
double IstimEnd = 50000.0f;
double IstimAmplitude = 0.2f;
double IstimPeriod = 500;
double IstimPulseDuration = 1.0f;
double tau_in = 0.3f;
double tau_open = 120.0f;
double tau_close = 150.0f;
double V_gate = 0.13f;
double tau_out = 6.0f;
int neq = 2;

// GLOBAL VARIABLES
double start_period = 400.0f;
double end_period = 150.0f;
double stim_start = 0.0f;
double stim_duration = 1.0f;
int nCycles = 5;
double period_step = 50.0f;

double istim (double time)
{
    double newPeriod; double newTime = 0.0f;
    double stim_period;
    
    for (newPeriod = start_period; newPeriod >= end_period; newPeriod-=period_step) 
    { 
        
        if (time >= newTime && (time < newTime + nCycles*newPeriod || newPeriod==120)) 
        {
            stim_period = newPeriod;
            time = time - newTime;
            break;
        }
        newTime += nCycles*newPeriod;
    }
    
    if( (time-floor(time/stim_period)*stim_period>=stim_start) && ( time - floor(time/stim_period)*stim_period <= stim_start + stim_duration ) )
    {
        //cout << time << endl;
        //return stim_amplitude;
        return IstimAmplitude;
    }
    else
    {
        return 0.0F;
    }
}

double j_in (double *y_old)
{
    return ( y_old[1]*( pow(y_old[0], 2.0f)*(1.0f - y_old[0])))/tau_in;
}

double j_out (double *y_old)
{
    return - (y_old[0]/tau_out);
}

void allocMemory (double **y_old, double **y_new)
{
    *y_old = new double[neq];
    *y_new = new double[neq];
}

void freeMemory (double *y_old, double *y_new)
{
    delete [] y_old;
    delete [] y_new;
}

void initialCondition (double *y_old)
{
    y_old[0] = v0;
    y_old[1] = h0;
}

void print_current_solution (double t, double *y)
{
    cout << fixed << setprecision(10) << "t = " << t << " || ";
    for (int i = 0; i < neq-1; i++)
        cout << fixed << setprecision(10) << "V = " << y[0] << " ";
    cout << fixed << setprecision(10) << "|| h = " << y[1] << endl;
}

void comp_rates (double t, double h, double *y_old, double *y_new)
{
    // Compute all the currents
    double Jstim = istim(t);
    double Jin = j_in(y_old);
    double Jout = j_out(y_old);

    // Compute rates
    double dvdt = Jstim + Jin + Jout;
    double dhdt = (y_old[0]<V_gate ? (1.0f - y_old[1])/tau_open : - y_old[1]/tau_close);
    
    // Forward Euler
    y_new[0] = y_old[0] + dvdt*h;
    y_new[1] = y_old[1] + dhdt*h;
}

int main ()
{
    double *y_old, *y_new;
    double tmax = 8250.0f;
    double h = 0.1f;
    int N = nearbyint(tmax/h);

    allocMemory(&y_old,&y_new);
    initialCondition(y_old);

    ofstream out("solution.dat");

    for (int i = 0; i < N; i++)
    {
        double t = i*h;

        // Write current solution to file
        out << fixed << setprecision(10) << t << " " << y_old[0] << " " << y_old[1] << endl;

        comp_rates(t,h,y_old,y_new);

        // Next timestep
        swap(y_old,y_new);
    }

    out.close();
    freeMemory(y_old,y_new);

    return 0;
}