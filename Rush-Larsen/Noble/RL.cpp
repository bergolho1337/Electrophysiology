// Rush-Larsen Method applied over the Noble model (1962)

#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

const int NEQ = 4;
const int PRINTRATE = 1;

const double Cm = 12.0;
const double g_Na_Max = 400000.0;
const double E_Na = 40.0;
const double g_L = 75.0;
const double E_L = -60.0;

void initialCondition (double *y)
{
    y[0] = -87.0;
    y[1] = 0.01;
    y[2] = 0.8;
    y[3] = 0.01;
}

void printSolution (const int k, const double t, const double y[])
{
    printf("%.10lf ",t);
    for (int i = 0; i < NEQ; i++)
        printf("%.10lf ",y[i]);
    printf("\n");
}
// ----------------------------------------------------------------------------------------
double I_Leak (const double t, const double y[])
{
    return g_L * (y[0] - E_L);
}

double g_K1 (const double t, const double y[])
{
    return 1200.0*exp((- y[0] - 90.0)/50.0)+ 15.0*exp((y[0]+90.0)/60.0);
}

double g_K2 (const double t, const double y[])
{
    return 1200.0*pow(y[3], 4.0);
} 

double I_K (const double t, const double y[])
{
    return ( g_K1(t,y) + g_K2(t,y) ) * ( y[0] + 100.0 );
}

double g_Na (const double t, const double y[])
{
    return pow(y[1], 3.0) * y[2] * g_Na_Max;
}

double I_Na (const double t, const double y[])
{
    return (g_Na(t,y) + 140.0) * (y[0] - E_Na);
}

double dvdt (const double t, const double y[])
{
    return -( I_Na(t,y) + I_K(t,y) + I_Leak(t,y) ) / Cm;
}

// ----------------------------------------------------------------------------------------
double alpha_m (const double t, const double y[])
{
    return ( 100.0*(- y[0] - 48.0))/(exp((- y[0] - 48.0)/15.0) - 1.0);
}

double beta_m (const double t, const double y[])
{
    return ( 120.0*(y[0]+8.0))/(exp((y[0]+8.0)/5.0) - 1.0);
}

double inf_m (const double t, const double y[])
{
    return alpha_m(t,y) / (alpha_m(t,y) + beta_m(t,y)); 
}

double tau_m (const double t, const double y[])
{
    return 1.0 / (alpha_m(t,y) + beta_m(t,y)); 
}
// ----------------------------------------------------------------------------------------
double alpha_h (const double t, const double y[])
{
    return 170.0*exp((- y[0] - 90.0)/20.0);
}

double beta_h (const double t, const double y[])
{
    return 1000.0/(1.0+exp((- y[0] - 42.0)/10.0));
}

double inf_h (const double t, const double y[])
{
    return alpha_h(t,y) / (alpha_h(t,y) + beta_h(t,y)); 
}

double tau_h (const double t, const double y[])
{
    return 1.0 / (alpha_h(t,y) + beta_h(t,y)); 
}
// ----------------------------------------------------------------------------------------
double alpha_n (const double t, const double y[])
{
    return ( 0.100000*(- y[0] - 50.0))/(exp((- y[0] - 50.0)/10.0) - 1.0);
}

double beta_n (const double t, const double y[])
{
    return 2.0*exp((- y[0] - 90.0)/80.0);
}

double inf_n (const double t, const double y[])
{
    return alpha_n(t,y) / (alpha_n(t,y) + beta_n(t,y)); 
}

double tau_n (const double t, const double y[])
{
    return 1.0 / (alpha_n(t,y) + beta_n(t,y)); 
}
// ----------------------------------------------------------------------------------------

int main ()
{
    double t;
    double dt = 1.0e-04;
    double tmax = 0.6;
    int N = nearbyint(tmax/dt);

    double *yOld = new double[NEQ];
    double *yNew = new double[NEQ];

    initialCondition(yOld);
    for (int k = 0; k < N; k++)
    {
        t = dt*k;
        if (k % PRINTRATE == 0)
            printSolution(k,t,yOld);
        
        // Forward Euler over the potential
        yNew[0] = yOld[0] + dvdt(t,yOld) * dt;
        // Rush-Larsen 
        yNew[1] = inf_m(t,yOld) + (yOld[1] - inf_m(t,yOld))*exp(-dt/tau_m(t,yOld));
        yNew[2] = inf_h(t,yOld) + (yOld[2] - inf_h(t,yOld))*exp(-dt/tau_h(t,yOld));
        yNew[3] = inf_n(t,yOld) + (yOld[3] - inf_n(t,yOld))*exp(-dt/tau_n(t,yOld));

        swap(yOld,yNew);

    }

    delete [] yOld;
    delete [] yNew;
}
