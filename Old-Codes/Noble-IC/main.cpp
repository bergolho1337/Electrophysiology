// Noble model for Purkinje cells

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>  

using namespace std;

// Parameters:
double g_Na_max = 4.0e+05;
double E_Na = 4.0e+01;
double g_L = 7.5e+01;
double E_L = -6.0e+01;
double Cm = 1.2e+01;

/* ********************************************************************************************************************** */
// dV/dt

double calc_i_Leak (double V_old)
{
	return ((g_L*(V_old-E_L)));
}

double calc_g_K2 (double n_old)
{
	return ((1.2e+03*pow(n_old,4.0e+00)));  
}

double calc_g_K1 (double V_old)
{
	return (((1.2e+03*exp((((-V_old)-9.0e+01)/5.0e+01)))+(1.5e+01*exp(((V_old+9.0e+01)/6.0e+01)))));
}

double calc_i_K (double V_old, double n_old)
{
	return (((calc_g_K1(V_old)+calc_g_K2(n_old))*(V_old+1.0e+02)));
}

double calc_g_Na (double m_old, double h_old)
{
	return ((pow(m_old,3.0e+00)*h_old*g_Na_max));
}

double calc_i_Na (double V_old, double m_old, double h_old)
{
	return ((calc_g_Na(m_old,h_old)+1.4e+02)*(V_old-E_Na));
}


double dvdt (double t,double V_old, double m_old, double h_old, double n_old)
{
	return ((-(calc_i_Na(V_old,m_old,h_old)+calc_i_K(V_old,n_old)+calc_i_Leak(V_old)))/Cm);
}

/* ********************************************************************************************************************** */
// dm/dt

double calc_beta_m (double V_old)
{
	return (((1.2e+02*(V_old+8.0e+00))/(exp(((V_old+8.0e+00)/5.0e+00))-1.0e+00)));
}

double calc_alpha_m (double V_old)
{
	return (((1.0e+02*((-V_old)-4.8e+01))/(exp((((-V_old)-4.8e+01)/1.5e+01))-1.0e+00)));
}

double dmdt (double t,double V_old, double m_old)
{
	return ((calc_alpha_m(V_old)*(1.0e+00-m_old))-(calc_beta_m(V_old)*m_old));
}

/* ********************************************************************************************************************** */
// dh/dt

double calc_beta_h (double V_old)
{
	return ((1.0e+03/(1.0e+00+exp((((-V_old)-4.2e+01)/1.0e+01)))));
}

double calc_alpha_h (double V_old)
{
	return ((1.7e+02*exp((((-V_old)-9.0e+01)/2.0e+01))));
}

double dhdt (double t,double V_old, double h_old)
{
	return ((calc_alpha_h(V_old)*(1.0e+00-h_old))-(calc_beta_h(V_old)*h_old));
}

/* ********************************************************************************************************************** */
// dn/dt

double calc_beta_n (double V_old)
{
	return ((2.0e+00*exp((((-V_old)-9.0e+01)/8.0e+01))));
}

double calc_alpha_n (double V_old)
{
	return (((1.0e-01*((-V_old)-5.0e+01))/(exp((((-V_old)-5.0e+01)/1.0e+01))-1.0e+00)));
}

double dndt (double t,double V_old, double n_old)
{
	return ((calc_alpha_n(V_old)*(1.0e+00-n_old))-(calc_beta_n(V_old)*n_old));
}


/* ********************************************************************************************************************** */

void makeGraph ()
{
	FILE *out = fopen("graph.plt","w+");
	fprintf(out,"set title \"Noble Model\"\n");
	fprintf(out,"set grid\n");
	fprintf(out,"set xlabel \"t (s)\"\n");
	fprintf(out,"set terminal png\n");
	fprintf(out,"set output \"noble.png\"\n");
	//fprintf(out,"plot \"data.dat\" using 1:2 title \"v\" w l, \"data.dat\" using 1:3 title \"m\" w l, \"data.dat\" using 1:4 title \"h\" w l, \"data.dat\" using 1:5 title \"n\" w l\n");
	fprintf(out,"plot \"data.dat\" using 1:2 title \"v\" w l");
	fclose(out);
	if (system("gnuplot graph.plt") == 0)
		cout << "[+] Graphic plotted !" << endl;
	else
		cout << "[-] ERROR ! Plotting graph !" << endl;
}

int main (int argc, char *argv[])
{
	FILE *out = fopen("data.dat","w+");
	double V_old, m_old, h_old, n_old;
	double V_new, m_new, h_new, n_new;
	
	double t_max, Dt, t;
	int k;
	// Initial conditions
	//V_old = -8.7e+01;
	V_old = -8.7e+01;
	m_old = 1.0e-02;
	h_old = 8.0e-01;
	n_old = 1.0e-02;
	
	if (argc-1 < 2)
	{
		cout << "Usage:> ./noble <t_max> <Dt>" << endl;
		exit(1);
	}	
	else
	{
		t_max = atof(argv[1]);
		Dt = atof(argv[2]);
		k = (int)nearbyint(t_max/Dt);
		cout << "Number of subintervals = " << k << endl;
		fprintf(out,"%e %e %e %e %e\n", 0.0, V_old, m_old, h_old, n_old);
		for (int i = 1; i < k; i++)
		{
			t = i*Dt;
			V_new = V_old + dvdt(t,V_old,m_old,h_old,n_old)*Dt;
			m_new = m_old + dmdt(t,V_old,m_old)*Dt;
			h_new = h_old + dhdt(t,V_old,h_old)*Dt;
			n_new = n_old + dndt(t,V_old,n_old)*Dt;
			fprintf(out,"%e %e %e %e %e\n", t, V_new, m_new, h_new, n_new);
			V_old = V_new;
			m_old = m_new;
			h_old = h_new;
			n_old = n_new;
		}
	}
	fclose(out);		
	makeGraph();
}
