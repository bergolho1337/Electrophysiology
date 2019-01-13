/* ================================================================================== */
// Versao 2 - Modelo de Controle Volume Celular
// Autor: Lucas Berg
// Last Update: 06/12/2016 - 23:38
// Simula o estado transiente do modelo Pump-Leak Model do livro do Keener
/* ================================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Bomba de sodio/potassio
const double  kk = 2.0;
const double  kna = 7.7;

/* Parametros do metodo de Euler explicito */
const double dt = 1.0e-05;                      // Tamanho do passo de tempo
const int num_EDO = 3;                          // Numero de equacoes diferencias
const int num_eq = 2;                           // Numero de equacoes algebricas

/* Parametros do modelo */
const double Na_e = 145.0;                      // Concentracao de Na+ no meio extracelular (mM=mmol/L)
//const double Na_e = 200.0;
const double K_e = 3.5;                         // Concentracao de K+ no meio extracelular (mM=mmol/L)
const double Na_i = 10.0;                       // Concentracao de Na+ no meio intracelular (mM=mmol/L)
const double K_i = 130.0;                       // Concentracao de K+ no meio intracelular (mM=mmol/L)
const double g_Na = 1.22258163e-02;             // Condutancia do canal de Na+ (mS/cm^2)
const double g_K = 6.28564508e-02;              // Condutancia do canal de K+  (mS/cm^2)
const double F = 9.64e+04;                      // Constante de Faraday (C/mol)
const double T = 370.0;                         // Temperatura (K)
const double Cm = 4.16121841e-03;               // Capacitancia total da membrana (s/ohm.cm^2)
const double Xi = 3.669e-06;                    // Numero de proteinas negativamente carregadas no meio intracelular (mmol)
const double RTF = 31.9;                        // (mV)
const double R = 8.314462;                      // (L.kPa/K.mol)
const double V0 = -70.0;                        // Potencial transmembranico inicial (mV)
const double Gamma = 6.4103e+03;                // Relacao area/volume (cm^-1)
const double nw_chap = 1.7e-05;                 // nw_chap = nw / RT --> (cm/s.Pa)
const double alpha_i_0 = 8.69565217e-01;        // Porcentagem de volume intracelular
//const double Imax = 13.0;                     // Corrente maxima (uA/cm^2) (Equilibrio)
const double Imax = 10000.0;
const double Pe = 266.64;                       // Pressao intertiscial (Pa)
const double S = 1.0;                           // Constante de rigidez
const double C1 = 8.5;                          // Constante 1 - Concentracao das proteinas internas

/* Funcoes do modelo */
double hNaK (double *yOld);
double I_Na (double *yOld);
double I_K (double *yOld);
double Naidt (double *yOld, double *yNew);
double Kidt (double *yOld, double *yNew);
double Vm (double *yOld);
double Pi (double *yOld);
double pi_e (double *yOld);
double pi_i (double *yOld);
double alpha_i (double *yOld);
double calcAlphaDt (double *yOld, double *yNew);

// yOld[0] = Vm  
// yOld[1] = Pi
// yOld[2] = alpha_i  
// yOld[3] = Nai  
// yOld[4] = Ki

double calcAlphaDt (double *yOld, double *yNew)
{
  return ( (yNew[2]-yOld[2])/(dt) );
}

double hNaK (double *yOld)
{
  return Imax*pow(1+(kk/K_e),-2)*pow(1+(kna/yOld[3]),-3);
}

double I_Na (double *yOld)
{
  return g_Na*( yOld[0] - ( RTF*log(Na_e/yOld[3]) ) );
}

double I_K (double *yOld)
{

  return g_K*( yOld[0] - ( RTF*log(K_e/yOld[4]) ) );
}

double pi_e (double *yOld)
{
  return R*T*(Na_e+K_e);
}

double pi_i (double *yOld)
{
  return R*T*(yOld[3]+yOld[4]+C1);
}

// 1
double Vm (double *yOld)
{
  return ( (-0.001*F*(1-yOld[2])*(Na_e+K_e)) / (Gamma*Cm) );
}

// 2
double Pi (double *yOld)
{
  return ( Pe + S*(yOld[2] - alpha_i_0) );
}

// 3
double alpha_i (double *yOld)
{
  return ( (-Gamma*nw_chap)*(yOld[1] - Pe + pi_e(yOld) - pi_i(yOld) ) );
}

// 4
double Naidt (double *yOld, double *yNew)
{
  return ( (1/yOld[2])*( ((-Gamma/F)*( I_Na(yOld) + (2*hNaK(yOld)) )) - (calcAlphaDt(yOld,yNew)*yOld[3]) ) );
}

// 5
double Kidt (double *yOld, double *yNew)
{
  return ( (1/yOld[2])*( ((-Gamma/F)*( I_K(yOld) - (3*hNaK(yOld)) )) - (calcAlphaDt(yOld,yNew)*yOld[4]) ) );
}

void allocMemory (double **v)
{
  *v = (double*)calloc(num_EDO+num_eq,sizeof(double));
}

void setInitialConditions (double *yOld)
{
  yOld[0] = V0;
  yOld[1] = Pe;
  yOld[2] = alpha_i_0;
  yOld[3] = Na_i;
  yOld[4] = K_i;
}

void writeSolution (double t, double *y, int n)
{
  int i;
  printf("%lf",t);
  for (i = 0; i < n; i++)
    printf(" %lf",y[i]);
  printf("\n");
}

int main ()
{
  int i, j, N;
  double t, f;
  double *yOld, *yNew;

  // Alocar memoria
  allocMemory(&yOld);
  allocMemory(&yNew);

  // Numero de subintervalos para o plot
  N = 1.0 / dt;
  // Inicializar as condicoes iniciais
  setInitialConditions(yOld);

  for (i = 0; i < N; i++)
  {
    t = i*dt;
    writeSolution(t,yOld,num_eq+num_EDO);
    
    // Resolver as equacoes algebricas
    yNew[0] = Vm(yOld);
    yNew[1] = Pi(yOld);

    // Resolver as EDOs
    for (j = 2; j < num_EDO+num_eq; j++)
    {
      switch (j)
      {
        case 2: {
                  f = alpha_i(yOld);
                  yNew[2] = yOld[2] + f*dt;
                  break;
                }

        case 3: {
                  f = Naidt(yOld,yNew);
                  yNew[3] = yOld[3] + f*dt;
                  break;
                }
        case 4: {
                  f = Kidt(yOld,yNew);
                  yNew[4] = yOld[4] + f*dt;
                  break;
                }
      }
    }
   
    // Passar para a proxima iteracao
    memcpy(yOld,yNew,sizeof(double)*(num_EDO+num_eq));
  }

  free(yOld);
  free(yNew);
  return 0;
}
