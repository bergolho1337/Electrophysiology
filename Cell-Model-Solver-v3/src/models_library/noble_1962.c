#include "noble_1962.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    // Original: -75.5344986658,0.0605467272,0.7259001355,0.4709239708
    // 300ms pacing steady state: -81.2427,0.0442248,0.852377,0.585062 
    // 280ms pacing steady-state: -81.4871,0.0436287,0.854551,0.597253
    // 260ms pacing steady-state: -81.582,0.043404,0.852091,0.608277
    // 240ms pacing steady-state: -81.4022,0.0438335,0.854737,0.590963
    // 230ms pacing steady-state: 
    // 220ms pacing steady-state: -79.0071 0.050064 0.809459 0.522273
    // 200ms pacing steady-state: -80.0059,0.0473721,0.830419,0.546311
    // 180ms pacing steady-state: -80.8008,0.0453256,0.845488,0.569096
    // 160ms pacing steady-state: -81.5934,0.0433682,0.857897,0.597179
    // 130ms pacing steady-state: -71.9473,0.0738602,0.429134,0.675337
    // 120ms pacing steady-state: -80.8949,0.0451545,0.794726,0.650358
    // 100ms pacing steady-state: -80.3817,0.0464561,0.791357,0.632244
    // 80ms pacing steady-state:   -74.8952,0.0629502,0.553061,0.656434
    
    // Original
    sv[0] = -75.5344986658;     // V millivolt 
    sv[1] = 0.0605467272;       // m dimensionless
    sv[2] = 0.7259001355;       // h dimensionless
    sv[3] = 0.4709239708;       // n dimensionless

    // 300 ms
    //sv[0] = -81.2427;     // V millivolt 
    //sv[1] = 0.0442248;       // m dimensionless
    //sv[2] = 0.852377;       // h dimensionless
    //sv[3] = 0.585062;       // n dimensionless

    // 280 ms
    //sv[0] = -81.4871;     // V millivolt 
    //sv[1] = 0.0436287;       // m dimensionless
    //sv[2] = 0.854551;       // h dimensionless
    //sv[3] = 0.597253;       // n dimensionless

    // 260 ms
    //sv[0] = -81.582;     // V millivolt 
    //sv[1] = 0.043404;       // m dimensionless
    //sv[2] = 0.852091;       // h dimensionless
    //sv[3] = 0.608277;       // n dimensionless

    // 240 ms
    //sv[0] = -81.4022;     // V millivolt 
    //sv[1] = 0.0438335;       // m dimensionless
    //sv[2] = 0.854737;       // h dimensionless
    //sv[3] = 0.590963;       // n dimensionless

    // 220 ms
    //sv[0] = -79.0071;     // V millivolt 
    //sv[1] = 0.050064;       // m dimensionless
    //sv[2] = 0.809459;       // h dimensionless
    //sv[3] = 0.522273;       // n dimensionless           

    // 200 ms
    //sv[0] = -80.0059;     // V millivolt 
    //sv[1] = 0.0473721;       // m dimensionless
    //sv[2] = 0.830419;       // h dimensionless
    //sv[3] = 0.546311;       // n dimensionless

    // 180 ms
    //sv[0] = -80.8008;     // V millivolt 
    //sv[1] = 0.0453256;       // m dimensionless
    //sv[2] = 0.845488;       // h dimensionless
    //sv[3] = 0.569096;       // n dimensionless

    // 160 ms
    //sv[0] = -81.5934;     // V millivolt 
    //sv[1] = 0.0433682;       // m dimensionless
    //sv[2] = 0.857897;       // h dimensionless
    //sv[3] = 0.597179;       // n dimensionless

    // 130 ms
    //sv[0] = -71.9473;     // V millivolt 
    //sv[1] = 0.0738602;       // m dimensionless
    //sv[2] = 0.429134;       // h dimensionless
    //sv[3] = 0.675337;       // n dimensionless

    // 120 ms
    //sv[0] = -80.8949;     // V millivolt 
    //sv[1] = 0.0451545;       // m dimensionless
    //sv[2] = 0.794726;       // h dimensionless
    //sv[3] = 0.650358;       // n dimensionless

    // 100 ms
    //sv[0] = -80.3817;     // V millivolt 
    //sv[1] = 0.0464561;       // m dimensionless
    //sv[2] = 0.791357;       // h dimensionless
    //sv[3] = 0.632244;       // n dimensionless

    // 80 ms
    //sv[0] = -74.8952;     // V millivolt 
    //sv[1] = 0.0629502;       // m dimensionless
    //sv[2] = 0.553061;       // h dimensionless
    //sv[3] = 0.656434;       // n dimensionless
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    for (int j = 0; j < num_steps; ++j) 
    {
        solve_model_ode_cpu(dt, sv, stim_currents);
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(dt, rY, rDY, stim_current);

    // Explicit Euler for solving the transmembrane potential ...
    sv[0] = dt*rDY[0] + rY[0];
    // Rush-Larsen over the gating variables ...
    for(int i = 1; i < NEQ; i++)
        sv[i] = rDY[i];
}

void RHS_cpu(const real dt, const real *sv, real *rDY_, real stim_current) 
{

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real n_old_ = sv[3]; 

    //___________________________________________________________________________
    //Parameters (miliseconds)
    const real Cm = 12.0f;                                 // (microF)
    const real g_na_max = 400.0f;                       // (microS)
    const real E_na = 40.0f;                               // (millivolt)
    const real g_L = 0.075f;                                // (microS)
    const real E_L = -60.0f;                               // (millivolt)

    real calc_I_stim = stim_current;

    // Algebraics
    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    real alpha_h =  ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    real i_na =  (g_na+1.4e-01)*(V_old_ - E_na);
    real beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    real beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    real beta_n =  ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    //real g_K1 =  1.3f*exp((- V_old_ - 90.0000)/50.0000)+ 0.015f*exp((V_old_+90.0000)/60.0000);
    real g_K1 = (((1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)))));
    real g_K2 =  1.2f*pow(n_old_, 4.00000);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    rDY_[0] = ( - (i_na + i_k + i_leak + calc_I_stim)) / Cm;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;

    // Rush-Larsen
    real m_inf = alpha_m / (alpha_m + beta_m);
    real tau_m = 1.0 / (alpha_m + beta_m);
    rDY_[1] = m_inf + ((m_old_ - m_inf)*exp(-dt/tau_m));

    real h_inf = alpha_h / (alpha_h + beta_h);
    real tau_h = 1.0 / (alpha_h + beta_h);
    rDY_[2] = h_inf + ((h_old_ - h_inf)*exp(-dt/tau_h));

    real n_inf = alpha_n / (alpha_n + beta_n);
    real tau_n = 1.0 / (alpha_n + beta_n);
    rDY_[3] = n_inf + ((n_old_ - n_inf)*exp(-dt/tau_n));
    
    //rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  (beta_m*m_old_) );
    //rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  (beta_h*h_old_) );
    //rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  (beta_n*n_old_) );

}
