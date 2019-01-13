#include "beeler_reuter_1977.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -84.624;        // V
    sv[1] = 0.011;          // m
    sv[2] = 0.988;          // h
    sv[3] = 0.975;          // j
    sv[4] = 1e-4;           // Cai
    sv[5] = 0.003;          // d
    sv[6] = 0.994;          // f
    sv[7] = 0.0001;         // x1
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    for (int j = 0; j < num_steps; ++j) 
    {
        solve_model_ode_cpu(dt, sv, stim_currents);
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    // Save old value of the state vector
    for(int i = 0; i < NEQ; i++)
    {
        rY[i] = sv[i];
    }
        
    // Compute Right-hand-side of the ODE's
    RHS_cpu(rY, rDY, stim_current);

    // Solve model using Forward Euler
    for(int i = 0; i < NEQ; i++)
    {
        sv[i] = dt*rDY[i] + rY[i];
    }
        
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    // State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real j_old_ = sv[3];
    const real Cai_old_ = sv[4];
    const real d_old_ = sv[5];
    const real f_old_ = sv[6];
    const real x1_old_ = sv[7];

    // Constants
    const real C = 0.01;
    const real g_na = 4e-2;
    const real E_na = 50;
    const real g_nac = 3e-5;
    const real g_s = 9e-4;

    // Algebraics
    real alpha_m = ( - 1.00000*(V_old_+47.0000))/(exp( - 0.100000*(V_old_+47.0000)) - 1.00000);
    real beta_m =  40.0000*exp( - 0.0560000*(V_old_+72.0000));
    real alpha_h =  0.126000*exp( - 0.250000*(V_old_+77.0000));
    real beta_h = 1.70000/(exp( - 0.0820000*(V_old_+22.5000))+1.00000);
    real alpha_j = ( 0.0550000*exp( - 0.250000*(V_old_+78.0000)))/(exp( - 0.200000*(V_old_+78.0000))+1.00000);
    real beta_j = 0.300000/(exp( - 0.100000*(V_old_+32.0000))+1.00000);
    real alpha_d = ( 0.0950000*exp(- (V_old_ - 5.00000)/100.000))/(1.00000+exp(- (V_old_ - 5.00000)/13.8900));
    real beta_d = ( 0.0700000*exp(- (V_old_+44.0000)/59.0000))/(1.00000+exp((V_old_+44.0000)/20.0000));
    real alpha_f = ( 0.0120000*exp(- (V_old_+28.0000)/125.000))/(1.00000+exp((V_old_+28.0000)/6.67000));
    real beta_f = ( 0.00650000*exp(- (V_old_+30.0000)/50.0000))/(1.00000+exp(- (V_old_+30.0000)/5.00000));
    real alpha_x1 = ( 0.000500000*exp((V_old_+50.0000)/12.1000))/(1.00000+exp((V_old_+50.0000)/17.5000));
    real beta_x1 = ( 0.00130000*exp(- (V_old_+20.0000)/16.6700))/(1.00000+exp(- (V_old_+20.0000)/25.0000));
    real E_s = - 82.3000 -  13.0287*log( Cai_old_*0.00100000);
    real i_s =  g_s*d_old_*f_old_*(V_old_ - E_s);
    real i_na =  ( g_na*pow(m_old_, 3.00000)*h_old_*j_old_+g_nac)*(V_old_ - E_na);
    real i_x1 = ( x1_old_*0.00800000*(exp( 0.0400000*(V_old_+77.0000)) - 1.00000))/exp( 0.0400000*(V_old_+35.0000));
    real i_k1 =  0.00350000*(( 4.00000*(exp( 0.0400000*(V_old_+85.0000)) - 1.00000))/(exp( 0.0800000*(V_old_+53.0000))+exp( 0.0400000*(V_old_+53.0000)))+( 0.200000*(V_old_+23.0000))/(1.00000 - exp( - 0.0400000*(V_old_+23.0000))));
    real i_stim = stim_current;

    // Rates
    rDY_[0] = (i_stim - (i_na+i_s+i_x1+i_k1))/C;
    rDY_[1] = alpha_m*(1.00000 - m_old_) -  beta_m*m_old_;
    rDY_[2] = alpha_h*(1.00000 - h_old_) -  beta_h*h_old_;
    rDY_[3] = alpha_j*(1.00000 - j_old_) -  beta_j*j_old_;
    rDY_[4] = ( - 0.0100000*i_s)/1.00000+ 0.0700000*(0.000100000 - Cai_old_);
    rDY_[5] = alpha_d*(1.00000 - d_old_) -  beta_d*d_old_;
    rDY_[6] = alpha_f*(1.00000 - f_old_) -  beta_f*f_old_;
    rDY_[7] = alpha_x1*(1.00000 - x1_old_) -  beta_x1*x1_old_;

}