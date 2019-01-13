#include "mitchell_shaeffer.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = 0.00000820413566106744f; //V millivolt 
    sv[1] = 0.8789655121804799f;     //h dimensionless 
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

    //State variables
    const real V_old_ = sv[0];
    const real h_old_ = sv[1];

    //Parameters
    const real tau_in = 0.3;
    const real tau_open = 120.0;
    const real tau_close = 150.0;
    const real V_gate = 0.13;
    const real tau_out = 6.0;

    real Jstim = stim_current;
    real Jin = ( h_old_*( pow(V_old_, 2.00000)*(1.00000 - V_old_)))/tau_in;
    real Jout = - (V_old_/tau_out);

    rDY_[0] = Jstim + Jin + Jout;
    rDY_[1] = (V_old_<V_gate ? (1.00000 - h_old_)/tau_open : - h_old_/tau_close);

}

