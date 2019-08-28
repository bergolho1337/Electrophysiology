#ifndef MONOALG3D_MODEL_ARPF_2009_H
#define MONOALG3D_MODEL_ARPF_2009_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 20
#define INITIAL_V (-69.1370441635924)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real stim_current);

#endif // MONOALG3D_MODEL_ARPF_2009_H

/*
   There are a total of 76 entries in the algebraic variable array.
   There are a total of 20 entries in each of the rate and state variable arrays.
   There are a total of 52 entries in the constant variable array.

 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (microF).
 * CONSTANTS[4] is V_c in component membrane (micrometre3).
 * ALGEBRAIC[51] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[58] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[60] is i_sus in component sustained_outward_current (picoA_per_picoF).
 * ALGEBRAIC[52] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[53] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[56] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[61] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[54] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[55] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[62] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[57] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[64] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[63] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[49] is i_f in component hyperpolarization_activated_current (picoA_per_picoF).
 * ALGEBRAIC[13] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[27] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[36] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[45] is E_Ca in component reversal_potentials (millivolt).
 * CONSTANTS[5] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[6] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[7] is Na_o in component sodium_dynamics (millimolar).
 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * CONSTANTS[8] is Ca_o in component calcium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * ALGEBRAIC[47] is i_f_Na in component hyperpolarization_activated_current (picoA_per_picoF).
 * ALGEBRAIC[48] is i_f_K in component hyperpolarization_activated_current (picoA_per_picoF).
 * CONSTANTS[9] is g_f_Na in component hyperpolarization_activated_current (nanoS_per_picoF).
 * CONSTANTS[10] is g_f_K in component hyperpolarization_activated_current (nanoS_per_picoF).
 * STATES[4] is y in component hyperpolarization_activated_current_y_gate (dimensionless).
 * ALGEBRAIC[0] is y_inf in component hyperpolarization_activated_current_y_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_y in component hyperpolarization_activated_current_y_gate (per_millisecond).
 * ALGEBRAIC[28] is beta_y in component hyperpolarization_activated_current_y_gate (per_millisecond).
 * ALGEBRAIC[37] is tau_y in component hyperpolarization_activated_current_y_gate (millisecond).
 * CONSTANTS[11] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * ALGEBRAIC[50] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * CONSTANTS[12] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[5] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[6] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[1] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[15] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[29] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[38] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[2] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[30] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[39] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * CONSTANTS[13] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[7] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[3] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[17] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[31] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[40] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * CONSTANTS[14] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[8] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[9] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[10] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[4] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[32] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[41] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[5] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[19] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[33] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[42] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[6] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[20] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[34] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[43] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * CONSTANTS[15] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[16] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * STATES[11] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[12] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[13] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[14] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[15] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[7] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[21] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[35] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[44] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[46] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[8] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[22] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[9] is f2_inf in component L_type_Ca_current_f2_gate (dimensionless).
 * ALGEBRAIC[23] is tau_f2 in component L_type_Ca_current_f2_gate (millisecond).
 * ALGEBRAIC[10] is fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[24] is tau_fCass in component L_type_Ca_current_fCass_gate (millisecond).
 * CONSTANTS[17] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[18] is g_to in component transient_outward_current (nanoS_per_picoF).
 * STATES[16] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[17] is r in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[11] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[25] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[12] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[26] is tau_r in component transient_outward_current_r_gate (millisecond).
 * CONSTANTS[19] is g_sus in component sustained_outward_current (nanoS_per_picoF).
 * ALGEBRAIC[59] is a in component sustained_outward_current (dimensionless).
 * CONSTANTS[20] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[21] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[22] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[23] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[24] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[25] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[26] is gamma in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[27] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[28] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[29] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[30] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[31] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * STATES[18] is Ca_SR in component calcium_dynamics (millimolar).
 * ALGEBRAIC[73] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[65] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[66] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[67] is i_xfer in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[72] is O in component calcium_dynamics (dimensionless).
 * STATES[19] is R_prime in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[70] is k1 in component calcium_dynamics (per_millimolar2_per_millisecond).
 * ALGEBRAIC[71] is k2 in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[32] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
 * CONSTANTS[33] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[34] is k3 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[35] is k4 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[36] is EC in component calcium_dynamics (millimolar).
 * CONSTANTS[37] is max_sr in component calcium_dynamics (dimensionless).
 * CONSTANTS[38] is min_sr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[68] is kcasr in component calcium_dynamics (dimensionless).
 * CONSTANTS[39] is V_rel in component calcium_dynamics (per_millisecond).
 * CONSTANTS[40] is V_xfer in component calcium_dynamics (per_millisecond).
 * CONSTANTS[41] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[42] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[43] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[69] is Ca_i_bufc in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[74] is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[75] is Ca_ss_bufss in component calcium_dynamics (dimensionless).
 * CONSTANTS[44] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[45] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[46] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[47] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[48] is Buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[49] is K_buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[50] is V_sr in component calcium_dynamics (micrometre3).
 * CONSTANTS[51] is V_ss in component calcium_dynamics (micrometre3).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt y in component hyperpolarization_activated_current_y_gate (dimensionless).
 * RATES[5] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[6] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[7] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[8] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[9] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[10] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[12] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[13] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[14] is d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * RATES[15] is d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * RATES[16] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[17] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[19] is d/dt R_prime in component calcium_dynamics (dimensionless).
 * RATES[3] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[18] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[11] is d/dt Ca_ss in component calcium_dynamics (millimolar).
 * RATES[2] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[1] is d/dt K_i in component potassium_dynamics (millimolar).
 * 
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = -69.1370441635924;
STATES[1] = 136.781894160227;
STATES[2] = 8.80420286531673;
STATES[3] = 0.000101878186157052;
STATES[4] = 0.0457562667986602;
STATES[5] = 0.00550281999719088;
STATES[6] = 0.313213286437995;
STATES[7] = 0.00953708522974789;
STATES[8] = 0.0417391656294997;
STATES[9] = 0.190678733735145;
STATES[10] = 0.238219836154029;
STATES[11] = 0.000446818714055411;
STATES[12] = 0.000287906256206415;
STATES[13] = 0.989328560287987;
STATES[14] = 0.995474890442185;
STATES[15] = 0.999955429598213;
STATES[16] = 0.96386101799501;
STATES[17] = 0.00103618091196912;
STATES[18] = 3.10836886659417;
STATES[19] = 0.991580051907845;


CONSTANTS[0] = 8314.472;
CONSTANTS[1] = 310;
CONSTANTS[2] = 96485.3415;
CONSTANTS[3] = 0.185;
CONSTANTS[4] = 0.016404;
CONSTANTS[5] = 0.03;
CONSTANTS[6] = 5.4;
CONSTANTS[7] = 140;
CONSTANTS[8] = 2;
CONSTANTS[9] = 0.0145654;
CONSTANTS[10] = 0.0234346;
CONSTANTS[11] = 0.065;
CONSTANTS[12] = 0.0918;
CONSTANTS[13] = 0.2352;
CONSTANTS[14] = 130.5744;
CONSTANTS[15] = 0.00029;
CONSTANTS[16] = 3.98e-5;
CONSTANTS[17] = 0.000592;
CONSTANTS[18] = 0.08184;
CONSTANTS[19] = 0.0227;
CONSTANTS[20] = 2.724;
CONSTANTS[21] = 1;
CONSTANTS[22] = 40;
CONSTANTS[23] = 1000;
CONSTANTS[24] = 0.1;
CONSTANTS[25] = 2.5;
CONSTANTS[26] = 0.35;
CONSTANTS[27] = 1.38;
CONSTANTS[28] = 87.5;
CONSTANTS[29] = 0.1238;
CONSTANTS[30] = 0.0005;
CONSTANTS[31] = 0.0146;
CONSTANTS[32] = 0.15;
CONSTANTS[33] = 0.045;
CONSTANTS[34] = 0.06;
CONSTANTS[35] = 0.005;
CONSTANTS[36] = 1.5;
CONSTANTS[37] = 2.5;
CONSTANTS[38] = 1;
CONSTANTS[39] = 0.102;
CONSTANTS[40] = 0.0038;
CONSTANTS[41] = 0.00025;
CONSTANTS[42] = 0.00036;
CONSTANTS[43] = 0.006375;
CONSTANTS[44] = 0.2;
CONSTANTS[45] = 0.001;
CONSTANTS[46] = 10;
CONSTANTS[47] = 0.3;
CONSTANTS[48] = 0.4;
CONSTANTS[49] = 0.00025;
CONSTANTS[50] = 0.001094;
CONSTANTS[51] = 5.468e-5;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[0] =  (- 1.00000/1.00000)*(ALGEBRAIC[51]+ALGEBRAIC[58]+ALGEBRAIC[60]+ALGEBRAIC[52]+ALGEBRAIC[53]+ALGEBRAIC[56]+ALGEBRAIC[61]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[62]+ALGEBRAIC[57]+ALGEBRAIC[64]+ALGEBRAIC[63]+ALGEBRAIC[49]);
RATES[1] =  (( - 1.00000*((ALGEBRAIC[51]+ALGEBRAIC[58]+ALGEBRAIC[48]+ALGEBRAIC[60]+ALGEBRAIC[52]+ALGEBRAIC[53]+ALGEBRAIC[64]) -  2.00000*ALGEBRAIC[61]))/( 1.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3];
RATES[2] =  (( - 1.00000*(ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[47]+ 3.00000*ALGEBRAIC[61]+ 3.00000*ALGEBRAIC[62]))/( 1.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3];
RATES[3] =  ALGEBRAIC[69]*((( (ALGEBRAIC[66] - ALGEBRAIC[65])*CONSTANTS[50])/CONSTANTS[4]+ALGEBRAIC[67]) - ( 1.00000*((ALGEBRAIC[57]+ALGEBRAIC[63]) -  2.00000*ALGEBRAIC[62])*CONSTANTS[3])/( 2.00000*1.00000*CONSTANTS[4]*CONSTANTS[2]));
RATES[4] = (ALGEBRAIC[0] - STATES[4])/ALGEBRAIC[37];
RATES[5] = (ALGEBRAIC[1] - STATES[5])/ALGEBRAIC[38];
RATES[6] = (ALGEBRAIC[2] - STATES[6])/ALGEBRAIC[39];
RATES[7] = (ALGEBRAIC[3] - STATES[7])/ALGEBRAIC[40];
RATES[8] = (ALGEBRAIC[4] - STATES[8])/ALGEBRAIC[41];
RATES[9] = (ALGEBRAIC[5] - STATES[9])/ALGEBRAIC[42];
RATES[10] = (ALGEBRAIC[6] - STATES[10])/ALGEBRAIC[43];
RATES[11] =  ALGEBRAIC[75]*((( - 1.00000*ALGEBRAIC[56]*CONSTANTS[3])/( 2.00000*1.00000*CONSTANTS[51]*CONSTANTS[2])+( ALGEBRAIC[73]*CONSTANTS[50])/CONSTANTS[51]) - ( ALGEBRAIC[67]*CONSTANTS[4])/CONSTANTS[51]);
RATES[12] = (ALGEBRAIC[7] - STATES[12])/ALGEBRAIC[46];
RATES[13] = (ALGEBRAIC[8] - STATES[13])/ALGEBRAIC[22];
RATES[14] = (ALGEBRAIC[9] - STATES[14])/ALGEBRAIC[23];
RATES[15] = (ALGEBRAIC[10] - STATES[15])/ALGEBRAIC[24];
RATES[16] = (ALGEBRAIC[11] - STATES[16])/ALGEBRAIC[25];
RATES[17] = (ALGEBRAIC[12] - STATES[17])/ALGEBRAIC[26];
RATES[18] =  ALGEBRAIC[74]*(ALGEBRAIC[65] - (ALGEBRAIC[73]+ALGEBRAIC[66]));
RATES[19] =  - ALGEBRAIC[71]*STATES[11]*STATES[19]+ CONSTANTS[35]*(1.00000 - STATES[19]);



ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[22] =  1102.50*exp(- pow(STATES[0]+27.0000, 2.00000)/225.000)+200.000/(1.00000+exp((13.0000 - STATES[0])/10.0000))+180.000/(1.00000+exp((STATES[0]+30.0000)/10.0000))+20.0000;
ALGEBRAIC[9] = 0.670000/(1.00000+exp((STATES[0]+35.0000)/7.00000))+0.330000;
ALGEBRAIC[23] =  562.000*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+31.0000/(1.00000+exp((25.0000 - STATES[0])/10.0000))+80.0000/(1.00000+exp((STATES[0]+30.0000)/10.0000));
ALGEBRAIC[10] = 0.600000/(1.00000+pow(STATES[11]/0.0500000, 2.00000))+0.400000;
ALGEBRAIC[24] = 80.0000/(1.00000+pow(STATES[11]/0.0500000, 2.00000))+2.00000;
ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+27.0000)/13.0000));
ALGEBRAIC[25] =  85.0000*exp(- pow(STATES[0]+25.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 40.0000)/5.00000))+42.0000;
ALGEBRAIC[12] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/13.0000));
ALGEBRAIC[26] =  10.4500*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+7.30000;
ALGEBRAIC[0] = 1.00000/(1.00000+exp((STATES[0]+80.6000)/6.80000));
ALGEBRAIC[14] =  1.00000*exp(- 2.90000 -  0.0400000*STATES[0]);
ALGEBRAIC[28] =  1.00000*exp(3.60000+ 0.110000*STATES[0]);
ALGEBRAIC[37] = 4000.00/(ALGEBRAIC[14]+ALGEBRAIC[28]);
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[15] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[29] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[29];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[16] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[30] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[30];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[17] = 1400.00/ pow((1.00000+exp((5.00000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[31] = 1.00000/(1.00000+exp((STATES[0] - 35.0000)/15.0000));
ALGEBRAIC[40] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[31]+80.0000;
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[18] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[32] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[41] =  1.00000*ALGEBRAIC[18]*ALGEBRAIC[32];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[33] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[42] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[33]);

ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[20] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[34] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[43] = 1.00000/(ALGEBRAIC[20]+ALGEBRAIC[34]);
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 8.00000 - STATES[0])/7.50000));
ALGEBRAIC[21] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[35] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[44] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[46] =  1.00000*ALGEBRAIC[21]*ALGEBRAIC[35]+ALGEBRAIC[44];
ALGEBRAIC[61] = (( (( CONSTANTS[20]*CONSTANTS[6])/(CONSTANTS[6]+CONSTANTS[21]))*STATES[2])/(STATES[2]+CONSTANTS[22]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[7]/STATES[2]);
ALGEBRAIC[54] =  CONSTANTS[14]*pow(STATES[8], 3.00000)*STATES[9]*STATES[10]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[55] =  CONSTANTS[15]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[62] = ( CONSTANTS[23]*( exp(( CONSTANTS[26]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[8] -  exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[7], 3.00000)*STATES[3]*CONSTANTS[25]))/( (pow(CONSTANTS[28], 3.00000)+pow(CONSTANTS[7], 3.00000))*(CONSTANTS[27]+CONSTANTS[8])*(1.00000+ CONSTANTS[24]*exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
ALGEBRAIC[47] =  STATES[4]*CONSTANTS[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[27] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[6]/STATES[1]);
ALGEBRAIC[50] = 1.00000/(1.00000+exp( 0.100000*(STATES[0]+75.4400)));
ALGEBRAIC[51] =  CONSTANTS[11]*ALGEBRAIC[50]*((STATES[0] - 8.00000) - ALGEBRAIC[27]);
ALGEBRAIC[58] =  CONSTANTS[18]*STATES[17]*STATES[16]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[59] = 1.00000/(1.00000+exp((5.00000 - STATES[0])/17.0000));
ALGEBRAIC[60] =  CONSTANTS[19]*ALGEBRAIC[59]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[52] =  CONSTANTS[12]* pow((CONSTANTS[6]/5.40000), 1.0 / 2)*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[36] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[6]+ CONSTANTS[5]*CONSTANTS[7])/(STATES[1]+ CONSTANTS[5]*STATES[2]));
ALGEBRAIC[53] =  CONSTANTS[13]*pow(STATES[7], 2.00000)*(STATES[0] - ALGEBRAIC[36]);
ALGEBRAIC[56] = ( (( CONSTANTS[16]*STATES[12]*STATES[13]*STATES[14]*STATES[15]*4.00000*(STATES[0] - 15.0000)*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( 0.250000*STATES[11]*exp(( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - CONSTANTS[8]))/(exp(( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[45] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[8]/STATES[3]);
ALGEBRAIC[57] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[45]);
ALGEBRAIC[64] = ( CONSTANTS[31]*(STATES[0] - ALGEBRAIC[27]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[63] = ( CONSTANTS[29]*STATES[3])/(STATES[3]+CONSTANTS[30]);
ALGEBRAIC[48] =  STATES[4]*CONSTANTS[10]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[49] = ALGEBRAIC[47]+ALGEBRAIC[48];
ALGEBRAIC[65] = CONSTANTS[43]/(1.00000+pow(CONSTANTS[41], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[66] =  CONSTANTS[42]*(STATES[18] - STATES[3]);
ALGEBRAIC[67] =  CONSTANTS[40]*(STATES[11] - STATES[3]);
ALGEBRAIC[69] = 1.00000/(1.00000+( CONSTANTS[44]*CONSTANTS[45])/pow(STATES[3]+CONSTANTS[45], 2.00000));
ALGEBRAIC[68] = CONSTANTS[37] - (CONSTANTS[37] - CONSTANTS[38])/(1.00000+pow(CONSTANTS[36]/STATES[18], 2.00000));
ALGEBRAIC[71] =  CONSTANTS[33]*ALGEBRAIC[68];
ALGEBRAIC[70] = CONSTANTS[32]/ALGEBRAIC[68];
ALGEBRAIC[72] = ( ALGEBRAIC[70]*pow(STATES[11], 2.00000)*STATES[19])/(CONSTANTS[34]+ ALGEBRAIC[70]*pow(STATES[11], 2.00000));
ALGEBRAIC[73] =  CONSTANTS[39]*ALGEBRAIC[72]*(STATES[18] - STATES[11]);
ALGEBRAIC[74] = 1.00000/(1.00000+( CONSTANTS[46]*CONSTANTS[47])/pow(STATES[18]+CONSTANTS[47], 2.00000));
ALGEBRAIC[75] = 1.00000/(1.00000+( CONSTANTS[48]*CONSTANTS[49])/pow(STATES[11]+CONSTANTS[49], 2.00000));

}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[22] =  1102.50*exp(- pow(STATES[0]+27.0000, 2.00000)/225.000)+200.000/(1.00000+exp((13.0000 - STATES[0])/10.0000))+180.000/(1.00000+exp((STATES[0]+30.0000)/10.0000))+20.0000;
ALGEBRAIC[9] = 0.670000/(1.00000+exp((STATES[0]+35.0000)/7.00000))+0.330000;
ALGEBRAIC[23] =  562.000*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+31.0000/(1.00000+exp((25.0000 - STATES[0])/10.0000))+80.0000/(1.00000+exp((STATES[0]+30.0000)/10.0000));
ALGEBRAIC[10] = 0.600000/(1.00000+pow(STATES[11]/0.0500000, 2.00000))+0.400000;
ALGEBRAIC[24] = 80.0000/(1.00000+pow(STATES[11]/0.0500000, 2.00000))+2.00000;
ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+27.0000)/13.0000));
ALGEBRAIC[25] =  85.0000*exp(- pow(STATES[0]+25.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 40.0000)/5.00000))+42.0000;
ALGEBRAIC[12] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/13.0000));
ALGEBRAIC[26] =  10.4500*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+7.30000;
ALGEBRAIC[0] = 1.00000/(1.00000+exp((STATES[0]+80.6000)/6.80000));
ALGEBRAIC[14] =  1.00000*exp(- 2.90000 -  0.0400000*STATES[0]);
ALGEBRAIC[28] =  1.00000*exp(3.60000+ 0.110000*STATES[0]);
ALGEBRAIC[37] = 4000.00/(ALGEBRAIC[14]+ALGEBRAIC[28]);
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[15] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[29] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[29];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[16] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[30] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[30];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[17] = 1400.00/ pow((1.00000+exp((5.00000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[31] = 1.00000/(1.00000+exp((STATES[0] - 35.0000)/15.0000));
ALGEBRAIC[40] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[31]+80.0000;
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[18] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[32] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[41] =  1.00000*ALGEBRAIC[18]*ALGEBRAIC[32];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[33] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[42] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[33]);
ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[20] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[34] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[43] = 1.00000/(ALGEBRAIC[20]+ALGEBRAIC[34]);
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 8.00000 - STATES[0])/7.50000));
ALGEBRAIC[21] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[35] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[44] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[46] =  1.00000*ALGEBRAIC[21]*ALGEBRAIC[35]+ALGEBRAIC[44];
ALGEBRAIC[61] = (( (( CONSTANTS[20]*CONSTANTS[6])/(CONSTANTS[6]+CONSTANTS[21]))*STATES[2])/(STATES[2]+CONSTANTS[22]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[7]/STATES[2]);
ALGEBRAIC[54] =  CONSTANTS[14]*pow(STATES[8], 3.00000)*STATES[9]*STATES[10]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[55] =  CONSTANTS[15]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[62] = ( CONSTANTS[23]*( exp(( CONSTANTS[26]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[8] -  exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[7], 3.00000)*STATES[3]*CONSTANTS[25]))/( (pow(CONSTANTS[28], 3.00000)+pow(CONSTANTS[7], 3.00000))*(CONSTANTS[27]+CONSTANTS[8])*(1.00000+ CONSTANTS[24]*exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
ALGEBRAIC[47] =  STATES[4]*CONSTANTS[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[27] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[6]/STATES[1]);
ALGEBRAIC[50] = 1.00000/(1.00000+exp( 0.100000*(STATES[0]+75.4400)));
ALGEBRAIC[51] =  CONSTANTS[11]*ALGEBRAIC[50]*((STATES[0] - 8.00000) - ALGEBRAIC[27]);
ALGEBRAIC[58] =  CONSTANTS[18]*STATES[17]*STATES[16]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[59] = 1.00000/(1.00000+exp((5.00000 - STATES[0])/17.0000));
ALGEBRAIC[60] =  CONSTANTS[19]*ALGEBRAIC[59]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[52] =  CONSTANTS[12]* pow((CONSTANTS[6]/5.40000), 1.0 / 2)*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[36] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[6]+ CONSTANTS[5]*CONSTANTS[7])/(STATES[1]+ CONSTANTS[5]*STATES[2]));
ALGEBRAIC[53] =  CONSTANTS[13]*pow(STATES[7], 2.00000)*(STATES[0] - ALGEBRAIC[36]);
ALGEBRAIC[56] = ( (( CONSTANTS[16]*STATES[12]*STATES[13]*STATES[14]*STATES[15]*4.00000*(STATES[0] - 15.0000)*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( 0.250000*STATES[11]*exp(( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - CONSTANTS[8]))/(exp(( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[45] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[8]/STATES[3]);
ALGEBRAIC[57] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[45]);
ALGEBRAIC[64] = ( CONSTANTS[31]*(STATES[0] - ALGEBRAIC[27]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[63] = ( CONSTANTS[29]*STATES[3])/(STATES[3]+CONSTANTS[30]);
ALGEBRAIC[48] =  STATES[4]*CONSTANTS[10]*(STATES[0] - ALGEBRAIC[27]);
ALGEBRAIC[49] = ALGEBRAIC[47]+ALGEBRAIC[48];
ALGEBRAIC[65] = CONSTANTS[43]/(1.00000+pow(CONSTANTS[41], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[66] =  CONSTANTS[42]*(STATES[18] - STATES[3]);
ALGEBRAIC[67] =  CONSTANTS[40]*(STATES[11] - STATES[3]);
ALGEBRAIC[69] = 1.00000/(1.00000+( CONSTANTS[44]*CONSTANTS[45])/pow(STATES[3]+CONSTANTS[45], 2.00000));
ALGEBRAIC[68] = CONSTANTS[37] - (CONSTANTS[37] - CONSTANTS[38])/(1.00000+pow(CONSTANTS[36]/STATES[18], 2.00000));
ALGEBRAIC[71] =  CONSTANTS[33]*ALGEBRAIC[68];
ALGEBRAIC[70] = CONSTANTS[32]/ALGEBRAIC[68];
ALGEBRAIC[72] = ( ALGEBRAIC[70]*pow(STATES[11], 2.00000)*STATES[19])/(CONSTANTS[34]+ ALGEBRAIC[70]*pow(STATES[11], 2.00000));
ALGEBRAIC[73] =  CONSTANTS[39]*ALGEBRAIC[72]*(STATES[18] - STATES[11]);
ALGEBRAIC[74] = 1.00000/(1.00000+( CONSTANTS[46]*CONSTANTS[47])/pow(STATES[18]+CONSTANTS[47], 2.00000));
ALGEBRAIC[75] = 1.00000/(1.00000+( CONSTANTS[48]*CONSTANTS[49])/pow(STATES[11]+CONSTANTS[49], 2.00000));
}
*/