# Size of variable arrays:
sizeAlgebraic = 76
sizeStates = 20
sizeConstants = 52
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "R in component membrane (joule_per_mole_kelvin)"
    legend_constants[1] = "T in component membrane (kelvin)"
    legend_constants[2] = "F in component membrane (coulomb_per_millimole)"
    legend_constants[3] = "Cm in component membrane (microF)"
    legend_constants[4] = "V_c in component membrane (micrometre3)"
    legend_algebraic[51] = "i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF)"
    legend_algebraic[58] = "i_to in component transient_outward_current (picoA_per_picoF)"
    legend_algebraic[60] = "i_sus in component sustained_outward_current (picoA_per_picoF)"
    legend_algebraic[52] = "i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF)"
    legend_algebraic[53] = "i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF)"
    legend_algebraic[56] = "i_CaL in component L_type_Ca_current (picoA_per_picoF)"
    legend_algebraic[61] = "i_NaK in component sodium_potassium_pump_current (picoA_per_picoF)"
    legend_algebraic[54] = "i_Na in component fast_sodium_current (picoA_per_picoF)"
    legend_algebraic[55] = "i_b_Na in component sodium_background_current (picoA_per_picoF)"
    legend_algebraic[62] = "i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF)"
    legend_algebraic[57] = "i_b_Ca in component calcium_background_current (picoA_per_picoF)"
    legend_algebraic[64] = "i_p_K in component potassium_pump_current (picoA_per_picoF)"
    legend_algebraic[63] = "i_p_Ca in component calcium_pump_current (picoA_per_picoF)"
    legend_algebraic[49] = "i_f in component hyperpolarization_activated_current (picoA_per_picoF)"
    legend_algebraic[13] = "E_Na in component reversal_potentials (millivolt)"
    legend_algebraic[27] = "E_K in component reversal_potentials (millivolt)"
    legend_algebraic[36] = "E_Ks in component reversal_potentials (millivolt)"
    legend_algebraic[45] = "E_Ca in component reversal_potentials (millivolt)"
    legend_constants[5] = "P_kna in component reversal_potentials (dimensionless)"
    legend_constants[6] = "K_o in component potassium_dynamics (millimolar)"
    legend_constants[7] = "Na_o in component sodium_dynamics (millimolar)"
    legend_states[1] = "K_i in component potassium_dynamics (millimolar)"
    legend_states[2] = "Na_i in component sodium_dynamics (millimolar)"
    legend_constants[8] = "Ca_o in component calcium_dynamics (millimolar)"
    legend_states[3] = "Ca_i in component calcium_dynamics (millimolar)"
    legend_algebraic[47] = "i_f_Na in component hyperpolarization_activated_current (picoA_per_picoF)"
    legend_algebraic[48] = "i_f_K in component hyperpolarization_activated_current (picoA_per_picoF)"
    legend_constants[9] = "g_f_Na in component hyperpolarization_activated_current (nanoS_per_picoF)"
    legend_constants[10] = "g_f_K in component hyperpolarization_activated_current (nanoS_per_picoF)"
    legend_states[4] = "y in component hyperpolarization_activated_current_y_gate (dimensionless)"
    legend_algebraic[0] = "y_inf in component hyperpolarization_activated_current_y_gate (dimensionless)"
    legend_algebraic[14] = "alpha_y in component hyperpolarization_activated_current_y_gate (per_millisecond)"
    legend_algebraic[28] = "beta_y in component hyperpolarization_activated_current_y_gate (per_millisecond)"
    legend_algebraic[37] = "tau_y in component hyperpolarization_activated_current_y_gate (millisecond)"
    legend_constants[11] = "g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF)"
    legend_algebraic[50] = "xK1_inf in component inward_rectifier_potassium_current (dimensionless)"
    legend_constants[12] = "g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF)"
    legend_states[5] = "Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)"
    legend_states[6] = "Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)"
    legend_algebraic[1] = "xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)"
    legend_algebraic[15] = "alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)"
    legend_algebraic[29] = "beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)"
    legend_algebraic[38] = "tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond)"
    legend_algebraic[2] = "xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)"
    legend_algebraic[16] = "alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)"
    legend_algebraic[30] = "beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)"
    legend_algebraic[39] = "tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond)"
    legend_constants[13] = "g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF)"
    legend_states[7] = "Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)"
    legend_algebraic[3] = "xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)"
    legend_algebraic[17] = "alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)"
    legend_algebraic[31] = "beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)"
    legend_algebraic[40] = "tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond)"
    legend_constants[14] = "g_Na in component fast_sodium_current (nanoS_per_picoF)"
    legend_states[8] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[9] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_states[10] = "j in component fast_sodium_current_j_gate (dimensionless)"
    legend_algebraic[4] = "m_inf in component fast_sodium_current_m_gate (dimensionless)"
    legend_algebraic[18] = "alpha_m in component fast_sodium_current_m_gate (dimensionless)"
    legend_algebraic[32] = "beta_m in component fast_sodium_current_m_gate (dimensionless)"
    legend_algebraic[41] = "tau_m in component fast_sodium_current_m_gate (millisecond)"
    legend_algebraic[5] = "h_inf in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[19] = "alpha_h in component fast_sodium_current_h_gate (per_millisecond)"
    legend_algebraic[33] = "beta_h in component fast_sodium_current_h_gate (per_millisecond)"
    legend_algebraic[42] = "tau_h in component fast_sodium_current_h_gate (millisecond)"
    legend_algebraic[6] = "j_inf in component fast_sodium_current_j_gate (dimensionless)"
    legend_algebraic[20] = "alpha_j in component fast_sodium_current_j_gate (per_millisecond)"
    legend_algebraic[34] = "beta_j in component fast_sodium_current_j_gate (per_millisecond)"
    legend_algebraic[43] = "tau_j in component fast_sodium_current_j_gate (millisecond)"
    legend_constants[15] = "g_bna in component sodium_background_current (nanoS_per_picoF)"
    legend_constants[16] = "g_CaL in component L_type_Ca_current (litre_per_farad_second)"
    legend_states[11] = "Ca_ss in component calcium_dynamics (millimolar)"
    legend_states[12] = "d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_states[13] = "f in component L_type_Ca_current_f_gate (dimensionless)"
    legend_states[14] = "f2 in component L_type_Ca_current_f2_gate (dimensionless)"
    legend_states[15] = "fCass in component L_type_Ca_current_fCass_gate (dimensionless)"
    legend_algebraic[7] = "d_inf in component L_type_Ca_current_d_gate (dimensionless)"
    legend_algebraic[21] = "alpha_d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_algebraic[35] = "beta_d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_algebraic[44] = "gamma_d in component L_type_Ca_current_d_gate (millisecond)"
    legend_algebraic[46] = "tau_d in component L_type_Ca_current_d_gate (millisecond)"
    legend_algebraic[8] = "f_inf in component L_type_Ca_current_f_gate (dimensionless)"
    legend_algebraic[22] = "tau_f in component L_type_Ca_current_f_gate (millisecond)"
    legend_algebraic[9] = "f2_inf in component L_type_Ca_current_f2_gate (dimensionless)"
    legend_algebraic[23] = "tau_f2 in component L_type_Ca_current_f2_gate (millisecond)"
    legend_algebraic[10] = "fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless)"
    legend_algebraic[24] = "tau_fCass in component L_type_Ca_current_fCass_gate (millisecond)"
    legend_constants[17] = "g_bca in component calcium_background_current (nanoS_per_picoF)"
    legend_constants[18] = "g_to in component transient_outward_current (nanoS_per_picoF)"
    legend_states[16] = "s in component transient_outward_current_s_gate (dimensionless)"
    legend_states[17] = "r in component transient_outward_current_r_gate (dimensionless)"
    legend_algebraic[11] = "s_inf in component transient_outward_current_s_gate (dimensionless)"
    legend_algebraic[25] = "tau_s in component transient_outward_current_s_gate (millisecond)"
    legend_algebraic[12] = "r_inf in component transient_outward_current_r_gate (dimensionless)"
    legend_algebraic[26] = "tau_r in component transient_outward_current_r_gate (millisecond)"
    legend_constants[19] = "g_sus in component sustained_outward_current (nanoS_per_picoF)"
    legend_algebraic[59] = "a in component sustained_outward_current (dimensionless)"
    legend_constants[20] = "P_NaK in component sodium_potassium_pump_current (picoA_per_picoF)"
    legend_constants[21] = "K_mk in component sodium_potassium_pump_current (millimolar)"
    legend_constants[22] = "K_mNa in component sodium_potassium_pump_current (millimolar)"
    legend_constants[23] = "K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF)"
    legend_constants[24] = "K_sat in component sodium_calcium_exchanger_current (dimensionless)"
    legend_constants[25] = "alpha in component sodium_calcium_exchanger_current (dimensionless)"
    legend_constants[26] = "gamma in component sodium_calcium_exchanger_current (dimensionless)"
    legend_constants[27] = "Km_Ca in component sodium_calcium_exchanger_current (millimolar)"
    legend_constants[28] = "Km_Nai in component sodium_calcium_exchanger_current (millimolar)"
    legend_constants[29] = "g_pCa in component calcium_pump_current (picoA_per_picoF)"
    legend_constants[30] = "K_pCa in component calcium_pump_current (millimolar)"
    legend_constants[31] = "g_pK in component potassium_pump_current (nanoS_per_picoF)"
    legend_states[18] = "Ca_SR in component calcium_dynamics (millimolar)"
    legend_algebraic[73] = "i_rel in component calcium_dynamics (millimolar_per_millisecond)"
    legend_algebraic[65] = "i_up in component calcium_dynamics (millimolar_per_millisecond)"
    legend_algebraic[66] = "i_leak in component calcium_dynamics (millimolar_per_millisecond)"
    legend_algebraic[67] = "i_xfer in component calcium_dynamics (millimolar_per_millisecond)"
    legend_algebraic[72] = "O in component calcium_dynamics (dimensionless)"
    legend_states[19] = "R_prime in component calcium_dynamics (dimensionless)"
    legend_algebraic[70] = "k1 in component calcium_dynamics (per_millimolar2_per_millisecond)"
    legend_algebraic[71] = "k2 in component calcium_dynamics (per_millimolar_per_millisecond)"
    legend_constants[32] = "k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond)"
    legend_constants[33] = "k2_prime in component calcium_dynamics (per_millimolar_per_millisecond)"
    legend_constants[34] = "k3 in component calcium_dynamics (per_millisecond)"
    legend_constants[35] = "k4 in component calcium_dynamics (per_millisecond)"
    legend_constants[36] = "EC in component calcium_dynamics (millimolar)"
    legend_constants[37] = "max_sr in component calcium_dynamics (dimensionless)"
    legend_constants[38] = "min_sr in component calcium_dynamics (dimensionless)"
    legend_algebraic[68] = "kcasr in component calcium_dynamics (dimensionless)"
    legend_constants[39] = "V_rel in component calcium_dynamics (per_millisecond)"
    legend_constants[40] = "V_xfer in component calcium_dynamics (per_millisecond)"
    legend_constants[41] = "K_up in component calcium_dynamics (millimolar)"
    legend_constants[42] = "V_leak in component calcium_dynamics (per_millisecond)"
    legend_constants[43] = "Vmax_up in component calcium_dynamics (millimolar_per_millisecond)"
    legend_algebraic[69] = "Ca_i_bufc in component calcium_dynamics (dimensionless)"
    legend_algebraic[74] = "Ca_sr_bufsr in component calcium_dynamics (dimensionless)"
    legend_algebraic[75] = "Ca_ss_bufss in component calcium_dynamics (dimensionless)"
    legend_constants[44] = "Buf_c in component calcium_dynamics (millimolar)"
    legend_constants[45] = "K_buf_c in component calcium_dynamics (millimolar)"
    legend_constants[46] = "Buf_sr in component calcium_dynamics (millimolar)"
    legend_constants[47] = "K_buf_sr in component calcium_dynamics (millimolar)"
    legend_constants[48] = "Buf_ss in component calcium_dynamics (millimolar)"
    legend_constants[49] = "K_buf_ss in component calcium_dynamics (millimolar)"
    legend_constants[50] = "V_sr in component calcium_dynamics (micrometre3)"
    legend_constants[51] = "V_ss in component calcium_dynamics (micrometre3)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[4] = "d/dt y in component hyperpolarization_activated_current_y_gate (dimensionless)"
    legend_rates[5] = "d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)"
    legend_rates[6] = "d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)"
    legend_rates[7] = "d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)"
    legend_rates[8] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[9] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[10] = "d/dt j in component fast_sodium_current_j_gate (dimensionless)"
    legend_rates[12] = "d/dt d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_rates[13] = "d/dt f in component L_type_Ca_current_f_gate (dimensionless)"
    legend_rates[14] = "d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless)"
    legend_rates[15] = "d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless)"
    legend_rates[16] = "d/dt s in component transient_outward_current_s_gate (dimensionless)"
    legend_rates[17] = "d/dt r in component transient_outward_current_r_gate (dimensionless)"
    legend_rates[19] = "d/dt R_prime in component calcium_dynamics (dimensionless)"
    legend_rates[3] = "d/dt Ca_i in component calcium_dynamics (millimolar)"
    legend_rates[18] = "d/dt Ca_SR in component calcium_dynamics (millimolar)"
    legend_rates[11] = "d/dt Ca_ss in component calcium_dynamics (millimolar)"
    legend_rates[2] = "d/dt Na_i in component sodium_dynamics (millimolar)"
    legend_rates[1] = "d/dt K_i in component potassium_dynamics (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -69.1370441635924
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 0.185
    constants[4] = 0.016404
    constants[5] = 0.03
    constants[6] = 5.4
    constants[7] = 140
    states[1] = 136.781894160227
    states[2] = 8.80420286531673
    constants[8] = 2
    states[3] = 0.000101878186157052
    constants[9] = 0.0145654
    constants[10] = 0.0234346
    states[4] = 0.0457562667986602
    constants[11] = 0.065
    constants[12] = 0.0918
    states[5] = 0.00550281999719088
    states[6] = 0.313213286437995
    constants[13] = 0.2352
    states[7] = 0.00953708522974789
    constants[14] = 130.5744
    states[8] = 0.0417391656294997
    states[9] = 0.190678733735145
    states[10] = 0.238219836154029
    constants[15] = 0.00029
    constants[16] = 3.98e-5
    states[11] = 0.000446818714055411
    states[12] = 0.000287906256206415
    states[13] = 0.989328560287987
    states[14] = 0.995474890442185
    states[15] = 0.999955429598213
    constants[17] = 0.000592
    constants[18] = 0.08184
    states[16] = 0.96386101799501
    states[17] = 0.00103618091196912
    constants[19] = 0.0227
    constants[20] = 2.724
    constants[21] = 1
    constants[22] = 40
    constants[23] = 1000
    constants[24] = 0.1
    constants[25] = 2.5
    constants[26] = 0.35
    constants[27] = 1.38
    constants[28] = 87.5
    constants[29] = 0.1238
    constants[30] = 0.0005
    constants[31] = 0.0146
    states[18] = 3.10836886659417
    states[19] = 0.991580051907845
    constants[32] = 0.15
    constants[33] = 0.045
    constants[34] = 0.06
    constants[35] = 0.005
    constants[36] = 1.5
    constants[37] = 2.5
    constants[38] = 1
    constants[39] = 0.102
    constants[40] = 0.0038
    constants[41] = 0.00025
    constants[42] = 0.00036
    constants[43] = 0.006375
    constants[44] = 0.2
    constants[45] = 0.001
    constants[46] = 10
    constants[47] = 0.3
    constants[48] = 0.4
    constants[49] = 0.00025
    constants[50] = 0.001094
    constants[51] = 5.468e-5
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+20.0000)/7.00000))
    algebraic[22] = 1102.50*exp(-(power(states[0]+27.0000, 2.00000))/225.000)+200.000/(1.00000+exp((13.0000-states[0])/10.0000))+180.000/(1.00000+exp((states[0]+30.0000)/10.0000))+20.0000
    rates[13] = (algebraic[8]-states[13])/algebraic[22]
    algebraic[9] = 0.670000/(1.00000+exp((states[0]+35.0000)/7.00000))+0.330000
    algebraic[23] = 562.000*exp(-(power(states[0]+27.0000, 2.00000))/240.000)+31.0000/(1.00000+exp((25.0000-states[0])/10.0000))+80.0000/(1.00000+exp((states[0]+30.0000)/10.0000))
    rates[14] = (algebraic[9]-states[14])/algebraic[23]
    algebraic[10] = 0.600000/(1.00000+power(states[11]/0.0500000, 2.00000))+0.400000
    algebraic[24] = 80.0000/(1.00000+power(states[11]/0.0500000, 2.00000))+2.00000
    rates[15] = (algebraic[10]-states[15])/algebraic[24]
    algebraic[11] = 1.00000/(1.00000+exp((states[0]+27.0000)/13.0000))
    algebraic[25] = 85.0000*exp(-(power(states[0]+25.0000, 2.00000))/320.000)+5.00000/(1.00000+exp((states[0]-40.0000)/5.00000))+42.0000
    rates[16] = (algebraic[11]-states[16])/algebraic[25]
    algebraic[12] = 1.00000/(1.00000+exp((20.0000-states[0])/13.0000))
    algebraic[26] = 10.4500*exp(-(power(states[0]+40.0000, 2.00000))/1800.00)+7.30000
    rates[17] = (algebraic[12]-states[17])/algebraic[26]
    algebraic[0] = 1.00000/(1.00000+exp((states[0]+80.6000)/6.80000))
    algebraic[14] = 1.00000*exp(-2.90000-0.0400000*states[0])
    algebraic[28] = 1.00000*exp(3.60000+0.110000*states[0])
    algebraic[37] = 4000.00/(algebraic[14]+algebraic[28])
    rates[4] = (algebraic[0]-states[4])/algebraic[37]
    algebraic[1] = 1.00000/(1.00000+exp((-26.0000-states[0])/7.00000))
    algebraic[15] = 450.000/(1.00000+exp((-45.0000-states[0])/10.0000))
    algebraic[29] = 6.00000/(1.00000+exp((states[0]+30.0000)/11.5000))
    algebraic[38] = 1.00000*algebraic[15]*algebraic[29]
    rates[5] = (algebraic[1]-states[5])/algebraic[38]
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+88.0000)/24.0000))
    algebraic[16] = 3.00000/(1.00000+exp((-60.0000-states[0])/20.0000))
    algebraic[30] = 1.12000/(1.00000+exp((states[0]-60.0000)/20.0000))
    algebraic[39] = 1.00000*algebraic[16]*algebraic[30]
    rates[6] = (algebraic[2]-states[6])/algebraic[39]
    algebraic[3] = 1.00000/(1.00000+exp((-5.00000-states[0])/14.0000))
    algebraic[17] = 1400.00/(power(1.00000+exp((5.00000-states[0])/6.00000), 1.0/2))
    algebraic[31] = 1.00000/(1.00000+exp((states[0]-35.0000)/15.0000))
    algebraic[40] = 1.00000*algebraic[17]*algebraic[31]+80.0000
    rates[7] = (algebraic[3]-states[7])/algebraic[40]
    algebraic[4] = 1.00000/(power(1.00000+exp((-56.8600-states[0])/9.03000), 2.00000))
    algebraic[18] = 1.00000/(1.00000+exp((-60.0000-states[0])/5.00000))
    algebraic[32] = 0.100000/(1.00000+exp((states[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((states[0]-50.0000)/200.000))
    algebraic[41] = 1.00000*algebraic[18]*algebraic[32]
    rates[8] = (algebraic[4]-states[8])/algebraic[41]
    algebraic[5] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[19] = custom_piecewise([less(states[0] , -40.0000), 0.0570000*exp(-(states[0]+80.0000)/6.80000) , True, 0.00000])
    algebraic[33] = custom_piecewise([less(states[0] , -40.0000), 2.70000*exp(0.0790000*states[0])+310000.*exp(0.348500*states[0]) , True, 0.770000/(0.130000*(1.00000+exp((states[0]+10.6600)/-11.1000)))])
    algebraic[42] = 1.00000/(algebraic[19]+algebraic[33])
    rates[9] = (algebraic[5]-states[9])/algebraic[42]
    algebraic[6] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[20] = custom_piecewise([less(states[0] , -40.0000), (((-25428.0*exp(0.244400*states[0])-6.94800e-06*exp(-0.0439100*states[0]))*(states[0]+37.7800))/1.00000)/(1.00000+exp(0.311000*(states[0]+79.2300))) , True, 0.00000])
    algebraic[34] = custom_piecewise([less(states[0] , -40.0000), (0.0242400*exp(-0.0105200*states[0]))/(1.00000+exp(-0.137800*(states[0]+40.1400))) , True, (0.600000*exp(0.0570000*states[0]))/(1.00000+exp(-0.100000*(states[0]+32.0000)))])
    algebraic[43] = 1.00000/(algebraic[20]+algebraic[34])
    rates[10] = (algebraic[6]-states[10])/algebraic[43]
    algebraic[7] = 1.00000/(1.00000+exp((-8.00000-states[0])/7.50000))
    algebraic[21] = 1.40000/(1.00000+exp((-35.0000-states[0])/13.0000))+0.250000
    algebraic[35] = 1.40000/(1.00000+exp((states[0]+5.00000)/5.00000))
    algebraic[44] = 1.00000/(1.00000+exp((50.0000-states[0])/20.0000))
    algebraic[46] = 1.00000*algebraic[21]*algebraic[35]+algebraic[44]
    rates[12] = (algebraic[7]-states[12])/algebraic[46]
    algebraic[61] = ((((constants[20]*constants[6])/(constants[6]+constants[21]))*states[2])/(states[2]+constants[22]))/(1.00000+0.124500*exp((-0.100000*states[0]*constants[2])/(constants[0]*constants[1]))+0.0353000*exp((-states[0]*constants[2])/(constants[0]*constants[1])))
    algebraic[13] = ((constants[0]*constants[1])/constants[2])*log(constants[7]/states[2])
    algebraic[54] = constants[14]*(power(states[8], 3.00000))*states[9]*states[10]*(states[0]-algebraic[13])
    algebraic[55] = constants[15]*(states[0]-algebraic[13])
    algebraic[62] = (constants[23]*(exp((constants[26]*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], 3.00000))*constants[8]-exp(((constants[26]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[7], 3.00000))*states[3]*constants[25]))/((power(constants[28], 3.00000)+power(constants[7], 3.00000))*(constants[27]+constants[8])*(1.00000+constants[24]*exp(((constants[26]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))))
    algebraic[47] = states[4]*constants[9]*(states[0]-algebraic[13])
    rates[2] = ((-1.00000*(algebraic[54]+algebraic[55]+algebraic[47]+3.00000*algebraic[61]+3.00000*algebraic[62]))/(1.00000*constants[4]*constants[2]))*constants[3]
    algebraic[27] = ((constants[0]*constants[1])/constants[2])*log(constants[6]/states[1])
    algebraic[50] = 1.00000/(1.00000+exp(0.100000*(states[0]+75.4400)))
    algebraic[51] = constants[11]*algebraic[50]*((states[0]-8.00000)-algebraic[27])
    algebraic[58] = constants[18]*states[17]*states[16]*(states[0]-algebraic[27])
    algebraic[59] = 1.00000/(1.00000+exp((5.00000-states[0])/17.0000))
    algebraic[60] = constants[19]*algebraic[59]*(states[0]-algebraic[27])
    algebraic[52] = constants[12]*(power(constants[6]/5.40000, 1.0/2))*states[5]*states[6]*(states[0]-algebraic[27])
    algebraic[36] = ((constants[0]*constants[1])/constants[2])*log((constants[6]+constants[5]*constants[7])/(states[1]+constants[5]*states[2]))
    algebraic[53] = constants[13]*(power(states[7], 2.00000))*(states[0]-algebraic[36])
    algebraic[56] = (((constants[16]*states[12]*states[13]*states[14]*states[15]*4.00000*(states[0]-15.0000)*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(0.250000*states[11]*exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-constants[8]))/(exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[45] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[8]/states[3])
    algebraic[57] = constants[17]*(states[0]-algebraic[45])
    algebraic[64] = (constants[31]*(states[0]-algebraic[27]))/(1.00000+exp((25.0000-states[0])/5.98000))
    algebraic[63] = (constants[29]*states[3])/(states[3]+constants[30])
    algebraic[48] = states[4]*constants[10]*(states[0]-algebraic[27])
    algebraic[49] = algebraic[47]+algebraic[48]
    rates[0] = (-1.00000/1.00000)*(algebraic[51]+algebraic[58]+algebraic[60]+algebraic[52]+algebraic[53]+algebraic[56]+algebraic[61]+algebraic[54]+algebraic[55]+algebraic[62]+algebraic[57]+algebraic[64]+algebraic[63]+algebraic[49])
    rates[1] = ((-1.00000*((algebraic[51]+algebraic[58]+algebraic[48]+algebraic[60]+algebraic[52]+algebraic[53]+algebraic[64])-2.00000*algebraic[61]))/(1.00000*constants[4]*constants[2]))*constants[3]
    algebraic[65] = constants[43]/(1.00000+(power(constants[41], 2.00000))/(power(states[3], 2.00000)))
    algebraic[66] = constants[42]*(states[18]-states[3])
    algebraic[67] = constants[40]*(states[11]-states[3])
    algebraic[69] = 1.00000/(1.00000+(constants[44]*constants[45])/(power(states[3]+constants[45], 2.00000)))
    rates[3] = algebraic[69]*((((algebraic[66]-algebraic[65])*constants[50])/constants[4]+algebraic[67])-(1.00000*((algebraic[57]+algebraic[63])-2.00000*algebraic[62])*constants[3])/(2.00000*1.00000*constants[4]*constants[2]))
    algebraic[68] = constants[37]-(constants[37]-constants[38])/(1.00000+power(constants[36]/states[18], 2.00000))
    algebraic[71] = constants[33]*algebraic[68]
    rates[19] = -algebraic[71]*states[11]*states[19]+constants[35]*(1.00000-states[19])
    algebraic[70] = constants[32]/algebraic[68]
    algebraic[72] = (algebraic[70]*(power(states[11], 2.00000))*states[19])/(constants[34]+algebraic[70]*(power(states[11], 2.00000)))
    algebraic[73] = constants[39]*algebraic[72]*(states[18]-states[11])
    algebraic[74] = 1.00000/(1.00000+(constants[46]*constants[47])/(power(states[18]+constants[47], 2.00000)))
    rates[18] = algebraic[74]*(algebraic[65]-(algebraic[73]+algebraic[66]))
    algebraic[75] = 1.00000/(1.00000+(constants[48]*constants[49])/(power(states[11]+constants[49], 2.00000)))
    rates[11] = algebraic[75]*(((-1.00000*algebraic[56]*constants[3])/(2.00000*1.00000*constants[51]*constants[2])+(algebraic[73]*constants[50])/constants[51])-(algebraic[67]*constants[4])/constants[51])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+20.0000)/7.00000))
    algebraic[22] = 1102.50*exp(-(power(states[0]+27.0000, 2.00000))/225.000)+200.000/(1.00000+exp((13.0000-states[0])/10.0000))+180.000/(1.00000+exp((states[0]+30.0000)/10.0000))+20.0000
    algebraic[9] = 0.670000/(1.00000+exp((states[0]+35.0000)/7.00000))+0.330000
    algebraic[23] = 562.000*exp(-(power(states[0]+27.0000, 2.00000))/240.000)+31.0000/(1.00000+exp((25.0000-states[0])/10.0000))+80.0000/(1.00000+exp((states[0]+30.0000)/10.0000))
    algebraic[10] = 0.600000/(1.00000+power(states[11]/0.0500000, 2.00000))+0.400000
    algebraic[24] = 80.0000/(1.00000+power(states[11]/0.0500000, 2.00000))+2.00000
    algebraic[11] = 1.00000/(1.00000+exp((states[0]+27.0000)/13.0000))
    algebraic[25] = 85.0000*exp(-(power(states[0]+25.0000, 2.00000))/320.000)+5.00000/(1.00000+exp((states[0]-40.0000)/5.00000))+42.0000
    algebraic[12] = 1.00000/(1.00000+exp((20.0000-states[0])/13.0000))
    algebraic[26] = 10.4500*exp(-(power(states[0]+40.0000, 2.00000))/1800.00)+7.30000
    algebraic[0] = 1.00000/(1.00000+exp((states[0]+80.6000)/6.80000))
    algebraic[14] = 1.00000*exp(-2.90000-0.0400000*states[0])
    algebraic[28] = 1.00000*exp(3.60000+0.110000*states[0])
    algebraic[37] = 4000.00/(algebraic[14]+algebraic[28])
    algebraic[1] = 1.00000/(1.00000+exp((-26.0000-states[0])/7.00000))
    algebraic[15] = 450.000/(1.00000+exp((-45.0000-states[0])/10.0000))
    algebraic[29] = 6.00000/(1.00000+exp((states[0]+30.0000)/11.5000))
    algebraic[38] = 1.00000*algebraic[15]*algebraic[29]
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+88.0000)/24.0000))
    algebraic[16] = 3.00000/(1.00000+exp((-60.0000-states[0])/20.0000))
    algebraic[30] = 1.12000/(1.00000+exp((states[0]-60.0000)/20.0000))
    algebraic[39] = 1.00000*algebraic[16]*algebraic[30]
    algebraic[3] = 1.00000/(1.00000+exp((-5.00000-states[0])/14.0000))
    algebraic[17] = 1400.00/(power(1.00000+exp((5.00000-states[0])/6.00000), 1.0/2))
    algebraic[31] = 1.00000/(1.00000+exp((states[0]-35.0000)/15.0000))
    algebraic[40] = 1.00000*algebraic[17]*algebraic[31]+80.0000
    algebraic[4] = 1.00000/(power(1.00000+exp((-56.8600-states[0])/9.03000), 2.00000))
    algebraic[18] = 1.00000/(1.00000+exp((-60.0000-states[0])/5.00000))
    algebraic[32] = 0.100000/(1.00000+exp((states[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((states[0]-50.0000)/200.000))
    algebraic[41] = 1.00000*algebraic[18]*algebraic[32]
    algebraic[5] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[19] = custom_piecewise([less(states[0] , -40.0000), 0.0570000*exp(-(states[0]+80.0000)/6.80000) , True, 0.00000])
    algebraic[33] = custom_piecewise([less(states[0] , -40.0000), 2.70000*exp(0.0790000*states[0])+310000.*exp(0.348500*states[0]) , True, 0.770000/(0.130000*(1.00000+exp((states[0]+10.6600)/-11.1000)))])
    algebraic[42] = 1.00000/(algebraic[19]+algebraic[33])
    algebraic[6] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[20] = custom_piecewise([less(states[0] , -40.0000), (((-25428.0*exp(0.244400*states[0])-6.94800e-06*exp(-0.0439100*states[0]))*(states[0]+37.7800))/1.00000)/(1.00000+exp(0.311000*(states[0]+79.2300))) , True, 0.00000])
    algebraic[34] = custom_piecewise([less(states[0] , -40.0000), (0.0242400*exp(-0.0105200*states[0]))/(1.00000+exp(-0.137800*(states[0]+40.1400))) , True, (0.600000*exp(0.0570000*states[0]))/(1.00000+exp(-0.100000*(states[0]+32.0000)))])
    algebraic[43] = 1.00000/(algebraic[20]+algebraic[34])
    algebraic[7] = 1.00000/(1.00000+exp((-8.00000-states[0])/7.50000))
    algebraic[21] = 1.40000/(1.00000+exp((-35.0000-states[0])/13.0000))+0.250000
    algebraic[35] = 1.40000/(1.00000+exp((states[0]+5.00000)/5.00000))
    algebraic[44] = 1.00000/(1.00000+exp((50.0000-states[0])/20.0000))
    algebraic[46] = 1.00000*algebraic[21]*algebraic[35]+algebraic[44]
    algebraic[61] = ((((constants[20]*constants[6])/(constants[6]+constants[21]))*states[2])/(states[2]+constants[22]))/(1.00000+0.124500*exp((-0.100000*states[0]*constants[2])/(constants[0]*constants[1]))+0.0353000*exp((-states[0]*constants[2])/(constants[0]*constants[1])))
    algebraic[13] = ((constants[0]*constants[1])/constants[2])*log(constants[7]/states[2])
    algebraic[54] = constants[14]*(power(states[8], 3.00000))*states[9]*states[10]*(states[0]-algebraic[13])
    algebraic[55] = constants[15]*(states[0]-algebraic[13])
    algebraic[62] = (constants[23]*(exp((constants[26]*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], 3.00000))*constants[8]-exp(((constants[26]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[7], 3.00000))*states[3]*constants[25]))/((power(constants[28], 3.00000)+power(constants[7], 3.00000))*(constants[27]+constants[8])*(1.00000+constants[24]*exp(((constants[26]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))))
    algebraic[47] = states[4]*constants[9]*(states[0]-algebraic[13])
    algebraic[27] = ((constants[0]*constants[1])/constants[2])*log(constants[6]/states[1])
    algebraic[50] = 1.00000/(1.00000+exp(0.100000*(states[0]+75.4400)))
    algebraic[51] = constants[11]*algebraic[50]*((states[0]-8.00000)-algebraic[27])
    algebraic[58] = constants[18]*states[17]*states[16]*(states[0]-algebraic[27])
    algebraic[59] = 1.00000/(1.00000+exp((5.00000-states[0])/17.0000))
    algebraic[60] = constants[19]*algebraic[59]*(states[0]-algebraic[27])
    algebraic[52] = constants[12]*(power(constants[6]/5.40000, 1.0/2))*states[5]*states[6]*(states[0]-algebraic[27])
    algebraic[36] = ((constants[0]*constants[1])/constants[2])*log((constants[6]+constants[5]*constants[7])/(states[1]+constants[5]*states[2]))
    algebraic[53] = constants[13]*(power(states[7], 2.00000))*(states[0]-algebraic[36])
    algebraic[56] = (((constants[16]*states[12]*states[13]*states[14]*states[15]*4.00000*(states[0]-15.0000)*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(0.250000*states[11]*exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-constants[8]))/(exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[45] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[8]/states[3])
    algebraic[57] = constants[17]*(states[0]-algebraic[45])
    algebraic[64] = (constants[31]*(states[0]-algebraic[27]))/(1.00000+exp((25.0000-states[0])/5.98000))
    algebraic[63] = (constants[29]*states[3])/(states[3]+constants[30])
    algebraic[48] = states[4]*constants[10]*(states[0]-algebraic[27])
    algebraic[49] = algebraic[47]+algebraic[48]
    algebraic[65] = constants[43]/(1.00000+(power(constants[41], 2.00000))/(power(states[3], 2.00000)))
    algebraic[66] = constants[42]*(states[18]-states[3])
    algebraic[67] = constants[40]*(states[11]-states[3])
    algebraic[69] = 1.00000/(1.00000+(constants[44]*constants[45])/(power(states[3]+constants[45], 2.00000)))
    algebraic[68] = constants[37]-(constants[37]-constants[38])/(1.00000+power(constants[36]/states[18], 2.00000))
    algebraic[71] = constants[33]*algebraic[68]
    algebraic[70] = constants[32]/algebraic[68]
    algebraic[72] = (algebraic[70]*(power(states[11], 2.00000))*states[19])/(constants[34]+algebraic[70]*(power(states[11], 2.00000)))
    algebraic[73] = constants[39]*algebraic[72]*(states[18]-states[11])
    algebraic[74] = 1.00000/(1.00000+(constants[46]*constants[47])/(power(states[18]+constants[47], 2.00000)))
    algebraic[75] = 1.00000/(1.00000+(constants[48]*constants[49])/(power(states[11]+constants[49], 2.00000)))
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)

