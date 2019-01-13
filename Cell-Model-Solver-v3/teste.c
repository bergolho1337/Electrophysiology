/*
   There are a total of 18 entries in the algebraic variable array.
   There are a total of 8 entries in each of the rate and state variable arrays.
   There are a total of 10 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (ms).
 * STATES[0] is V in component membrane (mV).
 * CONSTANTS[0] is C in component membrane (uF_per_mm2).
 * ALGEBRAIC[0] is i_Na in component sodium_current (uA_per_mm2).
 * ALGEBRAIC[14] is i_s in component slow_inward_current (uA_per_mm2).
 * ALGEBRAIC[15] is i_x1 in component time_dependent_outward_current (uA_per_mm2).
 * ALGEBRAIC[16] is i_K1 in component time_independent_outward_current (uA_per_mm2).
 * ALGEBRAIC[17] is Istim in component stimulus_protocol (uA_per_mm2).
 * CONSTANTS[1] is g_Na in component sodium_current (mS_per_mm2).
 * CONSTANTS[2] is E_Na in component sodium_current (mV).
 * CONSTANTS[3] is g_Nac in component sodium_current (mS_per_mm2).
 * STATES[1] is m in component sodium_current_m_gate (dimensionless).
 * STATES[2] is h in component sodium_current_h_gate (dimensionless).
 * STATES[3] is j in component sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[8] is beta_m in component sodium_current_m_gate (per_ms).
 * ALGEBRAIC[2] is alpha_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[9] is beta_h in component sodium_current_h_gate (per_ms).
 * ALGEBRAIC[3] is alpha_j in component sodium_current_j_gate (per_ms).
 * ALGEBRAIC[10] is beta_j in component sodium_current_j_gate (per_ms).
 * CONSTANTS[4] is g_s in component slow_inward_current (mS_per_mm2).
 * ALGEBRAIC[7] is E_s in component slow_inward_current (mV).
 * STATES[4] is Cai in component slow_inward_current (concentration_units).
 * STATES[5] is d in component slow_inward_current_d_gate (dimensionless).
 * STATES[6] is f in component slow_inward_current_f_gate (dimensionless).
 * ALGEBRAIC[4] is alpha_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[11] is beta_d in component slow_inward_current_d_gate (per_ms).
 * ALGEBRAIC[5] is alpha_f in component slow_inward_current_f_gate (per_ms).
 * ALGEBRAIC[12] is beta_f in component slow_inward_current_f_gate (per_ms).
 * STATES[7] is x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 * ALGEBRAIC[6] is alpha_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * ALGEBRAIC[13] is beta_x1 in component time_dependent_outward_current_x1_gate (per_ms).
 * CONSTANTS[5] is IstimStart in component stimulus_protocol (ms).
 * CONSTANTS[6] is IstimEnd in component stimulus_protocol (ms).
 * CONSTANTS[7] is IstimAmplitude in component stimulus_protocol (uA_per_mm2).
 * CONSTANTS[8] is IstimPeriod in component stimulus_protocol (ms).
 * CONSTANTS[9] is IstimPulseDuration in component stimulus_protocol (ms).
 * RATES[0] is d/dt V in component membrane (mV).
 * RATES[1] is d/dt m in component sodium_current_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_current_h_gate (dimensionless).
 * RATES[3] is d/dt j in component sodium_current_j_gate (dimensionless).
 * RATES[4] is d/dt Cai in component slow_inward_current (concentration_units).
 * RATES[5] is d/dt d in component slow_inward_current_d_gate (dimensionless).
 * RATES[6] is d/dt f in component slow_inward_current_f_gate (dimensionless).
 * RATES[7] is d/dt x1 in component time_dependent_outward_current_x1_gate (dimensionless).
 */