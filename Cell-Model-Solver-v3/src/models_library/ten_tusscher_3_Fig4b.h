#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_3_4B_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_3_4B_H

#include <stdio.h>
#include <stdint.h>
#include "model_common.h"

#define NEQ 12
#define INITIAL_V (-79.550919f)

#include "../utils/logfile_utils.h"

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real *extra_parameters);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, real *extra_parameters);

#endif //MONOALG3D_MODEL_TEN_TUSSCHER_3_COMMON_H
