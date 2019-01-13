#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_2006_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_2006_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 19
#define INITIAL_V (-85.23f)

#include "../utils/logfile_utils.h"

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
void solve_model_ode_cpu(real dt, real *sv, real stim_current);

#endif //MONOALG3D_MODEL_TEN_TUSSCHER_2006_H