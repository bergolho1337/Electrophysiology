#ifndef MONOALG3D_MODEL_BEELER_REUTER_1977_H
#define MONOALG3D_MODEL_BEELER_REUTER_1977_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 8
#define INITIAL_V (-84.624f)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real stim_current);

#endif // MONOALG3D_MODEL_BEELER_REUTER_1977_H

