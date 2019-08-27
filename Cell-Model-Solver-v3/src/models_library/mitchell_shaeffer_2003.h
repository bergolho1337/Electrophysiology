#ifndef MONOALG3D_MODEL_MITCHELL_SHAEFFER_2002_H
#define MONOALG3D_MODEL_MITCHELL_SHAEFFER_2002_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 2
#define INITIAL_V (0.0f)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real stim_current);

#endif // MONOALG3D_MODEL_MITCHELL_SHAEFFER_2002_H

