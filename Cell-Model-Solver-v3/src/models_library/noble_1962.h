#ifndef MONOALG3D_MODEL_NOBLE_1962_H
#define MONOALG3D_MODEL_NOBLE_1962_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 4
#define INITIAL_V (-75.5344986658f)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real dt, const real *sv, real *rDY_, real stim_current);

#endif //MONOALG3D_MODEL_NOBLE_1962_H
