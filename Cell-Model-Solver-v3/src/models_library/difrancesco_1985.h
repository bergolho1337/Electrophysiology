#ifndef MONOALG3D_MODEL_DIFRANCESCO_1985_H
#define MONOALG3D_MODEL_DIFRANCESCO_1985_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "model_common.h"

#define NEQ 16
#define INITIAL_V (-87.0)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(real *sv, real *rDY_, real stim_current);

#endif // MONOALG3D_MODEL_DIFRANCESCO_1985_H


