#ifndef MONOALG3D_MODEL_LIRUDY_2011_H
#define MONOALG3D_MODEL_LIRUDY_2011_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "model_common.h"

#define NEQ 36
#define INITIAL_V (-84.058830)

#include "../utils/logfile_utils.h"

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(real *sv, real *rDY_, real stim_current, real dt);

#endif // MONOALG3D_MODEL_LIRUDY_2011_H


