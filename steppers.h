#ifndef STEPPERS_H
#define STEPPERS_H

#include "gravity.h"

typedef struct NystromButcherTableau {
    int kappa;  // The order of the tableau
    const double *alpha;
    const double *gamma; // 2D array like a staircase disguised as 1D.
    const double *cdot;
} NBT_t;

void step_euler(Universe *uni, double h);

void step_rk4(Universe *uni, double h);

double step_rkn45(Universe *uni, double h);

double step_rkn45_tableau(Universe *uni, double h);

#endif /* STEPPERS_H */