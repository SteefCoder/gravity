#ifndef STEPPERS_H
#define STEPPERS_H

#include "gravity.h"

void step_euler(Universe *uni, double h);

void step_rk4(Universe *uni, double h);

#endif /* STEPPERS_H */