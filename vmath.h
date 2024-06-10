#ifndef VMATH_H
#define VMATH_H

#include <stdlib.h>

typedef struct Vector {
    double x;
    double y;
} Vector;

double length(Vector p);

double distance(Vector p1, Vector p2);

void vadd(Vector *out, const Vector *m, const Vector *n, size_t N);

void vmul(Vector *out, const Vector *m, double k, size_t N);

#endif /* VMATH_H */