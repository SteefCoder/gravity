#ifndef VMATH_H
#define VMATH_H

#include <stdlib.h>

#define distance(p1, p2) sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y))

#define inv_distance(p1, p2) 1.0 / sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
typedef struct Vector {
    double x;
    double y;
} Vector;

double uniform(double low, double high);

double length(Vector p);

void vadd(Vector *out, const Vector *m, const Vector *n, size_t N);

void vmul(Vector *out, const Vector *m, double k, size_t N);

#endif /* VMATH_H */