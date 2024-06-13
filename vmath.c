#include "vmath.h"
#include <math.h>
#include <stdlib.h>

double uniform(double low, double high) {
    return ((double)rand() / RAND_MAX) * (high - low) + low;
}

double length(Vector p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

void vmul(Vector *out, const Vector *m, double k, size_t N) {
    for (int i = 0; i < N; ++i) {
        out[i].x = m[i].x * k;
        out[i].y = m[i].y * k;
    }
}

void vadd(Vector *out, const Vector *m, const Vector *n, size_t N) {
    for (int i = 0; i < N; ++i) {
        out[i].x = m[i].x + n[i].x;
        out[i].y = m[i].y + n[i].y;
    }
}