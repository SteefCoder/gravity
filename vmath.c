#include "vmath.h"
#include <math.h>
#include <stdlib.h>

double length(Vector p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

double distance(Vector p1, Vector p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
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