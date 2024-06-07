#include "gravity.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void step_euler(Universe *uni, double h) {
    // v_i += a_i * h
    // x_i += v_i * h

    Vector a[uni->N];
    acc(uni, a);

    for (int i = 0; i < uni->N; ++i) {
        uni->v[i].x += a[i].x * h;
        uni->p[i].x += uni->v[i].x * h;

        uni->v[i].y += a[i].y * h;
        uni->p[i].y += uni->v[i].y * h;
    }
}

void step_rk4(Universe *uni, double h) {
    Vector k1v[uni->N], k2v[uni->N], k3v[uni->N], k4v[uni->N];
    Vector k1x[uni->N], k2x[uni->N], k3x[uni->N], k4x[uni->N];

    Vector p[uni->N];
    memcpy(p, uni->p, sizeof(Vector) * uni->N);

    // k1v = acc(uni)
    acc(uni, k1v);
    // k1x = v
    memcpy(k1x, uni->v, sizeof(Vector) * uni->N);

    // k2v = acc(uni + k1x * h / 2)
    vmul(k1x, h / 2, uni->N, uni->p);
    vadd(uni->p, p, uni->N, uni->p);
    acc(uni, k2v);

    // k2x = v + k1v * h / 2
    vmul(k1v, h / 2, uni->N, k2x);
    vadd(k2x, uni->v, uni->N, k2x);

    // k3v = acc(uni + k2x * h / 2)
    vmul(k2x, h / 2, uni->N, uni->p);
    vadd(uni->p, p, uni->N, uni->p);
    acc(uni, k3v);

    // k3x = v + k2v * h / 2
    vmul(k2v, h / 2, uni->N, k3x);
    vadd(k3x, uni->v, uni->N, k3x);

    // k4v = acc(uni + k3x * h)
    vmul(k3x, h, uni->N, uni->p);
    vadd(uni->p, p, uni->N, uni->p);
    acc(uni, k4v);

    // k4x = v + k3v * h
    vmul(k3v, h, uni->N, k4x);
    vadd(k4x, uni->v, uni->N, k4x);
    
    for (int i = 0; i < uni->N; ++i) {
        // x = x + h / 6 * (k1x + 2*k2x + 2*k3x + k4x)
        uni->p[i].x = p[i].x + h * (k1x[i].x + 2*k2x[i].x + 2*k3x[i].x + k4x[i].x) / 6;
        uni->p[i].y = p[i].y + h * (k1x[i].y + 2*k2x[i].y + 2*k3x[i].y + k4x[i].y) / 6;

        // v = v + h / 6 * (k1v + 2*k2v + 2*k3v + k4v)
        uni->v[i].x = uni->v[i].x + h * (k1v[i].x + 2*k2v[i].x + 2*k3v[i].x + k4v[i].x) / 6;
        uni->v[i].y = uni->v[i].y + h * (k1v[i].y + 2*k2v[i].y + 2*k3v[i].y + k4v[i].y) / 6;
    }
}