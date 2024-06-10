#include "gravity.h"
#include "vmath.h"

#include <string.h>
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
    vmul(uni->p, k1x, h / 2, uni->N);
    vadd(uni->p, uni->p, p, uni->N);
    acc(uni, k2v);

    // k2x = v + k1v * h / 2
    vmul(k2x, k1v, h / 2, uni->N);
    vadd(k2x, k2x, uni->v, uni->N);

    // k3v = acc(uni + k2x * h / 2)
    vmul(uni->p, k2x, h / 2, uni->N);
    vadd(uni->p, uni->p, p, uni->N);
    acc(uni, k3v);

    // k3x = v + k2v * h / 2
    vmul(k3x, k2v, h / 2, uni->N);
    vadd(k3x, k3x, uni->v, uni->N);

    // k4v = acc(uni + k3x * h)
    vmul(uni->p, k3x, h, uni->N);
    vadd(uni->p, uni->p, p, uni->N);
    acc(uni, k4v);

    // k4x = v + k3v * h
    vmul(k4x, k3v, h, uni->N);
    vadd(k4x, k4x, uni->v, uni->N);
    
    for (int i = 0; i < uni->N; ++i) {
        // x = x + h / 6 * (k1x + 2*k2x + 2*k3x + k4x)
        uni->p[i].x = p[i].x + h * (k1x[i].x + 2*k2x[i].x + 2*k3x[i].x + k4x[i].x) / 6;
        uni->p[i].y = p[i].y + h * (k1x[i].y + 2*k2x[i].y + 2*k3x[i].y + k4x[i].y) / 6;

        // v = v + h / 6 * (k1v + 2*k2v + 2*k3v + k4v)
        uni->v[i].x = uni->v[i].x + h * (k1v[i].x + 2*k2v[i].x + 2*k3v[i].x + k4v[i].x) / 6;
        uni->v[i].y = uni->v[i].y + h * (k1v[i].y + 2*k2v[i].y + 2*k3v[i].y + k4v[i].y) / 6;
    }
}

double step_rkn45(Universe *uni, double h) {
    double tol = 1e-9;
    double safety = 0.5; // 50%

    Vector k0[uni->N], k1[uni->N], k2[uni->N], k3[uni->N], k4[uni->N];

    Vector p[uni->N];
    memcpy(p, uni->p, sizeof(Vector) * uni->N);

    // k0 = acc(uni)
    acc(uni, k0);

    // k1 = acc(uni.p + uni.v*h/3 + h*h*k0/18)
    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h/3 + k0[i].x * h*h/18;
        uni->p[i].y = p[i].y + uni->v[i].y * h/3 + k0[i].y * h*h/18;
    }
    acc(uni, k1);

    // k2 = acc(uni.p + uni.v * 2/3h + k1 * h*h*2/9)
    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h*2/3 + k1[i].x * 2*h*h/9;
        uni->p[i].y = p[i].y + uni->v[i].y * h*2/3 + k1[i].y * 2*h*h/9;
    }
    acc(uni, k2);

    // k3 = acc(uni.p + h*uni.v + h*h*(k0/3 + k2/6))
    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h + (k0[i].x * 2 + k2[i].x) * h*h/6;
        uni->p[i].y = p[i].y + uni->v[i].y * h + (k0[i].y * 2 + k2[i].y) * h*h/6;
    }
    acc(uni, k3);

    // uni.p = uni.p + h*uni.v + h*h*(k0*13/120 + k1*3/10 + k2*3/40 + k3/60)
    // k4 = acc(uni.p)
    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h + (k0[i].x * 13 + k1[i].x * 36 + k2[i].x * 9 + k3[i].x * 2) * h*h/120;
        uni->p[i].y = p[i].y + uni->v[i].y * h + (k0[i].y * 13 + k1[i].y * 36 + k2[i].y * 9 + k3[i].y * 2) * h*h/120;
    }

    for (int i = 0; i < uni->N; ++i) {
        uni->v[i].x += (k0[i].x + k1[i].x * 3 + k2[i].x * 3 + k3[i].x) * h/8;
        uni->v[i].y += (k0[i].y + k1[i].y * 3 + k2[i].y * 3 + k3[i].y) * h/8;
    }

    acc(uni, k4);

    double error = 0;
    for (int i = 0; i < uni->N; ++i) {
        error += distance(k3[i], k4[i]);
    }
    if (0 == error) {
        // one third at max
        return 3 * h;
    } else {
        error *= h*h/60;
        return h * pow(tol * safety / error, 0.2);
    }
}