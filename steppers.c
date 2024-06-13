#include "steppers.h"
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
    double tol = 1e-12;
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
    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h + (k0[i].x * 13 + k1[i].x * 36 + k2[i].x * 9 + k3[i].x * 2) * h*h/120;
        uni->p[i].y = p[i].y + uni->v[i].y * h + (k0[i].y * 13 + k1[i].y * 36 + k2[i].y * 9 + k3[i].y * 2) * h*h/120;
    }
    // k4 = acc(uni.p)
    acc(uni, k4);

    // uni.v = uni.v + (k0 + 3k1 + 3k2 + k3) * h/8
    for (int i = 0; i < uni->N; ++i) {
        uni->v[i].x += (k0[i].x + k1[i].x * 3 + k2[i].x * 3 + k3[i].x) * h/8;
        uni->v[i].y += (k0[i].y + k1[i].y * 3 + k2[i].y * 3 + k3[i].y) * h/8;
    }

    double error = 0;
    for (int i = 0; i < uni->N; ++i) {
        error += distance(k3[i], k4[i]);
    }
    error *= h*h/60;
    return h * pow(tol * safety / error, 0.2);
}

double step_rkn_tableau(Universe *uni, double h, const NBT_t *tableau) {
    double tol = 1e-9;
    double safety = 0.5; // 50%

    int tk = tableau->kappa;

    Vector f[tk + 2][uni->N];

    Vector p[uni->N];
    memcpy(p, uni->p, sizeof(Vector) * uni->N);

    for (int kappa = 0; kappa < tk + 1; ++kappa) {
        int start_ind = kappa * (kappa - 1) / 2;
        for (int i = 0; i < uni->N; ++i) {
            double Tx = 0;
            double Ty = 0;
            for (int lambda = 0; lambda < kappa; ++lambda) {
                Tx += f[lambda][i].x * tableau->gamma[start_ind + lambda];
                Ty += f[lambda][i].y * tableau->gamma[start_ind + lambda];
            }
            uni->p[i].x = p[i].x + h * tableau->alpha[kappa] * uni->v[i].x + h*h * Tx;
            uni->p[i].y = p[i].y + h * tableau->alpha[kappa] * uni->v[i].y + h*h * Ty;
        }
        acc(uni, f[kappa]);
    }

    for (int i = 0; i < uni->N; ++i) {
        for (int kappa = 0; kappa < tk + 1; ++kappa) {
            uni->v[i].x += h * f[kappa][i].x * tableau->cdot[kappa];
            uni->v[i].y += h * f[kappa][i].y * tableau->cdot[kappa];
        }
    }

    double error = 0;
    for (int i = 0; i < uni->N; ++i) {
        error += distance(f[tk][i], f[tk + 1][i]);
    }
    error *= h*h * tableau->gamma[tk*(tk + 3) / 2];
    return h * pow(tol * safety / error, 1 / (tk + 1));
}

double step_rkn45_tableau(Universe *uni, double h) {
    static const double alpha[5] = { 0, 1.0/3, 2.0/3, 1, 1 };
    static const double gamma[10] = { 1.0/18, 0, 2.0/9, 1.0/3, 0, 1.0/6, 13.0/120, 3.0/10, 3.0/40, 1.0/60 };
    static const double cdot[4] = { 1.0/8, 3.0/8, 3.0/8, 1.0/8 };
    static const NBT_t tableau = { 4, alpha, gamma, cdot };
    return step_rkn_tableau(uni, h, &tableau);
}

double step_rkn67(Universe *uni, double h) {
    double tol = 1e-9;
    double safety = 0.5; // 50%

    Vector k0[uni->N], k1[uni->N], k2[uni->N], k3[uni->N];
    Vector k4[uni->N], k5[uni->N], k6[uni->N], k7[uni->N];

    Vector p[uni->N];
    memcpy(p, uni->p, sizeof(Vector) * uni->N);

    acc(uni, k0);

    for (int i = 0; i < uni->N; ++i) {
        uni->p[i].x = p[i].x + uni->v[i].x * h/10 + k0[i].x * h*h/200;
        uni->p[i].y = p[i].y + uni->v[i].y * h/10 + k0[i].y * h*h/200;
    }
    acc(uni, k1);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x + k1[i].x * 2) / 150;
        double Ty = (k0[i].y + k1[i].y * 2) / 150;
        uni->p[i].x = p[i].x + uni->v[i].x * h/5 + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h/5 + Ty * h*h;
    }
    acc(uni, k2);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x + k2[i].x*2) * 2/75;
        double Ty = (k0[i].y + k2[i].y*2) * 2/75;
        uni->p[i].x = p[i].x + uni->v[i].x * h*2/5 + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h*2/5 + Ty * h*h;
    }
    acc(uni, k3);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x + k2[i].x*2 + k3[i].x) * 9/200;
        double Ty = (k0[i].y + k2[i].y*2 + k3[i].y) * 9/200;
        uni->p[i].x = p[i].x + uni->v[i].x * h*3/5 + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h*3/5 + Ty * h*h;
    }
    acc(uni, k4);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x*199 - k1[i].x*456 + k2[i].x*1410 - k3[i].x*357 + k4[i].x*356) / 3600;
        double Ty = (k0[i].y*199 - k1[i].y*456 + k2[i].y*1410 - k3[i].y*357 + k4[i].y*356) / 3600;
        uni->p[i].x = p[i].x + uni->v[i].x * h*4/5 + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h*4/5 + Ty * h*h;
    }
    acc(uni, k5);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (-k0[i].x*179 + k1[i].x*816 - k3[i].x*444 + k4[i].x*876 - k5[i].x*157) / 1824;
        double Ty = (-k0[i].y*179 + k1[i].y*816 - k3[i].y*444 + k4[i].y*876 - k5[i].y*157) / 1824;

        uni->p[i].x = p[i].x + uni->v[i].x * h + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h + Ty * h*h;
    }
    acc(uni, k6);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x*122 + k2[i].x*475 + k3[i].x*100 + k4[i].x*250 + k5[i].x*50 + k6[i].x*11) / 2016;
        double Ty = (k0[i].y*122 + k2[i].y*475 + k3[i].y*100 + k4[i].y*250 + k5[i].y*50 + k6[i].y*11) / 2016;

        uni->p[i].x = p[i].x + uni->v[i].x * h + Tx * h*h;
        uni->p[i].y = p[i].y + uni->v[i].y * h + Ty * h*h;
    }
    acc(uni, k7);

    for (int i = 0; i < uni->N; ++i) {
        double Tx = (k0[i].x*19 + k2[i].x*75 + k3[i].x*50 + k4[i].x*50 + k5[i].x*75 + k6[i].x*19) / 288;
        double Ty = (k0[i].y*19 + k2[i].y*75 + k3[i].y*50 + k4[i].y*50 + k5[i].y*75 + k6[i].y*19) / 288;
        uni->v[i].x += h * Tx;
        uni->v[i].y += h * Ty;
    }

    double error = 0;
    for (int i = 0; i < uni->N; ++i) {
        error += distance(k6[i], k7[i]);
    }
    error *= h*h * 11/2016;
    return h * pow(safety * tol / error, 1/7);
}