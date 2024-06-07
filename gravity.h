#ifndef GRAVITY_H
#define GRAVITY_H

typedef struct Vector {
    double x;
    double y;
} Vector;

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a < _b ? _a : _b; })

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })

static const double G = 6.67430e-11;

typedef struct Universe {
    int N;
    Vector *p;
    Vector *v;
    double *m;
} Universe;

double distance(Vector p1, Vector p2);

void vmul(Vector *m, double k, int N, Vector *out);

void vadd(Vector *m, Vector *n, int N, Vector *out);

Vector center_of_gravity(const Universe *uni);

void acc(const Universe *uni, Vector *a);

double kinetic_energy(const Universe *uni);

double gravitational_energy(const Universe *uni);

#endif /* GRAVITY_H */