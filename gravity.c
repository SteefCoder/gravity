#include "gravity.h"
#include <math.h>
#include <time.h>

double length(Vector p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

double distance(Vector p1, Vector p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

double distance2(Vector p1, Vector p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

void vmul(Vector *m, double k, int N, Vector *out) {
    for (int i = 0; i < N; ++i) {
        out[i].x = m[i].x * k;
        out[i].y = m[i].y * k;
    }
}

void vadd(Vector *m, Vector *n, int N, Vector *out) {
    for (int i = 0; i < N; ++i) {
        out[i].x = m[i].x + n[i].x;
        out[i].y = m[i].y + n[i].y;
    }
}

Vector center_of_gravity(const Universe *uni) {
    double total_mass = 0;

    for (int i = 0; i < uni->N; ++i) {
        total_mass += uni->m[i];
    }

    Vector center = { 0, 0 };
    for (int i = 0; i < uni->N; ++i) {
        center.x += uni->p[i].x * uni->m[i] / total_mass;
        center.y += uni->p[i].y * uni->m[i] / total_mass;
    }

    return center;
}

void acc(const Universe *uni, Vector *a) {
    // Inverse cubes of Euclidian distances between objects.
    // Note that dist[i][j] == dist[j][i] and dist[i][i] == 1 (instead of 0)
    double dist3[uni->N][uni->N];
    for (int i = 1; i < uni->N; ++i) {
        for (int j = 0; j < i; ++j) {
            double d = distance(uni->p[i], uni->p[j]);
            dist3[i][j] = dist3[j][i] = 1 / (d * d * d);
            // printf("[%d][%d]\tdist3: %f", i, j, dist3[i][j]);
        }
    }

    // Fill the diagonal to avoid division by zero.
    for (int i = 0; i < uni->N; ++i) {
        dist3[i][i] = 1;
    }

    for (int i = 0; i < uni->N; ++i) {
        double sumx = 0;
        double sumy = 0;

        // Set the acceleration by summing over the influence of eath other object.
        for (int j = 0; j < uni->N; ++j) {
            sumx += dist3[i][j] * uni->m[j] * (uni->p[j].x - uni->p[i].x);
            sumy += dist3[i][j] * uni->m[j] * (uni->p[j].y - uni->p[i].y);
        }
        a[i] = (Vector) { G * sumx, G * sumy };
        // printf("[%d]\tacc: (%f, %f)\n", i, sumx, sumy);
    }
}

double kinetic_energy(const Universe *uni) {
    double energy = 0;

    for (int i = 0; i < uni->N; ++i) {
        double l = length(uni->v[i]);
        energy += l * l * uni->m[i];
    }

    return energy / 2;
}

double gravitational_energy(const Universe *uni) {
    double energy = 0;

    for (int i = 1; i < uni->N; ++i) {
        for (int j = 0; j < i; ++j) {
            energy -= uni->m[i] * uni->m[j] / distance(uni->p[i], uni->p[j]);
        }
    }

    return G * energy;
}