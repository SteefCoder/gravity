#include "gravity.h"
#include "vmath.h"

#include <math.h>
#include <time.h>


// Calculate the center of gravity of a universe.
// ### c = ∑ᵢ mᵢ ⋅ pᵢ / ∑ᵢ mᵢ
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

// Calculate the accelerations of the objects in a Universe.
// ### aᵢ = ∑ⱼ mⱼ ⋅ (pⱼ − pᵢ) / d(pⱼ, pᵢ)³
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
    }
}

// Calculate the total kinetic energy of the Universe.
// ### Eₖ = ∑ᵢ mᵢ ⋅ |vᵢ|² / 2
double kinetic_energy(const Universe *uni) {
    double energy = 0;

    for (int i = 0; i < uni->N; ++i) {
        double l = length(uni->v[i]);
        energy += l * l * uni->m[i];
    }

    return energy / 2;
}

// Calculate the total gravitational energy of the Universe.
// ### Eg = − ∑ᵢⱼ mᵢ ⋅ mⱼ / d(pᵢ, pⱼ)
double gravitational_energy(const Universe *uni) {
    double energy = 0;

    for (int i = 1; i < uni->N; ++i) {
        for (int j = 0; j < i; ++j) {
            energy -= uni->m[i] * uni->m[j] / distance(uni->p[i], uni->p[j]);
        }
    }

    return G * energy;
}

// Calculate the total energy of the Universe. By conservation of energy, this must remain equal at all times.
// Deviation from the initial value can be seen as the global error.
// ### Eₜ = Eₖ + Eg = constant
double total_energy(const Universe *uni) {
    return gravitational_energy(uni) + kinetic_energy(uni);
}