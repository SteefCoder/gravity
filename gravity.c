#include "gravity.h"
#include "vmath.h"

#include <math.h>
#include <time.h>
#include <string.h>


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
    memset(a, 0, sizeof(Vector) * uni->N);

    for (int i = 0; i < uni->N; ++i) {
        // Set the acceleration by summing over the influence of eath other object.
        for (int j = 0; j < i; ++j) {
            double dx = uni->p[j].x - uni->p[i].x;
            double dy = uni->p[j].y - uni->p[i].y;
            double d3 = 1.0 / sqrt(dx*dx + dy*dy);
            d3 = d3 * d3 * d3;

            a[i].x += d3 * uni->m[j] * dx;
            a[i].y += d3 * uni->m[j] * dy;
            
            a[j].x -= d3 * uni->m[i] * dx;
            a[j].y -= d3 * uni->m[i] * dy;
        }
    }

    for (int i = 0; i < uni->N; ++i) {
        a[i].x *= G;
        a[i].y *= G;
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
            energy -= uni->m[i] * uni->m[j] * inv_distance(uni->p[i], uni->p[j]);
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