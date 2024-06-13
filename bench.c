#include "steppers.h"

#include <time.h>
#include <stdio.h>


Universe* create_random_universe(int N) {
    Universe *uni = calloc(1, sizeof(Universe));
    Vector *p = calloc(N, sizeof(Vector));
    Vector *v = calloc(N, sizeof(Vector));
    double *m = calloc(N, sizeof(double));

    for (int i = 0; i < N; ++i) {
        p[i] = (Vector) { uniform(-1e+9, 1e+9), uniform(-1e+9, 1e+9) };
        v[i] = (Vector) { uniform(-3e+2, 3e+2), uniform(-3e+2, 3e+2) };
        m[i] = uniform(-1e+22, 1e+25);
    }

    uni->p = p;
    uni->v = v;
    uni->m = m;
    uni->N = N;
    return uni;
}

void destroy_universe(Universe *uni) {
    free(uni->p);
    free(uni->v);
    free(uni->m);
    free(uni);
}


int main() {
    // srand(time(NULL));
    
    int iters = 1000000;
    int N = 20;
    Universe *uni = create_random_universe(N);
    double start_energy = total_energy(uni);

    clock_t start = clock();
    
    double time_passed = 0;
    double h = 20.0;
    for (int i = 0; i < iters; ++i) {
        time_passed += h;
        h = step_rkn45(uni, h);
    }

    clock_t end = clock();

    double real_time = (double)(end - start) / CLOCKS_PER_SEC;
    double error = (total_energy(uni) - start_energy) / start_energy;
    time_passed /= 86400;

    printf("Simulated %f days in %f real seconds (error %E)\n", time_passed, real_time, error);
    printf("Average %f simulated seconds per second (@%f iters/sec).\n", time_passed / real_time, iters / real_time);

    destroy_universe(uni);
    return 0;
}