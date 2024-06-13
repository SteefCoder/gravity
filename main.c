#include "gravity.h"
#include "graphics.h"
#include "steppers.h"
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <time.h>


Universe create_random_universe2(int N) {
    Vector *p = calloc(N, sizeof(Vector));
    Vector *v = calloc(N, sizeof(Vector));
    double *m = calloc(N, sizeof(double));

    for (int i = 0; i < N; ++i) {
        p[i] = (Vector) { uniform(-1e+9, 1e+9), uniform(-1e+9, 1e+9) };
        v[i] = (Vector) { uniform(-3e+2, 3e+2), uniform(-3e+2, 3e+2) };
        m[i] = uniform(-1e+22, 1e+25);
    }

    return (Universe) { N, p, v, m };
}

Universe create_earth_moon(int N) {
    Vector *p = calloc(N, sizeof(Vector));
    Vector *v = calloc(N, sizeof(Vector));
    double *m = calloc(N, sizeof(double));

    p[0] = (Vector) { 0, 0 };
    v[0] = (Vector) { 0, 0 };
    m[0] = 5.9724e+24;

    p[1] = (Vector) { 0, 3.85e+8 };
    v[1] = (Vector) { 1.022e+3, 0 };
    m[1] = 0.07346e+24;

    return (Universe) { N, p, v, m };
}


int main2() {
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    if (graphics_init(&window, &renderer) < 0) {
        return -1;
    }
    
    // srand(time(NULL));

    const int N = 3;

    Universe uni = create_random_universe2(N);
    // Universe uni = create_earth_moon(N);

    SDL_Point coords[N];
    int radii[N];

    // main event handling loop
    SDL_Event e;
    bool quit = false;
    double smoothing = 0.9;
    int fps = 0;
    // double days = 0;
    double energy = kinetic_energy(&uni) + gravitational_energy(&uni);
    double h = 20.0;
    clock_t FRAME_CLOCKS = CLOCKS_PER_SEC / 120;
    while (!quit) {
        clock_t start = clock();

        // get all current events
        // it does not block
        while (SDL_PollEvent(&e)) {
            if (SDL_QUIT == e.type) {
                quit = true;
            }
        }

        // set white background
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        scale_to_screen(&uni, 2e+9, 2e+9, 30, coords, radii);

        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
        for (int i = 0; i < uni.N; ++i) {
            render_circle(renderer, coords[i], radii[i]);
        }
        
        double rts = 0;  // real-time seconds
        int i = 0;
        while ((clock() - start < FRAME_CLOCKS) && (i < 5000)) {
            rts += h;
            // step_rk4(&uni, h);
            h = step_rkn45(&uni, h);
            // h = step_rkn67(&uni, h);
            // step_euler(&uni, h);
            ++i;
        }

        // render to the screen
        SDL_RenderPresent(renderer);

        clock_t end = clock();
        double frame_time = (double)(end - start) / CLOCKS_PER_SEC;
        double error = (kinetic_energy(&uni) + gravitational_energy(&uni) - energy) / energy;
        fps = (int)(fps * smoothing + (1 - smoothing) / frame_time);
        rts = 120 * rts / 86400;
        printf("error: %E parts\t\r", error);
    }
    printf("\n");

    // destroy_universe(uni);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}