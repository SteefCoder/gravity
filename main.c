#include "gravity.h"
#include "graphics.h"
#include "steppers.h"
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <time.h>


int init(SDL_Window **window, SDL_Renderer **renderer) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("Error initializing SDL: %s\n", SDL_GetError());
        SDL_Quit();
        return -1;
    }

    *window = SDL_CreateWindow(
        "Test window",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        SCREEN_WIDTH,
        SCREEN_HEIGHT,
        0
    );
    if (NULL == *window) {
        printf("Error creating window: %s\n", SDL_GetError());
        SDL_Quit();
        return -1;
    }

    *renderer = SDL_CreateRenderer(*window, -1, 0);
    if (NULL == *renderer) {
        printf("Error creating renderer: %s\n", SDL_GetError());
        SDL_Quit();
        return -1;
    }

    return 0;
}


double uniform(double low, double high) {
    return ((double)rand() / RAND_MAX) * (high - low) + low;
}


Universe create_random_universe(int N) {
    Vector *p = calloc(N, sizeof(Vector));
    Vector *v = calloc(N, sizeof(Vector));
    double *m = calloc(N, sizeof(double));

    for (int i = 0; i < N; ++i) {
        p[i] = (Vector) { uniform(-1e+9, 1e+9), uniform(-1e+9, 1e+9) };
        v[i] = (Vector) { uniform(-1e+3, 1e+3), uniform(-1e+3, 1e+3) };
        m[i] = uniform(-1e+22, 1e+25);
    }

    return (Universe) { N, p, v, m };
}

void destroy_universe(Universe uni) {
    free(uni.p);
    free(uni.v);
    free(uni.m);
}


int main() {
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    if (init(&window, &renderer) < 0) {
        return -1;
    }
    
    srand(time(NULL));

    const int N = 4;

    Universe uni = create_random_universe(N);

    SDL_Point coords[N];
    int radii[N];

    // main event handling loop
    SDL_Event e;
    bool quit = false;
    double smoothing = 0.5;
    int fps = 0;
    double energy = kinetic_energy(&uni) + gravitational_energy(&uni);
    double h = 200.0;
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

        scale_to_screen(&uni, 2e+9, 2e+9, 20, coords, radii);

        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
        for (int i = 0; i < uni.N; ++i) {
            // printf("[%d]\tcoord: (%d, %d)\tradius: %d\treal coord: (%f, %f) \n", i, coords[i].x, coords[i].y, radii[i], uni.p[i].x, uni.p[i].y);
            render_circle(renderer, coords[i], radii[i]);
        }
        
        // step_euler(&uni, 200);
        step_rk4(&uni, h);

        // render to the screen
        SDL_RenderPresent(renderer);

        clock_t end = clock();
        double frame_time = (double)(end - start) / CLOCKS_PER_SEC;
        double error = (kinetic_energy(&uni) + gravitational_energy(&uni) - energy) / energy;
        fps = (int)(fps * smoothing + (1 - smoothing) / frame_time);

        printf("error: %f\r", error);
    }
    printf("\n");

    destroy_universe(uni);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}