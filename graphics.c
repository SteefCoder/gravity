#include "graphics.h"
#include "gravity.h"
#include <SDL2/SDL.h>

//const int SCREEN_WIDTH = 500;
//const int SCREEN_HEIGHT = 500;


int graphics_init(SDL_Window **window, SDL_Renderer **renderer) {
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
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

    *renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_PRESENTVSYNC | SDL_RENDERER_ACCELERATED);
    if (NULL == *renderer) {
        printf("Error creating renderer: %s\n", SDL_GetError());
        SDL_Quit();
        return -1;
    }

    return 0;
}


int render_circle(SDL_Renderer *renderer, SDL_Point centre, int radius) {
    int r2 = radius * radius;
    
    SDL_Point *points = calloc(4 * r2, sizeof(SDL_Point));

    int point_count = 0;
    for (int x = -radius+1; x < radius; ++x) {
        for (int y = -radius+1; y < radius; ++y) {
            if ((x*x + y*y) < r2) {
                points[point_count++] = (SDL_Point) { centre.x + x, centre.y + y };
            }
        }
    }

    int success = SDL_RenderDrawPoints(renderer, points, point_count);
    free(points);
    return success;
}


int scale_to_screen(const Universe *uni, double xscale, double yscale, int max_radius, SDL_Point *points, int *radii) {
    double maxm = uni->m[0];
    for (int i = 0; i < uni->N; ++i) {
        if (uni->m[i] > maxm) {
            maxm = uni->m[i];
        }
    }

    Vector center = center_of_gravity(uni);

    double k = (double)max_radius / sqrt(maxm);
    for (int i = 0; i < uni->N; ++i) {
        points[i] = (SDL_Point) {
            ((uni->p[i].x - center.x) / xscale + 1) * SCREEN_WIDTH / 2,
            ((uni->p[i].y - center.y) / yscale + 1) * SCREEN_HEIGHT / 2
        };

        radii[i] = k * sqrt(uni->m[i]);
    }
}