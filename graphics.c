#include "graphics.h"
#include "gravity.h"
#include <SDL2/SDL.h>


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

    // Enabling vsync gives a major speedup. FPS are capped at 120Hz for my laptop.
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
    
    // Maybe cache these points and shift them when needed?
    // Performance impact not yet clear.
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

// Scale the radii so that r ~ sqrt(m) and the heaviest object has a radius of max_radius (in pixels).
// Set the center of the screen at the center of gravity and scale down by xscale and yscale.
void scale_to_screen(const Universe *uni, double xscale, double yscale, int max_radius, SDL_Point *points, int *radii) {
    double maxm = uni->m[0];
    for (int i = 0; i < uni->N; ++i) {
        maxm = max(uni->m[i], maxm);
    }
    // r = k * sqrt(m)
    // k = r_max / sqrt(m_max)
    double k = (double)max_radius / sqrt(maxm);

    Vector center = center_of_gravity(uni);

    for (int i = 0; i < uni->N; ++i) {
        // Shift by the center of gravity, scale down by xscale and yscale and fit to screen.
        points[i] = (SDL_Point) {
            ((uni->p[i].x - center.x) / xscale + 1) * SCREEN_WIDTH / 2,
            ((uni->p[i].y - center.y) / yscale + 1) * SCREEN_HEIGHT / 2
        };

        radii[i] = k * sqrt(uni->m[i]);
    }
}