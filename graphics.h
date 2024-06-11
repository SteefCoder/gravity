#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "gravity.h"
#include <SDL2/SDL.h>

static const int SCREEN_WIDTH = 500;
static const int SCREEN_HEIGHT = 500;

int render_circle(SDL_Renderer *renderer, SDL_Point centre, int radius);

void scale_to_screen(const Universe *uni, double xscale, double yscale, int max_radius, SDL_Point *points, int *radii);

int graphics_init(SDL_Window **window, SDL_Renderer **renderer);

#endif /* GRAPHICS_H */