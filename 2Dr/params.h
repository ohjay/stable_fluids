#ifndef PARAMS_H
#define PARAMS_H

// relevant systemwide parameters should go here

/* GUI parameters */
#define WINDOW_HEIGHT 600
#define WINDOW_WIDTH  600
#define WINDOW_Y      100
#define WINDOW_X      400
#define DISPLAY_KEY     1

/* Grid parameters */
#define NDIM     2 // currently assumed to be 2
#define CELLS_Y 50
#define CELLS_X 50

/* Fluid parameters */
#define VISCOSITY   0.1
#define DIFFUSION   0.2
#define DISSIPATION 0.3

/* Simulation parameters */
#define NUM_ITER  30
#define DT       0.1

/* Computed */
#define num_cells_uy ((CELLS_Y + 1) * CELLS_X)
#define num_cells_ux (CELLS_Y * (CELLS_X + 1))
#define num_cells_s  (CELLS_Y * CELLS_X)

/* Functions */
#define idx2d(y, x) (y * CELLS_X + x)
#define idx2dx(y, x) (y * (CELLS_X + 1) + x) // for staggered horizontal grid

#endif
