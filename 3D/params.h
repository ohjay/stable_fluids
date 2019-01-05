#ifndef PARAMS_H
#define PARAMS_H

// relevant systemwide parameters should go here

/* GUI parameters */
#define WINDOW_HEIGHT         600
#define WINDOW_WIDTH          600
#define WINDOW_Y              100
#define WINDOW_X              400
#define DISPLAY_KEY             0
#define ADD_AMT_INIT         0.5f
#define FORCE_SCALE         10.0f
#define ALPHA_OPTION            3
#define COLOR_SCALE            20
#define RAINBOW_HOLD_NSTEPS    20

/* Colors */
#define RED         {1.0f, 0.0f, 0.0f}
#define GREEN       {0.0f, 1.0f, 0.0f}
#define BLUE        {0.0f, 0.0f, 1.0f}
#define YELLOW      {0.5f, 0.5f, 0.0f}
#define CYAN        {0.0f, 0.5f, 0.5f}
#define MAGENTA     {0.5f, 0.0f, 0.5f}
#define WHITE       {0.33f, 0.33f, 0.33f}
#define ALL_COLORS  {RED, GREEN, BLUE, YELLOW, CYAN, MAGENTA, WHITE}

/* Grid parameters */
#define NDIM         3
#define CELLS_Z     30
#define CELLS_Y     30
#define CELLS_X     30
#define NUM_FLUIDS   7

/* Fluid parameters */
#define DISSIPATION 0.01

/* Simulation parameters */
#define NUM_ITER    5
#define DT        0.5
#define CLEANUP false

/* Computed */
#define num_cells (CELLS_Z * CELLS_Y * CELLS_X)

/* Functions */
inline int idx3d(int z, int y, int x) { return z * CELLS_Y * CELLS_X + y * CELLS_X + x; }

#endif
