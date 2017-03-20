#ifndef PARAMS_H
#define PARAMS_H

// relevant systemwide parameters should go here

/* GUI parameters */
#define WINDOW_HEIGHT 600
#define WINDOW_WIDTH  600
#define WINDOW_Y      100
#define WINDOW_X      400
#define DISPLAY_KEY     0

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
#define NDIM     3 // currently assumed to be 2
#define CELLS_Z 8
#define CELLS_Y 8
#define CELLS_X 8
#define NUM_FLUIDS 1

/* Fluid parameters */
#define VISCOSITY    0.1
#define DIFFUSION    0.2
#define DISSIPATION 0.05

/* Simulation parameters */
#define NUM_ITER   4
#define DT       0.1
#define CLEANUP false

/* Computed */
#define num_cells_uz ((CELLS_Z + 1) * CELLS_Y * CELLS_X)
#define num_cells_uy (CELLS_Z * (CELLS_Y + 1) * CELLS_X)
#define num_cells_ux (CELLS_Z * CELLS_Y * (CELLS_X + 1))
#define num_cells_s  (CELLS_Z * CELLS_Y * CELLS_X)

/* Functions */
inline int idx2d(int y, int x) { return y * CELLS_X + x; }
inline int idx2dx(int y, int x) { return y * (CELLS_X + 1) + x; } // for staggered grid
inline int idx3d(int z, int y, int x) { return z * CELLS_Y * CELLS_X + y * CELLS_X + x; }
inline int idx3dx(int z, int y, int x) { return z * (CELLS_Y + 1) * (CELLS_X + 1) + y * (CELLS_X + 1) + x; }

#endif
