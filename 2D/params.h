#ifndef PARAMS_H
#define PARAMS_H

// relevant systemwide parameters should go here
#define NDIM           2                                 // number of dimensions
#define SIDE_LEN       20                                // for now, all of our grids are square (replaces L)
#define CELLS_PER_SIDE 20                                // # of cells in each coordinate (replaces N)
#define NUM_CELLS      (CELLS_PER_SIDE * CELLS_PER_SIDE) // total # of cells in each grid

#define VISC           0.1                               // viscosity parameter
#define KS             0.2                               // diffusion constant
#define AS             0.3                               // dissipation rate
#define DT             0.01                              // time step

#endif
