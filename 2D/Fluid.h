#include <iostream>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "solver.h"

class Fluid {
private:
    float **U0, **U1; // velocity grids (positive directions: down and right)
    // U0[0][i] is the y-direction at index i, U0[1][i] is the x-direction at index i
    // if 3D, U[2][i] is the z-direction at index i
    // purportedly row-major order within the second dimension
    
    float *S0, *S1; // scalar density grids
    float O[NDIM]; // origin
    int L[NDIM]; // length of each side
    int N[NDIM]; // number of cells in each coordinate
    float D[NDIM]; // size of each voxel
    int num_cells; // total number of cells in each grid
public:
    void init(float visc, float ks, float as, float dt);
    void step(float F[2], float Ssource);

    // properties of fluid
    float visc; // viscosity

    // properties of substance
    float ks; // diffusion constant
    float as; // dissipation rate

    // speed of interactivity
    float dt; // time step
    
    int num_cells_x(); // get number of cells in the x-direction
    int num_cells_y(); // get number of cells in the y-direction
    
    float grid_spacing(); // get maximum voxel size in any dimension
    
    float S_at(int y, int x); // scalar value at 2D location in grid
};
