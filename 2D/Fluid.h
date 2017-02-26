#include <iostream>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "solver.h"

#define NDIM 2

class Fluid {
private:
    float **U0, **U1; // velocity grids
    float **S0, **S1; // scalar grids
    float O[NDIM]; // origin
    int L[NDIM]; // length of each side
    int N[NDIM]; // number of cells in each coordinate
    float D[NDIM]; // size of each voxel
public:
    void init(float visc, float ks, float as, float dt);
    void step();

    // properties of fluid
    float visc; // viscosity

    // properties of substance
    float ks; // diffusion constant
    float as; // dissipation rate

    // speed of interactivity
    float dt; // time step
};
