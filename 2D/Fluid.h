#include <iostream>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "solver.h"

class Fluid {
private:
    double **U0, **U1; // velocity grids (positive directions: down and right)
    // U0[0][i] is the y-direction at index i, U0[1][i] is the x-direction at index i
    // if 3D, U[2][i] is the z-direction at index i
    // purportedly row-major order within the second dimension

    double *S0, *S1; // scalar density grids
    double O[NDIM]; // origin
    double D[NDIM]; // size of each voxel
public:
    void init(double visc, double ks, double as, double dt);
    void step(double F[2], double Ssource, int Fy, int Fx);

    // properties of fluid
    double visc; // viscosity

    // properties of substance
    double ks; // diffusion constant
    double as; // dissipation rate

    // speed of interactivity
    double dt; // time step

    double grid_spacing(); // get maximum voxel size in any dimension

    double S_at(int y, int x); // scalar value at 2D location in grid

    double U10_at(int y, int x);

    double U11_at(int y, int x);

    void add_S_at(int y, int x, double source); // add source amount of fluid at 2D location

    void cleanup(void);
};
