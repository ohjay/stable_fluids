#include "solver.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define NDIM 2

class Fluid {
private:
    float *U0, *U1; // velocity grids
    float *S0, *S1; // scalar grids
    float* O; // origin
    float* L; // length of each side
    float* N; // number of cells in each coordinate
    float* D; // size of each voxel
    
    // we should probably change some of these names
public:
    void init();
    void step();
};
