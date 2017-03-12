#include <math.h>
#include <eigen3/Eigen/IterativeLinearSolvers> // symlink or install into /usr/local/include/
#include <iostream>
#include <numeric>

#include "params.h"
#include "debug.h"

namespace solver {
    // functions we will actually call from outside
    
    int xyz_to_idx(int xyz[NDIM]);
    int idx2d(int y, int x);
    
    void v_step(float** U1, float** U0, float visc, float* F, float dt, float O[NDIM], float D[NDIM]);
    void s_step(float* S1, float* S0, float ks, float as, float** U, float source, float dt,
            float O[NDIM], float D[NDIM], int Fy, int Fx);
}
