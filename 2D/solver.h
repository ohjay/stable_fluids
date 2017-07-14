#include <math.h>
// #include <eigen3/Eigen/IterativeLinearSolvers> // symlink or install into /usr/local/include/
#include <iostream>
#include <numeric>

#include "params.h"
#include "debug.h"

namespace solver {
    // functions we will actually call from outside
    int idx2d(int y, int x);

    void v_step(double** U1, double** U0, double visc, double* F, double dt, double O[NDIM], double D[NDIM]);
    void s_step(double* S1, double* S0, double ks, double as, double** U, double source, double dt,
            double O[NDIM], double D[NDIM], int Fy, int Fx);
}
