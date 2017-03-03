#include <math.h>
#include <eigen3/Eigen/IterativeLinearSolvers> // symlink or install into /usr/local/include/

#include "params.h"

namespace solver {
    // functions we will actually call from outside
    void v_step(float** U1, float** U0, float visc, float* F, float dt, int num_cells,
            int N[NDIM], float O[NDIM], float D[NDIM]);
    void s_step(float* S1, float* S0, float ks, float as, float** U, float source, float dt,
            int num_cells, int N[NDIM], float O[NDIM], float D[NDIM]);
    void transport(float* S1, float* S0, float** U, float dt, int num_cells,
            int N[NDIM], float O[NDIM], float D[NDIM]);
}
