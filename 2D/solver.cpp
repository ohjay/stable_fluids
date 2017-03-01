#include "solver.h"

// given a 1D index and the number of cells in each coordinate, sets XYZ to be a 2D/3D position within the array
// N[NDIM] should be of the form (height, width) in 2D, and (height, width, depth) in 3D
static void idx_to_xyz(int idx, int N[NDIM], int xyz[NDIM]) {
    if (NDIM == 2) {
        xyz[0] = (int) idx / N[1]; // first dimension in a NumPy sense
        xyz[1] = idx % N[1];
    } else if (NDIM == 3) {
        xyz[0] = (int) idx / (N[0] * N[1]);
        xyz[1] = (int) idx / N[1];
        xyz[2] = idx % N[2];
    }
}

// given a 2D/3D position within the array, returns a 1D index
// N[NDIM] should be of the form (height, width) in 2D, and (height, width, depth) in 3D
static int xyz_to_idx(int xyz[NDIM], int N[NDIM]) {
    int idx = 0;
    int multiplier = 1;
    for (int j = NDIM; --j >= 0;) {
        idx += multiplier * xyz[j];
        multiplier *= N[j];
    }
    return idx;
}

// add the force field multiplied by the time step to each value of the field
static void add_force(float** field, float* force, float dt, int num_cells) {
    for (int i = 0; i < NDIM; ++i) {
        float scaled_force = force[i] * dt;
        for (int j = 0; j < num_cells; ++j) {
            field[i][j] += scaled_force;
        }
    }
}

// returns the distance between two 2D points
static float dist2(float point0_y, float point0_x, float point1_y, float point1_x) {
    return sqrt(pow(point0_y - point1_y, 2) + pow(point0_x - point1_x, 2));
}

// linearly interpolate value of scalar field S at the location X0
// this is annoying, and we should be using Foster's staggered grid for this
static float lin_interp(float* X0, float** S, int N[NDIM], int dim) {
    float result = 0.0f;
    if (NDIM == 2) {
        int y0 = (int) X0[0];
        int x0 = (int) X0[1];
        
        // change to bilinear; I don't know what I was thinking (TODO)
        
        // interpolation weights
        float weight_tl = dist2(X0[0], X0[1], (float) y0, (float) x0);
        float weight_tr = dist2(X0[0], X0[1], (float) y0, (float) x0 + 1);
        float weight_bl = dist2(X0[0], X0[1], (float) y0 + 1, (float) x0);
        float weight_br = dist2(X0[0], X0[1], (float) y0 + 1, (float) x0 + 1);
        
        // normalize
        float sum = weight_tl + weight_tr + weight_bl + weight_br;
        weight_tl /= sum; weight_tr /= sum; weight_bl /= sum; weight_br /= sum;
        
        int xyz[NDIM]; // NDIM should be 2
        xyz[0] = y0;
        xyz[1] = x0;     result += S[dim][xyz_to_idx(xyz, N)] * weight_tl;
        xyz[1] = x0 + 1; result += S[dim][xyz_to_idx(xyz, N)] * weight_tr;
        xyz[0] = y0 + 1; result += S[dim][xyz_to_idx(xyz, N)] * weight_br;
        xyz[1] = x0;     result += S[dim][xyz_to_idx(xyz, N)] * weight_bl;
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }
    return result;
}

// trace a path starting at X through the field U over a time -dt; store result in X0
static void trace_particle(float* X, float** U, float dt, float* X0, int N[NDIM]) {
    if (NDIM == 2) {
        float f_mid[NDIM];
        f_mid[0] = X[0] - dt / 2.0f * lin_interp(X, U, N, 0);; // U[0][idx] = y-dir @ index IDX
        f_mid[1] = X[1] - dt / 2.0f * lin_interp(X, U, N, 1);;
        
        // interpolate in order to evaluate U at the midpoint
        X0[0] = X[0] - dt * lin_interp(f_mid, U, N, 0);
        X0[1] = X[1] - dt * lin_interp(f_mid, U, N, 1);
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }
    
    // (TODO) add adaptive step size? (vary dt)
}

// solve for the diffusion
static void diffuse(float** S0, float** S1, float k, float dt) {
    // (TODO) requires linear solver
}

// perform the projection
static void project(float** U1, float** U0, float dt) {
    // (TODO) requires linear solver
}

// sparse linear solver for the equation Ax = b (method: conjugate gradient)
static float* solve_lin(float** A_, float* b_, int n) {
    // (TODO) switch everything to Eigen Vectors (?)
    // (TODO) write this yourself ??
    
    Eigen::VectorXf x(n);
    Eigen::Map<Eigen::VectorXf> b(b_, n);
    Eigen::SparseMatrix<float> A(n, n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A.insert(i, j) = A_[i][j]; // switch because A is column-major?
        }
    }
    
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    x = cg.solve(b);
    
    return x.data();
}

// divide each element of S0 by (1 + dt * as) and store in S1
static void dissipate(float** S1, float** S0, float as, float dt, int num_cells) {
    for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < num_cells; ++j) {
            S1[i][j] = S0[i][j] / (1.f + dt * as);
        }
    }
}

// considerations:
// should I move these methods to the Fluid class? (probably, I think) - oj

// velocity field solver
void solver::v_step(float** U1, float** U0, float visc, float* F, float dt, int num_cells,
        int N[NDIM], float O[NDIM], float D[NDIM]) {
    add_force(U0, F, dt, num_cells);
    transport(U1, U0, U0, dt, num_cells, N, O, D);
    diffuse(U0, U1, visc, dt);
    project(U1, U0, dt);
}

// scalar field solver
void solver::s_step(float** S1, float** S0, float ks, float as, float** U, float* source, float dt,
        int num_cells, int N[NDIM], float O[NDIM], float D[NDIM]) {
    add_force(S0, source, dt, num_cells);
    transport(S1, S0, U, dt, num_cells, N, O, D);
    diffuse(S0, S1, ks, dt);
    dissipate(S1, S0, as, dt, num_cells);
}

// accounts for movement of substance due to velocity field
void solver::transport(float** S1, float** S0, float** U, float dt, int num_cells,
        int N[NDIM], float O[NDIM], float D[NDIM]) {
    for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < num_cells; ++j) {
            int xyz[NDIM];
            idx_to_xyz(j, N, xyz);
            
            // add 0.5 to each coordinate in order to get to the center of the cell
            for (int k = 0; k < NDIM; ++k) {
                xyz[k] += 0.5f;
            }
            
            float X[NDIM];
            for (int m = 0; m < NDIM; ++m) {
                X[m] = O[m] + xyz[m] * D[m];
            }
            
            float X0[NDIM];
            trace_particle(X, U, -dt, X0, N);
            S1[i][j] = lin_interp(X0, S0, N, i);
        }
    }
}
