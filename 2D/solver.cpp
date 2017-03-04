#include "solver.h"

using namespace solver;
using namespace std;

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
int solver::xyz_to_idx(int xyz[NDIM], int N[NDIM]) {
    int idx = 0;
    int multiplier = 1;
    for (int j = NDIM; --j >= 0;) {
        idx += multiplier * xyz[j];
        multiplier *= N[j];
    }
    return idx;
}

// add the force field multiplied by the time step to each value of the field
static void add_force(float* field, float force, float dt, int num_cells) {
    for (int j = 0; j < num_cells; ++j) {
        field[j] += force;
    }
}

// returns the distance between two 2D points
static float dist2(float point0_y, float point0_x, float point1_y, float point1_x) {
    return sqrt(pow(point0_y - point1_y, 2) + pow(point0_x - point1_x, 2));
}

// linearly interpolate value of scalar field S at the location X0
// this is annoying, and we should be using Foster's staggered grid for this
static float lin_interp(float* X0, float* S, int N[NDIM], int num_cells) {
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
        
        int i, xyz[NDIM]; // NDIM should be 2
        xyz[0] = y0;
        xyz[1] = x0;     i = xyz_to_idx(xyz, N);
        if (i >= 0 && i < num_cells) {
            result += S[i] * weight_tl;
        }
        xyz[1] = x0 + 1; i = xyz_to_idx(xyz, N);
        if (i >= 0 && i < num_cells) {
            result += S[i] * weight_tr;
        }
        xyz[0] = y0 + 1; i = xyz_to_idx(xyz, N);
        if (i >= 0 && i < num_cells) {
            result += S[i] * weight_br;
        }
        xyz[1] = x0;     i = xyz_to_idx(xyz, N);
        if (i >= 0 && i < num_cells) {
            result += S[i] * weight_bl;
        }
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }
    return result;
}

// trace a path starting at X through the field U over a time -dt; store result in X0
static void trace_particle(float* X, float** U, float dt, float* X0, int N[NDIM], int num_cells) {
    if (NDIM == 2) {
        float f_mid[NDIM];
        f_mid[0] = X[0] - dt / 2.0f * lin_interp(X, U[0], N, num_cells); // U[0][idx] = y-dir @ index IDX
        f_mid[1] = X[1] - dt / 2.0f * lin_interp(X, U[1], N, num_cells);
        
        // interpolate in order to evaluate U at the midpoint
        X0[0] = X[0] - dt * lin_interp(f_mid, U[0], N, num_cells);
        X0[1] = X[1] - dt * lin_interp(f_mid, U[1], N, num_cells);
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }
    
    // (TODO) add adaptive step size? (vary dt)
}

// sparse linear solver for the equation Ax = b (method: conjugate gradient)
static float* solve_lin(Eigen::SparseMatrix<float, Eigen::RowMajor> A, float* b_, int n) {
    // (TODO) switch everything to Eigen Vectors (?)
    // (TODO) write this yourself ?? Can start with Jacobi or Gauss-Seidel
    
    Eigen::VectorXf x(n);
    Eigen::Map<Eigen::VectorXf> b(b_, n);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    x = cg.solve(b);
    
    return x.data();
}

// solve for the diffusion (currently only for 2D)
static void diffuse(float* S1, float* S0, float ks, float dt, int num_cells,
        int N[NDIM], float D[NDIM]) {
    // construct the A matrix (n^2 * n^2)
    // in each row, I believe there should be five nonzero entries
    
    float neighbor_multiplier = -dt * ks / (D[0] * D[0]);
    float center_multiplier = 1.0f - (4 * dt * ks) / (D[0] * D[0]);
    
    Eigen::SparseMatrix<float, Eigen::RowMajor> A(num_cells, num_cells);
    A.reserve(num_cells * 5);
    
    int xyz_intermed[NDIM];
    int curr_xyz[NDIM];
    int col;
    
    // handle boundary conditions somehow (currently neglected)
    for (int row = 1; row < num_cells - 1; ++row) {
        idx_to_xyz(row, N, curr_xyz);
        
        // 1. (i - 1, j) contribution
        xyz_intermed[0] = curr_xyz[0] - 1;
        xyz_intermed[1] = curr_xyz[1];
        if (xyz_intermed[0] >= 0) {
            col = xyz_to_idx(xyz_intermed, N);
            A.insert(row, col) = neighbor_multiplier;
        }
        
        // 2. (i + 1, j) contribution
        xyz_intermed[0] = curr_xyz[0] + 1;
        if (xyz_intermed[0] <= N[0] - 2) {
            col = xyz_to_idx(xyz_intermed, N);
            A.insert(row, col) = neighbor_multiplier;
        }
        
        // 3. (i, j - 1) contribution
        xyz_intermed[0] = curr_xyz[0];
        xyz_intermed[1] = curr_xyz[1] - 1;
        if (xyz_intermed[1] >= 0) {
            col = xyz_to_idx(xyz_intermed, N);
            A.insert(row, col) = neighbor_multiplier;
        }
        
        // 4. (i, j + 1) contribution
        xyz_intermed[1] = curr_xyz[1] + 1;
        if (xyz_intermed[1] <= N[1] - 2) {
            col = xyz_to_idx(xyz_intermed, N);
            A.insert(row, col) = neighbor_multiplier;
        }
        
        // 5. (i, j) contribution
        A.insert(row, row) = center_multiplier;
    }
    
    // store the result in S1
    float* diffusion_result = solve_lin(A, S0, num_cells);
    for (int i = 0; i < num_cells; ++i) {
        S1[i] = diffusion_result[i];
    }
}

// perform the projection
static void project(float** U1, float** U0, float dt, int num_cells, int N[NDIM], float D[NDIM]) {
    // again, assuming 2D here
    // again #2: in the P matrix, there should be five nonzero entries per row
    
    Eigen::SparseMatrix<float, Eigen::RowMajor> P(num_cells, num_cells);
    P.reserve(num_cells * 5);
    
    float b[num_cells];
    
    int xyz_intermed[NDIM];
    int curr_xyz[NDIM];
    int col;
    
    // we construct the P matrix and the b vector below
    // (TODO) look at boundaries
    for (int row = 1; row < num_cells - 1; ++row) {
        idx_to_xyz(row, N, curr_xyz);
        b[row] = 0.f;
        
        // 1. (i - 1, j)
        xyz_intermed[0] = curr_xyz[0] - 1;
        xyz_intermed[1] = curr_xyz[1];
        if (xyz_intermed[0] >= 0) {
            col = xyz_to_idx(xyz_intermed, N);
            P.insert(row, col) = 1.0f;
            b[row] -= U0[1][col]; // - U_{i-1,j}
        }
        
        // 2. (i + 1, j)
        xyz_intermed[0] = curr_xyz[0] + 1;
        if (xyz_intermed[0] <= N[0] - 2) {
            col = xyz_to_idx(xyz_intermed, N);
            P.insert(row, col) = 1.0f;
            b[row] += U0[1][col]; // + U_{i+1,j}
        }
        
        // 3. (i, j - 1)
        xyz_intermed[0] = curr_xyz[0];
        xyz_intermed[1] = curr_xyz[1] - 1;
        if (xyz_intermed[1] >= 0) {
            col = xyz_to_idx(xyz_intermed, N);
            P.insert(row, col) = 1.0f;
            b[row] -= U0[0][col]; // - U_{i,j-1}
        }
        
        // 4. (i, j + 1)
        xyz_intermed[1] = curr_xyz[1] + 1;
        if (xyz_intermed[1] <= N[1] - 2) {
            col = xyz_to_idx(xyz_intermed, N);
            P.insert(row, col) = 1.0f;
            b[row] += U0[0][col]; // + U_{i,j+1}
        }
        
        // 5. (i, j)
        P.insert(row, row) = -4.0f;
        
        b[row] *= D[0]; // multiply by h
    }
    
    float* S = solve_lin(P, b, num_cells);
    int ij[2], i1j[2], i_1j[2], ij1[2], ij_1[2]; // "i1j" is (i + 1, j)
    // likewise, "i_1j" is (i - 1, j)
    
    int idx_ij, idx_i1j, idx_i_1j, idx_ij1, idx_ij_1;
    
    // subtract the gradient from the previous solution
    for (int r = 1; r < N[0] - 1; ++r) { // boundaries
        for (int c = 1; c < N[1] - 1; ++c) {
            ij[0]   = r;     ij[1]   = c;
            i1j[0]  = r + 1; i1j[1]  = c;
            i_1j[0] = r - 1; i_1j[1] = c;
            ij1[0]  = r;     ij1[1]  = c + 1;
            ij_1[0] = r;     ij_1[1] = c - 1;
            
            idx_ij   = xyz_to_idx(ij,   N);
            idx_i1j  = xyz_to_idx(i1j,  N);
            idx_i_1j = xyz_to_idx(i_1j, N);
            idx_ij1  = xyz_to_idx(ij1,  N);
            idx_ij_1 = xyz_to_idx(ij_1, N);
            
            U1[1][idx_ij] = U0[1][idx_ij] - 0.5f * (S[idx_i1j] - S[idx_i_1j]) / D[1];
            U1[0][idx_ij] = U0[0][idx_ij] - 0.5f * (S[idx_ij1] - S[idx_ij_1]) / D[0];
        }
    }
}

// divide each element of S0 by (1 + dt * as) and store in S1
static void dissipate(float* S1, float* S0, float as, float dt, int num_cells) {
    for (int j = 0; j < num_cells; ++j) {
        S1[j] = S0[j] / (1.f + dt * as);
    }
}

// accounts for movement of substance due to velocity field
static void transport(float* S1, float* S0, float** U, float dt, int num_cells,
        int N[NDIM], float O[NDIM], float D[NDIM]) {
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
        trace_particle(X, U, -dt, X0, N, num_cells);
        S1[j] = lin_interp(X0, S0, N, num_cells);
    }
}

// considerations:
// should I move these methods to the Fluid class? (probably, I think) - oj

// velocity field solver
void solver::v_step(float** U1, float** U0, float visc, float* F, float dt, int num_cells,
        int N[NDIM], float O[NDIM], float D[NDIM]) {
    for (int i = 0; i < NDIM; ++i) {
        float scaled_force = F[i] * dt;
        add_force(U0[i], scaled_force, dt, num_cells);
    }
    for (int i = 0; i < NDIM; ++i) {
        transport(U1[i], U0[i], U0, dt, num_cells, N, O, D);
    }
    for (int i = 0; i < NDIM; ++i) {
        diffuse(U1[i], U0[i], visc, dt, num_cells, N, D); // notice that U0 and U1 switch
        // (TODO) resolve all of the U0 and U1 switches
    }
    project(U1, U0, dt, num_cells, N, D);
}

// scalar field solver
void solver::s_step(float* S1, float* S0, float ks, float as, float** U, float source, float dt,
        int num_cells, int N[NDIM], float O[NDIM], float D[NDIM]) {
    add_force(S0, source, dt, num_cells);
    transport(S1, S0, U, dt, num_cells, N, O, D);
    diffuse(S1, S0, ks, dt, num_cells, N, D);
    dissipate(S1, S0, as, dt, num_cells);
}
