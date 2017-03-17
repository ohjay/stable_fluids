#include "solver.h"

using namespace solver;
using namespace std;

// given a 1D index and the number of cells in each coordinate,
// sets XYZ to be a 2D/3D position within the array
static void idx_to_xyz(int idx, int xyz[NDIM]) {
    if (NDIM == 2) {
        xyz[0] = (int) idx / CELLS_PER_SIDE; // first dimension in a NumPy sense
        xyz[1] = idx % CELLS_PER_SIDE;
    } else if (NDIM == 3) {
        xyz[0] = (int) idx / NUM_CELLS;
        xyz[1] = (int) idx / CELLS_PER_SIDE;
        xyz[2] = idx % CELLS_PER_SIDE;
    }
}

// given a 2D/3D position within the array, returns a 1D index
int solver::xyz_to_idx(int xyz[NDIM]) {
    int idx = 0;
    int multiplier = 1;
    for (int j = NDIM; --j >= 0;) {
        idx += multiplier * xyz[j];
        multiplier *= CELLS_PER_SIDE;
    }
    return idx;
}

// get flattened index from 2D grid coordinates (y, x)
// it's kind of annoying to have to use an array every time
// y is like the column; x is like the row
int solver::idx2d(int y, int x) {
    return x + CELLS_PER_SIDE * y;
}

// add the force field multiplied by the time step to each value of the field
static void add_force(float* field, float force, float dt) {
    // for (int j = 1000; j < NUM_CELLS - 1000; ++j) {
    //     field[j] += force;
    // }

    for (int j = 0; j < NUM_CELLS; ++j) {
        field[j] += force * dt;
    }
}

// add a force at some specified position in the grid
static void add_force_at(float* field, float force, float dt, int y, int x) {
    for (int i = y - 5; i < y + 5; ++i) {
        for (int j = x - 5; j < x + 5; ++j) {
            field[idx2d(i, j)] += force * dt;
        }
    }
}

/**********************************************************************/
/* (TODO) REMOVE THIS IN FAVOR OF THE COMMENTED-OUT BILERP CODE BELOW */
/**********************************************************************/

// returns the distance between two 2D points
static float dist2(float point0_y, float point0_x, float point1_y, float point1_x) {
    return sqrt(pow(point0_y - point1_y, 2) + pow(point0_x - point1_x, 2));
}

// linearly interpolate value of scalar field S at the location X0
static float lin_interp1(float* X0, float* S) {
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
        xyz[1] = x0;     i = xyz_to_idx(xyz);
        if (i >= 0 && i < NUM_CELLS) {
            result += S[i] * weight_tl;
        }
        xyz[1] = x0 + 1; i = xyz_to_idx(xyz);
        if (i >= 0 && i < NUM_CELLS) {
            result += S[i] * weight_tr;
        }
        xyz[0] = y0 + 1; i = xyz_to_idx(xyz);
        if (i >= 0 && i < NUM_CELLS) {
            result += S[i] * weight_br;
        }
        xyz[1] = x0;     i = xyz_to_idx(xyz);
        if (i >= 0 && i < NUM_CELLS) {
            result += S[i] * weight_bl;
        }
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }
    return result;
}

// linearly interpolate value of scalar field S at the location X0
// we should maybe be using Foster's staggered grid for this
static float lin_interp(float* X0, float* S) {
    if (NDIM == 2) {
        X0[0] = fmin(CELLS_PER_SIDE * 1.0f, fmax(0.0f, X0[0]));
        X0[1] = fmin(CELLS_PER_SIDE * 1.0f, fmax(0.0f, X0[1]));
        int y0 = (int) X0[0];
        int x0 = (int) X0[1];
        if (X0[0] - y0 < 0.5f) y0--;
        if (X0[1] - x0 < 0.5f) x0--;
        int y1 = fmin(CELLS_PER_SIDE - 1, y0 + 1);
        int x1 = fmin(CELLS_PER_SIDE - 1, x0 + 1);
        y0 = max(0, y0); x0 = max(0, x0);
        // X0 is the actual point
        // x0 is the integer x-coordinate to the left
        // y0 is the integer y-coordinate to the top

        float top_left = S[idx2d(y0, x0)];
        float top_right = S[idx2d(y0, x1)];
        float bottom_left = S[idx2d(y1, x0)];
        float bottom_right = S[idx2d(y1, x1)];
        float lw = fabs(x1 + 0.5f - X0[1]);
        float rw = fabs(X0[1] - (x0 + 0.5f));
        float tw = fabs(y1 + 0.5f - X0[0]);
        float bw = fabs(X0[0] - (y0 + 0.5f));
        return tw * (lw * top_left + rw * top_right) + bw * (lw * bottom_left + rw * bottom_right);


        // float top_left = S[idx2d(y0, x0)];
        // float top_right = (x0 + 1 < CELLS_PER_SIDE) ? S[idx2d(y0, x0 + 1)] : 0.0f;
        // float bottom_left = (y0 + 1 < CELLS_PER_SIDE) ? S[idx2d(y0 + 1, x0)] : 0.0f;
        // float bottom_right = (y0 + 1 < CELLS_PER_SIDE && x0 + 1 < CELLS_PER_SIDE) ?
        //         S[idx2d(y0 + 1, x0 + 1)] : 0.0f;

        // float tl_weight = (x0 + 1 < CELLS_PER_SIDE) ? x0 + 1 - X0[1] : 1.0f;
        // float bl_weight = (y0 + 1 < CELLS_PER_SIDE && x0 + 1 < CELLS_PER_SIDE) ?
        //         x0 + 1 - X0[1] : 1.0f;
        // float br_weight = (y0 + 1 < CELLS_PER_SIDE) ? X0[1] - x0 : 1.0f;
        //
        // float x_result1, x_result2; // upper and lower portions, respectively
        // x_result1 = tl_weight * top_left + (X0[1] - x0) * top_right;
        // x_result2 = bl_weight * bottom_left + br_weight * bottom_right;
        //
        // return (X0[0] - y0) * x_result2 + (y0 + 1 - X0[0]) * x_result1;
    } else if (NDIM == 3) {
        // (TODO) add trilinear interpolation code here
        return 0.0f;
    } else {
        // currently we only support 2D and 3D animations
        return 0.0f;
    }
}


// trace a path starting at X through the field U over a time -dt; store result in X0
static void trace_particle(float* X, float** U, float dt, float* X0) {
    if (NDIM == 2) {
        float f_mid[NDIM];
        f_mid[0] = X[0] - dt / 2.0f * lin_interp(X, U[0]); // U[0][idx] = y-dir @ index IDX
        f_mid[1] = X[1] - dt / 2.0f * lin_interp(X, U[1]);
        // interpolate in order to evaluate U at the midpoint
        X0[0] = X[0] - dt * lin_interp(f_mid, U[0]);
        X0[1] = X[1] - dt * lin_interp(f_mid, U[1]);
    } else if (NDIM == 3) {
        // currently we don't support this (TODO)
    }

    // (TODO) add adaptive step size? (vary dt)
}

// adds two arrays ("vectors") together
template <class Type, size_t size>
static void add(const Type(&a)[size], const Type(&b)[size], Type(&result)[size]) {
    std::transform(a, a + size, b, result, std::plus<Type>());
}

// subtracts the second vector from the first
template <class Type, size_t size>
static void subtract(const Type(&a)[size], const Type(&b)[size], Type(&result)[size]) {
    std::transform(a, a + size, b, result, std::minus<Type>());
}

// for reference:
// to compute the dot product of two vectors A and B, use
// std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0);

// set the boundaries of the array ARR (with dimensions N) to VAL
static void set_boundaries2d(float* arr, float val) {
    // r = 0, r = N[0] - 1
    for (int c = 0; c < CELLS_PER_SIDE; ++c) {
        arr[idx2d(0, c)] = val;
        arr[idx2d(CELLS_PER_SIDE - 1, c)] = val;
    }

    // c = 0, c = N[1] - 1
    for (int r = 1; r < CELLS_PER_SIDE - 1; ++r) {
        arr[idx2d(r, 0)] = val;
        arr[idx2d(r, CELLS_PER_SIDE - 1)] = val;
    }
}

// reverses direction of field at boundaries; option specifies which dimension we're handling
// 0: vertical component, 1: horizontal component, 2: scalar field
static void boundary_reverse(float* arr, int option) {
    int nc = CELLS_PER_SIDE;
    for (int i = 1; i < nc - 1; ++i) {
        // vertical boundary reverse
        arr[idx2d(0, i)] = (option == 0) ? fabs(arr[idx2d(1, i)]) : arr[idx2d(1, i)];
        arr[idx2d(nc - 1, i)] = (option == 0) ? -fabs(arr[idx2d(nc - 2, i)]) : arr[idx2d(nc - 2, i)];

        // horizontal boundary reverse
        arr[idx2d(i, 0)] = (option == 1) ? fabs(arr[idx2d(i, 1)]) : arr[idx2d(i, 1)];
        arr[idx2d(i, nc - 1)] = (option == 1) ? -fabs(arr[idx2d(i, nc - 2)]) : arr[idx2d(i, nc - 2)];
    }

    // corners
    arr[idx2d(0, 0)] = 0.5f * (arr[idx2d(1, 0)] + arr[idx2d(0, 1)]);
    arr[idx2d(0, nc - 1)] = 0.5f * (arr[idx2d(1, nc - 1)] + arr[idx2d(0, nc - 2)]);
    arr[idx2d(nc - 1, 0)] = 0.5f * (arr[idx2d(nc - 2, 0)] + arr[idx2d(nc - 1, 1)]);
    arr[idx2d(nc - 1, nc - 1)] = 0.5f * (arr[idx2d(nc - 2, nc - 1)] + arr[idx2d(nc - 1, nc - 2)]);
}

// solves discretized 2d poisson equation using conjugate gradient
// finite difference version of eqn is as follows:
//
//   K1 (S1[i + 1, j] - 2 * S1[i, j] + S1[i - 1, j])
// + K2 (S1[i, j + 1] - 2 * S1[i, j] + S1[i, j - 1])
// + S1[i, j]
// = S0[i, j]
//
// (we are solving for S1 here; this fn will save its result in the provided array)
// (CG reference: https://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html)
static void poisson2d(float k1, float k2, float* S1, float* S0, int option, int num_iter=20) {
    // we will assume that S1 is already the initial solution guess (can theoretically be whatever)

    int i, j;

    // r = b - Ax
    // compute Ax, aka A * S1, and subtract it from b (aka S0) at the same time in order to arrive at r
    float r[NUM_CELLS];
    for (i = 0; i < CELLS_PER_SIDE; ++i) { // row
        for (j = 0; j < CELLS_PER_SIDE; ++j) { // column
            float Ax_ij;
            int idx_ij = idx2d(i, j);
            Ax_ij = (1 - 2 * k1 - 2 * k2) * S1[idx_ij];
            Ax_ij += k1 * S1[idx2d((i + 1) % CELLS_PER_SIDE, j)];
            Ax_ij += k2 * S1[idx2d((i - 1 + CELLS_PER_SIDE) % CELLS_PER_SIDE, j)];
            Ax_ij += k1 * S1[idx2d(i, (j + 1) % CELLS_PER_SIDE)];
            Ax_ij += k2 * S1[idx2d(i, (j - 1 + CELLS_PER_SIDE) % CELLS_PER_SIDE)];
            r[idx_ij] = S0[idx_ij] - Ax_ij;
        }
    }

    // p = r
    float p[NUM_CELLS];
    std::copy(std::begin(r), std::end(r), std::begin(p));

    // new_r = r
    float new_r[NUM_CELLS];
    std::copy(std::begin(r), std::end(r), std::begin(new_r));
    // for some # of iterations:
    for (int _ = 0; _ < num_iter; ++_) {
        // compute v = Ap
        float v[NUM_CELLS];
        for (i = 0; i < CELLS_PER_SIDE; ++i) { // row
            for (j = 0; j < CELLS_PER_SIDE; ++j) { // column
                int idx_ij = idx2d(i, j);
                v[idx_ij] = (1 - 2 * k1 - 2 * k2) * p[idx_ij];
                v[idx_ij] += k1 * p[idx2d((i + 1) % CELLS_PER_SIDE, j)];
                v[idx_ij] += k2 * p[idx2d((i - 1) % CELLS_PER_SIDE, j)];
                v[idx_ij] += k1 * p[idx2d(i, (j + 1) % CELLS_PER_SIDE)];
                v[idx_ij] += k2 * p[idx2d(i, (j - 1) % CELLS_PER_SIDE)];
            }
        }
        // compute a = dot(r, r) / dot(p, v)
        float rTr = std::inner_product(std::begin(r), std::end(r), std::begin(r), 0.0);
        float pTv = std::inner_product(std::begin(p), std::end(p), std::begin(v), 0.0);
        float a = (pTv != 0.0f) ? rTr / pTv : 0.0f;
        // x = x + a * p
        for (i = 0; i < NUM_CELLS; ++i) {
            S1[i] = S1[i] + a * p[i];
        }
        // new_r = new_r - av (compute the updated residual)
        for (i = 0; i < NUM_CELLS; ++i) {
            new_r[i] -= a * v[i];
        }
        // g = dot(new_r, new_r) / dot(r, r)
        float g = (rTr != 0) ? std::inner_product(std::begin(new_r), std::end(new_r), std::begin(new_r), 0.0) / rTr : 0.0f;

        // p = new_r + g * p
        for (i = 0; i < NUM_CELLS; ++i) {
            p[i] = new_r[i] + g * p[i];
        }

        // r = new_r
        memcpy(r, new_r, sizeof(r));

        // set the boundaries of our current solution to 0 (Neumann condition)
        // set_boundaries2d(S1, 0);

        if (option == -1) {
            set_boundaries2d(S1, 0);
        } else if (option == 2) {
            boundary_reverse(S1, 0);
            boundary_reverse(S1, 1);
        } else {
            boundary_reverse(S1, option);
        }
    }
}

// solve for the diffusion (currently only for 2D)
static void diffuse(float* S1, float* S0, float ks, float dt, float D[NDIM]) {
    float k1 = -dt * ks / (D[0] * D[0]);
    float k2 = -dt * ks / (D[1] * D[1]);
    poisson2d(k1, k2, S1, S0, -1);
}

// perform the projection (again, assuming 2D here)
static void project(float** U1, float** U0, float dt, float D[NDIM]) {
    int i, j, idx_ij;

    float k1 = 1.0f / (D[0] * D[0]);
    float k2 = 1.0f / (D[1] * D[1]);

    // construct initial guess for the solution (x) as a bunch of 0s
    float x[NUM_CELLS] = {};

    // compute the divergence of the velocity field, to be used as b
    float divergence[NUM_CELLS];
    for (i = 0; i < CELLS_PER_SIDE; ++i) { // row
        for (j = 0; j < CELLS_PER_SIDE; ++j) { // column
            idx_ij = idx2d(i, j);
            divergence[idx_ij] = ((U0[0][idx2d((i + 1) % CELLS_PER_SIDE, j)]
                    - U0[0][idx2d((i - 1) % CELLS_PER_SIDE, j)]) / D[0]
                    + (U0[1][idx2d(i, (j + 1) % CELLS_PER_SIDE)]
                    - U0[1][idx2d(i, (j - 1) % CELLS_PER_SIDE)]) / D[1]) * 0.5f;
        }
    }

    boundary_reverse(divergence, -1);
    // set_boundaries2d(divergence, 0);
    poisson2d(k1, k2, x, divergence, 2);

    int idx_i1j, idx_i_1j, idx_ij1, idx_ij_1;

    // subtract the gradient from the previous solution
    // x is the solution here; it's called S in the paper
    for (i = 0; i < CELLS_PER_SIDE; ++i) { // row
        for (j = 0; j < CELLS_PER_SIDE; ++j) { // column
            idx_ij   = idx2d(i, j);
            idx_i1j  = idx2d((i + 1) % CELLS_PER_SIDE, j);
            idx_i_1j = idx2d((i - 1) % CELLS_PER_SIDE, j);
            idx_ij1  = idx2d(i, (j + 1) % CELLS_PER_SIDE);
            idx_ij_1 = idx2d(i, (j - 1) % CELLS_PER_SIDE);

            U1[0][idx_ij] = U0[0][idx_ij] - 0.5f * (x[idx_i1j] - x[idx_i_1j]) / D[0];
            U1[1][idx_ij] = U0[1][idx_ij] - 0.5f * (x[idx_ij1] - x[idx_ij_1]) / D[1];
        }
    }

    // set the boundaries one final time
    boundary_reverse(U1[0], 0);
    boundary_reverse(U1[1], 1);

    // set_boundaries2d(U1[0], 0);
    // set_boundaries2d(U1[1], 0);
}

// divide each element of S0 by (1 + dt * as) and store in S1
static void dissipate(float* S1, float* S0, float as, float dt) {
    for (int j = 0; j < NUM_CELLS; ++j) {
        S0[j] = S1[j] / (1.f + dt * as);
    }
}

// accounts for movement of substance due to velocity field
static void transport(float* S1, float* S0, float** U, float dt, float O[NDIM], float D[NDIM], int option) {
    for (int j = 0; j < NUM_CELLS; ++j) {
        int xyz[NDIM];
        idx_to_xyz(j, xyz);
        // add 0.5 to each coordinate in order to get to the center of the cell
        // this didn't work because it's an int array :P
        // for (int k = 0; k < NDIM; ++k) {
        //     xyz[k] += 0.5f;
        // }

        float X[NDIM];
        for (int m = 0; m < NDIM; ++m) {
            X[m] = O[m] + (0.5f + xyz[m]) * D[m];
        }

        float X0[NDIM];
        trace_particle(X, U, -dt, X0);
        S1[j] = lin_interp(X0, S0);
    }
    boundary_reverse(S1, option);
}

// considerations:
// should I move these methods to the Fluid class? (probably, I think) - oj

// velocity field solver
void solver::v_step(float** U1, float** U0, float visc, float* F, float dt,
        float O[NDIM], float D[NDIM]) {
    for (int i = 0; i < NDIM; ++i) {
        add_force(U0[i], F[i], dt);
    }
    for (int i = 0; i < NDIM; ++i) {
        transport(U1[i], U0[i], U0, dt, O, D, i);
    }
    for (int i = 0; i < NDIM; ++i) {
        diffuse(U0[i], U1[i], visc, dt, D); // notice that U0 and U1 switch
        // (TODO) resolve all of the U0 and U1 switches
        // (TODO) shouldn't it be U1, U0?
    }
    project(U1, U0, dt, D);
}

// scalar field solver
void solver::s_step(float* S1, float* S0, float ks, float as, float** U, float source, float dt,
        float O[NDIM], float D[NDIM], int Fy, int Fx) {
    add_force(S0, source, dt);
    transport(S1, S0, U, dt, O, D, -1);
    diffuse(S0, S1, ks, dt, D);
    dissipate(S1, S0, as, dt);
}
