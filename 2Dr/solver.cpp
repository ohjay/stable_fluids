#include "solver.h"

using namespace solver;
using namespace std;

const int num_cells[3] = { num_cells_s, num_cells_uy, num_cells_ux };

static void add_force(float* field, float force, int key) {
    for (int i = 1; i < num_cells[key] - 1; ++i) {
        field[i] += force;
    }
}

// interpolate the value of field S at (y, x)
static float lin_interp(float y, float x, float* S, int key) {
    int y0 = (int) y, x0 = (int) x;
    float ydiff = y - (float) y0, xdiff = x - (float) x0;
    float vl, vr, tl, tr, bl, br;
    int tl_idx, bl_idx, tr_idx, br_idx;

    switch (key) {
        case 0:
            tl_idx = idx2d(y0,     x0);
            bl_idx = idx2d(y0 + 1, x0);
            tr_idx = idx2d(y0,     x0 + 1);
            br_idx = idx2d(y0 + 1, x0 + 1);
            break;
        case 1:
            if (ydiff > 0.5f) {
                ydiff -= 0.5f;
                tl_idx = idx2d(y0 + 1, x0);
                bl_idx = idx2d(y0 + 2, x0);
                tr_idx = idx2d(y0 + 1, x0 + 1);
                br_idx = idx2d(y0 + 2, x0 + 1);
            } else {
                ydiff += 0.5f;
                tl_idx = idx2d(y0,     x0);
                bl_idx = idx2d(y0 + 1, x0);
                tr_idx = idx2d(y0,     x0 + 1);
                br_idx = idx2d(y0 + 1, x0 + 1);
            }
            break;
        case 2:
            if (xdiff > 0.5f) {
                xdiff -= 0.5f;
                tl_idx = idx2dx(y0,     x0 + 1);
                bl_idx = idx2dx(y0 + 1, x0 + 1);
                tr_idx = idx2dx(y0,     x0 + 2);
                br_idx = idx2dx(y0 + 1, x0 + 2);
            } else {
                xdiff += 0.5f;
                tl_idx = idx2dx(y0,     x0);
                bl_idx = idx2dx(y0 + 1, x0);
                tr_idx = idx2dx(y0,     x0 + 1);
                br_idx = idx2dx(y0 + 1, x0 + 1);
            }
            break;
        default:
            return 0.0f;
    }

    tl = (tl_idx >= 0 && tl_idx < num_cells[key]) ? S[tl_idx] : 0.0f;
    bl = (bl_idx >= 0 && bl_idx < num_cells[key]) ? S[bl_idx] : 0.0f;
    tr = (tr_idx >= 0 && tr_idx < num_cells[key]) ? S[tr_idx] : 0.0f;
    br = (br_idx >= 0 && br_idx < num_cells[key]) ? S[br_idx] : 0.0f;

    vl = (1.0f - ydiff) * tl + ydiff * bl;
    vr = (1.0f - ydiff) * tr + ydiff * br;

    return (1.0f - xdiff) * vl + xdiff * vr;
}

static void trace_particle(float y, float x, float* U_y, float* U_x, float* X0) {
    X0[0] = y - DT * lin_interp(y, x, U_y, 1);
    X0[1] = x - DT * lin_interp(y, x, U_x, 2);
}

static void transport(float* S1, float* S0, float* U_y, float* U_x, int key) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;

    for (int y = 1; y < cells_y - 1; ++y) {
        for (int x = 1; x < cells_x - 1; ++x) {
            float X0[2];
            trace_particle(y, x, U_y, U_x, X0);
            int idx = (key == 2) ? idx2dx(y, x) : idx2d(y, x);
            S1[idx] = lin_interp(X0[0], X0[1], S0, 0);
        }
    }
}

static void set_boundaries(float* field, float val, int key) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;

    int idx;
    for (int y = 0; y < cells_y; ++y) {
        idx = (key == 2) ? idx2dx(y, 0) : idx2d(y, 0);
        field[idx] = val;
        idx = (key == 2) ? idx2dx(y, cells_x - 1) : idx2d(y, cells_x - 1);
        field[idx] = val;
    }

    for (int x = 0; x < cells_x; ++x) {
        idx = (key == 2) ? idx2dx(0, x) : idx2d(0, x);
        field[idx] = val;
        idx = (key == 2) ? idx2dx(cells_y - 1, x) : idx2d(cells_y - 1, x);
        field[idx] = val;
    }
}

static void boundary_reverse(float* field, int key) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;

    int (*idxf)(int, int);
    idxf = (key == 2) ? &idx2dx : &idx2d;

    for (int y = 1; y < cells_y - 1; ++y) {
        // horizontal boundary reverse
        field[idxf(y, 0)] = (key == 2) ? fabs(field[idxf(y, 1)]) : field[idxf(y, 1)];
        field[idxf(y, cells_x - 1)] = (key == 2) ? -fabs(field[idxf(y, cells_x - 2)])
                                                 : field[idxf(y, cells_x - 2)];
    }

    for (int x = 1; x < cells_x - 1; ++x) {
        // vertical boundary reverse
        field[idxf(0, x)] = (key == 1) ? fabs(field[idxf(1, x)]) : field[idxf(1, x)];
        field[idxf(cells_y - 1, x)] = (key == 1) ? -fabs(field[idxf(cells_y - 2, x)])
                                                 : field[idxf(cells_y - 2, x)];
    }

    // corners
    field[idxf(0, 0)] = (field[idxf(1, 0)] + field[idxf(0, 1)]) * 0.5f;
    field[idxf(0, cells_x - 1)] = (field[idxf(1, cells_x - 1)] +
                                   field[idxf(0, cells_x - 2)]) * 0.5f;
    field[idxf(cells_y - 1, 0)] = (field[idxf(cells_y - 2, 0)] +
                                   field[idxf(cells_y - 1, 1)]) * 0.5f;
    field[idxf(cells_y - 1, cells_x - 1)] = (field[idxf(cells_y - 2, cells_x - 1)]
                                           + field[idxf(cells_y - 1, cells_x - 2)]) * 0.5f;
}

static float dot(float* vec0, float* vec1, int size) {
    float result = 0.0f;
    for (int i = 0; i < size; ++i) {
        result += vec0[i] * vec1[i];
    }
    return result;
}

static float curl(float y, float x, float* U_y, float* U_x) {
    return (U_y[idx2d(y, x + 1)] - U_y[idx2d(y, x - 1)]
            - U_x[idx2dx(y + 1, x)] + U_x[idx2dx(y - 1, x)]) * 0.5f;
}

static void confine_vorticity(float* U_y, float* U_x) {
    // compute |w|, the curl, at each position in the velocity field
    float w[num_cells_s];
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            w[idx2d(y, x)] = fabs(curl(y, x, U_y, U_x));
        }
    }

    float dw_dy, dw_dx, norm, w_yx;
    float fy_conf[num_cells_s], fx_conf[num_cells_s];

    for (int y = 2; y < CELLS_Y - 2; ++y) {
        for (int x = 2; x < CELLS_X - 2; ++x) {
            // now compute the gradient of |w|, again using central differences
            dw_dy = (w[idx2d(y + 1, x)] - w[idx2d(y - 1, x)]) * 0.5f;
            dw_dx = (w[idx2d(y, x + 1)] - w[idx2d(y, x - 1)]) * 0.5f;

            // normalize to obtain N (unit vector pointing to center of rotation)
            norm = sqrtf(dw_dy * dw_dy + dw_dx * dw_dx) + 1e-5;
            dw_dy /= norm;
            dw_dx /= norm;

            // f_conf = N x w
            w_yx = curl(y, x, U_y, U_x);
            fy_conf[idx2d(y, x)] = VORTICITY * dw_dx * w_yx;
            fx_conf[idx2d(y, x)] = VORTICITY * dw_dy * -w_yx;
        }
    }

    int idx_yx;
    for (int y = 2; y < CELLS_Y - 2; ++y) {
        for (int x = 2; x < CELLS_X - 2; ++x) {
            idx_yx = idx2d(y, x);
            U_y[idx_yx] += fy_conf[idx_yx] * DT;
            U_x[idx2dx(y, x)] += fx_conf[idx_yx] * DT;
        }
    }
}

static void lin_solve(float* S1, float* S0, float a, float b, int key) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;

    int (*idxf)(int, int);
    idxf = (key == 2) ? &idx2dx : &idx2d;

    for (int _ = 0; _ < NUM_ITER; ++_) {
        for (int y = 1; y < cells_y - 1; ++y) {
            for (int x = 1; x < cells_x - 1; ++x) {
                S1[idxf(y, x)] = (S0[idxf(y, x)]
                        + a * (S1[idxf(y + 1, x)] + S1[idxf(y - 1, x)]
                             + S1[idxf(y, x + 1)] + S1[idxf(y, x - 1)])) / b;
            }
        }
        boundary_reverse(S1, key);
    }
}

static void diffuse(float* S1, float* S0, float diff, int key) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;

    float a = DT * diff * (cells_y - 2) * (cells_x - 2);
    lin_solve(S1, S0, a, 1 + 4 * a, key);
}

static void project(float* U1_y, float* U1_x, float* U0_y, float* U0_x) {
    int y, x, idx_yx;

    // construct initial guess for the solution
    float S[num_cells_s];
    memset(S, 0, sizeof(float) * num_cells_s);

    // compute the divergence of the velocity field
    float divergence[num_cells_s];
    memset(divergence, 0, sizeof(float) * num_cells_s);
    for (y = 1; y < CELLS_Y - 1; ++y) {
        for (x = 1; x < CELLS_X - 1; ++x) {
            idx_yx = idx2d(y, x);
            divergence[idx_yx] = U0_y[idx2d(y + 1, x)] - U0_y[idx_yx]
                               + U0_x[idx2dx(y, x + 1)] - U0_x[idx2dx(y, x)];
        }
    }

    boundary_reverse(divergence, 0);
    boundary_reverse(S, 0);
    lin_solve(S, divergence, 1, 4, 0);

    // subtract the gradient from the previous solution
    for (y = 1; y < CELLS_Y - 1; ++y) {
        for (x = 1; x < CELLS_X - 1; ++x) {
            idx_yx = idx2d(y, x);
            U1_y[idx_yx] = U0_y[idx_yx] - 0.5f * (S[idx2d(y + 1, x)] - S[idx2d(y, x)]);
            idx_yx = idx2dx(y, x);
            U1_x[idx_yx] = U0_x[idx_yx] - 0.5f * (S[idx2dx(y, x + 1)] - S[idx2dx(y, x)]);
        }
    }

    boundary_reverse(U1_y, 1);
    boundary_reverse(U1_x, 1);
}

static void dissipate(float* S1, float* S0) {
    for (int i = 0; i < num_cells_s; ++i) {
        S1[i] = S0[i] / (1.0f + DT * DISSIPATION);
    }
}

void solver::v_step(float* U1_y, float* U1_x, float* U0_y, float* U0_x, float force_y, float force_x) {
    // add forces
    add_force(U0_y, force_y * DT, 1);
    add_force(U0_x, force_x * DT, 2);

    // self-advect
    transport(U1_y, U0_y, U0_y, U0_x, 1);
    transport(U1_x, U0_x, U0_y, U0_x, 2);

    // diffuse
    diffuse(U0_y, U1_y, VISCOSITY, 1);
    diffuse(U0_x, U1_x, VISCOSITY, 2);

    // add vorticity
    confine_vorticity(U0_y, U0_x);

    // ensure incompressibility via pressure correction
    project(U1_y, U1_x, U0_y, U0_x);
}

void solver::s_step(float* S1, float* S0, float* U_y, float* U_x, float source) {
    add_force(S0, source * DT, 0);
    transport(S1, S0, U_y, U_x, 0);
    diffuse(S0, S1, DIFFUSION, 0);
    dissipate(S1, S0);
}
