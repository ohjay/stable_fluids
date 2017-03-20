#include "solver.h"

using namespace solver;
using namespace std;

const int num_cells[4] = { num_cells_s, num_cells_uz, num_cells_uy, num_cells_ux };

static void add_force(float* field, float force, int key) {
    for (int i = 1; i < num_cells[key] - 1; ++i) {
        field[i] += force;
    }
}

// interpolate the value of field S at (z, y, x)
static float lin_interp(float z, float y, float x, float* S, int key) {
    int z0 = (int) z, y0 = (int) y, x0 = (int) x;
    float zdiff = z - (float) z0, ydiff = y - (float) y0, xdiff = x - (float) x0;
    float ff, fb, vfl, vfr, vbl, vbr, ftl, ftr, fbl, fbr, btl, btr, bbl, bbr;
    int ftl_idx, fbl_idx, ftr_idx, fbr_idx, btl_idx, bbl_idx, btr_idx, bbr_idx;

    switch (key) {
        case 0:
            ftl_idx = idx3d(z0,     y0,     x0);
            fbl_idx = idx3d(z0,     y0 + 1, x0);
            ftr_idx = idx3d(z0,     y0,     x0 + 1);
            fbr_idx = idx3d(z0,     y0 + 1, x0 + 1);
            btl_idx = idx3d(z0 + 1, y0,     x0);
            bbl_idx = idx3d(z0 + 1, y0 + 1, x0);
            btr_idx = idx3d(z0 + 1, y0,     x0 + 1);
            bbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
            break;
        case 1:
            if (zdiff > 0.5f) {
                zdiff -= 0.5f;
                ftl_idx = idx3d(z0 + 1, y0,     x0);
                fbl_idx = idx3d(z0 + 1, y0 + 1, x0);
                ftr_idx = idx3d(z0 + 1, y0,     x0 + 1);
                fbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
                btl_idx = idx3d(z0 + 2, y0,     x0);
                bbl_idx = idx3d(z0 + 2, y0 + 1, x0);
                btr_idx = idx3d(z0 + 2, y0,     x0 + 1);
                bbr_idx = idx3d(z0 + 2, y0 + 1, x0 + 1);
            } else {
                zdiff += 0.5f;
                ftl_idx = idx3d(z0,     y0,     x0);
                fbl_idx = idx3d(z0,     y0 + 1, x0);
                ftr_idx = idx3d(z0,     y0,     x0 + 1);
                fbr_idx = idx3d(z0,     y0 + 1, x0 + 1);
                btl_idx = idx3d(z0 + 1, y0,     x0);
                bbl_idx = idx3d(z0 + 1, y0 + 1, x0);
                btr_idx = idx3d(z0 + 1, y0,     x0 + 1);
                bbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
            }
            break;
        case 2:
            if (ydiff > 0.5f) {
                ydiff -= 0.5f;
                ftl_idx = idx3d(z0,     y0 + 1, x0);
                fbl_idx = idx3d(z0,     y0 + 2, x0);
                ftr_idx = idx3d(z0,     y0 + 1, x0 + 1);
                fbr_idx = idx3d(z0,     y0 + 2, x0 + 1);
                btl_idx = idx3d(z0 + 1, y0 + 1, x0);
                bbl_idx = idx3d(z0 + 1, y0 + 2, x0);
                btr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
                bbr_idx = idx3d(z0 + 1, y0 + 2, x0 + 1);
            } else {
                ydiff += 0.5f;
                ftl_idx = idx3d(z0,     y0,     x0);
                fbl_idx = idx3d(z0,     y0 + 1, x0);
                ftr_idx = idx3d(z0,     y0,     x0 + 1);
                fbr_idx = idx3d(z0,     y0 + 1, x0 + 1);
                btl_idx = idx3d(z0 + 1, y0,     x0);
                bbl_idx = idx3d(z0 + 1, y0 + 1, x0);
                btr_idx = idx3d(z0 + 1, y0,     x0 + 1);
                bbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
            }
            break;
        case 3:
            if (xdiff > 0.5f) {
                xdiff -= 0.5f;
                ftl_idx = idx3d(z0,     y0,     x0 + 1);
                fbl_idx = idx3d(z0,     y0 + 1, x0 + 1);
                ftr_idx = idx3d(z0,     y0,     x0 + 2);
                fbr_idx = idx3d(z0,     y0 + 1, x0 + 2);
                btl_idx = idx3d(z0 + 1, y0,     x0 + 1);
                bbl_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
                btr_idx = idx3d(z0 + 1, y0,     x0 + 2);
                bbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 2);
            } else {
                xdiff += 0.5f;
                ftl_idx = idx3d(z0,     y0,     x0);
                fbl_idx = idx3d(z0,     y0 + 1, x0);
                ftr_idx = idx3d(z0,     y0,     x0 + 1);
                fbr_idx = idx3d(z0,     y0 + 1, x0 + 1);
                btl_idx = idx3d(z0 + 1, y0,     x0);
                bbl_idx = idx3d(z0 + 1, y0 + 1, x0);
                btr_idx = idx3d(z0 + 1, y0,     x0 + 1);
                bbr_idx = idx3d(z0 + 1, y0 + 1, x0 + 1);
            }
            break;
        default:
            return 0.0f;
    }

    ftl = (ftl_idx >= 0 && ftl_idx < num_cells[key]) ? S[ftl_idx] : 0.0f;
    fbl = (fbl_idx >= 0 && fbl_idx < num_cells[key]) ? S[fbl_idx] : 0.0f;
    ftr = (ftr_idx >= 0 && ftr_idx < num_cells[key]) ? S[ftr_idx] : 0.0f;
    fbr = (fbr_idx >= 0 && fbr_idx < num_cells[key]) ? S[fbr_idx] : 0.0f;
    btl = (btl_idx >= 0 && btl_idx < num_cells[key]) ? S[btl_idx] : 0.0f;
    bbl = (bbl_idx >= 0 && bbl_idx < num_cells[key]) ? S[bbl_idx] : 0.0f;
    btr = (btr_idx >= 0 && btr_idx < num_cells[key]) ? S[btr_idx] : 0.0f;
    bbr = (bbr_idx >= 0 && bbr_idx < num_cells[key]) ? S[bbr_idx] : 0.0f;


    vfl = (1.0f - ydiff) * ftl + ydiff * fbl;
    vfr = (1.0f - ydiff) * ftr + ydiff * fbr;
    vbl = (1.0f - ydiff) * btl + ydiff * bbl;
    vbr = (1.0f - ydiff) * btr + ydiff * bbr;
    ff = (1.0f - xdiff) * vfl + xdiff * vfr;
    fb = (1.0f - xdiff) * vbl + xdiff * vbr;

    return (1.0f - zdiff) * ff + zdiff * fb;
}

static void trace_particle(float z, float y, float x, float* U_z, float* U_y, float* U_x, float* X0) {
    X0[0] = z - DT * lin_interp(z, y, x, U_z, 1);
    X0[1] = y - DT * lin_interp(z, y, x, U_y, 2);
    X0[2] = x - DT * lin_interp(z, y, x, U_x, 3);
}

static void transport(float* S1, float* S0, float* U_z, float* U_y, float* U_x, int key) {
    int cells_z = (key == 1) ? CELLS_Z + 1 : CELLS_Z;
    int cells_y = (key == 2) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 3) ? CELLS_X + 1 : CELLS_X;

    for (int z = 1; z < cells_z - 1; ++z) {
        for (int y = 1; y < cells_y - 1; ++y) {
            for (int x = 1; x < cells_x - 1; ++x) {
                float X0[3];
                trace_particle(z, y, x, U_z, U_y, U_x, X0);
                int idx = (key == 0) ? idx3d(z, y, x) : idx3dx(z, y, x);
                S1[idx] = lin_interp(X0[0], X0[1], X0[2], S0, 0);
            }
        }
    }
}

// static void set_boundaries(float* field, float val, int key) {
//     int cells_z = (key == 1) ? CELLS_Z + 1 : CELLS_Z;
//     int cells_y = (key == 2) ? CELLS_Y + 1 : CELLS_Y;
//     int cells_x = (key == 3) ? CELLS_X + 1 : CELLS_X;
//
//     int idx;
//     for (int z = 0; z < cells_z; ++z) {
//         idx = (key == 2) ? idx2dx(y, 0) : idx2d(y, 0);
//         field[idx] = val;
//         idx = (key == 2) ? idx2dx(y, cells_x - 1) : idx2d(y, cells_x - 1);
//         field[idx] = val;
//     }
//
//     for (int y = 0; y < cells_y; ++y) {
//         idx = (key == 2) ? idx2dx(y, 0) : idx2d(y, 0);
//         field[idx] = val;
//         idx = (key == 2) ? idx2dx(y, cells_x - 1) : idx2d(y, cells_x - 1);
//         field[idx] = val;
//     }
//
//     for (int x = 0; x < cells_x; ++x) {
//         idx = (key == 2) ? idx2dx(0, x) : idx2d(0, x);
//         field[idx] = val;
//         idx = (key == 2) ? idx2dx(cells_y - 1, x) : idx2d(cells_y - 1, x);
//         field[idx] = val;
//     }
// }

static void boundary_reverse(float* field, int key) {
    int cells_z = (key == 1) ? CELLS_Z + 1 : CELLS_Z;
    int cells_y = (key == 2) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 3) ? CELLS_X + 1 : CELLS_X;

    int (*idxf)(int, int, int);
    idxf = (key == 0) ? &idx3d : &idx3dx;

    for (int z = 1; z < cells_z - 1; ++z) {
        // horizontal boundary reverse
        field[idxf(z, 0, 0)] = (key != 1) ? fabs(field[idxf(z, 1, 1)]) : field[idxf(z, 1, 1)];
        field[idxf(z, 0, cells_x - 1)] = (key != 1) ? -fabs(field[idxf(z, 1, cells_x - 2)])
                                                 : field[idxf(z, 1, cells_x - 2)];
        field[idxf(z, cells_y - 1, 0)] = (key != 1) ? fabs(field[idxf(z, cells_y - 2, 1)])
                                                 : field[idxf(z, cells_y - 2, 1)];
        field[idxf(z, cells_y - 1, cells_x - 1)] = (key != 1) ? -fabs(field[idxf(z, cells_y - 2, cells_x - 2)])
                                                 : field[idxf(z, cells_y - 2, cells_x - 2)];
    }

    for (int y = 1; y < cells_y - 1; ++y) {
        // horizontal boundary reverse
        field[idxf(0, y, 0)] = (key != 2) ? fabs(field[idxf(1, y, 1)]) : field[idxf(1, y, 1)];
        field[idxf(0, y, cells_x - 1)] = (key != 2) ? -fabs(field[idxf(1, y, cells_x - 2)])
                                                 : field[idxf(1, y, cells_x - 2)];
        field[idxf(cells_z - 1, y, 0)] = (key != 2) ? fabs(field[idxf(cells_z - 2, y, 1)])
                                                 : field[idxf(cells_z - 2, y, 1)];
        field[idxf(cells_z - 1, y, cells_x - 1)] = (key != 2) ? -fabs(field[idxf(cells_z - 2, y, cells_x - 2)])
                                                 : field[idxf(cells_z - 2, y, cells_x - 2)];
    }

    for (int x = 1; x < cells_x - 1; ++x) {
        // vertical boundary reverse
        field[idxf(0, 0, 0)] = (key != 3) ? fabs(field[idxf(1, 1, x)]) : field[idxf(1, 1, x)];
        field[idxf(cells_z - 1, 0, cells_x - 1)] = (key != 3) ? -fabs(field[idxf(cells_z - 2, 1, x)])
                                                 : field[idxf(cells_z - 2, 1, x)];
        field[idxf(0, cells_y - 1, 0)] = (key != 3) ? fabs(field[idxf(1, cells_y - 2, x)])
                                                 : field[idxf(1, cells_y - 2, x)];
        field[idxf(cells_z - 1, cells_y - 1, cells_x - 1)] = (key != 3) ? -fabs(field[idxf(cells_z - 2, cells_y - 2, x)])
                                                 : field[idxf(cells_z - 2, cells_y - 2, x)];
    }

    // corners
    field[idxf(0, 0, 0)] = (field[idxf(0, 1, 0)] + field[idxf(0, 0, 1)] + field[idxf(1, 0, 0)]) / 3.0f;
    field[idxf(0, 0, cells_x - 1)] = (field[idxf(0, 1, cells_x - 1)] +
                                      field[idxf(0, 0, cells_x - 2)] +
                                      field[idxf(1, 0, cells_x - 1)]) / 3.0f;
    field[idxf(0, cells_y - 1, 0)] = (field[idxf(0, cells_y - 2, 0)] +
                                      field[idxf(0, cells_y - 1, 1)] +
                                      field[idxf(1, cells_y - 1, 0)]) / 3.0f;
    field[idxf(0, cells_y - 1, cells_x - 1)] = (field[idxf(0, cells_y - 2, cells_x - 1)] +
                                                field[idxf(0, cells_y - 1, cells_x - 2)] +
                                                field[idxf(1, cells_y - 1, cells_x - 1)]) / 3.0f;
    field[idxf(cells_z - 1, 0, 0)] = (field[idxf(cells_z - 1, 1, 0)] +
                                      field[idxf(cells_z - 1, 0, 1)] +
                                      field[idxf(cells_z - 2, 0, 0)]) / 3.0f;
    field[idxf(cells_z - 1, 0, cells_x - 1)] = (field[idxf(cells_z - 1, 1, cells_x - 1)] +
                                                field[idxf(cells_z - 1, 0, cells_x - 2)] +
                                                field[idxf(cells_z - 2, 0, cells_x - 1)]) / 3.0f;
    field[idxf(cells_z - 1, cells_y - 1, 0)] = (field[idxf(cells_z - 1, cells_y - 2, 0)] +
                                                field[idxf(cells_z - 1, cells_y - 1, 1)] +
                                                field[idxf(cells_z - 2, cells_y - 1, 0)]) / 3.0f;
    field[idxf(cells_z - 1, cells_y - 1, cells_x - 1)] = (field[idxf(cells_z - 1, cells_y - 2, cells_x - 1)] +
                                                          field[idxf(cells_z - 1, cells_y - 1, cells_x - 2)] +
                                                          field[idxf(cells_z - 2, cells_y - 1, cells_x - 1)]) / 3.0f;
}

static float dot(float* vec0, float* vec1, int size) {
    float result = 0.0f;
    for (int i = 0; i < size; ++i) {
        result += vec0[i] * vec1[i];
    }
    return result;
}

static void lin_solve(float* S1, float* S0, float a, float b, int key) {
    int cells_z = (key == 1) ? CELLS_Z + 1 : CELLS_Z;
    int cells_y = (key == 2) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 3) ? CELLS_X + 1 : CELLS_X;

    int (*idxf)(int, int, int);
    idxf = (key == 0) ? &idx3d : &idx3dx;

    for (int _ = 0; _ < NUM_ITER; ++_) {
        for (int z = 1; z < cells_z - 1; ++z) {
            for (int y = 1; y < cells_y - 1; ++y) {
                for (int x = 1; x < cells_x - 1; ++x) {
                    S1[idxf(z, y, x)] = (S0[idxf(z, y, x)]
                            + a * (S1[idxf(z + 1, y, x)] + S1[idxf(z - 1, y, x)]
                                 + S1[idxf(z, y + 1, x)] + S1[idxf(z, y - 1, x)]
                                 + S1[idxf(z, y, x + 1)] + S1[idxf(z, y, x - 1)])) / b;
                }
            }
        }
        boundary_reverse(S1, key);
    }
}

static void diffuse(float* S1, float* S0, float diff, int key) {
    int cells_z = (key == 1) ? CELLS_Z + 1 : CELLS_Z;
    int cells_y = (key == 2) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 3) ? CELLS_X + 1 : CELLS_X;

    float a = DT * diff * (cells_z - 2) * (cells_y - 2) * (cells_x - 2);
    lin_solve(S1, S0, a, 1 + 4 * a, key);
}

static void project(float* U1_z, float* U1_y, float* U1_x, float* U0_z, float* U0_y, float* U0_x) {
    int z, y, x, idx_zyx;

    // construct initial guess for the solution
    float S[num_cells_s];
    memset(S, 0, sizeof(float) * num_cells_s);

    // compute the divergence of the velocity field
    float divergence[num_cells_s];
    memset(divergence, 0, sizeof(float) * num_cells_s);
    for (z = 1; z < CELLS_Z - 1; ++z) {
        for (y = 1; y < CELLS_Y - 1; ++y) {
            for (x = 1; x < CELLS_X - 1; ++x) {
                idx_zyx = idx3d(z, y, x);
                divergence[idx_zyx] = U0_z[idx3d(z + 1, y, x)] - U0_z[idx_zyx]
                                   + U0_y[idx3d(z, y + 1, x)] - U0_y[idx_zyx]
                                   + U0_x[idx3dx(z, y, x + 1)] - U0_x[idx_zyx];
            }
        }
    }

    boundary_reverse(divergence, 0);
    boundary_reverse(S, 0);
    lin_solve(S, divergence, 1, 4, 0);

    // subtract the gradient from the previous solution
    for (z = 1; z < CELLS_Z; ++z) {
        for (y = 1; y < CELLS_Y - 1; ++y) {
            for (x = 1; x < CELLS_X - 1; ++x) {
                idx_zyx = idx3d(z, y, x);
                U1_y[idx_zyx] = U0_y[idx_zyx] - 0.5f * (S[idx3d(z + 1, y, x)] - S[idx_zyx]);
                U1_y[idx_zyx] = U0_y[idx_zyx] - 0.5f * (S[idx3d(z, y + 1, x)] - S[idx_zyx]);
                U1_x[idx_zyx] = U0_x[idx_zyx] - 0.5f * (S[idx3dx(z, y, x + 1)] - S[idx_zyx]);
            }
        }
    }

    boundary_reverse(U1_z, 1);
    boundary_reverse(U1_y, 1);
    boundary_reverse(U1_x, 1);
}

static void dissipate(float* S1, float* S0) {
    for (int i = 0; i < num_cells_s; ++i) {
        S1[i] = S0[i] / (1.0f + DT * DISSIPATION);
    }
}

void solver::v_step(float* U1_z, float* U1_y, float* U1_x, float* U0_z, float* U0_y, float* U0_x, float force_z, float force_y, float force_x) {
    // add forces
    add_force(U0_z, force_z * DT, 1);
    add_force(U0_y, force_y * DT, 2);
    add_force(U0_x, force_x * DT, 3);

    // self-advect
    transport(U1_z, U0_z, U0_z, U0_y, U0_x, 1);
    transport(U1_y, U0_y, U0_z, U0_y, U0_x, 2);
    transport(U1_x, U0_x, U0_z, U0_y, U0_x, 3);

    // diffuse
    diffuse(U0_z, U1_z, VISCOSITY, 1);
    diffuse(U0_y, U1_y, VISCOSITY, 2);
    diffuse(U0_x, U1_x, VISCOSITY, 3);

    // ensure incompressibility via pressure correction
    project(U1_z, U1_y, U1_x, U0_z, U0_y, U0_x);
}

void solver::s_step(float* S1, float* S0, float* U_z, float* U_y, float* U_x, float source) {
    add_force(S0, source * DT, 0);
    transport(S1, S0, U_z, U_y, U_x, 0);
    diffuse(S0, S1, DIFFUSION, 0);
    dissipate(S1, S0);
}
