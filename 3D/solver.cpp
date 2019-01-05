#include "solver.h"

using namespace solver;
using namespace std;

static void negate_field(float* field) {
    for (int i = 0; i < num_cells; ++i) {
        field[i] = -field[i];
    }
}

static void set_boundary_values(float* field, int key) {
    switch (key) {
        case 1:
            // z-velocity
            for (int z = 1; z < CELLS_Z - 1; ++z) {
                field[idx3d(z, 0, 0)] = field[idx3d(z, 1, 1)];
                field[idx3d(z, 0, CELLS_X - 1)] = field[idx3d(z, 1, CELLS_X - 2)];
                field[idx3d(z, CELLS_Y - 1, 0)] = field[idx3d(z, CELLS_Y - 2, 1)];
                field[idx3d(z, CELLS_Y - 1, CELLS_X - 1)] = field[idx3d(z, CELLS_Y - 2, CELLS_X - 2)];
            }
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx3d(0, y, 0)] = -field[idx3d(1, y, 1)];
                field[idx3d(0, y, CELLS_X - 1)] = -field[idx3d(1, y, CELLS_X - 2)];
                field[idx3d(CELLS_Z - 1, y, 0)] = -field[idx3d(CELLS_Z - 2, y, 1)];
                field[idx3d(CELLS_Z - 1, y, CELLS_X - 1)] = -field[idx3d(CELLS_Z - 2, y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx3d(0, 0, x)] = -field[idx3d(1, 1, x)];
                field[idx3d(0, CELLS_Y - 1, x)] = -field[idx3d(1, CELLS_Y - 2, x)];
                field[idx3d(CELLS_Z - 1, 0, x)] = -field[idx3d(CELLS_Z - 2, 1, x)];
                field[idx3d(CELLS_Z - 1, CELLS_Y - 1, x)] = -field[idx3d(CELLS_Z - 2, CELLS_Y - 2, x)];
            }
            break;
        case 2:
            // y-velocity
            for (int z = 1; z < CELLS_Z - 1; ++z) {
                field[idx3d(z, 0, 0)] = -field[idx3d(z, 1, 1)];
                field[idx3d(z, 0, CELLS_X - 1)] = -field[idx3d(z, 1, CELLS_X - 2)];
                field[idx3d(z, CELLS_Y - 1, 0)] = -field[idx3d(z, CELLS_Y - 2, 1)];
                field[idx3d(z, CELLS_Y - 1, CELLS_X - 1)] = -field[idx3d(z, CELLS_Y - 2, CELLS_X - 2)];
            }
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx3d(0, y, 0)] = field[idx3d(1, y, 1)];
                field[idx3d(0, y, CELLS_X - 1)] = field[idx3d(1, y, CELLS_X - 2)];
                field[idx3d(CELLS_Z - 1, y, 0)] = field[idx3d(CELLS_Z - 2, y, 1)];
                field[idx3d(CELLS_Z - 1, y, CELLS_X - 1)] = field[idx3d(CELLS_Z - 2, y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx3d(0, 0, x)] = -field[idx3d(1, 1, x)];
                field[idx3d(0, CELLS_Y - 1, x)] = -field[idx3d(1, CELLS_Y - 2, x)];
                field[idx3d(CELLS_Z - 1, 0, x)] = -field[idx3d(CELLS_Z - 2, 1, x)];
                field[idx3d(CELLS_Z - 1, CELLS_Y - 1, x)] = -field[idx3d(CELLS_Z - 2, CELLS_Y - 2, x)];
            }
            break;
        case 3:
            // x-velocity
            for (int z = 1; z < CELLS_Z - 1; ++z) {
                field[idx3d(z, 0, 0)] = -field[idx3d(z, 1, 1)];
                field[idx3d(z, 0, CELLS_X - 1)] = -field[idx3d(z, 1, CELLS_X - 2)];
                field[idx3d(z, CELLS_Y - 1, 0)] = -field[idx3d(z, CELLS_Y - 2, 1)];
                field[idx3d(z, CELLS_Y - 1, CELLS_X - 1)] = -field[idx3d(z, CELLS_Y - 2, CELLS_X - 2)];
            }
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx3d(0, y, 0)] = -field[idx3d(1, y, 1)];
                field[idx3d(0, y, CELLS_X - 1)] = -field[idx3d(1, y, CELLS_X - 2)];
                field[idx3d(CELLS_Z - 1, y, 0)] = -field[idx3d(CELLS_Z - 2, y, 1)];
                field[idx3d(CELLS_Z - 1, y, CELLS_X - 1)] = -field[idx3d(CELLS_Z - 2, y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx3d(0, 0, x)] = field[idx3d(1, 1, x)];
                field[idx3d(0, CELLS_Y - 1, x)] = field[idx3d(1, CELLS_Y - 2, x)];
                field[idx3d(CELLS_Z - 1, 0, x)] = field[idx3d(CELLS_Z - 2, 1, x)];
                field[idx3d(CELLS_Z - 1, CELLS_Y - 1, x)] = field[idx3d(CELLS_Z - 2, CELLS_Y - 2, x)];
            }
            break;
        default:
            // scalar
            for (int z = 1; z < CELLS_Z - 1; ++z) {
                field[idx3d(z, 0, 0)] = field[idx3d(z, 1, 1)];
                field[idx3d(z, 0, CELLS_X - 1)] = field[idx3d(z, 1, CELLS_X - 2)];
                field[idx3d(z, CELLS_Y - 1, 0)] = field[idx3d(z, CELLS_Y - 2, 1)];
                field[idx3d(z, CELLS_Y - 1, CELLS_X - 1)] = field[idx3d(z, CELLS_Y - 2, CELLS_X - 2)];
            }
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx3d(0, y, 0)] = field[idx3d(1, y, 1)];
                field[idx3d(0, y, CELLS_X - 1)] = field[idx3d(1, y, CELLS_X - 2)];
                field[idx3d(CELLS_Z - 1, y, 0)] = field[idx3d(CELLS_Z - 2, y, 1)];
                field[idx3d(CELLS_Z - 1, y, CELLS_X - 1)] = field[idx3d(CELLS_Z - 2, y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx3d(0, 0, x)] = field[idx3d(1, 1, x)];
                field[idx3d(0, CELLS_Y - 1, x)] = field[idx3d(1, CELLS_Y - 2, x)];
                field[idx3d(CELLS_Z - 1, 0, x)] = field[idx3d(CELLS_Z - 2, 1, x)];
                field[idx3d(CELLS_Z - 1, CELLS_Y - 1, x)] = field[idx3d(CELLS_Z - 2, CELLS_Y - 2, x)];
            }
            break;
    }
    // corner values
    field[idx3d(0, 0, 0)] = (field[idx3d(0, 1, 0)] +
                             field[idx3d(0, 0, 1)] +
                             field[idx3d(1, 0, 0)]) / 3.0f;
    field[idx3d(0, 0, CELLS_X - 1)] = (field[idx3d(0, 1, CELLS_X - 1)] +
                                       field[idx3d(0, 0, CELLS_X - 2)] +
                                       field[idx3d(1, 0, CELLS_X - 1)]) / 3.0f;
    field[idx3d(0, CELLS_Y - 1, 0)] = (field[idx3d(0, CELLS_Y - 2, 0)] +
                                       field[idx3d(0, CELLS_Y - 1, 1)] +
                                       field[idx3d(1, CELLS_Y - 1, 0)]) / 3.0f;
    field[idx3d(0, CELLS_Y - 1, CELLS_X - 1)] = (field[idx3d(0, CELLS_Y - 2, CELLS_X - 1)] +
                                                 field[idx3d(0, CELLS_Y - 1, CELLS_X - 2)] +
                                                 field[idx3d(1, CELLS_Y - 1, CELLS_X - 1)]) / 3.0f;
    field[idx3d(CELLS_Z - 1, 0, 0)] = (field[idx3d(CELLS_Z - 1, 1, 0)] +
                                       field[idx3d(CELLS_Z - 1, 0, 1)] +
                                       field[idx3d(CELLS_Z - 2, 0, 0)]) / 3.0f;
    field[idx3d(CELLS_Z - 1, 0, CELLS_X - 1)] = (field[idx3d(CELLS_Z - 1, 1, CELLS_X - 1)] +
                                                 field[idx3d(CELLS_Z - 1, 0, CELLS_X - 2)] +
                                                 field[idx3d(CELLS_Z - 2, 0, CELLS_X - 1)]) / 3.0f;
    field[idx3d(CELLS_Z - 1, CELLS_Y - 1, 0)] = (field[idx3d(CELLS_Z - 1, CELLS_Y - 2, 0)] +
                                                 field[idx3d(CELLS_Z - 1, CELLS_Y - 1, 1)] +
                                                 field[idx3d(CELLS_Z - 2, CELLS_Y - 1, 0)]) / 3.0f;
    field[idx3d(CELLS_Z - 1, CELLS_Y - 1, CELLS_X - 1)] = (field[idx3d(CELLS_Z - 1, CELLS_Y - 2, CELLS_X - 1)] +
                                                           field[idx3d(CELLS_Z - 1, CELLS_Y - 1, CELLS_X - 2)] +
                                                           field[idx3d(CELLS_Z - 2, CELLS_Y - 1, CELLS_X - 1)]) / 3.0f;
}

static float lin_interp(float z, float y, float x, float* field) {
    int zfloor = (int) (z - 0.5f);
    int yfloor = (int) (y - 0.5f);
    int xfloor = (int) (x - 0.5f);

    float zdiff = (z - 0.5f) - (float) zfloor;
    float ydiff = (y - 0.5f) - (float) yfloor;
    float xdiff = (x - 0.5f) - (float) xfloor;

    float ftl = field[idx3d(zfloor, yfloor, xfloor)];
    float fbl = field[idx3d(zfloor, yfloor + 1, xfloor)];
    float ftr = field[idx3d(zfloor, yfloor, xfloor + 1)];
    float fbr = field[idx3d(zfloor, yfloor + 1, xfloor + 1)];
    float btl = field[idx3d(zfloor + 1, yfloor, xfloor)];
    float bbl = field[idx3d(zfloor + 1, yfloor + 1, xfloor)];
    float btr = field[idx3d(zfloor + 1, yfloor, xfloor + 1)];
    float bbr = field[idx3d(zfloor + 1, yfloor + 1, xfloor + 1)];

    float vfl = (1.0f - ydiff) * ftl + ydiff * fbl;
    float vfr = (1.0f - ydiff) * ftr + ydiff * fbr;
    float vbl = (1.0f - ydiff) * btl + ydiff * bbl;
    float vbr = (1.0f - ydiff) * btr + ydiff * bbr;

    float ff  = (1.0f - xdiff) * vfl + xdiff * vfr;
    float fb  = (1.0f - xdiff) * vbl + xdiff * vbr;

    return (1.0f - zdiff) * ff + zdiff * fb;
}

static void transport(float* S1, float* S0, float* U_z, float* U_y, float* U_x, int key) {
    for (int z = 1; z < CELLS_Z - 1; ++z) {
        for (int y = 1; y < CELLS_Y - 1; ++y) {
            for (int x = 1; x < CELLS_X - 1; ++x) {
                float z0 = ((float) z + 0.5f) - DT * U_z[idx3d(z, y, x)];
                float y0 = ((float) y + 0.5f) - DT * U_y[idx3d(z, y, x)];
                float x0 = ((float) x + 0.5f) - DT * U_x[idx3d(z, y, x)];

                z0 = fmax(1.0f, fmin(((float) CELLS_Z) - 2.0f, z0));
                y0 = fmax(1.0f, fmin(((float) CELLS_Y) - 2.0f, y0));
                x0 = fmax(1.0f, fmin(((float) CELLS_X) - 2.0f, x0));

                S1[idx3d(z, y, x)] = lin_interp(z0, y0, x0, S0);
            }
        }
    }
    set_boundary_values(S1, key);
}

static void lin_solve(float* S1, float* S0, float a, float b, int key) {
    for (int _ = 0; _ < NUM_ITER; ++_) {
        for (int z = 1; z < CELLS_Z - 1; ++z) {
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                for (int x = 1; x < CELLS_X - 1; ++x) {
                    S1[idx3d(z, y, x)] = (S0[idx3d(z, y, x)]
                            + a * (S1[idx3d(z + 1, y, x)] + S1[idx3d(z - 1, y, x)]
                                 + S1[idx3d(z, y + 1, x)] + S1[idx3d(z, y - 1, x)]
                                 + S1[idx3d(z, y, x + 1)] + S1[idx3d(z, y, x - 1)])) / b;
                }
            }
        }
        set_boundary_values(S1, key);
    }
}

static void project(float* U1_z, float* U1_y, float* U1_x, float* U0_z, float* U0_y, float* U0_x) {
    // construct initial guess for the solution
    float S[num_cells];
    memset(S, 0, sizeof(float) * num_cells);

    // compute the divergence of the velocity field
    float divergence[num_cells];
    for (int z = 1; z < CELLS_Z - 1; ++z) {
        for (int y = 1; y < CELLS_Y - 1; ++y) {
            for (int x = 1; x < CELLS_X - 1; ++x) {
                divergence[idx3d(z, y, x)] = U0_z[idx3d(z + 1, y, x)] - U0_z[idx3d(z - 1, y, x)]
                                           + U0_y[idx3d(z, y + 1, x)] - U0_y[idx3d(z, y - 1, x)]
                                           + U0_x[idx3d(z, y, x + 1)] - U0_x[idx3d(z, y, x - 1)];
            }
        }
    }
    set_boundary_values(divergence, 0);

    // solve the Poisson equation
    negate_field(divergence);
    lin_solve(S, divergence, 1.0f, 8.0f, 0);

    // subtract the gradient from the previous solution
    for (int z = 1; z < CELLS_Z - 1; ++z) {
        for (int y = 1; y < CELLS_Y - 1; ++y) {
            for (int x = 1; x < CELLS_X - 1; ++x) {
                U1_z[idx3d(z, y, x)] = U0_z[idx3d(z, y, x)] - (S[idx3d(z + 1, y, x)] - S[idx3d(z - 1, y, x)]) / 2.0f;
                U1_y[idx3d(z, y, x)] = U0_y[idx3d(z, y, x)] - (S[idx3d(z, y + 1, x)] - S[idx3d(z, y - 1, x)]) / 2.0f;
                U1_x[idx3d(z, y, x)] = U0_x[idx3d(z, y, x)] - (S[idx3d(z, y, x + 1)] - S[idx3d(z, y, x - 1)]) / 2.0f;
            }
        }
    }
    set_boundary_values(U1_z, 1);
    set_boundary_values(U1_y, 2);
    set_boundary_values(U1_x, 3);
}

static void dissipate(float* S1, float* S0) {
    for (int i = 0; i < num_cells; ++i) {
        S1[i] = S0[i] / (1.0f + DT * DISSIPATION);
    }
}

void solver::v_step(float* U1_z, float* U1_y, float* U1_x, float* U0_z, float* U0_y, float* U0_x) {
    // aftermath of adding forces
    set_boundary_values(U1_z, 1);
    set_boundary_values(U1_y, 2);
    set_boundary_values(U1_x, 3);

    // ensure incompressibility via pressure correction (1)
    project(U0_z, U0_y, U0_x, U1_z, U1_y, U1_x);

    // self-advect
    transport(U1_z, U0_z, U0_z, U0_y, U0_x, 1);
    transport(U1_y, U0_y, U0_z, U0_y, U0_x, 2);
    transport(U1_x, U0_x, U0_z, U0_y, U0_x, 3);

    // ensure incompressibility via pressure correction (2)
    project(U0_z, U0_y, U0_x, U1_z, U1_y, U1_x);
}

void solver::s_step(float* S1, float* S0, float* U_z, float* U_y, float* U_x) {
    // aftermath of adding source
    set_boundary_values(S1, 0);

    // advect according to velocity field
    transport(S0, S1, U_z, U_y, U_x, 0);

    // dissipate
    std::swap(S1, S0);
    dissipate(S0, S1);
}
