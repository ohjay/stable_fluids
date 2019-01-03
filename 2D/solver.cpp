#include "solver.h"

using namespace solver;
using namespace std;

// -------------
// STABLE FLUIDS
// -------------

// negate all values in a field
static void negate_field(float* field) {
    for (int i = 0; i < CELLS_Y * CELLS_X; ++i) {
        field[i] = -field[i];
    }
}

// assign boundary values which aren't handled otherwise
static void set_boundary_values(float* field, int key) {
    switch (key) {
        case 1:
            // vertical velocity
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx2d(y, 0)] = field[idx2d(y, 1)];
                field[idx2d(y, CELLS_X - 1)] = field[idx2d(y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx2d(0, x)] = -field[idx2d(1, x)];
                field[idx2d(CELLS_Y - 1, x)] = -field[idx2d(CELLS_Y - 2, x)];
            }
            break;
        case 2:
            // horizontal velocity
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx2d(y, 0)] = -field[idx2d(y, 1)];
                field[idx2d(y, CELLS_X - 1)] = -field[idx2d(y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx2d(0, x)] = field[idx2d(1, x)];
                field[idx2d(CELLS_Y - 1, x)] = field[idx2d(CELLS_Y - 2, x)];
            }
            break;
        default:
            // scalar
            for (int y = 1; y < CELLS_Y - 1; ++y) {
                field[idx2d(y, 0)] = field[idx2d(y, 1)];
                field[idx2d(y, CELLS_X - 1)] = field[idx2d(y, CELLS_X - 2)];
            }
            for (int x = 1; x < CELLS_X - 1; ++x) {
                field[idx2d(0, x)] = field[idx2d(1, x)];
                field[idx2d(CELLS_Y - 1, x)] = field[idx2d(CELLS_Y - 2, x)];
            }
            break;
    }
    // corner values
    field[idx2d(0, 0)] = (field[idx2d(0, 1)] + field[idx2d(1, 0)]) / 2.0f;
    field[idx2d(0, CELLS_X - 1)] = (field[idx2d(0, CELLS_X - 2)] + field[idx2d(1, CELLS_X - 1)]) / 2.0f;
    field[idx2d(CELLS_Y - 1, 0)] = (field[idx2d(CELLS_Y - 1, 1)] + field[idx2d(CELLS_Y - 2, 0)]) / 2.0f;
    field[idx2d(CELLS_Y - 1, CELLS_X - 1)] = (field[idx2d(CELLS_Y - 1, CELLS_X - 2)] + field[idx2d(CELLS_Y - 2, CELLS_X - 1)]) / 2.0f;
}

// add values to a field
static void add_force(float* field, float* force, int key) {
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            field[idx2d(y, x)] += force[idx2d(y, x)];
        }
    }
    set_boundary_values(field, key);
}

// interpolate the value of the field at (y, x)
static float lin_interp(float y, float x, float* field) {
    int yfloor = (int) (y - 0.5f);
    int xfloor = (int) (x - 0.5f);

    float ydiff = (y - 0.5f) - (float) yfloor;
    float xdiff = (x - 0.5f) - (float) xfloor;

    float tl = field[idx2d(yfloor, xfloor)];
    float bl = field[idx2d(yfloor + 1, xfloor)];
    float tr = field[idx2d(yfloor, xfloor + 1)];
    float br = field[idx2d(yfloor + 1, xfloor + 1)];

    float vl = (1.0f - ydiff) * tl + ydiff * bl;
    float vr = (1.0f - ydiff) * tr + ydiff * br;

    return (1.0f - xdiff) * vl + xdiff * vr;
}

// advect by tracing backward in time
static void transport(float* S1, float* S0, float* U_y, float* U_x, int key) {
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            // trace particle
            float y0 = ((float) y + 0.5f) - DT * U_y[idx2d(y, x)];
            float x0 = ((float) x + 0.5f) - DT * U_x[idx2d(y, x)];

            y0 = fmax(1.0f, fmin(((float) CELLS_Y) - 2.0f, y0));
            x0 = fmax(1.0f, fmin(((float) CELLS_X) - 2.0f, x0));

            S1[idx2d(y, x)] = lin_interp(y0, x0, S0);
        }
    }
    set_boundary_values(S1, key);
}

static void lin_solve(float* S1, float* S0, float a, float b, int key) {
    for (int _ = 0; _ < NUM_ITER; ++_) {
        for (int y = 1; y < CELLS_Y - 1; ++y) {
            for (int x = 1; x < CELLS_X - 1; ++x) {
                S1[idx2d(y, x)] = (S0[idx2d(y, x)]
                        + a * (S1[idx2d(y + 1, x)] + S1[idx2d(y - 1, x)]
                             + S1[idx2d(y, x + 1)] + S1[idx2d(y, x - 1)])) / b;
            }
        }
        set_boundary_values(S1, key);
    }
}

static void diffuse(float* S1, float* S0, float diff, int key) {
    float a = DT * diff * CELLS_Y * CELLS_X;
    lin_solve(S1, S0, a, 1.0f + 4.0f * a, key);
}

static void project(float* U1_y, float* U1_x, float* U0_y, float* U0_x) {
    // construct initial guess for the solution
    float S[num_cells_s];
    memset(S, 0, sizeof(float) * num_cells_s);

    // compute the divergence of the velocity field
    float divergence[num_cells_s];
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            divergence[idx2d(y, x)] = U0_y[idx2d(y + 1, x)] - U0_y[idx2d(y - 1, x)]
                                    + U0_x[idx2d(y, x + 1)] - U0_x[idx2d(y, x - 1)];
        }
    }
    set_boundary_values(divergence, 0);

    // solve the Poisson equation
    negate_field(divergence);
    lin_solve(S, divergence, 1.0f, 4.0f, 0);

    // subtract the gradient from the previous solution
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            U1_y[idx2d(y, x)] = U0_y[idx2d(y, x)] - (S[idx2d(y + 1, x)] - S[idx2d(y - 1, x)]) / 2.0f;
            U1_x[idx2d(y, x)] = U0_x[idx2d(y, x)] - (S[idx2d(y, x + 1)] - S[idx2d(y, x - 1)]) / 2.0f;
        }
    }
    set_boundary_values(U1_y, 1);
    set_boundary_values(U1_x, 2);
}

static void dissipate(float* S1, float* S0) {
    for (int i = 0; i < num_cells_s; ++i) {
        S1[i] = S0[i] / (1.0f + DT * DISSIPATION);
    }
}

static float curl(int y, int x, float* U_y, float* U_x) {
    return (U_y[idx2d(y, x + 1)] - U_y[idx2d(y, x - 1)]
            - U_x[idx2d(y + 1, x)] + U_x[idx2d(y - 1, x)]) / 2.0f;
}

static void confine_vorticity(float* U_y, float* U_x) {
    // compute |w|, the curl, at each position in the velocity field
    float w[num_cells_s];
    float abs_w[num_cells_s];
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            w[idx2d(y, x)] = curl(y, x, U_y, U_x);
            abs_w[idx2d(y, x)] = fabs(w[idx2d(y, x)]);
        }
    }
    set_boundary_values(w, 0);
    set_boundary_values(abs_w, 0);

    float dw_dy, dw_dx, norm, w_yx;
    float fy_conf[num_cells_s], fx_conf[num_cells_s];

    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            // now compute the gradient of |w|, again using central differences
            dw_dy = (abs_w[idx2d(y + 1, x)] - abs_w[idx2d(y - 1, x)]) / 2.0f;
            dw_dx = (abs_w[idx2d(y, x + 1)] - abs_w[idx2d(y, x - 1)]) / 2.0f;

            // normalize to obtain N (unit vector pointing to center of rotation)
            norm = sqrtf(dw_dy * dw_dy + dw_dx * dw_dx);
            dw_dy /= norm + 1e-5;
            dw_dx /= norm + 1e-5;

            // f_conf = N x w
            fy_conf[idx2d(y, x)] = VORTICITY * dw_dx * w[idx2d(y, x)];
            fx_conf[idx2d(y, x)] = VORTICITY * dw_dy * -w[idx2d(y, x)];
        }
    }
    set_boundary_values(fy_conf, 0);
    set_boundary_values(fx_conf, 0);

    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            U_y[idx2d(y, x)] += fy_conf[idx2d(y, x)] * DT;
            U_x[idx2d(y, x)] += fx_conf[idx2d(y, x)] * DT;
        }
    }
    set_boundary_values(U_y, 1);
    set_boundary_values(U_x, 2);
}

void solver::v_step(float* U1_y, float* U1_x, float* U0_y, float* U0_x) {
    // aftermath of adding forces
    set_boundary_values(U1_y, 1);
    set_boundary_values(U1_x, 2);

    // add vorticity
    confine_vorticity(U1_y, U1_x);

    // diffuse
    if (VISCOSITY > 0.0f) {
        diffuse(U0_y, U1_y, VISCOSITY, 1);
        diffuse(U0_x, U1_x, VISCOSITY, 2);
        std::swap(U1_y, U0_y);
        std::swap(U1_x, U0_x);
    }

    // ensure incompressibility via pressure correction (1)
    project(U0_y, U0_x, U1_y, U1_x);

    // self-advect
    transport(U1_y, U0_y, U0_y, U0_x, 1);
    transport(U1_x, U0_x, U0_y, U0_x, 2);

    // ensure incompressibility via pressure correction (2)
    project(U0_y, U0_x, U1_y, U1_x);
}

void solver::s_step(float* S1, float* S0, float* U_y, float* U_x) {
    // advect according to velocity field
    transport(S0, S1, U_y, U_x, 0);

    // diffuse
    if (DIFFUSION > 0.0f) {
        diffuse(S1, S0, DIFFUSION, 0);
        std::swap(S0, S1);
    }

    // dissipate
    dissipate(S1, S0);
}

// ------------------------
// TARGET-DRIVEN SIMULATION
// ------------------------

static void drive_force(float* U1_y, float* U1_x, float* U0_y, float* U0_x, float* target_p, float* p) {
    float hf_p, hf_p_star, dp_star;
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            // add the vertical driving force
            hf_p = (p[idx2d(y, x)] + p[idx2d(y + 1, x)]) / 2.0f;
            hf_p_star = (target_p[idx2d(y, x)] + target_p[idx2d(y + 1, x)]) / 2.0f;
            dp_star = target_p[idx2d(y + 1, x)] - target_p[idx2d(y, x)];
            if (hf_p_star == 0.0f) {
                U1_y[idx2d(y, x)] = U0_y[idx2d(y, x)];
            } else {
                U1_y[idx2d(y, x)] = U0_y[idx2d(y, x)] + DT * DRIVING_FORCE * (hf_p * dp_star / hf_p_star);
            }

            // add the horizontal driving force
            hf_p = (p[idx2d(y, x)] + p[idx2d(y, x + 1)]) / 2.0f;
            hf_p_star = (target_p[idx2d(y, x)] + target_p[idx2d(y, x + 1)]) / 2.0f;
            dp_star = target_p[idx2d(y, x + 1)] - target_p[idx2d(y, x)];
            if (hf_p_star == 0.0f) {
                U1_x[idx2d(y, x)] = U0_x[idx2d(y, x)];
            } else {
                U1_x[idx2d(y, x)] = U0_x[idx2d(y, x)] + DT * DRIVING_FORCE * (hf_p * dp_star / hf_p_star);
            }
        }
    }
    set_boundary_values(U1_y, 1);
    set_boundary_values(U1_x, 2);
}

static void attenuate(float* field, int key) {
    for (int y = 1; y < CELLS_Y - 1; ++y) {
        for (int x = 1; x < CELLS_X - 1; ++x) {
            field[idx2d(y, x)] /= DT * ATTENUATION + 1.0f;
        }
    }
    set_boundary_values(field, key);
}

// TODO - causes NaNs
// make sure to pass in the target densities here
static void gather(float* S1, float* S0, float* ps, float* p_sblur) {
    memcpy(S1, S0, sizeof(float) * num_cells_s);
    float update0, update1;
    for (int _ = 0; _ < NUM_ITER; ++_) {
        for (int y = 1; y < CELLS_Y - 1; ++y) {
            for (int x = 1; x < CELLS_X - 1; ++x) {
                update0 = (S1[idx2d(y + 1, x)] + ps[idx2d(y, x)] - ps[idx2d(y + 1, x)]) * p_sblur[idx2d(y + 1, x)] * S1[idx2d(y + 1, x)]
                        + (S1[idx2d(y - 1, x)] + ps[idx2d(y, x)] - ps[idx2d(y - 1, x)]) * p_sblur[idx2d(y - 1, x)] * S1[idx2d(y - 1, x)]
                        + (S1[idx2d(y, x + 1)] + ps[idx2d(y, x)] - ps[idx2d(y, x + 1)]) * p_sblur[idx2d(y, x + 1)] * S1[idx2d(y, x + 1)]
                        + (S1[idx2d(y, x - 1)] + ps[idx2d(y, x)] - ps[idx2d(y, x - 1)]) * p_sblur[idx2d(y, x - 1)] * S1[idx2d(y, x - 1)];

                update1 = p_sblur[idx2d(y + 1, x)] * S1[idx2d(y + 1, x)] +
                          p_sblur[idx2d(y, x + 1)] * S1[idx2d(y, x + 1)] +
                          p_sblur[idx2d(y, x)] * 2.0f;

                S1[idx2d(y, x)] = (S0[idx2d(y, x)] + DT * GATHER_RATE * update0)
                        / (DT * GATHER_RATE * update1 + 1.0f);
            }
        }
        set_boundary_values(S1, 0);
    }
}

// convolves the field with a 3x3 Gaussian kernel (zero/"same" padding)
void solver::gaussian_blur(float* outfield, float* infield, int key) {
    for (int y = 0; y < CELLS_Y; ++y) {
        for (int x = 0; x < CELLS_X; ++x) {
            if (y > 0 && x > 0)
                outfield[idx2d(y, x)] += 1.0f * infield[idx2d(y - 1, x - 1)];
            if (y > 0)
                outfield[idx2d(y, x)] += 2.0f * infield[idx2d(y - 1, x)];
            if (y > 0 && x < CELLS_X - 1)
                outfield[idx2d(y, x)] += 1.0f * infield[idx2d(y - 1, x + 1)];
            if (x > 0)
                outfield[idx2d(y, x)] += 2.0f * infield[idx2d(y, x - 1)];
            outfield[idx2d(y, x)]     += 4.0f * infield[idx2d(y, x)];
            if (x < CELLS_X - 1)
                outfield[idx2d(y, x)] += 2.0f * infield[idx2d(y, x + 1)];
            if (y < CELLS_Y - 1 && x > 0)
                outfield[idx2d(y, x)] += 1.0f * infield[idx2d(y + 1, x - 1)];
            if (y < CELLS_Y - 1)
                outfield[idx2d(y, x)] += 2.0f * infield[idx2d(y + 1, x)];
            if (y < CELLS_Y - 1 && x < CELLS_X - 1)
                outfield[idx2d(y, x)] += 1.0f * infield[idx2d(y + 1, x + 1)];
            outfield[idx2d(y, x)] /= 16.0f;
        }
    }
}

void solver::v_step_td(float* U1_y, float* U1_x, float* U0_y, float* U0_x,
                       float* target_p, float* target_p_blur, float* S, float* S_blur) {
    // apply the driving force
    drive_force(U0_y, U0_x, U1_y, U1_x, target_p_blur, S_blur);

    // attenuate momentum
    if (ATTENUATION > 0.0f) {
        attenuate(U0_y, 1);
        attenuate(U0_x, 2);
    }

    // advect the field through itself
    transport(U1_y, U0_y, U0_y, U0_x, 1);
    transport(U1_x, U0_x, U0_y, U0_x, 2);

    // add vorticity
    confine_vorticity(U1_y, U1_x);

    // ensure incompressibility via pressure correction
    project(U0_y, U0_x, U1_y, U1_x);
}

void solver::s_step_td(float* S1, float* S0, float* U_y, float* U_x,
                       float* target_p, float* target_p_blur) {
    transport(S0, S1, U_y, U_x, 0);
    // gather(S1, S0, target_p, target_p_blur);
    // std::swap(S0, S1);
}
