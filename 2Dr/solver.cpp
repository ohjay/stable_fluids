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

    int idx_outer, idx_inner;
    for (int y = 1; y < cells_y - 1; ++y) {
        // horizontal boundary reverse
        idx_outer = (key == 2) ? idx2dx(y, 0) : idx2d(y, 0);
        idx_inner = (key == 2) ? idx2dx(y, 1) : idx2d(y, 1);
        field[idx_outer] = (key == 2) ? fabs(field[idx_inner]) : field[idx_inner];

        idx_outer = (key == 2) ? idx2dx(y, cells_x - 1) : idx2d(y, cells_x - 1);
        idx_inner = (key == 2) ? idx2dx(y, cells_x - 2) : idx2d(y, cells_x - 2);
        field[idx_outer] = (key == 2) ? -fabs(field[idx_inner]) : field[idx_inner];
    }

    for (int x = 1; x < cells_x - 1; ++x) {
        // vertical boundary reverse
        idx_outer = (key == 2) ? idx2dx(0, x) : idx2d(0, x);
        idx_inner = (key == 2) ? idx2dx(1, x) : idx2d(1, x);
        field[idx_outer] = (key == 1) ? fabs(field[idx_inner]) : field[idx_inner];

        idx_outer = (key == 2) ? idx2dx(cells_y - 1, x) : idx2d(cells_y - 1, x);
        idx_inner = (key == 2) ? idx2dx(cells_y - 2, x) : idx2d(cells_y - 2, x);
        field[idx_outer] = (key == 1) ? -fabs(field[idx_inner]) : field[idx_inner];
    }

    // corners
    field[idx2d(0, 0)] = (field[idx2d(1, 0)] + field[idx2d(0, 1)]) * 0.5f;
    field[idx2d(0, cells_x - 1)] = (field[idx2d(1, cells_x - 1)] +
                                    field[idx2d(0, cells_x - 2)]) * 0.5f;
    field[idx2d(cells_y - 1, 0)] = (field[idx2d(cells_y - 2, 0)] +
                                    field[idx2d(cells_y - 1, 1)]) * 0.5f;
    field[idx2d(cells_y - 1, cells_x - 1)] = (field[idx2d(cells_y - 2, cells_x - 1)]
                                            + field[idx2d(cells_y - 1, cells_x - 2)]) * 0.5f;
}

static float dot(float* vec0, float* vec1, int size) {
    float result = 0.0f;
    for (int i = 0; i < size; ++i) {
        result += vec0[i] * vec1[i];
    }
    return result;
}

static bool done = true;
static bool done2 = false;
static int i = 0;

static void poisson2d(float k1, float k2, float* S1, float* S0, int key, int b) {
    int cells_y = (key == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (key == 2) ? CELLS_X + 1 : CELLS_X;
    int y, x, i;
    float pTv, rTr, a, nrTnr, g;

    // r = b - Ax
    float r[num_cells[key]];
    for (y = 0; y < cells_y; ++y) {
        for (x = 0; x < cells_x; ++x) {
            float Ax_yx;
            int idx_yx = idx2d(y, x);
            Ax_yx = (1 - 2 * k1 - 2 * k2) * S1[idx_yx]
                  + k1 * ((y < cells_y - 2) ? S1[idx2d(y + 1, x)] : 0.0f)
                  + k2 * ((y > 1) ?           S1[idx2d(y - 1, x)] : 0.0f)
                  + k1 * ((x < cells_x - 2) ? S1[idx2d(y, x + 1)] : 0.0f)
                  + k2 * ((x > 1) ?           S1[idx2d(y, x - 1)] : 0.0f);

            r[idx_yx] = S0[idx_yx] - Ax_yx;
        }
    }

    // p = r
    float p[num_cells[key]];
    memcpy(p, r, sizeof(p));

    // new_r = r
    float new_r[num_cells[key]];
    memcpy(new_r, r, sizeof(new_r));

    float v[num_cells[key]];
    for (int _ = 0; _ < NUM_ITER; ++_) {
        // v = Ap
        for (y = 0; y < cells_y; ++y) {
            for (x = 0; x < cells_x; ++x) {
                int idx_yx = idx2d(y, x);
                v[idx_yx] = (1 - 2 * k1 - 2 * k2) * p[idx_yx]
                          + k1 * ((y < cells_y - 2) ? p[idx2d(y + 1, x)] : 0.0f)
                          + k2 * ((y > 1) ?           p[idx2d(y - 1, x)] : 0.0f)
                          + k1 * ((x < cells_x - 2) ? p[idx2d(y, x + 1)] : 0.0f)
                          + k2 * ((x > 1) ?           p[idx2d(y, x - 1)] : 0.0f);
            }
        }

        if (!done) { cout << _ << "\nv[40]: " << v[40] << endl; }

        // a = dot(r, r) / dot(p, v)
        pTv = dot(p, v, num_cells[key]);
        if (!done) { cout << "pTv: " << pTv << endl; }
        rTr = dot(r, r, num_cells[key]);
        if (!done) { cout << "rTr: " << rTr << endl; }
        a = (pTv != 0.0f) ? rTr / pTv : 0.0f;
        if (!done) { cout << "a: " << a << endl; }

        // x += a * p; new_r -= av
        for (i = 0; i < num_cells[key]; ++i) {
            S1[i] += a * p[i];
            if (fabs(v[i]) < 21) {
                new_r[i] -= a * v[i];
            }

            // if (!done && fabs(new_r[i]) > 21) { cout << "i: " << i << ": " << v[i] << " / of " << num_cells[key] - 1 << endl; }
        }

        if (!done) { cout << "S1[40]: " << S1[40] << endl; }
        if (!done) { cout << "new_r[40]: " << new_r[40] << endl; }

        // g = dot(new_r, new_r) / dot(r, r)
        nrTnr = dot(new_r, new_r, num_cells[key]);
        if (!done) { cout << "nrTnr: " << nrTnr << endl; }
        g = (rTr != 0.0f) ? nrTnr / rTr : 0.0f;
        if (!done) { cout << "g: " << g << endl; }

        // p = new_r + g * p
        for (i = 0; i < num_cells[key]; ++i) {
            p[i] = new_r[i] + g * p[i];
        }

        if (!done) { cout << "p[40]: " << p[40] << endl; }

        // r = new_r
        memcpy(r, new_r, sizeof(r));

        if (!done) { cout << "r[40]: " << r[40] << "\n\n ---\n\n" << endl; }
        if (_ >= 10) { done = true; }

        if (b == -1) {
            set_boundaries(S1, 0, key);
        } else if (b > 2) {
            boundary_reverse(S1, 1);
            boundary_reverse(S1, 2);
        }
    }
}

static void diffuse(float* S1, float* S0, int key) {
    memset(S1, 0, sizeof(float) * num_cells[key]);
    float k1 = -DT * DIFFUSION;
    float k2 = -DT * DIFFUSION;
    poisson2d(k1, k2, S1, S0, key, -1);
}

static void project(float* U1_y, float* U1_x, float* U0_y, float* U0_x) {
    int y, x, idx_yx;

    float k1 = 1.0f;
    float k2 = 1.0f;

    // construct initial guess for the solution
    float S[num_cells_s];
    memset(S, 0, sizeof(float) * num_cells_s);

    // compute the divergence of the velocity field
    float divergence[num_cells[0]];
    for (y = 1; y < CELLS_Y - 1; ++y) {
        for (x = 1; x < CELLS_X - 1; ++x) {
            idx_yx = idx2d(y, x);
            divergence[idx_yx] = U0_y[idx2d(y + 1, x)] - U0_y[idx_yx]
                               + U0_x[idx2dx(y, x + 1)] - U0_x[idx2dx(y, x)];
        }
    }

    boundary_reverse(divergence, -1);
    poisson2d(k1, k2, S, divergence, 0, 3);

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
    boundary_reverse(U1_x, 2);
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
    diffuse(U0_y, U1_y, 1);
    // if (!done2) { done = false; done2 = true; }
    diffuse(U0_x, U1_x, 2);

    // ensure incompressibility via pressure correction
    project(U1_y, U1_x, U0_y, U0_x);
}

void solver::s_step(float* S1, float* S0, float* U_y, float* U_x, float source) {
    add_force(S0, source * DT, 0);
    transport(S1, S0, U_y, U_x, 0);
    diffuse(S0, S1, 0);
    // dissipate(S1, S0);
}
