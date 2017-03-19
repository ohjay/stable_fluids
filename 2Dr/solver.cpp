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

void print_fl_array(float* arr, int n, string label="") {
    label = (label.empty()) ? "" : "[" + label + "] ";

    int idx1 = (int) n / 4, idx2 = (int) n / 2, idx3 = (int) 3 * n / 4;
    float val1 = arr[idx1], val2 = arr[idx2], val3 = arr[idx3];
    cout << "---" << endl;
    cout << label << "idx " << idx1 << ": " << val1 << "; idx "
            << idx2 << ": " << val2 << "; idx " << idx3 << ": " << val3 << endl;
}

void solver::v_step(float* U1_y, float* U1_x, float* U0_y, float* U0_x, float force_y, float force_x) {
    // add forces
    add_force(U0_y, force_y * DT, 1);
    add_force(U0_x, force_x * DT, 2);

    // self-advect
    transport(U1_y, U0_y, U0_y, U0_x, 1);
    transport(U1_x, U0_x, U0_y, U0_x, 2);

    /*
    // diffuse
    diffuse(U0_y, U1_y, 1);
    diffuse(U0_x, U1_x, 2);

    // ensure incompressibility via pressure correction
    project(U1_y, U1_x, U0_y, U0_x);
    */
}

void solver::s_step(float* S1, float* S0, float* U_y, float* U_x, float source) {
    add_force(S0, source * DT, 0);

    // print_fl_array(S0, CELLS_X * CELLS_Y, "1");
    transport(S1, S0, U_y, U_x, 0);
    // print_fl_array(S1, CELLS_X * CELLS_Y, "2");
    
    /*
    diffuse(S0, S1, 0);
    dissipate(S1, S0);
    */
}
