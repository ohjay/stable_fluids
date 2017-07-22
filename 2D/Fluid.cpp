#include "Fluid.h"

void Fluid::swap_grids(void) {
    float* temp;
    temp = U0_y; U0_y = U1_y; U1_y = temp;
    temp = U0_x; U0_x = U1_x; U1_x = temp;
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S0[i] = S1[i];
    }
}

void Fluid::init(void) {
    U0_y = new float[num_cells_uy]();
    U0_x = new float[num_cells_ux]();
    U1_y = new float[num_cells_uy]();
    U1_x = new float[num_cells_ux]();

    S0 = new float*[NUM_FLUIDS]();
    S1 = new float*[NUM_FLUIDS]();
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S0[i] = new float[num_cells_s]();
        S1[i] = new float[num_cells_s]();
    }

    target_driven = false;
    memset(S_blur, 0, sizeof(float) * num_cells_s);
}

void Fluid::step(float force_y, float force_x, float source) {
    if (target_driven) {
        solver::gaussian_blur(S_blur, S0[target_fluid], 0);
        solver::v_step_td(U1_y, U1_x, U0_y, U0_x, target_p, target_p_blur, S0[target_fluid], S_blur);
        solver::s_step_td(S1[target_fluid], S0[target_fluid], U1_y, U1_x, source, target_p, target_p_blur);
    } else {
        solver::v_step(U1_y, U1_x, U0_y, U0_x, force_y, force_x);
        for (int i = 0; i < NUM_FLUIDS; i++) {
            solver::s_step(S1[i], S0[i], U1_y, U1_x, source);
        }
    }
    swap_grids();
}

void Fluid::cleanup(void) {
    delete[] U0_y;
    delete[] U0_x;
    delete[] U1_y;
    delete[] U1_x;

    for (int i = 0; i < NUM_FLUIDS; i++) {
        delete[] S0[i];
        delete[] S1[i];
    }
    delete[] S0;
    delete[] S1;
}

void Fluid::add_source_at(int y, int x, int i, float source) {
    if (y > 0 && y < CELLS_Y - 1 && x > 0 && x < CELLS_X - 1) {
        S1[i][idx2d(y, x)] += source;
    }
}

float Fluid::Uy_at(int y, int x) {
    return U1_y[idx2d(y, x)];
}

float Fluid::Ux_at(int y, int x) {
    return U1_x[idx2dx(y, x)];
}

float Fluid::S_at(int y, int x, int i) {
    return S1[i][idx2d(y, x)];
}

void Fluid::save_density(int i) {
    memcpy(target_p, S1[i], sizeof(float) * num_cells_s);
    solver::gaussian_blur(target_p_blur, target_p, 0);
    target_fluid = i;
    memset(S1[i], 0, sizeof(float) * num_cells_s);
}

void Fluid::toggle_target_driven(void) {
    target_driven = !target_driven;
}
