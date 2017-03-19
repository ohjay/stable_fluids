#include "Fluid.h"

void Fluid::swap_grids(void) {
    float* temp;
    temp = U0_y; U0_y = U1_y; U1_y = temp;
    temp = U0_x; U0_x = U1_x; U1_x = temp;
    temp = S0;   S0 = S1;     S1 = S0;
}

void Fluid::init(void) {
    U0_y = new float[num_cells_uy]();
    U0_x = new float[num_cells_ux]();
    U1_y = new float[num_cells_uy]();
    U1_x = new float[num_cells_ux]();

    S0 = new float[num_cells_s]();
    S1 = new float[num_cells_s]();
}

void Fluid::step(float force_y, float force_x, float source) {
    solver::v_step(U1_y, U1_x, U0_y, U0_x, force_y, force_x);
    solver::s_step(S1, S0, U1_y, U1_x, source);
    swap_grids();
}

void Fluid::cleanup(void) {
    delete[] U0_y;
    delete[] U0_x;
    delete[] U1_y;
    delete[] U1_x;

    delete[] S0;
    delete[] S1;
}

void Fluid::add_source_at(int y, int x, float source) {
    if (y <= 0) { y = 1; }
    else if (y >= CELLS_Y - 1) { y = CELLS_Y - 2; }
    if (x <= 0) { x = 1; }
    else if (x >= CELLS_X - 1) { x = CELLS_X - 2; }

    S1[idx2d(y, x)] += source;
}

float Fluid::Uy_at(int y, int x) {
    return U1_y[idx2d(y, x)];
}

float Fluid::Ux_at(int y, int x) {
    return U1_x[idx2dx(y, x)];
}

float Fluid::S_at(int y, int x) {
    return S1[idx2d(y, x)];
}
