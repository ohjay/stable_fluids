#include "Fluid.h"

void Fluid::swap_grids(void) {
    float* temp;
    temp = U0_z; U0_z = U1_z; U1_z = temp;
    temp = U0_y; U0_y = U1_y; U1_y = temp;
    temp = U0_x; U0_x = U1_x; U1_x = temp;
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S0[i] = S1[i];
    }
}

void Fluid::init(void) {
    U0_z = new float[num_cells_uz]();
    U0_y = new float[num_cells_uy]();
    U0_x = new float[num_cells_ux]();
    U1_z = new float[num_cells_uz]();
    U1_y = new float[num_cells_uy]();
    U1_x = new float[num_cells_ux]();

    S0 = new float*[NUM_FLUIDS]();
    S1 = new float*[NUM_FLUIDS]();
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S0[i] = new float[num_cells_s]();
        S1[i] = new float[num_cells_s]();
    }
}

void Fluid::step(float force_z, float force_y, float force_x, float source) {
    solver::v_step(U1_z, U1_y, U1_x, U0_z, U0_y, U0_x, force_z, force_y, force_x);
    for (int i = 0; i < NUM_FLUIDS; i++) {
        solver::s_step(S1[i], S0[i], U1_z, U1_y, U1_x, source);
    }
    swap_grids();
}

void Fluid::cleanup(void) {
    delete[] U0_z;
    delete[] U0_y;
    delete[] U0_x;
    delete[] U1_z;
    delete[] U1_y;
    delete[] U1_x;

    for (int i = 0; i < NUM_FLUIDS; i++) {
        delete[] S0[i];
        delete[] S1[i];
    }
    delete[] S0;
    delete[] S1;
}

void Fluid::add_source_at(int z, int y, int x, int i, float source) {
    if (z > 0 && z < CELLS_Z - 1 && y > 0 && y < CELLS_Y - 1 && x > 0 && x < CELLS_X - 1) {
        S1[i][idx3d(z, y, x)] += source;
    }
}

float Fluid::Uz_at(int z, int y, int x) {
    return U1_z[idx3d(z, y, x)];
}

float Fluid::Uy_at(int z, int y, int x) {
    return U1_y[idx3d(z, y, x)];
}

float Fluid::Ux_at(int z, int y, int x) {
    return U1_x[idx3dx(z, y, x)];
}

float Fluid::S_at(int z, int y, int x, int i) {
    return S1[i][idx3d(z, y, x)];
}
