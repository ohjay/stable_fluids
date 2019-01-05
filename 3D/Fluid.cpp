#include "Fluid.h"

void Fluid::swap_grids(void) {
    float* temp;
    temp = U0_z; U0_z = U1_z; U1_z = temp;
    temp = U0_y; U0_y = U1_y; U1_y = temp;
    temp = U0_x; U0_x = U1_x; U1_x = temp;
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S1[i] = S0[i];
    }
}

void Fluid::init(void) {
    U0_z = new float[num_cells]();
    U0_y = new float[num_cells]();
    U0_x = new float[num_cells]();
    U1_z = new float[num_cells]();
    U1_y = new float[num_cells]();
    U1_x = new float[num_cells]();
    memset(U0_z, 0, sizeof(float) * num_cells);
    memset(U0_y, 0, sizeof(float) * num_cells);
    memset(U0_x, 0, sizeof(float) * num_cells);
    memset(U1_z, 0, sizeof(float) * num_cells);
    memset(U1_y, 0, sizeof(float) * num_cells);
    memset(U1_x, 0, sizeof(float) * num_cells);

    S0 = new float*[NUM_FLUIDS]();
    S1 = new float*[NUM_FLUIDS]();
    for (int i = 0; i < NUM_FLUIDS; i++) {
        S0[i] = new float[num_cells]();
        S1[i] = new float[num_cells]();
        memset(S0[i], 0, sizeof(float) * num_cells);
        memset(S1[i], 0, sizeof(float) * num_cells);
    }
}

void Fluid::step(void) {
    solver::v_step(U1_z, U1_y, U1_x, U0_z, U0_y, U0_x);
    for (int i = 0; i < NUM_FLUIDS; i++) {
        solver::s_step(S1[i], S0[i], U0_z, U0_y, U0_x);
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

void Fluid::add_U_z_force_at(int z, int y, int x, float force) {
    if (z > 0 && z < CELLS_Z - 1 && y > 0 && y < CELLS_Y - 1 && x > 0 && x < CELLS_X - 1) {
        U1_z[idx3d(z, y, x)] += force;
    }
}

void Fluid::add_U_y_force_at(int z, int y, int x, float force) {
    if (z > 0 && z < CELLS_Z - 1 && y > 0 && y < CELLS_Y - 1 && x > 0 && x < CELLS_X - 1) {
        U1_y[idx3d(z, y, x)] += force;
    }
}

void Fluid::add_U_x_force_at(int z, int y, int x, float force) {
    if (z > 0 && z < CELLS_Z - 1 && y > 0 && y < CELLS_Y - 1 && x > 0 && x < CELLS_X - 1) {
        U1_x[idx3d(z, y, x)] += force;
    }
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
    return U1_x[idx3d(z, y, x)];
}

float Fluid::S_at(int z, int y, int x, int i) {
    return S1[i][idx3d(z, y, x)];
}
