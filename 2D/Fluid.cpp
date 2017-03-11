#include "Fluid.h"

static void swap2d(float*** arr0, float*** arr1) {
    float** temp = *arr0;
    *arr0 = *arr1;
    *arr1 = temp;
}

static void swap1d(float** arr0, float** arr1) {
    float* temp = *arr0;
    *arr0 = *arr1;
    *arr1 = temp;
}

void Fluid::init(float visc, float ks, float as, float dt) {
    (*this).visc = visc;
    (*this).ks = ks;
    (*this).as = as;
    (*this).dt = dt;

    // set grid attributes
    int num_cells = 1;
    for (int i = 0; i < NDIM; ++i) {
        O[i] = 0.f;
        L[i] = 50;
        N[i] = 50;
        D[i] = L[i] / N[i];
        num_cells *= N[i];
    }

    // initialize grids [(TODO) don't forget to delete this memory later]
    U0 = new float*[NDIM]; U1 = new float*[NDIM];
    for (int i = 0; i < NDIM; ++i) {
        // for velocity fields (fluid)
        U0[i] = new float[num_cells];
        U1[i] = new float[num_cells];
    }
    
    // for density fields (substance)
    S0 = new float[num_cells];
    S1 = new float[num_cells];
    
    (*this).num_cells = num_cells;
}

void Fluid::step(float F[2], float Ssource) {
    // handle display and user interaction
    // get forces F and sources Ssource from UI

    // swap U1 and U0, swap S1 and S0
    swap2d(&U1, &U0); swap1d(&S1, &S0);

    // perform a velocity step (using U1, U0, visc, F, and dt)
    solver::v_step(U1, U0, visc, F, dt, num_cells, N, O, D);

    // perform a scalar step (using S1, S0, kS, aS, U1, Ssource, and dt)
    solver::s_step(S1, S0, ks, as, U1, Ssource, dt, num_cells, N, O, D);
}

int Fluid::num_cells_x() {
    return N[1];
}

int Fluid::num_cells_y() {
    return N[0];
}

float Fluid::grid_spacing() {
    return (D[0] > D[1]) ? D[0] : D[1];
}

float Fluid::S_at(int y, int x) {
    int xyz[2];
    xyz[0] = y; xyz[1] = x;
    
    return S1[solver::xyz_to_idx(xyz, N)];
}
