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
    for (int i = 0; i < NDIM; ++i) {
        O[i] = 0.f;
        D[i] = SIDE_LEN / CELLS_PER_SIDE;
    }

    // initialize grids [(TODO) don't forget to delete this memory later]
    U0 = new float*[NDIM]; U1 = new float*[NDIM];
    for (int i = 0; i < NDIM; ++i) {
        // for velocity fields (fluid)
        U0[i] = new float[NUM_CELLS];
        U1[i] = new float[NUM_CELLS];
    }

    // for density fields (substance)
    S0 = new float[NUM_CELLS];
    S1 = new float[NUM_CELLS];
}

void Fluid::step(float F[2], float Ssource, int Fy, int Fx) {
    // handle display and user interaction
    // get forces F and sources Ssource from UI
    print_fl_array_perc(U0[0], NUM_CELLS, 0.02f, "U00");
    print_fl_array_perc(U0[1], NUM_CELLS, 0.02f, "U01");

    // swap U1 and U0, swap S1 and S0
    swap2d(&U1, &U0); swap1d(&S1, &S0);

    // perform a velocity step (using U1, U0, visc, F, and dt)
    solver::v_step(U1, U0, visc, F, dt, O, D);

    // perform a scalar step (using S1, S0, kS, aS, U1, Ssource, and dt)
    solver::s_step(S1, S0, ks, as, U1, Ssource, dt, O, D, Fy, Fx);
}

float Fluid::grid_spacing() {
    return (D[0] > D[1]) ? D[0] : D[1];
}

float Fluid::S_at(int y, int x) {
    return S1[solver::idx2d(y, x)];
}

void Fluid::add_S_at(int y, int x, float source) {
    S1[solver::idx2d(y, x)] += source;
}

void Fluid::cleanup(void) {
    for (int i = 0; i < NDIM; ++i) {
        delete[] U0[i];
        delete[] U1[i];
    }
    delete[] U0;
    delete[] U1;
    delete[] S0;
    delete[] S1;
}
