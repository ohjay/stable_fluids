#include "Fluid.h"

static void swap2d(double*** arr0, double*** arr1) {
    double** temp = *arr0;
    *arr0 = *arr1;
    *arr1 = temp;
}

static void swap1d(double** arr0, double** arr1) {
    double* temp = *arr0;
    *arr0 = *arr1;
    *arr1 = temp;
}

void Fluid::init(double visc, double ks, double as, double dt) {
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
    U0 = new double*[NDIM]; U1 = new double*[NDIM];
    for (int i = 0; i < NDIM; ++i) {
        // for velocity fields (fluid)
        U0[i] = new double[NUM_CELLS]();
        U1[i] = new double[NUM_CELLS]();
    }

    // for density fields (substance)
    S0 = new double[NUM_CELLS]();
    S1 = new double[NUM_CELLS]();
}

void Fluid::step(double F[2], double Ssource, int Fy, int Fx) {
    // handle display and user interaction
    // get forces F and sources Ssource from UI

    // swap U1 and U0, swap S1 and S0
    swap2d(&U1, &U0); swap1d(&S1, &S0);

    // perform a velocity step (using U1, U0, visc, F, and dt)
    solver::v_step(U1, U0, visc, F, dt, O, D);

    // perform a scalar step (using S1, S0, kS, aS, U1, Ssource, and dt)
    solver::s_step(S1, S0, ks, as, U1, Ssource, dt, O, D, Fy, Fx);
}

double Fluid::grid_spacing() {
    return (D[0] > D[1]) ? D[0] : D[1];
}

double Fluid::S_at(int y, int x) {
    return S1[solver::idx2d(y, x)];
}

double Fluid::U10_at(int y, int x) {
    return U1[0][solver::idx2d(y, x)] * 1000;
}

double Fluid::U11_at(int y, int x) {
    return U1[1][solver::idx2d(y, x)] * 1000;
}

void Fluid::add_S_at(int y, int x, double source) {
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
