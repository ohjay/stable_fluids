#include "Fluid.h"

static void swap(float*** arr0, float*** arr1) {
    float** temp = *arr0;
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

    // initialize grids (don't forget to delete this memory later)
    U0 = new float*[NDIM]; U1 = new float*[NDIM];
    S0 = new float*[NDIM]; S1 = new float*[NDIM];
    for (int i = 0; i < NDIM; ++i) {
        // for velocity fields (fluid)
        U0[i] = new float[num_cells];
        U1[i] = new float[num_cells];

        // for scalar fields (substance)
        S0[i] = new float[num_cells];
        S1[i] = new float[num_cells];
    }
}

void Fluid::step() {
    // handle display and user interaction

    // get forces F and sources Ssource from UI
    float F[3] = {0.f, 0.f, 0.f};

    // swap U1 and U0, swap S1 and S0
    swap(&U1, &U0); swap(&S1, &S0);

    // perform a velocity step (using U1, U0, visc, F, and dt)
    solver::v_step(U1, U0, visc, F, dt);

    // perform a scalar step (using S1, S0, kS, aS, U1, Ssource, and dt)
}
