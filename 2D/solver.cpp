#include "solver.h"

static float trace_particle(float* X, float* U, float dt, float* X0) {
    // trace a path starting at X through the field U over a time -dt
    return 0.0f;
}

static float lin_interp(float* X0, float* S) {
    // linearly interpolates value of scalar field S at the location X0
    return 0.0f;
}

static float diffuse(float* S0, float* S1, float k, float dt) {
    // solve for the diffusion
    return 0.0f;
}

static float project(float* U1, float* U0, float dt) {
    // perform the projection
    return 0.0f;
}

static float solve_lin(float** A, float* b) {
    // sparse linear solver for the equation Ax = b (recommended: conjugate gradient)
    return 0.0f;
}

static float dissipate(float* S1, float* S0, float a, float dt) {
    // divide each element of first array by (1 + dt * a) and store in new array
    return 0.0f;
}

float solver::v_step(float* U1, float* U0, float visc, float* F, float dt) {
    // velocity solver
    
    /* Vstep pseudocode from paper:
     * 
     * for (i = 0; i < NDIM; i++)
     *   addForce(U0[i], F[i], dt);
     * for (i = 0; i < NDIM; i++)
     *   Transport(U1[i], U0[i], U0, dt);
     * for (i = 0; i < NDIM; i++)
     *   Diffuse(U0[i], U1[i], visc, dt);
     * Project(U1, U0, dt);
     */
    return 0.0f;
}

float solver::s_step(float* S1, float* S0, float a, float* U, float* source, float dt) {
    // scalar field solver
    
   /* Sstep pseudocode from paper:
    * 
    * addForce(S0, source, dt);
    * Transport(S1, S0, U, dt);
    * Diffuse(S0, S1, k, dt);
    * Dissipate(S1, S0, a, dt);
    */
    return 0.0f;
}

float solver::transport(float* S1, float* S0, float* U, float dt) {
    // accounts for movement of substance due to velocity field
    
   /* Transport pseudocode from paper:
    * 
    * for each cell (i, j, k) do
    *   X = O + (i + 0.5, j + 0.5, k + 0.5) * D
    *   TraceParticle(X, U, -dt, X0);
    *   S1[i, j, k] = LinInterp(X0, S0);
    * end
    */
    return 0.0f;
}
