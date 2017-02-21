#include <math.h>

namespace solver {
    // functions we will actually call from outside
    float v_step(float* U1, float* U0, float visc, float* F, float dt);
    float s_step(float* S1, float* S0, float a, float* U, float* source, float dt);
    float transport(float* S1, float* S0, float* U, float dt);
}
