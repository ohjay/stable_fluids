#include <math.h>
#include <iostream>
#include <numeric>

#include "params.h"

namespace solver {
    void v_step(float* U1_z, float* U1_y, float* U1_x, float* U0_z, float* U0_y, float* U0_x, float force_z, float force_y, float force_x);
    void s_step(float* S1, float* S0, float* U1_z, float* U1_y, float* U1_x, float source);
}
