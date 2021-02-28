#include <math.h>
#include <cstring>
#include <iostream>
#include <numeric>
#include <utility>

#include "params.h"

namespace solver {
    void v_step(float* U1_y, float* U1_x, float* U0_y, float* U0_x);
    void s_step(float* S1, float* S0, float* U1_y, float* U1_x);

    // target-driven variants
    void v_step_td(float* U1_y, float* U1_x, float* U0_y, float* U0_x,
                   float* target_p, float* target_p_blur, float* S, float* S_blur);
    void s_step_td(float* S1, float* S0, float* U1_y, float* U1_x,
                   float* target_p, float* target_p_blur);

    void gaussian_blur(float* outfield, float* infield);
}
