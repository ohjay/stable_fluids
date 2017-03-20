#include <iostream>
#include "solver.h"

class Fluid {
private:
    float *U0_y, *U0_x, *U1_y, *U1_x; // staggered velocity grids
    float **S0, **S1; // scalar grids
    void swap_grids(void);
public:
    void init(void);
    void step(float force_y, float force_x, float source);
    void cleanup(void);

    // setters, essentially
    void add_source_at(int y, int x, int i, float source);

    // getters
    float Uy_at(int y, int x);
    float Ux_at(int y, int x);
    float S_at(int y, int x, int i);
};
