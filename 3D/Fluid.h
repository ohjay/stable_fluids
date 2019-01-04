#include <iostream>
#include "solver.h"

class Fluid {
private:
    float *U0_z, *U0_y, *U0_x, *U1_z, *U1_y, *U1_x; // velocity grids
    float **S0, **S1; // scalar grids

    void swap_grids(void);
public:
    void init(void);
    void step(void);
    void cleanup(void);

    // setters, essentially
    void add_U_z_force_at(int z, int y, int x, float force);
    void add_U_y_force_at(int z, int y, int x, float force);
    void add_U_x_force_at(int z, int y, int x, float force);
    void add_source_at(int z, int y, int x, int i, float source);

    // getters
    float Uz_at(int z, int y, int x);
    float Uy_at(int z, int y, int x);
    float Ux_at(int z, int y, int x);
    float S_at(int z, int y, int x, int i);
};
