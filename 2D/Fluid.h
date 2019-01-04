#include <iostream>
#include "solver.h"

class Fluid {
private:
    float *U0_y, *U0_x, *U1_y, *U1_x; // velocity grids
    float **S0, **S1; // scalar grids

    // target-driven attributes
    bool target_driven;
    float target_p[num_cells];
    float target_p_blur[num_cells];
    int target_fluid;
    float S_blur[num_cells];

    void swap_grids(void);
public:
    void init(void);
    void step(void);
    void cleanup(void);

    // setters, essentially
    void add_U_y_force_at(int y, int x, float force);
    void add_U_x_force_at(int y, int x, float force);
    void add_source_at(int y, int x, int i, float source);

    // getters
    float Uy_at(int y, int x);
    float Ux_at(int y, int x);
    float S_at(int y, int x, int i);

    // target-driven methods
    void save_density(int i);
    void toggle_target_driven(void);
};
