#include "Fluid.h"

void Fluid::init() {
    // initialize here
}

void Fluid::step() {
    // simulate a step here
    
   /* Simulation pseudocode from paper:
    * 
    * handle display and user interaction
    * get forces F and sources Ssource from UI
    * Swap(U1, U0);
    * Swap(S1, S0);
    * Vstep(U1, U0, visc, F, dt);
    * Sstep(S1, S0, kS, aS, U1, Ssource, dt);
    */
}
