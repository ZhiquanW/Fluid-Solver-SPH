//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_FLUIDSOLVER_H
#define SPH_FLUID_SOLVER_FLUIDSOLVER_H

#include "FluidDatabase.h"
class FluidSolver {
private:
    FluidDatabase fluid_database;
    vector<Particle> current_particle_list;

};


#endif //SPH_FLUID_SOLVER_FLUIDSOLVER_H
