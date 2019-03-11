#include <utility>

#include <utility>

//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_FLUIDSOLVER_H
#define SPH_FLUID_SOLVER_FLUIDSOLVER_H

#include <math.h>

#include "FluidParameter.h"
#include "FluidDatabase.h"
#include "RestrictionBox.h"

#define vec3 Vector3

class FluidSolver {
private:
    FluidParameter fluid_parameter;
    RestrictionBox restriction_box;
    FluidDatabase fluid_database;
    vector<Particle> realtime_particle_list;

public:
    FluidSolver() = default;

    FluidSolver(const FluidParameter &_fluid_parameter, const RestrictionBox &_restriction_box,
                const FluidDatabase &_fluid_database)
            : fluid_parameter(_fluid_parameter),
              restriction_box(_restriction_box),
              fluid_database(_fluid_database) {
    }

private:
    double_t compute_kernel_poly6(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 315 / (64 * M_PI * pow(_core_radius, 9)) * pow(pow(_core_radius, 2) - squared_distance(_offset_vec), 3);
    }

    const vec3 compute_hamiltonian_kernel_spiky(const vec3 &_offset_vec, const double_t &_core_radius){
        double_t  tmp_dis = _offset_vec.length();
        if(tmp_dis > _core_radius){
            return vec3();
        }
        return -45/(M_PI*pow(_core_radius,6))*pow(_core_radius-tmp_dis,2)*_offset_vec.normalize();
    }

    const vec3 compute_Laplacian_kernel_viscosity(const vec3 &_offset_vec,const double_t & _core_radius){
        double_t  tmp_dis = _offset_vec.length();
    }

};


#endif //SPH_FLUID_SOLVER_FLUIDSOLVER_H
