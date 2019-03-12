#include <utility>

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
                FluidDatabase _fluid_database)
            : fluid_parameter(_fluid_parameter),
              restriction_box(_restriction_box),
              fluid_database(std::move(_fluid_database)) {
    }

private:
    const double_t compute_kernel_poly6(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 315.0 / (64.0 * M_PI * pow(_core_radius, 9)) *
               pow(pow(_core_radius, 2) - _offset_vec.squared_distance(), 3);
    }

    const double_t compte_laplacian_kernel_poly6(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        double_t squared_distance = _offset_vec.squared_distance();
        double_t squared_radius = pow(_core_radius, 2);
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 945 / (8 * M_PI * pow(_core_radius, 9)) * pow(squared_radius - squared_distance, 2) *
               (squared_distance - (3.0 / 4.0) * (squared_distance - squared_radius));
    }

    const vec3 compute_hamiltonian_kernel_spiky(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return vec3();
        }
        return -45.0 / (M_PI * pow(_core_radius, 6)) * pow(_core_radius - tmp_dis, 2) * _offset_vec.normalize();
    }

    const double_t compute_laplacian_kernel_viscosity(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 45.0 / (M_PI * pow(_core_radius, 6)) * (_core_radius - tmp_dis);
    }


};


#endif //SPH_FLUID_SOLVER_FLUIDSOLVER_H
