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
    FluidParameter fluid_parameter{};
    RestrictionBox restriction_box{};
    FluidDatabase fluid_database{};
    vector<Particle> realtime_particle_list{};

public:

    FluidSolver(const FluidParameter &_fluid_parameter,
                const RestrictionBox &_restriction_box,
                FluidDatabase &_fluid_database)
            : fluid_parameter(_fluid_parameter),
              restriction_box(_restriction_box),
              fluid_database(_fluid_database) {
    }

    void test_kernel_functions() {
        cout << compute_kernel_poly6(vec3(), 2) << endl;
        cout << compute_hamiltonian_kernel_spiky(vec3(0, 1, 0), 2) << endl;
    }

    void initialize_particles(const vec3 &center_pos, double _cube_edge_len) {
        realtime_particle_list.reserve(fluid_parameter.get_particle_num());
        double_t particle_interval = _cube_edge_len / pow((double_t) fluid_parameter.get_particle_num(), 1 / 3.0);
        auto initial_pos = vec3(center_pos.x() - _cube_edge_len / 2, center_pos.y() - _cube_edge_len / 2,
                                center_pos.z() - _cube_edge_len / 2);
        double_t particle_num_per_edge = _cube_edge_len / particle_interval;
        vec3 tmp_vec;
        size_t cur_num = 0;
        for (size_t i = 0; i < particle_num_per_edge; ++i) {
            for (size_t j = 0; j < particle_num_per_edge; ++j) {
                for (size_t k = 0; k < particle_num_per_edge; ++k) {
                    ++cur_num;
                    if (cur_num > fluid_parameter.get_particle_num()) {
                        return;
                    }
                    tmp_vec = initial_pos + vec3(i * particle_interval, j * particle_interval, k * particle_interval);
                    Particle tmp_particle(cur_num);
                    tmp_particle.set_position(tmp_vec);
                    realtime_particle_list.emplace_back(tmp_particle);
                }
            }
        }
    }

    void compute_density() {
        auto tmp_vector = realtime_particle_list;
        //cout << "mass:" << fluid_parameter.get_particle_mass() << endl;
        int counter = 0;
        for (auto &p_i:realtime_particle_list) {
            double_t tmp_density = 0;
            for (auto &p_j:tmp_vector) {
                if (p_i.get_index() == p_j.get_index()) {
                    continue;
                }
                tmp_density += fluid_parameter.get_particle_mass() *
                               compute_kernel_poly6(p_i.get_position() - p_j.get_position(),
                                                    fluid_parameter.get_core_radius());
            }
            p_i.set_density(tmp_density);
            cout << ++counter << " " << tmp_density << endl;

        }
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

public:
    const vector<Particle> &get_realtime_particle_list() const {
        return realtime_particle_list;
    }

    void set_realtime_particle_list(const vector<Particle> &realtime_particle_list) {
        FluidSolver::realtime_particle_list = realtime_particle_list;
    }
};


#endif //SPH_FLUID_SOLVER_FLUIDSOLVER_H
