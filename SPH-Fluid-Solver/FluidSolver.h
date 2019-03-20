#include <utility>

#include <utility>

#include <utility>

//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_FLUIDSOLVER_H
#define SPH_FLUID_SOLVER_FLUIDSOLVER_H

#include <cmath>
#include <chrono>

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

    void initialize_particles(const vec3 &center_pos, double _cube_edge_len) {
        cout << "Start Initializing Particles : " << fluid_parameter.get_particle_num() << endl;
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
                        cout << "Particles Initialized" << endl;
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

    void simulate_particles() {
        //Record Start Time
        auto start = std::chrono::system_clock::now();
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        cout << "Start Simulation : " << fluid_database.get_frame_num() << " frames" << endl;
        for (size_t i = 0; i < fluid_database.get_frame_num(); ++i) {
            cout << flush << '\r';
            printf("%.2lf%%", (i + 1.0) * 100.0 / fluid_database.get_frame_num());
            end = std::chrono::system_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            printf("\t%.2f min", (double) duration.count() / std::chrono::microseconds::period::den / 60);
            printf("\tEstimated Time : %.2f min",
                   (double) (fluid_database.get_frame_num() - i) / i * (double) duration.count() /
                   std::chrono::microseconds::period::den / 60);
            compute_particle_acceleration();
            if (false) {
                cout << endl;
                cout << "Index : " << realtime_particle_list[6].get_index() << endl;
                cout << "Density : " << realtime_particle_list[6].get_density() << endl;
                cout << "Pressure : " << realtime_particle_list[6].get_pressure() << endl;
                cout << "Pressure_Acceleration : " << realtime_particle_list[6].get_pressure_acceleration() << endl;
                cout << "Viscosity_Acceleration : " << realtime_particle_list[6].get_viscosity_acceleration() << endl;
                cout << "Surface_Tension_Acceleration : "
                     << realtime_particle_list[6].get_surface_tension_acceleration() << endl;
                cout << "Total_Acceleration : " << realtime_particle_list[6].get_acceleration() << endl;
                cout << "Velocity : " << realtime_particle_list[6].get_velocity() << endl;

            }
            restrict_particles();
//            cout << "AFTER RESTRICT : " << realtime_particle_list[6].get_acceleration() << endl;
//            cout << "Velocity : " << realtime_particle_list[6].get_velocity() << endl;
//            cout << "Position : " << realtime_particle_list[6].get_position() << endl;
//
//            cout << endl;

            update_position();
            fluid_database.append_particle_list(i, realtime_particle_list);
        }
    }

    void output_data() {
        fluid_database.export_database();
    }


private:
    void compute_particle_acceleration() {
        compute_density();
        compute_pressure();
        compute_pressure_acceleration();
        compute_viscosity_acceleration();
        compute_surface_tension_acceleration();
        for (auto &p:realtime_particle_list) {
            p.set_acceleration(p.get_pressure_acceleration() +
                               p.get_viscosity_acceleration() +
                               p.get_surface_tension_acceleration() +
                               fluid_parameter.get_gravity_acceleration_coefficient() * vec3(0, -1, 0));
        }
    }

    void restrict_particles() {
        for (auto &p:realtime_particle_list) {
            restriction_box.restrict_particle(p);
        }
    }

    void update_position() {
        for (auto &p:realtime_particle_list) {
            p.set_velocity(p.get_velocity() + p.get_acceleration() * fluid_database.get_frame_interval());
            p.set_position(p.get_position() + p.get_velocity() * fluid_database.get_frame_interval());
        }
    }

    void compute_density() {
        auto tmp_vector = realtime_particle_list;
        for (auto &p_i:realtime_particle_list) {
            double_t tmp_density = 0;
            for (auto &p_j:tmp_vector) {
//                if (p_i.get_index() == p_j.get_index()) {
//                    continue;
//                }
                if ((p_i.get_position() - p_j.get_position()).length() <= fluid_parameter.get_core_radius()) {
                    tmp_density += fluid_parameter.get_particle_mass() *
                                   compute_kernel_poly6(p_i.get_position() - p_j.get_position(),
                                                        fluid_parameter.get_core_radius());
                }
            }
            p_i.set_density(tmp_density);

        }
    }

    void compute_pressure() {
        for (auto &p_i:realtime_particle_list) {
            p_i.set_pressure(
                    fluid_parameter.get_gas_constant() * (p_i.get_density() - fluid_parameter.get_rest_density()));
        }
    }

    void compute_pressure_acceleration() {
        auto tmp_vector = realtime_particle_list;
        for (auto &p_i:realtime_particle_list) {
            auto tmp_force = vec3();
            for (auto &p_j:tmp_vector) {
                if ((p_i.get_position() - p_j.get_position()).length() <= fluid_parameter.get_core_radius()) {
                    tmp_force += (p_i.get_pressure() + p_j.get_pressure()) / p_j.get_density() *
                                 compute_kernel_spiky_hamiltonian(p_i.get_position() - p_j.get_position(),
                                                                  fluid_parameter.get_core_radius());
                }

            }
            p_i.set_pressure_acceleration(
                    -1.0 * fluid_parameter.get_particle_mass() / (p_i.get_density() * 2) * tmp_force);
        }
    }

    void compute_viscosity_acceleration() {
        auto tmp_vector = realtime_particle_list;
        for (auto &p_i:realtime_particle_list) {
            auto tmp_force = vec3();
            for (auto &p_j:tmp_vector) {
                if ((p_i.get_position() - p_j.get_position()).length() <= fluid_parameter.get_core_radius()) {

                    tmp_force += (p_i.get_velocity() - p_j.get_velocity()) * p_j.get_density() *
                                 compute_kernel_viscosity_laplacian(p_i.get_position() - p_j.get_position(),
                                                                    fluid_parameter.get_core_radius());
                }
            }
            p_i.set_viscosity_acceleration(
                    fluid_parameter.get_viscosity_coefficient() * fluid_parameter.get_particle_mass() /
                    p_i.get_density() * tmp_force);
        }
    }

    void compute_surface_tension_acceleration() {
        auto tmp_vector = realtime_particle_list;
        for (auto &p_i:realtime_particle_list) {
            auto tmp_force = vec3();
            auto tmp_gradient = vec3();
            auto tmp_curvature = 0.0;
            for (auto &p_j:tmp_vector) {
                if ((p_i.get_position() - p_j.get_position()).length() <= fluid_parameter.get_core_radius()) {

                    tmp_gradient += (1.0 / p_j.get_density()) *
                                    compute_kernel_poly6_hamiltonian(p_i.get_position() - p_j.get_position(),
                                                                     fluid_parameter.get_core_radius());
                    tmp_curvature += (1.0 / p_j.get_density()) *
                                     compute_kernel_poly6_laplacian(p_i.get_position() - p_j.get_position(),
                                                                    fluid_parameter.get_core_radius());
                }
            }
            tmp_gradient *= -1.0 * fluid_parameter.get_particle_mass();
            tmp_curvature *= fluid_parameter.get_particle_mass();
            p_i.set_surface_tension_acceleration(
                    -1.0 * fluid_parameter.get_tension_coefficient() / p_i.get_density() * tmp_curvature *
                    tmp_gradient.normalize());
        }
    }

    const double_t compute_kernel_poly6(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 315.0 / (64.0 * M_PI * pow(_core_radius, 9)) *
               pow(pow(_core_radius, 2) - _offset_vec.squared_distance(), 3);
    }

    const vec3 compute_kernel_poly6_hamiltonian(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return vec3();
        }
        return -945.0 / (32.0 * M_PI * pow(_core_radius, 9)) *
               pow((pow(_core_radius, 2) - _offset_vec.squared_distance()), 2) *
               _offset_vec;
    }

    const double_t compute_kernel_poly6_laplacian(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        double_t squared_distance = _offset_vec.squared_distance();
        double_t squared_radius = pow(_core_radius, 2);
        if (tmp_dis > _core_radius) {
            return 0.0;
        }
        return 945.0 / (8.0 * M_PI * pow(_core_radius, 9)) * pow(squared_radius - squared_distance, 2) *
               (squared_distance - (3.0 / 4.0) * (squared_distance - squared_radius));
    }

    const vec3 compute_kernel_spiky_hamiltonian(const vec3 &_offset_vec, const double_t &_core_radius) {
        double_t tmp_dis = _offset_vec.length();
        if (tmp_dis > _core_radius) {
            return vec3();
        }
        return -45.0 / (M_PI * pow(_core_radius, 6)) * pow(_core_radius - tmp_dis, 2) * _offset_vec.normalize();
    }

    const double_t compute_kernel_viscosity_laplacian(const vec3 &_offset_vec, const double_t &_core_radius) {
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
