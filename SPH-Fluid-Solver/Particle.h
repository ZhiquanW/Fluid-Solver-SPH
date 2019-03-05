//
// Created by zhiquan on 3/5/19.
//

#ifndef SPH_FLUID_SOLVER_PARTICLE_H
#define SPH_FLUID_SOLVER_PARTICLE_H

#include "Vector3.h"
#define vec3 Vector3
class Particle {
private:
    size_t index;
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    vec3 gravity_acceleration;
    vec3 pressure_acceleration;
    vec3 viscosity_acceleration;
    vec3 surface_tension_acceleration;
public:
    Particle() = default;

    explicit Particle(const size_t &index):index(index){
    }

    const Vector3 &get_position() const {
        return position;
    }

    size_t get_index() const {
        return index;
    }

    void set_index(size_t index) {
        Particle::index = index;
    }


    void set_position(const Vector3 &position) {
        Particle::position = position;
    }

    const Vector3 &get_velocity() const {
        return velocity;
    }

    void set_velocity(const Vector3 &velocity) {
        Particle::velocity = velocity;
    }

    const Vector3 &get_acceleration() const {
        return acceleration;
    }

    void set_acceleration(const Vector3 &acceleration) {
        Particle::acceleration = acceleration;
    }

    const Vector3 &get_gravity_acceleration() const {
        return gravity_acceleration;
    }

    void set_gravity_acceleration(const Vector3 &gravity_acceleration) {
        Particle::gravity_acceleration = gravity_acceleration;
    }

    const Vector3 &get_pressure_acceleration() const {
        return pressure_acceleration;
    }

    void set_pressure_acceleration(const Vector3 &pressure_acceleration) {
        Particle::pressure_acceleration = pressure_acceleration;
    }

    const Vector3 &get_viscosity_acceleration() const {
        return viscosity_acceleration;
    }

    void set_viscosity_acceleration(const Vector3 &viscosity_acceleration) {
        Particle::viscosity_acceleration = viscosity_acceleration;
    }

    const Vector3 &get_surface_tension_acceleration() const {
        return surface_tension_acceleration;
    }

    void set_surface_tension_acceleration(const Vector3 &surface_tension_acceleration) {
        Particle::surface_tension_acceleration = surface_tension_acceleration;
    }



};


#endif //SPH_FLUID_SOLVER_PARTICLE_H