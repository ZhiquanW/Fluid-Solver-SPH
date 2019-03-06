//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_RESTRICTIONBOX_H
#define SPH_FLUID_SOLVER_RESTRICTIONBOX_H

#include "Vector3.h"
#include "Particle.h"

class RestrictionBox {
    double left_right_detection[2]{};
    double front_back_detection[2]{};
    double bottom_top_detection[2]{};

    RestrictionBox(const Vector3 &vertex0, const Vector3 &vertex1) {
        left_right_detection[0] = vertex0.x();
        left_right_detection[1] = vertex1.x();
        bottom_top_detection[0] = vertex0.y();
        bottom_top_detection[1] = vertex1.y();
        front_back_detection[0] = vertex0.z();
        front_back_detection[1] = vertex1.z();
    }

    bool detect_restriction(const Particle &in_particle) {
        bool is_outside = false;
        if (left_right_detection[0] > in_particle.get_position().x() ||
            left_right_detection[1] < in_particle.get_position().x() ||
            bottom_top_detection[0] > in_particle.get_position().y() ||
            bottom_top_detection[1] < in_particle.get_position().y() ||
            front_back_detection[0] > in_particle.get_position().z() ||
            front_back_detection[1] < in_particle.get_position().z()) {
            is_outside = true;
        }else{
            return false;
        }
        Vector3 normal_arr[6] = {Vector3(-1, 0, 0), Vector3(1, 0, 0),
                                 Vector3(0, -1, 0), Vector3(0, 1, 0),
                                 Vector3(0, 0, -1), Vector3(0, 0, 1)};

    }
};


#endif //SPH_FLUID_SOLVER_RESTRICTIONBOX_H
