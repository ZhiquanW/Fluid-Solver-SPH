//
// Created by zhiquan on 3/6/19.
//

#ifndef SPH_FLUID_SOLVER_RESTRICTIONBOX_H
#define SPH_FLUID_SOLVER_RESTRICTIONBOX_H

#include <cmath>
#include "Vector3.h"
#include "Particle.h"

#define vec3 Vector3

class RestrictionBox {
private:
    double left_right_detection[2]{};
    double front_back_detection[2]{};
    double bottom_top_detection[2]{};
public:
    RestrictionBox(const vec3 &vertex0, const vec3 &vertex1) {
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
        } else {
            return false;
        }
        vec3 normal_arr[6] = {vec3(-1, 0, 0), vec3(1, 0, 0),
                              vec3(0, -1, 0), vec3(0, 1, 0),
                              vec3(0, 0, -1), vec3(0, 0, 1)};
        
        return true;

    }

private:
    double intersection_detection(const vec3 &in_vec_pos,const vec3 &in_vec_dir,
                                  const vec3 &in_plane_point,const vec3 &in_plane_dir) {
        double vec_dot_plane_n = dot(in_vec_dir,in_plane_dir);
        if(fabs(vec_dot_plane_n) <= 0.001f){
            return 0;
        }
        if(in_vec_pos.x() == in_plane_point.x() &&
           in_vec_pos.y() == in_plane_point.y() &&
           in_vec_pos.z() == in_plane_point.z()){
            return 0;
        }
        double tmp_value = dot(in_vec_dir,in_plane_dir-in_vec_pos);
        return tmp_value/vec_dot_plane_n;
    }

};


#endif //SPH_FLUID_SOLVER_RESTRICTIONBOX_H
