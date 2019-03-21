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
    RestrictionBox() = default;

    RestrictionBox(const vec3 &vertex0, const vec3 &vertex1) {
//        cout << vertex0 << " " << vertex1 << endl;
        left_right_detection[0] = vertex0.x();
        left_right_detection[1] = vertex1.x();
        bottom_top_detection[0] = vertex0.y();
        bottom_top_detection[1] = vertex1.y();
        front_back_detection[0] = vertex0.z();
        front_back_detection[1] = vertex1.z();
    }

    void restrict_particle(Particle &in_particle) {
        double_t boundary_stiffness = 100.0;
        double_t boundary_dampening = 200.0;
        vec3 normal_arr[6] = {vec3(1, 0, 0), vec3(-1, 0, 0),
                              vec3(0, 1, 0), vec3(0, -1, 0),
                              vec3(0, 0, 1), vec3(0, 0, -1)};
        const vec3 &tmp_pos = in_particle.get_position();
        const vec3 &tmp_vel = in_particle.get_velocity();
        vec3 new_vel = tmp_vel;
        if (tmp_pos.x() < left_right_detection[0]) {
            auto offset = left_right_detection[0] - tmp_pos.x();
            new_vel -= 2 * dot(normal_arr[0], tmp_vel) * normal_arr[0];
        }
        if (tmp_pos.x() > left_right_detection[1]) {
            auto offset = tmp_pos.x() - left_right_detection[1];
            new_vel -= 2 * dot(normal_arr[1], tmp_vel) * normal_arr[1];
        }
        if (tmp_pos.y() < bottom_top_detection[0]) {
            auto offset = bottom_top_detection[0] - tmp_pos.y();

            new_vel -= 2 * dot(normal_arr[2], tmp_vel) * normal_arr[2];
        }
        if (tmp_pos.y() > bottom_top_detection[1]) {
            auto offset = tmp_pos.y() - bottom_top_detection[1];
            new_vel -= 2 * dot(normal_arr[3], tmp_vel) * normal_arr[3];
        }

        if (tmp_pos.z() < front_back_detection[0]) {
            auto offset = front_back_detection[0] - tmp_pos.z();
            new_vel -= 2 * dot(normal_arr[4], tmp_vel) * normal_arr[4];
        }
        if (tmp_pos.z() > front_back_detection[1]) {
            auto offset = tmp_pos.z() - front_back_detection[1];
            new_vel -= 2 * dot(normal_arr[5], tmp_vel) * normal_arr[5];
        }
        in_particle.set_velocity(new_vel);
    }

    vec3 reflect_vector(const vec3 &in_vec, const vec3 &in_normal) {
        return in_vec - 2 * dot(in_vec, in_normal) * in_normal;
    }

private:

    bool detect_restriction(const vec3 &in_pos) {
        return left_right_detection[0] > in_pos.x() ||
               left_right_detection[1] < in_pos.x() ||
               bottom_top_detection[0] > in_pos.y() ||
               bottom_top_detection[1] < in_pos.y() ||
               front_back_detection[0] > in_pos.z() ||
               front_back_detection[1] < in_pos.z();
    }

    bool intersect_detection(const vec3 &in_vec_pos, const vec3 &in_vec_dir,
                             const vec3 &in_plane_point, const vec3 &in_plane_dir,
                             double &in_dis) {
        double vec_dot_plane_n = dot(in_vec_dir, in_plane_dir);
        if (vec_dot_plane_n >= 0) {
            in_dis = MAXFLOAT;
            return false;
        }

        double tmp_value = dot(in_plane_dir, in_plane_point - in_vec_pos);
        in_dis = tmp_value / vec_dot_plane_n;
//        cout << in_vec_pos + in_dis * in_vec_dir << endl;
        return in_dis <= 0;
    }


};


#endif //SPH_FLUID_SOLVER_RESTRICTIONBOX_H
