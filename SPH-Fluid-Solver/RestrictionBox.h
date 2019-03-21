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
        int flag[6] = {0};
        if (tmp_pos.x() < left_right_detection[0]) {
            auto offset = left_right_detection[0] - tmp_pos.x();
            new_vel -=  2 * dot(normal_arr[0], tmp_vel) * normal_arr[0];
            flag[0] = 1;
//            cout << "AAAAAAAAAAAAA: " << tmp_pos <<" == " << left_right_detection[0] << endl;
        }
        if (tmp_pos.x() > left_right_detection[1]) {
            auto offset = tmp_pos.x() - left_right_detection[1];
            new_vel -=  2 * dot(normal_arr[1], tmp_vel) * normal_arr[1];
            flag[1] = 1;

        }
        if (tmp_pos.y() < bottom_top_detection[0]) {
            auto offset = bottom_top_detection[0] - tmp_pos.y();

            new_vel -=  2 * dot(normal_arr[2], tmp_vel) * normal_arr[2];
            cout << endl;
            cout << normal_arr[2] << endl;
            cout << tmp_vel << endl;
            cout << dot(normal_arr[2], tmp_vel) << endl;
            cout << 2 * dot(normal_arr[2], tmp_vel) * normal_arr[2] << endl;
            cout << new_vel << endl;
            flag[2] = 1;

        }
        if (tmp_pos.y() > bottom_top_detection[1]) {
            auto offset = tmp_pos.y() - bottom_top_detection[1];
            new_vel -= 2 * dot(normal_arr[3], tmp_vel) * normal_arr[3];
            flag[3] = 1;

        }

        if (tmp_pos.z() < front_back_detection[0]) {
            auto offset = front_back_detection[0] - tmp_pos.z();
            new_vel -=  2 * dot(normal_arr[4], tmp_vel) * normal_arr[4];
            flag[4] = 1;

        }
        if (tmp_pos.z() > front_back_detection[1]) {
            auto offset = tmp_pos.z() - front_back_detection[1];
            new_vel -=  2 * dot(normal_arr[5], tmp_vel) * normal_arr[5];
            flag[5] = 1;

        }
//        cout << in_particle.get_index() << " : " << tmp_pos << endl;
//        for(int i = 0;i < 6;++ i){
//            cout << flag[i] <<" ";
//        }cout << endl;
//        cout << in_particle.get_velocity() << " - " << new_vel << endl;
        in_particle.set_velocity(new_vel);
//        if (in_particle.get_index() == 450) {
//            cout << in_particle.get_index() << endl << in_particle.get_velocity() << endl
//                 << in_particle.get_position() << endl;
//        }
//        if (in_particle.get_velocity().length() == 0) {
//            return;
//        }
//        if (detect_restriction(in_particle.get_position())) {
//            vec3 normal_arr[6] = {vec3(1, 0, 0), vec3(-1, 0, 0),
//                                  vec3(0, 1, 0), vec3(0, -1, 0),
//                                  vec3(0, 0, 1), vec3(0, 0, -1)};
//            vec3 points_arr[6] = {vec3(left_right_detection[0], bottom_top_detection[0], front_back_detection[0]),
//                                  vec3(left_right_detection[1], bottom_top_detection[1], front_back_detection[1]),
//                                  vec3(left_right_detection[0], bottom_top_detection[0], front_back_detection[0]),
//                                  vec3(left_right_detection[1], bottom_top_detection[1], front_back_detection[1]),
//                                  vec3(left_right_detection[0], bottom_top_detection[0], front_back_detection[0]),
//                                  vec3(left_right_detection[1], bottom_top_detection[1], front_back_detection[1])};;
//            double nearest_dis = -MAXFLOAT;//non-positive
//            int nearest_index = -1;
//            for (int i = 0; i < 6; ++i) {
//                double tmp_dis = 0;
//                bool is_intersected = intersect_detection(in_particle.get_position(),
//                                                          in_particle.get_velocity().normalize(),
//                                                          points_arr[i], normal_arr[i], tmp_dis);
////                cout << "tmp dis" << tmp_dis << endl;
//                if (is_intersected && tmp_dis > nearest_dis&&
//                    !detect_restriction(in_particle.get_position() + tmp_dis * in_particle.get_velocity().normalize())) {
//                    nearest_dis = tmp_dis;
//                    nearest_index = i;
//                }
//            }
//
//            if (nearest_index == -1) {
//                cout << in_particle.get_index() << endl << in_particle.get_velocity() << endl
//                     << in_particle.get_position() << endl;
//            }
//            in_particle.set_position(in_particle.get_position() + nearest_dis * in_particle.get_velocity().normalize());
//            in_particle.set_velocity(reflect_vector(in_particle.get_velocity(), normal_arr[nearest_index]));
//
//        }
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
