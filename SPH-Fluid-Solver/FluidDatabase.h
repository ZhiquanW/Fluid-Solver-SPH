#include <utility>

#include <utility>

//
// Created by zhiquan on 3/5/19.
//

#ifndef SPH_FLUID_SOLVER_FLUIDDATABASE_H
#define SPH_FLUID_SOLVER_FLUIDDATABASE_H

#include <vector>
#include <fstream>
#include <algorithm>
#include "Vector3.h"
#include "Particle.h"

using namespace std;

class FluidDatabase {
private:
    string file_name;
    size_t particle_num{};
    size_t frame_num{};
    size_t frame_interval{};
    size_t animation_duration{};
    vector<vector<Particle>> particle_matrix;
public:

    FluidDatabase() = default;

    explicit FluidDatabase(string file_name, size_t particle_num, size_t frame_num, size_t frame_interval)
            : file_name(std::move(file_name)), particle_num(particle_num), frame_num(frame_num),
              frame_interval(frame_interval) {
        this->animation_duration = frame_num*frame_interval;

        particle_matrix.resize(frame_num);
        for(auto &list:this->particle_matrix){
            list.resize(particle_num);
        }
    }

    void export_database() {
        fstream fout;
        fout.open(this->file_name.c_str(), ios::out);
        fout << file_name << endl;
        fout << particle_num << " " << frame_num << " " << frame_interval << " " << frame_interval << " "
             << animation_duration << endl;
        for (size_t i = 0; i < frame_num; ++i) {
            fout << i << endl;
            for (const auto &list:this->particle_matrix) {
                for (const auto &p:list) {
                    fout << p.get_position() << endl;
                }
            }
        }
    }

    const string &get_file_name() const {
        return file_name;
    }

    void set_file_name(const string &file_name) {
        FluidDatabase::file_name = file_name;
    }

    size_t get_particle_num() const {
        return particle_num;
    }

    void set_particle_num(size_t particle_num) {
        FluidDatabase::particle_num = particle_num;
    }

    size_t get_frame_num() const {
        return frame_num;
    }

    void set_frame_num(size_t frame_num) {
        FluidDatabase::frame_num = frame_num;
    }

    size_t get_frame_interval() const {
        return frame_interval;
    }

    void set_frame_interval(size_t frame_interval) {
        FluidDatabase::frame_interval = frame_interval;
    }

    size_t get_animation_duration() const {
        return animation_duration;
    }

    void set_animation_duration(size_t animation_duration) {
        FluidDatabase::animation_duration = animation_duration;
    }

    const vector<vector<Particle>> &get_particle_matrix() const {
        return particle_matrix;
    }

    void set_particle_matrix(const vector<vector<Particle>> &particle_matrix) {
        FluidDatabase::particle_matrix = particle_matrix;
    }
};


#endif //SPH_FLUID_SOLVER_FLUIDDATABASE_H