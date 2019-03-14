#include <iostream>
#include <vector>
#include <random>
#include "FluidSolver.h"

using namespace std;
#define vec3 Vector3

void test_class_RestrictionBox_Particle();

void test_RestrictionBox();

void test_fluid_solver();

void test_iterator();

void test_initialize_particles();
int main() {
    std::cout << "Hello, World!" << std::endl;
    test_initialize_particles();
//    cout  << pow(8,1/3.0)<< endl;
    // Test
//    test_iterator();
    return 0;
}

void test_iterator() {
    vector<Particle> tmp_vector;
    tmp_vector.reserve(10);
    for (int i = 0; i < 10; ++i) {
        tmp_vector.emplace_back(Particle(static_cast<size_t>(i)));
    }
    cout << tmp_vector.size()<< endl;
    for (auto &p :tmp_vector) {
        cout << p.get_index()<< endl;
    }
    for(auto&p:tmp_vector){
        p.set_index(100);
    }
    for (auto &p :tmp_vector) {
        cout << p.get_index()<< endl;
    }

}

void test_initialize_particles() {
    const FluidParameter tmp_paras(1000, 0.01, 1, 1, 1, 1, 1);
    const RestrictionBox tmp_box(vec3(0, 0, 0), vec3(100, 100, 100));
    FluidDatabase tmp_database("test0", tmp_paras.get_particle_num(), 60, 0.01);
    FluidSolver tmp_solver(tmp_paras, tmp_box, tmp_database);

    tmp_solver.initialize_particles(vec3(100, 100, 100), 100.0);
    for (auto &p :tmp_solver.get_realtime_particle_list()) {
        cout << p.get_index() <<" " << p.get_position() << endl;
    }
}
void test_fluid_solver() {
    const FluidParameter tmp_paras(1000, 0.01, 1, 1, 1, 1, 1);
    const RestrictionBox tmp_box(vec3(0, 0, 0), vec3(100, 100, 100));
    FluidDatabase tmp_database("test0", tmp_paras.get_particle_num(), 60, 0.01);
    FluidSolver tmp_solver(tmp_paras, tmp_box, tmp_database);
    tmp_solver.test_kernel_functions();
}

void test_RestrictionBox() {
    auto tmp_box = RestrictionBox(vec3(10, 10, 10), vec3(100, 100, 100));
    auto tmp_particle = Particle(1);
    tmp_particle.set_position(vec3(60, 0, 0));
    tmp_particle.set_velocity(vec3(0, 2, -2));
    tmp_box.restrict_particle(tmp_particle);
    //cout << tmp_particle.get_velocity()<<endl;
    cout << tmp_box.reflect_vector(vec3(2, 2, 2), vec3(0, -1, 0));

}

void test_class_RestrictionBox_Particle() {
    // Test RestrictionBox and Particle
    auto fluidDatabase = FluidDatabase("test3.fb", 100, 100, 0.1);
    for (size_t i = 0; i < 100; ++i) {
        vector<Particle> tmp_list;
        for (size_t j = 0; j < 100; ++j) {
            random_device rd;
            uniform_real_distribution<double> dist(0, 100);

            auto tmp_particle = Particle(j);
            tmp_particle.set_position(Vector3(dist(rd), dist(rd), dist(rd)));
            tmp_list.emplace_back(tmp_particle);
        }
        fluidDatabase.append_particle_list(i, tmp_list);
    }
    fluidDatabase.export_database();
}
