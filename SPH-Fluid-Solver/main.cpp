#include <iostream>
#include <vector>
#include <random>
#include "FluidDatabase.h"

using namespace std;

int main() {
    std::cout << "Hello, World!" << std::endl;
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
    return 0;
}