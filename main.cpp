#include "Simulation.h"
#include "input_json.h"
#include <thread>

using namespace std;

int main(int argc, char** argv) {

    std::vector<Simulation*> simulations;
    std::vector<std::thread> threads;
    auto hw_concurrency = 1;//std::thread::hardware_concurrency();
    for (int i = 0; i < hw_concurrency; ++i) {
        threads.emplace_back( [&]() {
            auto simulation = read_simulation_file("../simulation.json");
            simulations.push_back( new Simulation(simulation) );
            {
                auto& simulation = *simulations.back();
                simulation.read_input_file("../input.json");
                simulation.run();
            }
        });
    }

    for( auto& t : threads ){
        t.join();
    }

    for( auto& simulation : simulations ){
        simulation->output();
    }
}
