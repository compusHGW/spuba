#ifndef PROJECT_INPUT_JSON_H
#define PROJECT_INPUT_JSON_H

#include <string>
#include "Simulation.h"
#include "yield.h"


Simulation read_simulation_file(std::string filename);
ParticleSource read_source_file( std::string filename = "source.json", string type = "beam");
Yield read_yield(string filename, string wall_type, string feedgas_type);

#endif //PROJECT_INPUT_JSON_H
