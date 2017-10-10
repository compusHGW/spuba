#ifndef PROJECT_SIMULATION_H
#define PROJECT_SIMULATION_H

#include <string>
#include "intersection.h"
#include "Thruster.hpp"
#include "Chamber.hpp"
#include <chrono>
#include <random>


class Simulation {

public:

    Simulation() :
            dis(0,1)
    {


    }
    const string &getName() const {
        return name;
    }

    void setName(const string &name) {
        Simulation::name = name;
    }

    vector<Thruster> &getThrusters() {
        return thrusters;
    }

    Chamber &getChamber() {
        return chamber;
    }

    void setChamber(const Chamber &chamber) {
        Simulation::chamber = chamber;
    }

    //void read_input_from_file( std::string filename = "input.txt" );
    void read_input_file(std::string filename);
    void run();
    void output();

private:

    void calculate_reflection(Particle3D& p, Particle3D& reflect);
    string name;
    //int job_id;

    Chamber chamber;
    std::vector<Thruster> thrusters;

    // input.json
    int nParticles;
    int job_id;
    string comment;
    string wall;
    string feedgas;
    bool baffles_flag;
    string source_type;

    // high precision timing
    std::chrono::high_resolution_clock::time_point begin_time;

    // random number generator
    std::minstd_rand rng;
    std::uniform_real_distribution<> dis;


};


#endif //PROJECT_SIMULATION_H
