#include "yield.h"

void Yield::setYield(vector<double> &_yield) {
    Yield::yield = _yield;
}

void Yield::setEnergy(vector<double> &_energy) {
    Yield::energy = _energy;
}

double Yield::interpolate_yield(double particle_energy){
    particle_energy;
    if ( particle_energy > energy.back() ) return yield.back();
    if ( particle_energy < energy.front() ) return yield.front();


    std::vector<double>::iterator low = lower_bound(energy.begin(),energy.end(), particle_energy);
    long int index = low - energy.begin();
    double interpolated_yield = yield[index]+(particle_energy-energy[index])*(yield[index+1]-yield[index])/(energy[index+1]-energy[index]);

    return interpolated_yield;
}