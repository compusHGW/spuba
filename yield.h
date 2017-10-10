#ifndef YIELD_INC
#define YIELD_INC
#include <vector>
#include <string>
#include <algorithm>


using namespace std;

class Yield
{  public:
    Yield(){};

    void   setYield( vector<double> &_yield);
    void   setEnergy( vector<double> &_energy);

    double  interpolate_yield(double particle_energy);

  private:
    vector<double> energy;
    vector<double> yield;

};

#endif