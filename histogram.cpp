#include <vector>
#include "circle.hpp"
#include <iostream>
#include "shared.hpp"

using namespace std;

class HistogramCircle : public R3::Circle {
public:
    HistogramCircle(double* pos , double* dir, double radius ):
            R3::Circle(pos,dir,radius)
    {

    }

    void fillHistogram( const Particle3D& reflect, double* hit, double weight = 1.0 ){

        double point_exit[3];
        double phi_exit, r_exit;

        point_exit[0]=hit[0]-mPosition[0];
        point_exit[1]=hit[1]-mPosition[1];
        point_exit[2]=hit[2];

        // phi in x-y-plane from 0 to 2*M_PI (0 to 360)
        phi_exit=atan2(point_exit[1],point_exit[0]);
        if (point_exit[1]<0) phi_exit+=2*M_PI;

        r_exit=sqrt(point_exit[0]*point_exit[0]+point_exit[1]*point_exit[1]);

        auto r_addr = r_exit / bin_width_r;
        auto phi_addr = phi_exit / bin_width_phi;

        if ( reflect.GetType() == TWOPLUS_IONS ) {
            bins[r_addr][phi_addr].xe2pp += weight;
        }else{
            bins[r_addr][phi_addr].xe2p += weight;
        }
        bins[r_addr][phi_addr].n++;

    }

private:

    double bin_width_r;
    double bin_width_phi;

    struct Bin {
        double xe2pp;
        double xe2p;
        int n;
    };
    /// @brief bins in r and phi
    std::vector<std::vector<Bin>> bins;

};