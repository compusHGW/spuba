#ifndef PROJECT_THRUSTER_HPP
#define PROJECT_THRUSTER_HPP

#include <string>
#include "intersection.h"
#include "Geometry.h"
#include "bins.h"
#include "shared.hpp"
#include "intersect_objects.h"
#include "ParticleSource.h"

class Thruster {

    enum FEEDGAS_TYPE{
        XE,
        AR,
        KR
    };
public:
    Thruster(){
        //source.init("0");
        init_all_bins();
    }

    void init_all_bins(){
        init_bins();
        init_bins_2d();
    }

    const string &getName() const {
        return name;
    }

    void setName(const string &name) {
        Thruster::name = name;
    }

    FEEDGAS_TYPE getFeedgas() const {
        return feedgas;
    }

    void setFeedgas(FEEDGAS_TYPE feedgas) {
        Thruster::feedgas = feedgas;
    }

    const string &getEmmision_file() const {
        return emmision_file;
    }

    void setEmmision_file(const string &emmision_file) {
        Thruster::emmision_file = emmision_file;
    }

    Geometry &getGeometry() {
        return geometry;
    }

    friend auto& operator << ( ostream& o, Thruster& t );

    bool Intersect(Particle3D& reflect, double weight, double dthES){

        double hit[3] = {-1.0, -1.0, -1.0}; // some random values
        double t;
        // TODO for the time beeing this is ok since we don't support something else but
        // TODO cylindrical thrusters
        auto& cylinder = this->getGeometry().getGeometry_objects().front().object.cyl;
        auto& exit_circle = cylinder->getTop();
        // TODO potentially wrong
        //      you need to intersect with all circles
        //      and based on t you have to decide which
        //      was the nearest. then go an execute the
        //      if statement body for the nearest

        if (exit_circle->Intersect(reflect.GetPosition(), reflect.GetDirection(), hit, &t))  // thruster exit itself
        {

            fillAbsorbBins1D(nBins,
                             bins.d1.absorbtionBins,
                             &counter.circleCounter,
                             weight,
                             reflect
            );

            if (reflect.GetDirection()[2] > 0)
                printf("direction.z=%f  > 0 !!\n", reflect.GetDirection()[2]);

            //intersction with thruster channel (cylinder+cap)
            intersectThru(cylinder->getBottom().get(),
                          cylinder->getWall().get(),
                          &reflect, dthES
            );

            fillAbsorbBins2D(nBins2Dphi, nBins2Dz, nBins2Dr,
                             bins.d2.absorbtion.exit,
                             bins.d2.absorbtion.chan,
                             bins.d2.absorbtion.bot,
                             &counter.exitCounter,
                             &counter.channelCounter,
                             &counter.bottomCounter,
                             weight,
                             reflect,
                             cylinder->getBottom()->getPosition(),
                             hit
            );
            return true;
        }
        return false;
    }

private:

    void init_bins(){

        auto& emissionBins = bins.d1.emissionBins;
        auto& absorbtionBins = bins.d1.absorbtionBins;
        double dth = 90. / (nBins - 1);
        absorbtionBins[0].start = emissionBins[0].start = 0 ;
        absorbtionBins[0].end = emissionBins[0].end = 0.5*dth;

        for (int i = 1; i < nBins-1; i++) {
            absorbtionBins[i].start = emissionBins[i].start = emissionBins[i-1].end ;
            absorbtionBins[i].end = emissionBins[i].end = emissionBins[i].start + dth;
            if ( emissionBins[i].start < 0 ) emissionBins[i].start = 0 ;
            if ( absorbtionBins[i].start < 0 ) absorbtionBins[i].start = 0 ;
        }
        absorbtionBins[nBins-1].start = emissionBins[nBins-1].start = emissionBins[nBins-2].end;
        absorbtionBins[nBins-1].end = emissionBins[nBins-1].end = emissionBins[nBins-1].start+0.5*dth;
    }

    void init_bins_2d(){

        if ( this->getGeometry().getGeometry_objects().empty() ) return;
        auto TLength = this->getGeometry().getGeometry_objects().front().object.cyl->GetLength();
        auto Rthr = this->getGeometry().getGeometry_objects().front().object.cyl->GetRadius();
        double dphi = 360. / nBins2Dphi;
        double dz = TLength / nBins2Dz;
        double dr = Rthr / nBins2Dr;

        auto& absorbtionBins2DExit = bins.d2.absorbtion.exit;
        auto& absorbtionBins2DBot = bins.d2.absorbtion.bot;
        auto& absorbtionBins2DChan = bins.d2.absorbtion.chan;

        for (int i=0; i<nBins2Dphi; ++i)
        {
            for (int j=0; j<nBins2Dr; ++j)
            { // exit, bottom
                int k=i*nBins2Dr+j;
                // r direction
                absorbtionBins2DExit[k].startJ=absorbtionBins2DBot[k].startJ=j*dr;
                absorbtionBins2DExit[k].endJ=absorbtionBins2DBot[k].endJ=(j+1)*dr;
                // phi direction
                absorbtionBins2DExit[k].startI=absorbtionBins2DBot[k].startI=(i)*2*M_PI*dphi/360;
                absorbtionBins2DExit[k].endI=absorbtionBins2DBot[k].endI=(i+1)*2*M_PI*dphi/360;

            }

            for (int m=0; m<nBins2Dz; ++m)
            { // channel
                int n=i*nBins2Dz+m;
                // z direction
                absorbtionBins2DChan[n].startJ=m*dz;
                absorbtionBins2DChan[n].endJ=(m+1)*dz;
                // phi direction
                absorbtionBins2DChan[n].startI=(i)*2*M_PI*dphi/360;
                absorbtionBins2DChan[n].endI=(i+1)*2*M_PI*dphi/360;
            }
        }

    }

    std::string name;
    FEEDGAS_TYPE feedgas;
    std::string emmision_file;
    Geometry geometry;

public:
    ParticleSource source;
    Bins bins;
    Counter counter;
};

inline auto& operator << ( ostream& o, Thruster& t ){
    o << "\"name\" : " << t.name << endl;
    o << "\"feedgas\" : " << t.feedgas << endl;
    o << "\"emmision_file\" : " << t.emmision_file << endl;
    o << "\"geometry\" : [ " << endl;
    o << t.geometry << endl;
    o << "]" << endl;
    return o;
}


#endif //PROJECT_THRUSTER_HPP
