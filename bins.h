#pragma once

#include "particle3d.hpp"

struct Bin{
    double start;
    double end;
    double nParticles;
    double xe2p;
    double xe2pp;
    double current;
    double current2;
    void operator+=( Bin& b ){
        this->nParticles += b.nParticles;
        this->xe2p += b.xe2p;
        this->xe2pp += b.xe2pp;
    }
};


struct Bin2D{ 	double startI,startJ;
    double endI,endJ;
    double nParticles;
    double xe2p;
    double xe2pp;
    double current;
    double current2;
    void operator+=( Bin2D& b ){
        this->nParticles += b.nParticles;
        this->xe2p += b.xe2p;
        this->xe2pp += b.xe2pp;
    }
};

void construct4Bins1D(Bin *emissionBins, Bin *absorbtionBins, Bin *absorbtionBins1, Bin *absorbtionBins2,
                      Bin *absorbtionBins3, int nBins, double dth);

void constructBins2D(Bin2D *absorbtionBins2DExit, Bin2D *absorbtionBins2DChan, Bin2D *absorbtionBins2DBot,
                     double nBins2Dphi, double dphi, int nBins2Dz, double dz, int nBins2Dr, double dr);

void fillEmissionBins1D(int nBins, Bin *emissionBins, Particle3D p);

void fillAbsorbBins1D(int nBins, Bin *absorbtionBins, double *circleCounter, double weight, Particle3D reflect);

void fillAbsorbBins2D(int nBins2Dphi, int nBins2Dz, int nBins2Dr, Bin2D *absorbtionBins2DExit, Bin2D *absorbtionBins2DChan,
                      Bin2D *absorbtionBins2DBot, double *exitCounter, double *channelCounter, double *bottomCounter,
                      double weight, Particle3D reflect, const double *ThrBottPos, double *hit);
	   
