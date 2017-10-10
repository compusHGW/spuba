#pragma once

#include "bins.h"

void outputAbsorbtion(std::string filename, double circleCounter, int nParticles, int nBins, Bin absorbtionBins[]);

void outputEmission(std::string filename, int nParticles, int nBins, Bin emissionBins[]);

void outputAbsorbtion2DChan(std::string filename, double circleCounter, int nParticles, int nBinsphi, int nBinsz, Bin2D absorbtionBins[], double r);

void outputAbsorbtion2DBot(std::string filename, double circleCounter, int nParticles, int nBinsphi, int nBinsz, Bin2D absorbtionBins[]);

