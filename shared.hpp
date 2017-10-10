#ifndef PROJECT_SHARED_HPP
#define PROJECT_SHARED_HPP

#include "bins.h"
#include "intersection.h"
#include <vector>

#define OUT_FOLDER "./out"s

enum {IONS=0,TWOPLUS_IONS,DIED};

constexpr int NUM_EM_SP_SEGS = 20;
constexpr int NUM_CYL_SEGS = 500;
constexpr int n = 100;

const double Rc = 1.2;
const double Rs = 1.825;

const double OrigRs = 0.275;

constexpr int number_of_cones_in_z = 13;
constexpr int number_of_cones_in_r = 11;

constexpr int nBins2Dphi = 360;
constexpr int nBins2Dz = 10;
constexpr int nBins2Dr = 10;
constexpr int nBins = 91;


struct T {
    T():
            l_v_nPLandCylSeg(NUM_CYL_SEGS, 0),
            l_v_nReDepCylSeg(NUM_CYL_SEGS, 0),
            l_v_nPLandSpSeg(NUM_EM_SP_SEGS, 0),
            l_v_nReDepSpSeg(NUM_EM_SP_SEGS, 0)
    {}
    std::vector<double> l_v_nReDepSpSeg;
    std::vector<double> l_v_nReDepCylSeg;
    std::vector<double> l_v_nPLandCylSeg;
    std::vector<double> l_v_nPLandSpSeg;
};

struct Counter {
    double circleCounter = 0.0;
    double channelCounter = 0.0;
    double exitCounter = 0.0;
    double bottomCounter = 0.0;
};

struct Bins {
    struct {
        Bin emissionBins[nBins];
        Bin absorbtionBins[nBins];
    } d1;
    struct {
        struct {
            Bin2D chan[nBins2Dphi * nBins2Dz];
            Bin2D bot[nBins2Dphi * nBins2Dr];
            Bin2D exit[nBins2Dphi * nBins2Dr];
        } absorbtion;
    }d2;
};

#endif //PROJECT_SHARED_HPP
