#include <vector>
#include <algorithm>
#include <array>
#include "shared.hpp"
#include "ParticleSource.h"
#include <fstream>
#include <sstream>

using namespace std;


#ifdef DEBUG
#define ASSERT(a) assert(a)
#else
#define ASSERT(a)
#endif

inline double get_pdf3d(double xim1, double xi, double fim1, double Fim1) {
    return 2. * Fim1 / (xi - xim1) - fim1;
}

void ParticleSource::calc_cum() {
    cum.push_back(0.0);
    partial_sum(pt.begin(), pt.end(), back_inserter(cum));
}

void ParticleSource::calc_cum3D() {
    cum3d.push_back(0.0);
    partial_sum(pt3d.begin(), pt3d.end(), back_inserter(cum3d));
}



int ParticleSource::lind3d(double p) {
    long r;
    return (r = lower_bound(cum3d.begin(), cum3d.end(), p) - cum3d.begin() - 1) > 0 ? r : 0;
}

int ParticleSource::lind(double p) {
    int r;
    return (r = lower_bound(cum.begin(), cum.end(), p) - cum.begin() - 1) > 0 ? r : 0;
}

void ParticleSource::Construct_Pdf() {
    vector<double> p;

    for (auto e : particle_density) {
        p.push_back(e / totFld);
    }

    pdf.push_back(0);
    double a = DTH * .5;
    pdf.push_back(get_pdf3d(0., a, 0., p[0]));
    for (int i = 2; i < pdf.size() - 1; ++i, a += DTH)
        pdf.push_back(get_pdf3d(a, a + DTH, pdf[i - 1], p[i - 1]));
    pdf.push_back(get_pdf3d(a, a + DTH * .5, pdf[pdf.size() - 2], p[pdf.size() - 2]));
}

double ParticleSource::GetFrac(double _th) {
    int i = _th / DTH;
    return i < particle_density.size() - 1 ? fraction_of_double_charged_ions[i] + (_th - DTH * i) * (fraction_of_double_charged_ions[i + 1] - fraction_of_double_charged_ions[i]) / DTH : fraction_of_double_charged_ions[i];
}

void ParticleSource::GetIonFluxDensE(double _th, double &_Fl, double &_E) {
    int i = _th / DTH;
    _Fl = i < particle_density.size() - 1 ? particle_density[i] + (_th - DTH * i) * (particle_density[i + 1] - particle_density[i]) / DTH : particle_density[i];
    _E = i < particle_density.size() - 1 ? mean_energy[i] + (_th - DTH * i) * (mean_energy[i + 1] - mean_energy[i]) / DTH : mean_energy[i];
}

void ParticleSource::init(string source_number) {
    std::string filename;

    filename = OUT_FOLDER + "/source_" + source_number + "_current.txt";
    ofstream ofc(filename);
    filename = OUT_FOLDER + "/source_" + source_number  + "_energy.txt";
    ofstream ofe(filename);
    filename  = OUT_FOLDER + "/source_" + source_number  + "_fraction.txt";
    ofstream off(filename);

    for (unsigned i = 0; i < n_angle_bins; ++i) {
        ofc << (double) i * delta_angle_theta << " " << current_density[i] << std::endl;
        ofe << (double) i * delta_angle_theta << " " << mean_energy[i] << std::endl;
        off << (double) i * delta_angle_theta << " " << fraction_of_double_charged_ions[i] << std::endl;
    }


    for (int i=0; i < n_angle_bins; i++){
        particle_density.push_back(current_density[i]);
        pt.push_back(current_density[i]);
    }


    //First we will obtain total current for each ring and also the half sphere
    double a = DTH * .5;

    pt[0] *= surf_sph_ring(0, a) / (1. + fraction_of_double_charged_ions[0]);
    particle_density[0] /= (1. + fraction_of_double_charged_ions[0]) * Qe;// Transform current density into particle density
    totc = pt[0];
    totFl = particle_density[0] * surf_sph_ring(0, a);
    totFld = particle_density[0];

    for (int i = 1; i < n_angle_bins - 1; ++i, a += DTH) {
        pt[i] *= surf_sph_ring(a, a + DTH) / (1. + fraction_of_double_charged_ions[i]);
        particle_density[i] /= (1. + fraction_of_double_charged_ions[i]) * Qe;
        totc += pt[i];
        totFl += surf_sph_ring(a, a + DTH) * particle_density[i];
        totFld += particle_density[i];
    }

    pt[n_angle_bins - 1] *= surf_sph_ring(a, a + DTH * .5) / (1. + fraction_of_double_charged_ions[n_angle_bins - 1]);
    particle_density[n_angle_bins - 1] /= (1. + fraction_of_double_charged_ions[n_angle_bins - 1]) * Qe;
    totc += pt[n_angle_bins - 1];
    totFl += surf_sph_ring(a, a + DTH * .5) * particle_density[n_angle_bins - 1];
    totFld += particle_density[n_angle_bins - 1];


    for (int i = 0; i < n_angle_bins; ++i) {
        pt[i] /= totc;
    }

    a = DTH_3D * .5;

    double ss = current_density[0] * surf_sph_ring(0, a) / (1. + fraction_of_double_charged_ions[0]);
    pt3d.push_back(ss);
    pfrac3d.push_back(fraction_of_double_charged_ions[0]);
    Ei3d.push_back(mean_energy[0]);
    mVi3d.push_back(mean_energy[0]);

    int nn = M_PI * .5 / DTH_3D + 1;
    double sc = DTH_3D / DTH;

    for (int i = 1; i < nn - 1; ++i, a += DTH_3D) {
        int j = (double) i * sc;
        double v = fraction_of_double_charged_ions[j] + ((double) i * sc - (double) j) * (fraction_of_double_charged_ions[j + 1] - fraction_of_double_charged_ions[j]);
        pfrac3d.push_back(v);
        v = current_density[j] + ((double) i * sc - (double) j) * (current_density[j + 1] - current_density[j]);
        pt3d.push_back(v * surf_sph_ring(a, a + DTH_3D) / (1. + pfrac3d[pfrac3d.size() - 1]));
        v = mean_energy[j] + ((double) i * sc - (double) j) * (mean_energy[j + 1] - mean_energy[j]);
        Ei3d.push_back(v);
        mVi3d.push_back(v);
        ss += pt3d[pt3d.size() - 1];
    }

    pfrac3d.push_back(fraction_of_double_charged_ions[fraction_of_double_charged_ions.size() - 1]);
    pt3d.push_back(current_density[current_density.size() - 1] * surf_sph_ring(a, a + DTH_3D * .5) / (1. + pfrac3d[pfrac3d.size() - 1]));
    Ei3d.push_back(mean_energy[mean_energy.size() - 1]);
    mVi3d.push_back(mean_energy[mean_energy.size() - 1]);
    ss += pt3d[pt3d.size() - 1];

    for (auto &e : pt3d) {
        e /= ss;
    }

    double ini = 0.;
    ini = 0.;
    ASSERT(abs(accumulate(pt.begin(), pt.end(), ini) - 1.) < EPSILON);
    // Now calculate cumulative probabilities
    calc_cum();
    calc_cum3D();
    ASSERT(abs(cum3d[cum3d.size() - 1] - 1.) < EPSILON);
    Construct_Pdf();

    // Init velocities
    for (int i = 0; i < mean_velocity.size(); ++i) {
        // V_Xep = sqrt(2*E_pot[mesuared]*e/[atom unit mass]/m_Xe[in aum]/(1+fraction_of_double_charged_ions))
        // for V_Xepp it have to be multiplyed by sqrt(2.)
        mean_velocity[i] = sqrt(2. * mean_energy[i] * 96485333.6 / 131 / (1. + fraction_of_double_charged_ions[i]));
    }
    // Init velocities 3d
    for (int i = 0; i < Ei3d.size(); ++i){
        mVi3d[i] = sqrt(2. * Ei3d[i] * 96485333.6 / 131 / (1. + pfrac3d[i]));
    }


}


Particle3D ParticleSource::gen_ion3D(std::minstd_rand& rng) {
    std::uniform_real_distribution<> dis(0,1);
    auto RNG = [&](){
        auto ret = dis(rng);
        return ret;
    };

    double r, theta, phi, E;
    double f = 1.0;
    int i;

    r = RNG();
    i = lind3d(r);
    theta = fabs(DTH_3D * (i - .5) + DTH_3D * (r - cum3d[i]) / (cum3d[i + 1] - cum3d[i]));
    if (theta > 0.5 * M_PI) theta = M_PI - theta;


    if ( isnan(constant_phi) ) { // is NaN comparison
        phi = RNG() * 2 * M_PI; // particles scattered on rings
    } else {
        phi = constant_phi * M_PI / 180.; //for a beam!!
    }


    Particle3D p(&getPosition()[0], theta, phi, 0);

    r = RNG();
    p.SetType(IONS);
    if (r < pfrac3d[i]) {
        p.SetType(TWOPLUS_IONS);
        f = sqrt(2);
    }

    if (p.GetType() == IONS) E = Ei3d[i] / (1. + pfrac3d[i]);
    else E = 2. * Ei3d[i] / (1. + pfrac3d[i]);

    p.SetVelocity(f * mVi3d[i ? --i : i]);
    p.SetEnergy(E);

    return p;
}

void ParticleSource::setType(string type) {
    ParticleSource::type = type;
}

string ParticleSource::getType(){
    return type;
}


void ParticleSource::setNumberOfAngleBins(int number_of_angle_bins) {
    ParticleSource::n_angle_bins = number_of_angle_bins;
}


void ParticleSource::setDeltaAngle(double delta_angle) {
    ParticleSource::delta_angle_theta = delta_angle;

}

double ParticleSource::getDeltaAngle(){
    return delta_angle_theta;
}

void ParticleSource::setDeltaTheta(double delta_angle){
    ParticleSource::DTH = delta_angle * M_PI / 180;
}

void ParticleSource::setCurrentDensity(vector<double> &current_density_vector) {
    ParticleSource::current_density = current_density_vector;

}

void ParticleSource::setEnergy(vector<double> &energy_vector) {
    ParticleSource::mean_energy = energy_vector;
}

void ParticleSource::setFraction(vector<double> &fraction_vector) {
    ParticleSource::fraction_of_double_charged_ions = fraction_vector;

}

void ParticleSource::setPosition(vector<double> &position_vector) {
    ParticleSource::position = position_vector;

}

vector<double> ParticleSource::getPosition() {
    return position;
}

void ParticleSource::setConstantPhi(double phi) {
    ParticleSource::constant_phi = phi;
}





