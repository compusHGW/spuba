#pragma once

 
#include <stdio.h>
#include <cmath>
#include "particle3d.hpp"
#include "const.h"
#include <random>

class ParticleSource {
public:
	void init(string source_number);
	Particle3D gen_ion3D(std::minstd_rand &rng);

	double GetTotIonFl() {
		return totFl;
	}

	void setType(string type);

	string getType();

	void setNumberOfAngleBins(int number_of_angle_bins);

	void setDeltaAngle(double delta_angle);

	double getDeltaAngle();

	void setDeltaTheta(double delta_angle);

	void setCurrentDensity(vector<double> &current_density_vector);

	void setEnergy(vector<double> &energy_vector);

	void setFraction(vector<double> &fraction_vector);

	void setPosition(vector<double> &position_vector);

	vector<double> getPosition();


	void setConstantPhi(double phi);

private:

	void calc_cum3D();
	int lind(double p);
	int lind3d(double p);
	void Construct_Pdf();
	void calc_cum();
	double GetFrac(double _th);
	void GetIonFluxDensE(double _th, double &_Fl, double &_E);

	int n_angle_bins;
	double DTH;
	const double DTH_3D = M_PI / 360.; // = 0.5 degree

	double totc;               // total current of particles
	double totFl;              // Total flux
	double totFld;             // Total flux  density

	vector<double> pt3d;
	vector<double> pfrac3d;
	vector<double> Ei3d;
	vector<double> mVi3d;

	vector<double> cum;
	vector<double> cum3d;
	vector<double> pdf;                // probability density function

	string type;
	double delta_angle_theta;
	vector<double> position;
	double constant_phi;

	vector<double> current_density;
	vector<double> particle_density;
	vector<double> pt;

	vector<double> mean_energy;
	vector<double> mean_velocity;

	vector<double> fraction_of_double_charged_ions; // fraction = X++ / X+

};


inline double surf_sph_ring(double t1, double t2) {
	return TWOPI * (cos(t1) - cos(t2));
};