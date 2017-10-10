#include "bins.h"
#include "shared.hpp"

double _doti ( const double* a , const double* b )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Anglei ( const double* a, const double* b ) 
{
  return acos(_doti( a  , b ) / sqrt( _doti ( a, a ) ) * sqrt( _doti ( b,b )));
}


void construct4Bins1D(Bin *emissionBins, Bin *absorbtionBins, Bin *absorbtionBins1, Bin *absorbtionBins2,
					  Bin *absorbtionBins3, int nBins, double dth)
{

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

// additional thruster (>=1) absorptionBins 

	for (int i = 0; i < nBins; i++) {
		absorbtionBins1[i].start = absorbtionBins2[i].start = absorbtionBins3[i].start = absorbtionBins[i].start;
		absorbtionBins1[i].end = absorbtionBins2[i].end = absorbtionBins3[i].end = absorbtionBins[i].end;

	}

}


void constructBins2D(Bin2D *absorbtionBins2DExit, Bin2D *absorbtionBins2DChan, Bin2D *absorbtionBins2DBot,
					 double nBins2Dphi, double dphi, int nBins2Dz, double dz, int nBins2Dr, double dr)
{

	int k,n;

	for (int i=0; i<nBins2Dphi; ++i)
	{
		for (int j=0; j<nBins2Dr; ++j)
		{ // exit, bottom
			k=i*nBins2Dr+j;
			// r direction
			absorbtionBins2DExit[k].startJ=absorbtionBins2DBot[k].startJ=j*dr;
			absorbtionBins2DExit[k].endJ=absorbtionBins2DBot[k].endJ=(j+1)*dr;
			// phi direction
			absorbtionBins2DExit[k].startI=absorbtionBins2DBot[k].startI=(i)*2*M_PI*dphi/360;
			absorbtionBins2DExit[k].endI=absorbtionBins2DBot[k].endI=(i+1)*2*M_PI*dphi/360;

		}

		for (int m=0; m<nBins2Dz; ++m)
		{ // channel
			n=i*nBins2Dz+m;
			// z direction
			absorbtionBins2DChan[n].startJ=m*dz;
			absorbtionBins2DChan[n].endJ=(m+1)*dz;
			// phi direction
			absorbtionBins2DChan[n].startI=(i)*2*M_PI*dphi/360;
			absorbtionBins2DChan[n].endI=(i+1)*2*M_PI*dphi/360;
		}
	}

}

void fillEmissionBins1D(int nBins, Bin *emissionBins, Particle3D p)
{

	double theta = p.GetDirectionPolar()[0];
    	theta *= 180.0/M_PI;
    	for (int j = 0; j < nBins; j++) {
      		if ( theta >= emissionBins[j].start && theta < emissionBins[j].end ) {
        		emissionBins[j].nParticles++;
        		if ( p.GetType() == TWOPLUS_IONS ) {
          			emissionBins[j].xe2pp++;
        		}else{
          			emissionBins[j].xe2p++;
        		}
        		break;
      		}
    	}

}

void fillAbsorbBins1D(int nBins, Bin *absorbtionBins, double *circleCounter, double weight, Particle3D reflect)
{

	double normalcircle[3]={ 0.0,0.0,-1.0 };
	double vect[3];
	// check conversion from reflect to vect
	for (int i = 0; i < 3; i++) {
		vect[i] = reflect.GetDirection()[i];
	}

	double angle = Anglei(normalcircle,vect);

	*circleCounter+=1;
	double theta = angle;
	theta *= 180.0/M_PI;
	for (int j = 0; j < nBins; j++)
	{
		if ( theta >= absorbtionBins[j].start && theta < absorbtionBins[j].end ) {
			absorbtionBins[j].nParticles+=weight;
			if ( reflect.GetType() == TWOPLUS_IONS ) {
				absorbtionBins[j].xe2pp+=weight;
			}else{
				absorbtionBins[j].xe2p+=weight;
			}
			break;
		}
	}
}

void fillAbsorbBins2D(int nBins2Dphi, int nBins2Dz, int nBins2Dr, Bin2D *absorbtionBins2DExit, Bin2D *absorbtionBins2DChan,
					  Bin2D *absorbtionBins2DBot, double *exitCounter, double *channelCounter, double *bottomCounter,
					  double weight, Particle3D reflect, const double *ThrBottPos, double *hit)
{
	int k;

	if (reflect.GetENumS()>-1 || reflect.GetENum()>-1)
	{
		//
		// THRUSTER EXIT - POSITION OF PARTICLE ENTRANCE
		//

		double point_exit[3];
		double phi_exit, r_exit;

		point_exit[0]=hit[0]-ThrBottPos[0];
		point_exit[1]=hit[1]-ThrBottPos[1];
		point_exit[2]=hit[2];

		// phi in x-y-plane from 0 to 2*M_PI (0 to 360)
		phi_exit=atan2(point_exit[1],point_exit[0]);
		if (point_exit[1]<0) phi_exit+=2*M_PI;

		r_exit=sqrt(point_exit[0]*point_exit[0]+point_exit[1]*point_exit[1]);
		//cout << phi_exit << " " << r_exit << endl;
		for (int i = 0; i < nBins2Dphi; i++)
		{
			for (int j = 0; j < nBins2Dr; j++)
			{	k= i*nBins2Dr+j;

				// find bin for local coordinates
				if(phi_exit >= absorbtionBins2DExit[k].startI && phi_exit <= absorbtionBins2DExit[k].endI
				   && r_exit >= absorbtionBins2DExit[k].startJ && r_exit <= absorbtionBins2DExit[k].endJ)
				{
					if ( reflect.GetType() == TWOPLUS_IONS ) {
						absorbtionBins2DExit[k].xe2pp+=weight;
						absorbtionBins2DExit[k].nParticles+=1;
						*exitCounter+=1;
					}else{
						absorbtionBins2DExit[k].xe2p+=weight;
						absorbtionBins2DExit[k].nParticles+=1;
						*exitCounter+=1;
					}
					break;
				}

			}
		}
	}

	//
	// THRUSTER BOTTOM and CHANNEL
	//
	double point[3];
	double phi_local, r_local, z_local;

	for (int i = 0; i < 3; ++i) {
		// use thrusterbottom position to get local values
		point[i]=reflect.GetHitpoint()[i]-ThrBottPos[i];
	}

	// phi in x-y-plane from 0 to 2*M_PI (0 to 360)
	phi_local=atan2(point[1],point[0]);
	if (point[1]<0) phi_local+=2*M_PI;

	if (reflect.GetENumS() > -1 && reflect.GetENum() == -1){

		//
		// BOTTOM: ignore z, create r from x,y
		//
		r_local=sqrt(point[0]*point[0]+point[1]*point[1]);

		for (int i = 0; i < nBins2Dphi; i++)
		{
			for (int j = 0; j < nBins2Dr; j++)
			{	k= i*nBins2Dr+j;

				// find bin for local coordinates
				if(phi_local >= absorbtionBins2DBot[k].startI && phi_local <= absorbtionBins2DBot[k].endI
				   && r_local >= absorbtionBins2DBot[k].startJ && r_local <= absorbtionBins2DBot[k].endJ)
				{
					if ( reflect.GetType() == TWOPLUS_IONS ) {
						absorbtionBins2DBot[k].xe2pp+=weight;
						absorbtionBins2DBot[k].nParticles+=1;
						*bottomCounter+=weight;
					}else{
						absorbtionBins2DBot[k].xe2p+=weight;
						absorbtionBins2DBot[k].nParticles+=1;
						*bottomCounter+=weight;
					}
					break;
				}

			}
		}

	} else if((reflect.GetENumS() == -1 && reflect.GetENum() > -1) ) {

		//
		// CHANNEL: z remains
		//
		z_local=point[2];

//	  	cout << "phi: " << phi_local << ",z: " << z_local << endl;
		for (int i = 0; i < nBins2Dphi; i++)
		{
			for (int j = 0; j < nBins2Dz; j++)
			{	k= i*nBins2Dz+j;

				// find bin for local coordinates
				if(phi_local >= absorbtionBins2DChan[k].startI && phi_local <= absorbtionBins2DChan[k].endI
				   && z_local >= absorbtionBins2DChan[k].startJ && z_local <= absorbtionBins2DChan[k].endJ)
				{
					if ( reflect.GetType() == TWOPLUS_IONS ) {
						absorbtionBins2DChan[k].xe2pp+=weight;
						absorbtionBins2DChan[k].nParticles+=1;
						*channelCounter+=weight;
					}else{
						absorbtionBins2DChan[k].xe2p+=weight;
						absorbtionBins2DChan[k].nParticles+=1;
						*channelCounter+=weight;
					}
					break;
				}

			}
		}

	} else {

		double point_exit[3];
		point_exit[0]=hit[0]-ThrBottPos[0];
		point_exit[1]=hit[1]-ThrBottPos[1];
		point_exit[2]=hit[2];
		cout << point_exit[0] << " " << point_exit[1] << " " << point_exit[2] << " " << endl;
		cout << hit[0] << " " << hit[1] << " " << hit[2] << " " << endl;
		cout << "particle vanished inside thruster" << endl;
		exit(-1);
	}
}	
 
