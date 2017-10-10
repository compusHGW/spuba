#include "tube.hpp"
#include "const.h"

Tube::Tube(double rad,double* pos,double _length) {
	//printf("Creating Tube\n");
	mRadius = rad;
	mPosition[0] = pos[0];
	mPosition[1] = pos[1];
	mPosition[2] = pos[2];
	mLength = _length;
	//printf("done creating tube\n");
}

bool
Tube::Intersect(const double* source, const double* direction, double* hit1, double* hit2, double* pT1, double* pT2) {
	//printf("** Tube intersect:  ");
	double a = source[0]-mPosition[0]; // in corrdinate system of the cylinder
	double b = source[1]-mPosition[1];
	double c = source[2]-mPosition[2];
	double u = direction[0];
	double v = direction[1];
	double w = direction[2];
	double r = mRadius;

	double asq = a*a;
	double bsq = b*b;
	double csq = c*c;
	double usq = u*u;
	double vsq = v*v;
	double wsq = w*w;
	double rsq = r*r;

	double temp;
	temp =  -asq*usq*vsq + 2.0*a*b*usq*u*v - bsq*usq*usq + rsq*usq*usq + rsq*usq*vsq; // = Diskriminante/(u*u)
	if ( temp <= 0 ){ 
		//printf("temp < 0 \n");
		return false;
	}
	temp = sqrt (temp);
	double x1 = (-temp + a*vsq - b*u*v) / (usq+vsq);
	double x2 = (+temp + a*vsq - b*u*v) / (usq+vsq);

	double x1sq = x1*x1;
	double x2sq = x2*x2;

	/*  Calculate the multiplyer for ray vector to the hit points */
	double l1,l2;
	if (u > -EPSILON && u < EPSILON ) {
		l1  = ( + sqrt(rsq-x1sq) - b) / v;
		l2  = ( - sqrt(rsq-x1sq) - b) / v;
	} else {
		l1 = (x1 - a) / u;
		l2 = (x2 - a) / u;
	}

	if (l1 < 0.0 && l2 < 0.0) {
		//printf("l1 and l2 < 0\n");
		return false;
	}
	  
	if (l2 < l1) {
		double temp = l2;
		l2 = l1;
		l1 = temp;
	}
	
	if ( l1 < EPSILON && l1 > - EPSILON) { return false;}
	if ( l2 < EPSILON && l2 > - EPSILON) { return false;}
	/* Calculate distance from tube beginning */
	double m1 = l1 * w;
	double m2 = l2 * w;
//  if ( l1 < EPSILON ) { // means pos is in the tube
//  //	D("inside");
//    if (m2 > mPosition[2] - mLength && m2 < mPosition[2]) {		
//      //D("can see inner wall");
//      t1  = l2;
//      t2 = l1;
//    }else{
//      if (m1 < mPosition[2] - mLength || m1 > mPosition[2]) return false;
//      //D("cant see inner wall");
//      t1 = l1;
//      t2 = l2;
//    }

//  }else{ 		// means pos is outside 
//    //D("outside");
//    if (m1 > mPosition[2] - mLength && m1 < mPosition[2]){
//      //D("can see outer wall");
//      t1  = l1;
//      t2  = l2;
//    }else{
//      if (m2 < mPosition[2] - mLength || m2 > mPosition[2]) return false;		
//      //D("can see inner wall");
//      t1 = l2;
//      t2 = l1;
//    }
//  }
        /*if( (m1<mPosition[2]||m1>(mPosition[2]+mLength)) && (m2<mPosition[2]||m2>(mPosition[2]+mLength)) ) 
	{ 
	    printf("** Tube intersect: not inside the tube  m1=%f  m2=%f\n",m1,m2); 
	    return false;
	}
        if( (m1<0||m1>fabs(mLength)) && (m2<0||m2>fabs(mLength)) ) 
	{ 
	    printf("** Tube intersect: not inside the tube  m1=%f  m2=%f\n",m1,m2); 
	    return false;
	}*/

	/* Test if particle hit tube*/
	double z1 = source[2] + direction[2] * l1;
	double z2 = source[2] + direction[2] * l2;
        if( (z1<mPosition[2]||z1>(mPosition[2]+mLength)) && (z2<mPosition[2]||z2>(mPosition[2]+mLength)) ) 
	{ 
	    //printf("** Tube intersect: not inside the tube  z1=%f  z2=%f\n",m1,m2); 
	    return false;
	}
	double t1,t2;
	t1 = l1;
	t2 = l2;
	//printf("ray->t %f\n",ray->t);

	 hit1[0] = source[0] + direction[0] * t1;
	 hit1[1] = source[1] + direction[1] * t1;
	 hit1[2] = source[2] + direction[2] * t1;

	 hit2[0] = source[0] + direction[0] * t2;
	 hit2[1] = source[1] + direction[1] * t2;
	 hit2[2] = source[2] + direction[2] * t2;

	 *pT1 = m1/mLength;
	 *pT2 = m2/mLength;

	 //printf(" end -> \thit1=(%f, %f, %f) pT1=%f \n",hit1[0],hit1[1],hit1[2],*pT1);
	 //printf(" \t\t\t\thit2=(%f, %f, %F) pT2=%f  \n",hit2[0],hit2[1],hit2[2],*pT2);
	 return true;
}

