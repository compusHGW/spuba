#include "sphere.hpp"

double Dot( const double* a, const double* b ){
	return a[0]* b[0] + a[1] * b[1] + a[2] * b[2];
}

// Code for sphere calculation from wiki.delphigl.com/index.php/Tutorial_Raytracing_-_Grundlagen_I

bool
Sphere::CalculateHitPoints(double* hit1,double* hit2, const double* source, const double* direction, double* pT1, double* pT2) {
	// Take the calculation for spheres
	double boundingSquare = mRadius * mRadius;

	double tempRayPos[3];

	tempRayPos[0] = source[0] - mPosition[0];
	tempRayPos[1] = source[1] - mPosition[1];
	tempRayPos[2] = source[2] - mPosition[2];

	double a = 0.0,b = 0.0 ,c = 0.0;

	//printf("s %f %f %f\n",source[0],source[1],source[2]);
	//printf("d %f %f %f\n",direction[0],direction[1],direction[2]);

	a = Dot(direction,direction);

	b = 2 * Dot(tempRayPos,direction);

	c = Dot(tempRayPos,tempRayPos) - boundingSquare;

	double t1,t2;

	int roots = CalcQuadricRoots(a,b,c,t1,t2);
	//printf("roots %d\n",roots);
	if (roots > 0) {
//    printf("roots >0!!! t1 %f t2 %f\n",t1,t2);
		//if (t1 < 0 || t2 < 0) {
		//	return false;
		//}
		double temp;
		// sort by size
		if (t1 > t2) {
			temp = t2;
			t2 = t1;
			t1 = temp;
		}
		hit1[0] = t1 * direction[0] + source[0];
		hit1[1] = t1 * direction[1] + source[1];
		hit1[2] = t1 * direction[2] + source[2];

		hit2[0] = t2 * direction[0] + source[0];
		hit2[1] = t2 * direction[1] + source[1];
		hit2[2] = t2 * direction[2] + source[2];

		*pT1 = t1;
		*pT2 = t2;

		return true;
	} else {
		return false;
	}
}

int
Sphere::CalcQuadricRoots(double& a, double& b, double& c , double& t1, double& t2) {
	double determinant = b * b - 4 * a * c;
	//printf("det %f\n",determinant);
	if (determinant < 0) {
		t1 = 0.0;
		t2 = 0.0;
		return 0;
	}
	determinant = sqrt(determinant);

	double sign;
	if (b < 0) {
		sign = -1;
	} else {
		sign = +1;
	}
	double q = -0.5 * (b + sign * determinant);

	t1 = q / a;
	t2 = c / q;

	if (t1 > t2) {
		q = t2;
		t2 = t1;
		t1 = q;
	}
	return t1 == t2 ? 1 :2;
}

bool
Sphere::Intersect(const double* source, const double* direction, double* hit1, double* hit2, double* pT1, double* pT2) {
	return CalculateHitPoints(hit1,hit2,source,direction,pT1,pT2);
	//printf("%d\n",hit);
//  if (hit == false) return false;
}

Sphere::Sphere(double _radius, const double *_position) {
	mRadius = _radius;
	mPosition[0] = _position[0];
	mPosition[1] = _position[1];
	mPosition[2] = _position[2];
}

