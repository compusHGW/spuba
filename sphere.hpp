#ifndef __SPHERE_H_
#define __SPHERE_H_


/* #####   HEADER FILE INCLUDES   ################################################### */

#include <math.h>
#include <stdio.h>
/*
 * =====================================================================================
 *        Class:  Sphere
 *  Description:  A spherical object based on Obj
 * =====================================================================================
 */
class Sphere {
public:

	/* ====================  LIFECYCLE     ======================================= */
	Sphere(double, const double *);                             /* constructor      */

	/* ====================  ACCESSORS     ======================================= */
	const double*		GetPosition();

	/* ====================  MUTATORS      ======================================= */
	bool 			Intersect( const double*, const double*, double*, double*, double*, double*);

	/* ====================  OPERATORS     ======================================= */


protected:
	/* ====================  DATA MEMBERS  ======================================= */

private:
	/* ====================  DATA MEMBERS  ======================================= */
	double	 		mRadius;
	double 			mPosition[3];

	/* ====================  MUTATORS      ======================================= */
	bool 			CalculateHitPoints(double*,double*,const double*,const double*, double*,double*);
	int  			CalcQuadricRoots(double&,double&,double&,double&,double&);

}; /* -----  end of class Sphere  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Sphere
 *      Method:  GetPosition
 *     Returns:  const double*
 *--------------------------------------------------------------------------------------
 */
inline const double*
Sphere::GetPosition (  )
{
	return mPosition;
}		/* -----  end of method Sphere::GetPosition  ----- */

#endif
