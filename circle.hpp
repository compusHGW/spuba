#include "const.h"
#include <stdio.h>
#include <string>
#include <cmath>
#include <ostream>

#ifndef CIRCLE_LDO6QD4A
#define CIRCLE_LDO6QD4A

namespace R3{
/*
 * =====================================================================================
 *        Class:  Circle
 *  Description:  
 * =====================================================================================
 */
class Circle
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		Circle (const double* , const double* direction, double radius);            /* constructor */

		/* ====================  ACCESSORS     ======================================= */

		/* ====================  METHODS       ======================================= */
		bool			Intersect(const double* , const double*, double* , double*);
                //double                  GetPosition(int i);
                //double                  GetRadius();

	    friend auto& operator << ( std::ostream& o, Circle& c );

		/* ====================  OPERATORS     ======================================= */

	protected:
		/* ====================  DATA MEMBERS  ======================================= */
		double			mPosition[3];
public:
	    const double *getPosition() const;
		double GetRadius();

private:
		/* ====================  DATA MEMBERS  ======================================= */

		double			mDirection[3];
		double			mRadius;
		std::string     mName;
public:
	const std::string &getMName() const;

	void setMName(const std::string &mName);

}; /* -----  end of class Circle  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Circle
 *      Method:  Circle
 * Description:  constructor
 *--------------------------------------------------------------------------------------
 */
inline
	Circle::Circle (const double* pos,const double* dir, double radius)
{
	mPosition[0] = pos[0];
	mPosition[1] = pos[1];
	mPosition[2] = pos[2];

	mRadius = radius;
}  /* -----  end of method Circle::Circle  (constructor)  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Circle
 *      Method:  Intersect
 *     Returns:  bool
 * Description: 
 *--------------------------------------------------------------------------------------
 */
inline bool
Circle::Intersect ( const double* source, const double* direction, double* hit, double* t )
{
	if ( direction[2] > - EPSILON && direction[2] < EPSILON) return false;
	
	double lambda 	= (mPosition[2] - source[2]) / direction[2];
	double x 	= lambda * direction[0] + source[0];
	double y        = lambda * direction[1] + source[1];
        double r2        = mRadius * mRadius;

	double temp = x-mPosition[0]; temp=temp*temp;
	if ( temp > r2 ) return false;
        
	x = x-mPosition[0]; y=y-mPosition[1];
	temp = +sqrt( pow(x,2)+pow(y,2));
        if(temp > mRadius) return false;

	if (lambda > EPSILON && lambda < EPSILON) return false;
	
	*t = lambda;
	hit[0] = source[0] + direction[0] * lambda;
	hit[1] = source[1] + direction[1] * lambda;
	hit[2] = source[2] + direction[2] * lambda;
	return true; 

}

inline const std::string &Circle::getMName() const {
	return mName;
}

inline void Circle::setMName(const std::string &mName) {
	Circle::mName = mName;
}

inline const double *
Circle::getPosition() const {
	return mPosition;
}


inline auto& operator << ( std::ostream& o, Circle& c ){
		using namespace std;
		o << "\"position\" : {" << c.mPosition[0] << ", " << c.mPosition[1] << ", " << c.mPosition[2] << "} "<< endl;
		o << "\"direction\" : {" << c.mDirection[0] << ", " << c.mDirection[1] << ", " << c.mDirection[2] << "} "<< endl;
		o << "\"radius\" : " << c.mRadius << endl;
		return o;
}

inline double R3::Circle::GetRadius(){
	return mRadius;
}

}

#endif
