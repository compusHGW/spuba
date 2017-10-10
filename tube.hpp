#ifndef TUBE_K3MZD58W
#define TUBE_K3MZD58W

/* #####   HEADER FILE INCLUDES   ################################################### */

#include <stdio.h>
#include <math.h>
#include <ostream>

/*
 * =====================================================================================
 *        Class:  Tube
 *  Description:  The Tube class
 * =====================================================================================
 */
class Tube
{
	public:

		/* ====================  LIFECYCLE     ======================================= */
		Tube(double,double*,double);

		/* ====================  ACCESSORS     ======================================= */
		const double*		GetPosition();
	    double GetLength(){
			return mLength;
		}

		/* ====================  MUTATORS      ======================================= */
		bool 			Intersect(const double*, const double*, double*,double*,double*,double*);

		/* ====================  OPERATORS     ======================================= */

		friend auto& operator << ( std::ostream& o, Tube& t );

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */
		double 		mLength;
		double 		mRadius;
		double 		mPosition[3];
		double 		mDirection[3];

}; /* -----  end of class Tube  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Tube
 *      Method:  GetPosition
 *     Returns:  const double*
 *--------------------------------------------------------------------------------------
 */
	inline const double*
Tube::GetPosition ( )
{
	return mPosition;
}		/* -----  end of method Tube::GetPosition  ----- */

inline auto& operator << ( std::ostream& o, Tube& t ){
	using namespace std;
	o << "\"position\" : {" << t.mPosition[0] << ", " << t.mPosition[1] << ", " << t.mPosition[2] << "} "<< endl;
	o << "\"direction\" : {" << t.mDirection[0] << ", " << t.mDirection[1] << ", " << t.mDirection[2] << "} "<< endl;
	o << "\"radius\" : " << t.mRadius << endl;
	o << "\"length\" : " << t.mLength << endl;
	return o;
}

#endif /* end of include guard: TUBE_K3MZD58W */
