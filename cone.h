#ifndef  CONE_INC
#define  CONE_INC

/* #####   HEADER FILE INCLUDES   ################################################### */
#include	<cmath>
#include	"bgeom.hpp"
/*
 * =====================================================================================
 *        Class:  Cone3D
 *  Description:  Cone class with ray-cone intersection routine
 * =====================================================================================
 */
class Cone3D
{
  public:
    /* ====================  LIFECYCLE     ======================================= */
    Cone3D (vector3d& _pos,double _height, double _Rl,double _a);   /* constructor */

    /* ====================  ACCESSORS     ======================================= */
    vector3d& GetIntP() {return mIntP;}
    vector3d& GetNormal() {return mNormal;}

    /* ====================  MUTATORS      ======================================= */
    bool Intersect (const vector3d& _pos, const vector3d& _dir);

    /* ====================  OPERATORS     ======================================= */

  protected:
    /* ====================  DATA MEMBERS  ======================================= */

  private:
    /* ====================  DATA MEMBERS  ======================================= */
    vector3d mPosition;
    vector3d mIntP;
    vector3d mNormal;
    double mHeight;
    double mRl;
    double mAngle;
    double mRs;
}; /* -----  end of class Cone3D  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Cone3D
 *      Method:  Cone3D
 * Description:  constructor
 *--------------------------------------------------------------------------------------
 */
inline Cone3D::Cone3D (vector3d& _pos,double _height, double _Rl,double _a)
{
  mPosition = _pos;
  mHeight = _height;
  mRl = _Rl;
  mAngle = _a;
  mRs = mRl-mHeight*tan(mAngle);
}  /* -----  end of method Cone3D::Cone3D  (constructor)  ----- */


#endif   /* ----- #ifndef CONE_INC  ----- */
