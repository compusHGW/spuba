#ifndef PARTICLE3D_NVKITAZG
#define PARTICLE3D_NVKITAZG

#include <stdio.h>
#include <math.h>
/*
 * =====================================================================================
 *        Class:  Particle3D
 *  Description:  
 * =====================================================================================
 */
class Particle3D
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		Particle3D ( const double* s,double t,double p,double v){
			SetParticle(s,t,p,v);
		};                           /* constructor */

		/* ====================  ACCESSORS     ======================================= */
		void			SetENumS(int);   // number of the hitten cylinder segment
		void			SetENum(int);    // number of the hitten sphere segment
		void                    SetENumB(int);   // shows if a baffle is hidden (-1,1), jd 11.09.13

		void			SetType(int);
		void			SetHitpoint( double* );
		void			SetHitNormal( double*);
		void			SetVelocity( double );
		void			SetEnergy(double );
		void			SetW(double );
		void			SetParticle( const double* source, double theta, double phi, double v);
		void			SetParticleCartesian( const double* , const double*, double);

		int			GetType() const;
		int			GetENum();  
		int			GetENumS();
		int                     GetENumB();          // jd 11.09.13
		double*			GetDirectionPolar();
		const double*		GetDirection();
		const double*		GetPosition();
		const double*		GetHitpoint();
		const double*		GetHitNormal();
 		double			GetVelocity();
		double			GetEnergy();
		double			GetW();

		/* ====================  METHODS       ======================================= */
		void			Print();

		/* ====================  OPERATORS     ======================================= */

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */
		int			mType;
		double			mDirection[3];
		double			mSource[3];
		double			mDirectionPolar[3];
		double			mHitpoint[3];
		double			mHitNormal[3];
		double 			mV;
		double 			mE;
		double 			mW;            // Weight
		int         		mN, mNS, mNB;  // element number of cylinder, sphere, baffel

}; /* -----  end of class Particle3D  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetParticle
 *     Returns:  void
 * Description: 
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetParticle ( const double* source, double theta, double phi , double v)
{
#if 1
	mV = v;
	mSource[0] = source[0];  //x
	mSource[1] = source[1];  //y
	mSource[2] = source[2];  //z
        // ! theta = 'attac angle' 
	mDirection[0] = sin(theta) * cos(phi) ;
	mDirection[1] = sin(theta) * sin(phi) ;
	mDirection[2] = cos(theta) ;

	mDirectionPolar[0] = theta;
	mDirectionPolar[1] = phi;
	mDirectionPolar[2] = 0;
#endif
}		/* -----  end of method Particle3D::SetParticle  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetParticleCartesian
 *     Returns:  void
 * Description: 
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetParticleCartesian ( const double* source, const double* direction, double velocity)
{
	mSource[0] = source[0];
	mSource[1] = source[1];
	mSource[2] = source[2];

	mDirection[0] = direction[0];
	mDirection[1] = direction[1];
	mDirection[2] = direction[2];

	mV = velocity;
}		/* -----  end of method Particle3D::SetParticleCartesian  ----- */

/*--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetENumS
 *     Returns:  void
 *--------------------------------------------------------------------------------------*/
inline int
Particle3D::GetENumS (  )
{
	return mNS;
}		/* -----  end of method Particle3D::GetENumS  ----- */
/*--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetENumS
 *-------------------------------------------------------------------------------------*/
inline void
Particle3D::SetENumS ( int value )
{
	mNS	= value;
}		/* -----  end of method Particle3D::SetENumS  ----- */


/*--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetENum
 *     Returns:  void
 *-------------------------------------------------------------------------------------*/
inline int
Particle3D::GetENum (  )
{
	return mN;
}		/* -----  end of method Particle3D::GetENum  ----- */
/*--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetENum
 *-------------------------------------------------------------------------------------*/
inline void
Particle3D::SetENum ( int value )
{
	mN	= value;
}		/* -----  end of method Particle3D::SetENum  ----- */


/*--------------------------------------------------------------------------------------
 * Class:  Particle3D
 * Method:  GetENumB                ,jd 11.09.13
 * -------------------------------------------------------------------------------------*/
inline int
Particle3D::GetENumB (  )
{
        return mNB;
}               /* -----  end of method Particle3D::GetENumBC  ----- */
inline void
Particle3D::SetENumB ( int value )
{
        mNB      = value;
}               /* -----  end of method Particle3D::SetENum  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetType
 *     Returns:  void
 *--------------------------------------------------------------------------------------
 */
inline int
Particle3D::GetType(  ) const
{
	return mType;
}		/* -----  end of method Particle3D::GetType  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetType
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetType ( int value )
{
	mType	= value;
}		/* -----  end of method Particle3D::SetType  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  Print
 *     Returns:  void
 * Description: 
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::Print ( )
{
	printf("%f %f %f %f %f %f\n",mSource[0],
					mSource[1],
					mSource[2],
					mDirection[0],
					mDirection[1],
					mDirection[2]
			);
}		/* -----  end of method Particle3D::Print  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetDirectionPolar
 *     Returns:  void
 *--------------------------------------------------------------------------------------
 */
inline double*
Particle3D::GetDirectionPolar (  )
{
	return mDirectionPolar;
}		/* -----  end of method Particle3D::GetDirectionPolar  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetDirection
 *     Returns:  double*
 *--------------------------------------------------------------------------------------
 */
inline const double*
Particle3D::GetDirection (  )
{
	return mDirection;
}		/* -----  end of method Particle3D::GetDirection  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetPosition
 *     Returns:  const double*
 *--------------------------------------------------------------------------------------
 */
inline const double*
Particle3D::GetPosition (  )
{
	return mSource;
}		/* -----  end of method Particle3D::GetPosition  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetHitpoint
 *     Returns:  double
 *--------------------------------------------------------------------------------------
 */
inline const double*
Particle3D::GetHitpoint (  )
{
	return mHitpoint;
}		/* -----  end of method Particle3D::GetHitpoint  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  Set_Hitpoint
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetHitpoint ( double* value )
{
	mHitpoint[0]	= value[0];
	mHitpoint[1]	= value[1];
	mHitpoint[2]	= value[2];
}		/* -----  end of method Particle3D::SetHitpoint  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetHitNormal
 *     Returns:  const double*
 *--------------------------------------------------------------------------------------
 */
inline const double*
Particle3D::GetHitNormal (  )
{
	return mHitNormal;
}		/* -----  end of method Particle3D::GetHitNormal  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  Set_HitNormal
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetHitNormal ( double* value )
{
	mHitNormal[0]	= value[0];
	mHitNormal[1]	= value[1];
	mHitNormal[2]	= value[2];
}		/* -----  end of method Particle3D::SetHitNormal  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetVelocity
 *     Returns:  const double
 *--------------------------------------------------------------------------------------
 */
inline double
Particle3D::GetVelocity (  )
{
	return mV;
}		/* -----  end of method Particle3D::GetVelocity  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  Set_Velocity
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetVelocity ( double value )
{
	mV	= value;
}		/* -----  end of method Particle3D::SetVelocity  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetEnergy
 *     Returns:  const double
 *--------------------------------------------------------------------------------------
 */
inline double
Particle3D::GetEnergy (  )
{
	return mE;
}		/* -----  end of method Particle3D::GetVelocity  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetEnergy
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetEnergy ( double value )
{
	mE	= value;
}		/* -----  end of method Particle3D::SetEnergy  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  GetW
 *     Returns:  const double
 *--------------------------------------------------------------------------------------
 */
inline double
Particle3D::GetW (  )
{
	return mW;
}		/* -----  end of method Particle3D::GetW----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Particle3D
 *      Method:  SetW
 *--------------------------------------------------------------------------------------
 */
inline void
Particle3D::SetW ( double value )
{
	mW	= value;
}		/* -----  end of method Particle3D::SetW  ----- */

#endif /* end of include guard: PARTICLE3D_NVKITAZG */

