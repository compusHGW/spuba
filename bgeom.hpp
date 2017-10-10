#ifndef  BGEOM_INC
#define  BGEOM_INC


#include <array>
#include	"intersection.h"
#include	"const.h"
#include	<cmath>

typedef std::array<double,3> vector3d;

inline vector3d operator +(const vector3d& a, const vector3d& b)
{
  vector3d v;
  v[0] = a[0]+b[0];
  v[1] = a[1]+b[1];
  v[2] = a[2]+b[2];
  return v;
}

inline vector3d operator -(const vector3d& a)
{
  vector3d v;
  v[0] = -a[0];
  v[1] = -a[1];
  v[2] = -a[2];
  return v;
}

inline vector3d operator -(const vector3d& a, const vector3d& b)
{
  vector3d v;
  v[0] = a[0]-b[0];
  v[1] = a[1]-b[1];
  v[2] = a[2]-b[2];
  return v;
}

inline double operator *(const vector3d& a, const vector3d& b)
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline vector3d operator *(const vector3d& v, const double c )
{
  vector3d r;
  r[0]=v[0]*c;
  r[1]=v[1]*c;
  r[2]=v[2]*c;
  return r;
}

inline vector3d operator *(const double c, const vector3d& v)
{
  vector3d r;
  r[0]=v[0]*c;
  r[1]=v[1]*c;
  r[2]=v[2]*c;
  return r;
}

inline vector3d operator *(const SqMatrix3& m, const vector3d& v)
{
  vector3d r;
  r[0] = m.m[0][0]*v[0]+m.m[0][1]*v[1]+m.m[0][2]*v[2];
  r[1] = m.m[1][0]*v[0]+m.m[1][1]*v[1]+m.m[1][2]*v[2];
  r[2] = m.m[2][0]*v[0]+m.m[2][1]*v[1]+m.m[2][2]*v[2];
  return r;
}

inline void MRotZ( SqMatrix3& m, double alpha)
{
  m.m[0][0] =  cos( alpha ) ;
  m.m[0][1] =  - sin( alpha ) ;
  m.m[0][2] =  0.0 ;

  m.m[1][0] = sin (alpha ) ;
  m.m[1][1] = cos ( alpha ) ;
  m.m[1][2] = 0.0 ;

  m.m[2][0] = 0.0 ;
  m.m[2][1] = 0.0 ;
  m.m[2][2] = 1.0 ;
}

inline void MRotX( SqMatrix3& m, double alpha)
{
  m.m[0][0] =  1.0;
  m.m[0][1] =  0.0 ;
  m.m[0][2] =  0.0 ;

  m.m[1][0] = 0.0 ;
  m.m[1][1] = cos ( alpha ) ;
  m.m[1][2] = -sin( alpha ) ;

  m.m[2][0] = 0.0 ;
  m.m[2][1] = sin( alpha ) ;
  m.m[2][2] = cos( alpha ) ;
}

inline void MRotY( SqMatrix3& m, double alpha)
{
  m.m[0][0] =  cos ( alpha );
  m.m[0][1] =  0.0 ;
  m.m[0][2] =  sin ( alpha ) ;

  m.m[1][0] = 0.0 ;
  m.m[1][1] = 1.0 ;
  m.m[1][2] = 0.0 ;

  m.m[2][0] = -sin( alpha) ;
  m.m[2][1] = 0.0 ;
  m.m[2][2] = cos( alpha ) ;
}

inline void Norm(vector3d& v)
{
  register double l = sqrt(SQR(v[0])+SQR(v[1])+SQR(v[2]));
  v[0] /= l;
  v[1] /= l;
  v[2] /= l;
}

inline vector3d Vec(const vector3d& a, const vector3d& b)
{
  vector3d r;
  r[0] = a[1]*b[2]-a[2]*b[1];
  r[1] = a[2]*b[0]-a[0]*b[2];
  r[2] = a[0]*b[1]-a[1]*b[0];
  return r;
}
#endif
