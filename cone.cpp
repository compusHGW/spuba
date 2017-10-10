#include	"cone.h"
 
bool Cone3D::Intersect (const vector3d& _pos, const vector3d& _dir)
{
  double    y1 = mPosition[2];
  double    x1 = sqrt(SQR(mPosition[0])+SQR(mPosition[1]))+mRl;
  double    y2 = mPosition[2]+mHeight;
  double    x2 = x1+mRs-mRl;
  Segment2D s(x1,y1,x2,y2);
  double    xr1 = sqrt(SQR(_pos[0])+SQR(_pos[1]));
  double    xr2 = xr1 + sqrt(SQR(_dir[0])+SQR(_dir[1]));
  Ray2D     r(xr1,_pos[2],xr2,_pos[2]+_dir[2]);
  bool      ret = r.FindIntersection(s);

  if(!ret) return false;
  
  mIntP = _pos+_dir*r.GetPar(0);
  vector3d vo = /* vt */ mIntP;
  unsigned nx = abs(_dir[0])>abs(_dir[1])?0:1;

  if( !nx )
  {
    double phi = atan2( vo[1] , vo[0] );

    mNormal[0] = 1;
    mNormal[1] = 0;
    mNormal[2] = 0;
    SqMatrix3 m(3);

    MRotY(m,mAngle+M_PI);
    mNormal = m*mNormal;

    MRotZ(m,phi);
    mNormal = m*mNormal;
  }
  else
  {
    double phi = atan2( vo[1] , vo[0] );
    mNormal[0] = 0;
    mNormal[1] = 1;
    mNormal[2] = 0;
    SqMatrix3 m(3);

    MRotX(m,mAngle+M_PI);
    mNormal = m*mNormal;

    MRotZ(m,phi-0.5*M_PI);
    mNormal = m*mNormal;
  }

  return true;
}

