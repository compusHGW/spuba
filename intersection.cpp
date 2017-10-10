#include	"intersection.h"

SqMatrix3::SqMatrix3(unsigned s)
{
}

Vector2D operator *(const SqMatrix3& m, const Vector2D& v)
{
  Vector2D r;
  r.x[0] = m.m[0][0]*v.x[0]+m.m[0][1]*v.x[1];
  r.x[1] = m.m[1][0]*v.x[0]+m.m[1][1]*v.x[1];
  return r;
}

Vector2D operator +(const Vector2D& a, const Vector2D& b)
{
  Vector2D v;
  v.x[0] = a.x[0]+b.x[0];
  v.x[1] = a.x[1]+b.x[1];
  return v;
}

Vector2D operator -(const Vector2D& a, const Vector2D& b)
{
  Vector2D v;
  v.x[0] = a.x[0]-b.x[0];
  v.x[1] = a.x[1]-b.x[1];
  return v;
}

double operator *(const Vector2D& a, const Vector2D& b)
{
  return a.x[0]*b.x[0]+a.x[1]*b.x[1];
}

Vector2D operator *(const double c, const Vector2D& v)
{
  Vector2D r;
  r.x[0]=v.x[0]*c;
  r.x[1]=v.x[1]*c;
  return r;
}

std::ostream & operator<<(std::ostream &_os, const Vector2D& _p)
{
  _os << _p.x[0] << " "<< _p.x[1];
  return _os; 
}

Ray2D::Ray2D(double _p1x1,double _p1x2, double _p2x1, double _p2x2)
{
  m_vO.x[0]=_p1x1;
  m_vO.x[1]=_p1x2;
  m_vD.x[0]=_p2x1-_p1x1;
  m_vD.x[1]=_p2x2-_p1x2;
}

bool Ray2D::FindIntersection(Segment2D &_s)
{
  m_vIntP.clear();
  Vector2D& o = _s.GetOrigin();
  Vector2D& d = _s.GetDirection();
  Vector2D    w = o - m_vO;
  double    D = perp(d,m_vD);

  // test if they are parallel (includes either being a point)
  if (fabs(D) < SMALL_NUM)           // S1 and S2 are parallel
  {
    if (perp(d,w) != 0 || perp(m_vD,w) != 0)
    {
      return false;                   // they are NOT collinear
    }
    //
    // they are collinear segments - get overlap (or not)
    double t0, t1;                   // endpoints of S1 in eqn for S2
    Vector2D w2 = (o + d) - m_vO;
    if (m_vD.x[0] != 0) {
      t0 = w.x[0] / m_vD.x[0];
      t1 = w2.x[0] / m_vD.x[0];
    }
    else {
      t0 = w.x[1] / m_vD.x[1];
      t1 = w2.x[1] / m_vD.x[1];
    }
    if (t0 > t1) {                  // must have t0 smaller than t1
      double t=t0; t0=t1; t1=t;    // swap if not
    }
    if (t1 < 0) {
      return false;     // NO overlap
    }
    t0 = t0<0? 0 : t0;              // clip to min 0
    if (t0 == t1) {                 // intersect is a point
      Point2D p0;
      p0 = m_vO + t0 * m_vD;
      m_vIntP.push_back(p0);
      return true;
    }

    // they overlap in a valid subsegment
    Point2D p0;
    p0 = m_vO + t0 * m_vD;
    m_vIntP.push_back(p0);
    p0 = m_vO + t1 * m_vD;
    m_vIntP.push_back(p0);
    return true;
  }

  // the segments are skew and may intersect in a point
  // get the intersect parameter for S1
  double     sI = perp(m_vD,w) / D;
  if (sI < 0 || sI > 1)               // no intersect with S1
    return false;

  // get the intersect parameter for S2
  double     tI = perp(d,w) / D;
  if (tI < 0)               // no intersect with ray
    return false;

  Point2D p0;
  p0 = o + sI * d;               // compute S1 intersect point
  m_vIntP.push_back(p0);
  m_vPar.push_back(tI);
  return true;
}

Segment2D::Segment2D(double _p1x1,double _p1x2, double _p2x1, double _p2x2)
{
  m_vO.x[0]=_p1x1;
  m_vO.x[1]=_p1x2;
  m_vD.x[0]=_p2x1-_p1x1;
  m_vD.x[1]=_p2x2-_p1x2;
}

bool Segment2D::FindInt(Segment2D& _s)
{
  m_vIntP.clear();
  Vector2D& o = _s.GetOrigin();
  Vector2D& d = _s.GetDirection();
  Vector2D  w = o - m_vO;
  double    D = perp(d,m_vD);

  // test if they are parallel (includes either being a point)
  if (fabs(D) < SMALL_NUM)           // S1 and S2 are parallel
  {
    if (perp(d,w) != 0 || perp(m_vD,w) != 0)
    {
      return false;                   // they are NOT collinear
    }
    //
    // they are collinear segments - get overlap (or not)
    double t0, t1;                   // endpoints of S1 in eqn for S2
    Vector2D w2 = (o + d) - m_vO;
    if (m_vD.x[0] != 0) {
      t0 = w.x[0] / m_vD.x[0];
      t1 = w2.x[0] / m_vD.x[0];
    }
    else {
      t0 = w.x[1] / m_vD.x[1];
      t1 = w2.x[1] / m_vD.x[1];
    }
    if (t0 > t1) {                  // must have t0 smaller than t1
      double t=t0; t0=t1; t1=t;    // swap if not
    }
    if (t0 > 1 || t1 < 0) {
      return false;     // NO overlap
    }
    t0 = t0<0? 0 : t0;              // clip to min 0
    t1 = t1>1? 1 : t1;              // clip to max 1
    if (t0 == t1) {                 // intersect is a point
      Point2D p0;
      p0 = m_vO + t0 * m_vD;
      m_vIntP.push_back(p0);
      return true;
    }

    // they overlap in a valid subsegment
    Point2D p0;
    p0 = m_vO + t0 * m_vD;
    m_vIntP.push_back(p0);
    p0 = m_vO + t1 * m_vD;
    m_vIntP.push_back(p0);
    return true;
  }

  // the segments are skew and may intersect in a point
  // get the intersect parameter for S1
  double     sI = perp(m_vD,w) / D;
  if (sI < 0 || sI > 1)               // no intersect with S1
    return false;

  // get the intersect parameter for S2
  double     tI = perp(d,w) / D;
  if (tI < 0 || tI > 1)               // no intersect with S2
    return false;

  Point2D p0;
  p0 = o + sI * d;               // compute S1 intersect point
  m_vIntP.push_back(p0);
  return true;
}
#if 1
Circle::Circle(double _x1,double _x2, double _r)
{
  m_vO.x[0]=_x1;
  m_vO.x[1]=_x2;
  m_dR = _r;
}

bool Circle::FindInt(Segment2D& _s)
{
  m_vIntP.clear();
  Vector2D& o = _s.GetOrigin();
  Vector2D& d = _s.GetDirection();
  Vector2D w = o - m_vO;
  double a = dot(d,d);
  double b = 2.*dot(d,w);
  double c = dot(w,w)-m_dR*m_dR;
  double D = b*b-4.*a*c;
  if(abs(D)<SMALL_NUM) // Line touches the circle
  {
    double t = -0.5*b / a;
    if(t < 0 || t > 1)
      return false;
    Point2D p;
    p = o + t*d;
    m_vIntP.push_back(p);
    return true;
  }
  double t1,t2;
  int ind = 2;
  t1 = -.5*(b+sqrt(D))/a;
  t2 = .5*(sqrt(D)-b)/a;
  Point2D p;
  if(t1 < -SMALL_NUM || t1 > 1+SMALL_NUM) --ind;
  else
  {
    p = o + t1*d;
    m_vIntP.push_back(p);
  }
  if(t2 < -SMALL_NUM || t2 > 1+SMALL_NUM) --ind;
  else
  {
    p = o + t2*d;
    m_vIntP.push_back(p);
  }
  return ind != 0;
}

bool Circle::FindInt(Circle& _c)
{
  m_vIntP.clear();
  double dx = _c.m_vO.x[0] - m_vO.x[0];
  double dy = _c.m_vO.x[1] - m_vO.x[1];
  double d = sqrt((dy*dy)+(dx*dx));
  if(d>m_dR+_c.m_dR) return false;//Circles do not intersect
  if(d<abs(m_dR-_c.m_dR)) return false;//One circle is contained into another
  // point 2' is the point where the line through the circle
  // intersection points crosses the line between the circle
  // centers.  

  // Determine the distance from point 0 to point 2.
  double a = ((m_dR*m_dR) - (_c.m_dR*_c.m_dR) + (d*d)) / (2.0 * d) ;

  // Determine the coordinates of point 2.
  double x2 = m_vO.x[0] + (dx * a/d);
  double y2 = m_vO.x[1] + (dy * a/d);

  // Determine the distance from point 2 to either of the
  // intersection points.
  //
  double h = sqrt((m_dR*m_dR) - (a*a));

  // Now determine the offsets of the intersection points from
  // point 2.
  double rx = -dy * (h/d);
  double ry = dx * (h/d);

  // Determine the absolute intersection points.
  Point2D p;
  // point 1
  p.x[0] = x2 + rx;
  p.x[1] = y2 + ry;
  m_vIntP.push_back(p);
  // point 2
  p.x[0] = x2 - rx;
  p.x[1] = y2 - ry;
  m_vIntP.push_back(p);
  return true;
}
#endif