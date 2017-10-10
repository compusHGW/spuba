#ifndef  INTERSECTION_INC
#define  INTERSECTION_INC
#include <array>

#include	<iostream>
#include	<vector>
#include	<cmath>
#include	<iostream>

using namespace std;

#define SMALL_NUM  0.0000000001

#define dot(u,v)   ((u).x[0] * (v).x[0] + (u).x[1] * (v).x[1] )
#define perp(u,v)  ((u).x[0] * (v).x[1] - (u).x[1] * (v).x[0] )  // perp product (2D)

enum {_1D=1,_2D,_3D};

template<int DIM> class Point
{
  public: 
    std::array<double,DIM> x;
};

#if 1
class SqMatrix3
{
  public:
    std::array<std::array<double,3>,3> m;
    //boost::multi_array<double,2> m;
    SqMatrix3(unsigned _s);
};
#endif
typedef Point<_2D> Point2D;
typedef Point<_2D> Vector2D;

Vector2D operator +(const Vector2D& a, const Vector2D& b);

Vector2D operator -(const Vector2D& a, const Vector2D& b);

Vector2D operator *(const double c, const Vector2D& v);

double operator *(const Vector2D& a, const Vector2D& b);

std::ostream & operator<<(std::ostream &_os, const Vector2D& _p);

class Segment2D;

template<int DIM> class Shape
{
  protected:
    std::vector< Point<DIM> > m_vIntP;
  public:

    int GetQuantity() const {  return m_vIntP.size();  };
    const Point<DIM>& GetPoint(int _i) const { return m_vIntP[_i]; }; //return iteraction point (for circle?)
};

typedef Shape<_2D> Shape2D;

#if 1
class Ray2D : public Shape2D
{
  private:
    Vector2D m_vO;
    Vector2D m_vD;
    vector<double> m_vPar;
  public:
    Vector2D& GetOrigin() {return m_vO;}
    Vector2D& GetDirection() {return m_vD;}
    double GetPar(int i){return m_vPar[i];}
    double GetTheta() {return atan2(m_vD.x[1],m_vD.x[0]);}
    double GetLength() {return sqrt(m_vD.x[1]*m_vD.x[1]+m_vD.x[0]*m_vD.x[0]);}
    Ray2D(const Vector2D& _vO,const Vector2D& _vD) : m_vO(_vO),m_vD(_vD){};
    Ray2D(double _p1x1,double _p1x2, double _p2x1, double _p2x2);
    Ray2D(){};
    virtual bool FindIntersection(Segment2D &);
};
#endif

#if 1
class Segment2D : public Shape2D
{
  private:
    Vector2D m_vO;
    Vector2D m_vD;
  public:
    Vector2D& GetOrigin() {return m_vO;}
    Vector2D& GetDirection() {return m_vD;}
    double GetTheta() {return atan2(m_vD.x[1],m_vD.x[0]);}
    double GetLength() {return sqrt(m_vD.x[1]*m_vD.x[1]+m_vD.x[0]*m_vD.x[0]);}
    Segment2D(const Vector2D& _vO,const Vector2D& _vD) : m_vO(_vO),m_vD(_vD){};
    Segment2D(double _p1x1,double _p1x2, double _p2x1, double _p2x2);
    Segment2D(){};
    virtual bool FindInt(Segment2D& );
};
#endif

#if 1
class Circle : public Shape2D
{
  private:
    Vector2D m_vO;
    double   m_dR;
  public:
    Vector2D& GetOrigin() {return m_vO;}
    double    GetRadius() {return m_dR;}
    Circle(const Vector2D& _vO,double _r) : m_vO(_vO),m_dR(_r){};
    Circle(double _x1,double _x2, double _r);
    virtual bool FindInt(Segment2D& );
    bool         FindInt(Circle& _c);
};

#endif
#if 0
inline double GetVecTheta(const Vector2D& v) {return atan2(v.x[1],v.x[0]);};
inline double GetLength(const Vector2D& v) {return sqrt(v.x[1]*v.x[1]+v.x[0]*v.x[0]);};
#endif

#endif   /* ----- #ifndef INTERSECTION_INC  ----- */
