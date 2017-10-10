#ifndef PROJECT_CYLINDER_H
#define PROJECT_CYLINDER_H


#include <memory>
#include "intersection.h"
#include "tube.hpp"
#include "circle.hpp"

class Cylinder {
public:
    Cylinder(
            double* pos,
            double* dir,
            double rad,
            double _length
    ) {
        bottom.reset(new R3::Circle(&pos[0],&dir[0],rad));
        double pos_top[3];
        pos_top[0] = pos[0]+dir[0]*_length;
        pos_top[1] = pos[1]+dir[1]*_length;
        pos_top[2] = pos[2]+dir[2]*_length;
        top.reset(new R3::Circle(&pos_top[0],&dir[0],rad));
        wall.reset(new Tube( rad, pos, _length));
    }

    double GetLength(){
        return wall->GetLength();
    }
    double GetRadius(){
        return bottom->GetRadius();
    }

    friend auto& operator << ( ostream& o, Cylinder& c );

public:
    unique_ptr<R3::Circle> &getBottom(){
        return bottom;
    }

    unique_ptr<Tube> &getWall() {
        return wall;
    }

    unique_ptr<R3::Circle> &getTop() {
        return top;
    }


private:
    unique_ptr<R3::Circle> top; // aka exit
    unique_ptr<R3::Circle> bottom;
    unique_ptr<Tube> wall;
};

inline auto& operator << ( ostream& o, Cylinder& c ){
    o << "\"top\"" << endl << *c.top << endl;
    o << "\"wall\"" << endl << *c.wall << endl;
    o << "\"bottom\"" << endl << *c.bottom << endl;
    return o;
}

#endif
