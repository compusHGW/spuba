#ifndef PROJECT_GEOMETRY_H
#define PROJECT_GEOMETRY_H


#include <string>
#include "sphere.hpp"
#include "intersection.h"
#include "tube.hpp"
#include "Cylinder.hpp"

class Geometry {
public:
    enum OBJECT_TYPE {
        CIRCLE,
        TUBE,
        SPHERE,
        CYLINDER
    };
    union ObjectsUnion {
        R3::Circle *c;
        Tube *t;
        Sphere *s;
        Cylinder* cyl;
    };
    struct Objects{
        ObjectsUnion object;

        ObjectsUnion &getObject(){
            return object;
        }

        OBJECT_TYPE type;

        OBJECT_TYPE getType() const {
            return type;
        }

        void setType(OBJECT_TYPE _type) {
            type = _type;
        }
        friend auto& operator << ( ostream& o, Objects& obj );
    };

    vector<Objects> &getGeometry_objects(){
        return geometry_objects;
    }

    void setGeometry_objects(const vector<Objects> &geometry_objects) {
        Geometry::geometry_objects = geometry_objects;
    }

    friend auto& operator << ( ostream& o, Geometry& g );

private:

    vector<Objects> geometry_objects;
};

inline auto& operator << ( ostream& o, Geometry::Objects& obj ){
    switch (obj.getType()) {
        case Geometry::CYLINDER :
            o << *obj.getObject().cyl << endl;
            break;
#if 0
        case Geometry::CIRCLE :
            o << *obj.getObject().c << endl;
            break;
        case Geometry::SPHERE :
            o << *obj.getObject().s << endl;
            break;
        case Geometry::TUBE :
            o << *obj.getObject().t << endl;
            break;
#endif
    }
    return o;
}

inline auto& operator << ( ostream& o, Geometry& g ){
    for( auto& obj : g.geometry_objects ){
        o << obj << endl;
    }

    return o;
}

#endif //PROJECT_GEOMETRY_H
