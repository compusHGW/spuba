#ifndef PROJECT_CHAMBER_H
#define PROJECT_CHAMBER_H

#include "Geometry.h"
#include "intersect_objects.h"

// TODO add baffels to the chamber
class Chamber {
public:
    enum MATERIAL_TYPE {
        A,
        C
    };
    MATERIAL_TYPE getMaterial_type() const {
        return material_type;
    }

    void setMaterial_type(MATERIAL_TYPE _material_type) {
        material_type = _material_type;
    }

    Geometry &getGeometry() {
        return geometry;
    }

    void setGeometry(const Geometry &_geometry) {
        geometry = _geometry;
    }

    std::vector<Cone3D>& getBaffels(){
        return baffels;
    }

    bool intersect_baffels(Particle3D& reflect, char baf){

        if (baf == 'Y') {
            for (auto& cone : baffels) {
                if (reflect.GetENumB() == -1) {
                    ::intersect(cone, reflect);
                } else {
                    return true;
                }
            }
            if(reflect.GetENumB() == 1){
                return true;
            }
        }
        return false;
    }

    bool intersect(Particle3D& p, char baf, T& t,double dthES){

        if (baf == 'Y') {
            for (auto& cone: baffels) {
                if (p.GetENumB() == -1) { ::intersect(cone, p); }
            }
        }

        bool hit = false;
        for( auto& obj : geometry.getGeometry_objects()) {
            if (p.GetENumB() == -1) {
                switch( obj.getType() ) {
                    case Geometry::SPHERE :
                        hit |= ::intersect(obj.object.s,&p,t,dthES);
                        break;
                    case Geometry::CIRCLE:
                        hit |= ::intersect(obj.object.c,&p,t,dthES);
                        break;
                    case Geometry::TUBE :
                        hit |= ::intersect(obj.object.t,&p,t);
                        break;
                }

            }
        }

        if ( !hit ){
            p.SetType(DIED);
        }
    }

private:
    MATERIAL_TYPE material_type;
    Geometry geometry;

    std::vector<Cone3D> baffels;
};


#endif
