#ifndef PROJECT_INTERSECT_OBJECTS_H
#define PROJECT_INTERSECT_OBJECTS_H

#include "shared.hpp"
#include "cone.h"

void intersectThru(R3::Circle *circle, Tube *tube, Particle3D *p, double dthES);
void intersect(Sphere *sphere, Tube *tube, Particle3D *p, T &t, double dthES);
void intersect(R3::Circle *circle, Tube *tube, Particle3D *p, T &t, double dthES);
void intersect(Cone3D &cone, Particle3D &p);

bool intersect( Sphere* sphere, Particle3D* p, T& t, double dthES);
bool intersect( Tube* tube, Particle3D* p, T&t );
bool intersect( R3::Circle* circle, Particle3D* p, T& t , double dthES);

#endif //PROJECT_INTERSECT_OBJECTS_H
