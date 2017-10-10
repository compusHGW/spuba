#include <cassert>
#include "cone.h"
#include "circle.hpp"
#include "tube.hpp"
#include "sphere.hpp"
#include "ParticleSource.h"
#include "intersect_objects.h"
#include "utility.hpp"

void intersect(Cone3D &cone, Particle3D &p) {

    /*-----------------------------------------------------------------------------
     *  intersect with the cone
     *-----------------------------------------------------------------------------*/
    vector3d pos, dir;
    pos[0] = p.GetPosition()[0];
    pos[1] = p.GetPosition()[1];
    pos[2] = p.GetPosition()[2];
    dir[0] = p.GetDirection()[0];
    dir[1] = p.GetDirection()[1];
    dir[2] = p.GetDirection()[2];
    if (cone.Intersect(pos, dir)) {
        p.SetHitpoint(cone.GetIntP().data());
        p.SetHitNormal(cone.GetNormal().data());
        p.SetENumB(1);
        return;
    }
    p.SetENumB(-1);

}

void intersect(Sphere *sphere, Tube *tube, Particle3D *p, T &t, double dthES) {

    double hit1[3] = {0.0, 0.0, 0.0};
    double hit2[3] = {0.0, 0.0, 0.0};
    double thit1[3] = {0.0, 0.0, 0.0};
    double thit2[3] = {0.0, 0.0, 0.0};
    double st1 = 0, st2 = 0;


    /*-----------------------------------------------------------------------------
     *  intersect with the tube
     *-----------------------------------------------------------------------------*/
    double tt1 = 0, tt2 = 0;
    if (tube->Intersect(p->GetPosition(), p->GetDirection(), thit1, thit2, &tt1, &tt2)) {
        double tt = tt2;
        if (tt1 > tt2) {
            tt = tt1;
            double change[3] = {thit1[0], thit1[1], thit1[2]};
            thit1[0] = thit2[0];
            thit1[1] = thit2[1];
            thit1[2] = thit2[2];
            thit2[0] = change[0];
            thit2[1] = change[1];
            thit2[2] = change[2];
        }
        double normal[3];
        int cyl_enum = int(tt * double(NUM_CYL_SEGS));
        if (cyl_enum >= 0 && cyl_enum < t.l_v_nPLandCylSeg.size())
            t.l_v_nPLandCylSeg[cyl_enum]++;
        else
            assert(0);
        p->SetENum(cyl_enum);

        normal[0] = -thit2[0];
        normal[1] = -thit2[1];
        normal[2] = 0;
        norm(normal);

        p->SetHitpoint(thit2);
        p->SetHitNormal(normal);
        p->SetENumS(-1);
        return;
    }
    p->SetENum(-1);

    /*-----------------------------------------------------------------------------
     *  intersect with the sphere
     *-----------------------------------------------------------------------------*/
    if (sphere->Intersect(p->GetPosition(), p->GetDirection(), hit1, hit2, &st1, &st2)) {
        double vv[3] = {0.0, 0.0, 1.0};
        double vv1[3];
        vv1[0] = hit2[0];
        vv1[1] = hit2[1];
        vv1[2] = hit2[2] - sphere->GetPosition()[2];
        double ang = Angle(vv1, vv);
        int seg = ang / dthES;
        if (seg >= 0 && seg < t.l_v_nPLandSpSeg.size())
            t.l_v_nPLandSpSeg[seg]++;
        else
                ASSERT(0);
        p->SetENumS(seg);
        p->SetHitpoint(hit2);
        double normal[3] = {sphere->GetPosition()[0] - hit2[0],
                            sphere->GetPosition()[1] - hit2[1],
                            sphere->GetPosition()[2] - hit2[2]
        };
        norm(normal);

        p->SetHitNormal(normal);
        return;
    }
    p->SetType(DIED);

}

bool intersect( Tube* tube, Particle3D* p, T&t ){

    double thit1[3] = {0.0, 0.0, 0.0};
    double thit2[3] = {0.0, 0.0, 0.0};
    double st1 = 0, st2 = 0;

    /*-----------------------------------------------------------------------------
     *  intersect with the tube
     *-----------------------------------------------------------------------------*/
    double tt1 = 0, tt2 = 0;
    if (tube->Intersect(p->GetPosition(), p->GetDirection(), thit1, thit2, &tt1, &tt2)) {
        double tt = tt2;
        if (tt1 > tt2) {
            tt = tt1;
            double change[3] = {thit1[0], thit1[1], thit1[2]};
            thit1[0] = thit2[0];
            thit1[1] = thit2[1];
            thit1[2] = thit2[2];
            thit2[0] = change[0];
            thit2[1] = change[1];
            thit2[2] = change[2];
        }
        double normal[3];
        int cyl_enum = int(tt * double(NUM_CYL_SEGS));
        if (cyl_enum >= 0 && cyl_enum < t.l_v_nPLandCylSeg.size())
            t.l_v_nPLandCylSeg[cyl_enum]++;
        else
            assert(0);
        p->SetENum(cyl_enum);

        normal[0] = -thit2[0];
        normal[1] = -thit2[1];
        normal[2] = 0;
        norm(normal);

        p->SetHitpoint(thit2);
        p->SetHitNormal(normal);
        p->SetENumS(-1);
        return true;
    }
    p->SetENum(-1);
    return false;
}

bool intersect( Sphere* sphere, Particle3D* p, T& t, double dthES) {


    double hit1[3] = {0.0, 0.0, 0.0};
    double hit2[3] = {0.0, 0.0, 0.0};
    double st1 = 0, st2 = 0;
    /*-----------------------------------------------------------------------------
     *  intersect with the sphere
     *-----------------------------------------------------------------------------*/
    if (sphere->Intersect(p->GetPosition(), p->GetDirection(), hit1, hit2, &st1, &st2)) {
        double vv[3] = {0.0, 0.0, 1.0};
        double vv1[3];
        vv1[0] = hit2[0];
        vv1[1] = hit2[1];
        vv1[2] = hit2[2] - sphere->GetPosition()[2];
        double ang = Angle(vv1, vv);
        int seg = ang / dthES;
        if (seg >= 0 && seg < t.l_v_nPLandSpSeg.size())
            t.l_v_nPLandSpSeg[seg]++;
        else
                ASSERT(0);
        p->SetENumS(seg);
        p->SetHitpoint(hit2);
        double normal[3] = {sphere->GetPosition()[0] - hit2[0],
                            sphere->GetPosition()[1] - hit2[1],
                            sphere->GetPosition()[2] - hit2[2]
        };

        norm(normal);
        p->SetHitNormal(normal);
        return true;
    }
    p->SetType(DIED);
    return false;
}

bool intersect( R3::Circle* circle, Particle3D* p, T& t , double dthES){
    /*-----------------------------------------------------------------------------
     *  intersect with the circle
     *-----------------------------------------------------------------------------*/
    double chit[3] = {0.0, 0.0, 0.0};
    double ct = 0;
    if (circle->Intersect(p->GetPosition(), p->GetDirection(), chit, &ct)) {
        double vv[3] = {0.0, 0.0, 1.0};
        double vv1[3];
        vv1[0] = chit[0];
        vv1[1] = chit[1];
        vv1[2] = chit[2];
        double ang = Angle(vv1, vv);
        int seg = ang / dthES;
        if (seg >= 0 && seg < t.l_v_nPLandSpSeg.size())
            t.l_v_nPLandSpSeg[seg]++;
        else
                ASSERT(0);
        p->SetENumS(seg);
        p->SetHitpoint(chit);
        double normal[3] = {0, 0, -chit[2]};
        norm(normal);

        p->SetHitNormal(normal);
        return true;
    }

    return false;
}

/* intersection of particle with thruster channel (tube + planar bottom cap) jd, 05.12'13*/
void intersectThru(R3::Circle *circle, Tube *tube, Particle3D *p, double dthES) {
    p->SetENum(-99999);
    p->SetENumS(-99999);

    double tt1 = 0, tt2 = 0;
    double thit1[3] = {0.0, 0.0, 0.0};
    double thit2[3] = {0.0, 0.0, 0.0};
    double chit[3] = {0.0, 0.0, 0.0};
    double ct = 0;


    /*-----------------------------------------------------------------------------
     *  intersect with the tube
     *-----------------------------------------------------------------------------*/
    if (tube->Intersect(p->GetPosition(), p->GetDirection(), thit1, thit2, &tt1, &tt2)) {
        double tt = tt2;
        if (thit1[2] < thit2[2]) {
            tt = tt1;
            double change[3] = {thit1[0], thit1[1], thit1[2]};
            thit1[0] = thit2[0];
            thit1[1] = thit2[1];
            thit1[2] = thit2[2];
            thit2[0] = change[0];
            thit2[1] = change[1];
            thit2[2] = change[2];
        }
        double normal[3];
        int cyl_enum = int(abs(tt) * double(NUM_CYL_SEGS));
        p->SetENum(cyl_enum);
        normal[0] = -thit2[0];
        normal[1] = -thit2[1];
        normal[2] = 0;
        norm(normal);
        p->SetHitpoint(thit2);
        p->SetHitNormal(normal);
        p->SetENumS(-1);
        return;
    }

        /*-----------------------------------------------------------------------------
         *  intersect with the circle
         *-----------------------------------------------------------------------------*/

    else if (circle->Intersect(p->GetPosition(), p->GetDirection(), chit, &ct)) {
        p->SetENum(-1);
        double vv[3] = {0.0, 0.0, 1.0};
        double vv1[3];
        vv1[0] = chit[0];
        vv1[1] = chit[1];
        vv1[2] = chit[2];
        double ang = Angle(vv1, vv);
        int seg = ang / dthES;
        p->SetENumS(seg);
        p->SetHitpoint(chit);
        double normal[3] = {0, 0, -chit[2]};
        norm(normal);
        p->SetHitNormal(normal);
        return;
    } else {
        p->SetENum(-1);
        p->SetENumS(-1);
        printf("!! Intersection with thruster exit plane but not with channel nor bottom !!\n");
        printf("@intersectThru: \tPos=(%f, %f, %f)\t direction = (%f, %f, %f)\n", p->GetPosition()[0],
               p->GetPosition()[1], p->GetPosition()[2], p->GetDirection()[0], p->GetDirection()[1],
               p->GetDirection()[2]);
        p->SetType(DIED);
    }

}
