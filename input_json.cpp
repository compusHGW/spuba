#include <iomanip>
#include "input_json.h"
#include <json/value.h>
#include <json/reader.h>
#include <fstream>
#include <cassert>
#include "yield.h"

using namespace std;

void parse_pos_dir_rad_from_json(const Json::Value &object, double* pos, double* dir, double& rad) {

    pos[0] = object["position"]["x"].asDouble();
    pos[1] = object["position"]["y"].asDouble();
    pos[2] = object["position"]["z"].asDouble();

    dir[0] = object["orientation"]["x"].asDouble();
    dir[1] = object["orientation"]["y"].asDouble();
    dir[2] = object["orientation"]["z"].asDouble();

    rad = object["radius"].asDouble();
}

void parse_vector_from_json(const Json::Value &object, vector<double>& vector) {

    for (int i = 0; i < object.size(); ++i) {
        vector.push_back(object[i].asDouble());
    }

}

void parse_circle_from_json(vector<Geometry::Objects> &objects, const Json::Value &object) {
    double pos[3];
    double dir[3];
    double rad;
    parse_pos_dir_rad_from_json(object, pos, dir, rad);

    auto obj = new R3::Circle(&pos[0],&dir[0],rad);
    Geometry::Objects o;
    o.getObject().c = obj;
    o.setType(Geometry::CIRCLE);
    objects.push_back(o);
}


void parse_tube_from_json(vector<Geometry::Objects> &objects, const Json::Value &object) {
    double pos[3];
    double dir[3];
    double rad;
    parse_pos_dir_rad_from_json(object, pos, dir, rad);

    double length = object["length"].asDouble();

    auto obj = new Tube(rad,&pos[0],length);
    Geometry::Objects o;
    o.getObject().t = obj;
    o.setType(Geometry::TUBE);
    objects.push_back(o);
}

void parse_cylinder_from_json(vector<Geometry::Objects> &objects, const Json::Value &object) {
    double pos[3];
    double dir[3];
    double rad;
    parse_pos_dir_rad_from_json(object, pos, dir, rad);

    double length = object["length"].asDouble();


    auto obj = new Cylinder(&pos[0],&dir[0], rad,length);
    Geometry::Objects o;
    o.getObject().cyl = obj;
    o.setType(Geometry::CYLINDER);
    objects.push_back(o);
}

void parse_sphere_from_json(vector<Geometry::Objects> &objects, const Json::Value &object) {
    double pos[3];
    double dir[3];
    double rad;
    parse_pos_dir_rad_from_json(object, pos, dir, rad);

    auto obj = new Sphere(rad,&pos[0]);
    Geometry::Objects o;
    o.getObject().s = obj;
    o.setType(Geometry::SPHERE);
    objects.push_back(o);
}

void parse_geometry_objects_from_json(const Json::Value &chamber_geometry, vector<Geometry::Objects> &geometry_objects) {
    for (auto &json_object : chamber_geometry) {
        if( json_object["type"] == "CIRCLE" ) {
            parse_circle_from_json(geometry_objects, json_object);
        }
        if( json_object["type"] == "TUBE" ) {
            parse_tube_from_json(geometry_objects, json_object);
        }
        if( json_object["type"] == "CYLINDER" ) {
            parse_cylinder_from_json(geometry_objects, json_object);
        }
        if( json_object["type"] == "SPHERE" ) {
            parse_sphere_from_json(geometry_objects, json_object);
        }
    }
}

void parse_test_chamber_from_json(Simulation &s, const Json::Value &simulation) {
    auto& json_chamber = simulation["TESTCHAMBER"];
    auto& chamber = s.getChamber();
    auto& chamber_geometry = json_chamber["geometry"];
    auto& geometry = chamber.getGeometry();
    auto& geometry_objects = geometry.getGeometry_objects();

    parse_geometry_objects_from_json(chamber_geometry, geometry_objects);
}

void parse_thrusters_from_json( Simulation &s, const Json::Value &simulation ) {
    for (const auto &json_thruster : simulation["THRUSTERS"]) {
        double pos[3];
        double dir[3];
        double rad;

        parse_pos_dir_rad_from_json(json_thruster,pos,dir,rad);
        Thruster thruster;
        thruster.setName(json_thruster["name"].asString());
        auto& geometry = json_thruster["geometry"];
        auto& geometry_objects = thruster.getGeometry().getGeometry_objects();

        parse_geometry_objects_from_json(geometry, geometry_objects);
        s.getThrusters().push_back(thruster);
    }
}

Simulation read_simulation_file(string filename){
    Json::Value root;
    ifstream in(filename, ifstream::binary);
    if ( !in.good() ){
        cerr << "could not open file: " << filename << endl;
        exit( EXIT_FAILURE );
    }

    in >> root;

    if ( !root["simulation"].isNull() ) {
        Simulation s;
        auto& simulation = root["simulation"];
        s.setName(simulation["name"].asString());

        parse_test_chamber_from_json(s, simulation);
        parse_thrusters_from_json(s, simulation );
        cout << "number of thrusters " << s.getThrusters().size() << endl;
        return s;
    }

    exit(1);
}

ParticleSource read_source_file(string filename, string source_type){
    Json::Value root;
    ifstream in(filename, ifstream::binary);
    if ( !in.good() ){
        cerr << "could not open file: " << filename <<endl;
        exit( EXIT_FAILURE );
    }

    in >> root;

    if ( !root["sources"].isNull() ) {
        ParticleSource source;
        auto& sources = root["sources"];

        for (auto &json_source : sources) {
            vector<double> current_density_vector;
            vector<double> energy_vector;
            vector<double> fraction_vector;

            if ( json_source["type"] == source_type) {
                source.setType(json_source["type"].asString());
                source.setNumberOfAngleBins(json_source["number_of_angle_bins"].asInt());
                source.setDeltaAngle(json_source["delta_angle_theta"].asDouble());
                source.setDeltaTheta(source.getDeltaAngle());

                parse_vector_from_json(json_source["current_density"], current_density_vector);
                source.setCurrentDensity(current_density_vector);

                parse_vector_from_json(json_source["energy"], energy_vector);
                source.setEnergy(energy_vector);

                parse_vector_from_json(json_source["fraction"], fraction_vector);
                source.setFraction(fraction_vector);

                if( current_density_vector.size() != energy_vector.size() ||
                    current_density_vector.size() != fraction_vector.size()

                    ){
                    assert( 0 == "vector sizes in" + filename + "(" + source_type + "don't match");
                }
                if ( json_source.isMember("constant_phi") ) {
                    source.setConstantPhi(json_source["constant_phi"].asDouble());
                    cout << "found a constant phi -> beam like behaviour!" << endl;
                } else {
                    source.setConstantPhi(NAN);
                }


            } else {
                continue;
            }
        }

        return source;

    }

    exit(1);

}


Yield read_yield(string filename, string wall_type, string feedgas_type) {
    Json::Value root;
    ifstream in(filename, ifstream::binary);
    if ( !in.good() ){
        cerr << "could not open file: " << filename <<endl;
        exit( EXIT_FAILURE );
    }

    in >> root;

    if ( !root["yields"].isNull() ) {
        Yield yield;
        auto& yields = root["yields"];

        for (auto &json_yield : yields) {
            vector<double> energy_vector;
            vector<double> yield_vector;

            if ( json_yield["wall"] == wall_type && json_yield["feedgas"] == feedgas_type ) {
                parse_vector_from_json(json_yield["yield"], yield_vector);
                yield.setYield(yield_vector);

                parse_vector_from_json(json_yield["energy"], energy_vector);
                yield.setEnergy(energy_vector);

                if( yield_vector.size() != energy_vector.size() ){
                    assert( 0 == "vector sizes in" + filename + "(" + wall_type + "," feedgas_type ") don't match");
                }
            } else {
                continue;
            }
        }

        return yield;

    }

    exit(1);

}

