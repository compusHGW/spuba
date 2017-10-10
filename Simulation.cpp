#include <chrono>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <json/value.h>
#include <json/reader.h>

#include "Simulation.h"
#include "input_json.h"
#include "output.h"
#include "utility.hpp"

namespace {
    auto setup_cones(std::vector<Cone3D> &cones, double tubeh, double tubel) {

        double coneh = tubeh / 14.;
        double alf = 10. * M_PI / 180.;       // shevron tilding

        cout << "Baffles for planar bottom cap implemented (alpha = "<< alf * 180 / M_PI << "° , coneh = " << coneh << "m )" << endl;

        double coneh_z = coneh; // cone length (on z-axis)
        double Rl = Rc - coneh * tan(alf);    // radius of cone (<Rc)

        for (int i = 0; i < number_of_cones_in_z; ++i) {
            vector3d cone_pos = {{0, 0, (i + 1) * coneh_z}};
            cones.emplace_back(cone_pos, coneh, Rl, alf);
        }

        for (int i = 0; i < number_of_cones_in_r; ++i) {
            vector3d cone_pos = {{0, 0, tubel}};
            cones.emplace_back(cone_pos, coneh, i * coneh * tan(alf), alf);
        }

        return cones;
    }
}
void Simulation::calculate_reflection(Particle3D& p, Particle3D& reflect) {

    auto RNG = [&](){
        auto ret = dis(rng);
        return ret;
    };

    reflect.SetParticle(p.GetHitpoint(), asin(sqrt(RNG())), RNG() * 2 * M_PI, p.GetVelocity());

    double tmpMatrix[3][3];
    double identity[3][3];
    /*-----------------------------------------------------------------------------
     *  Turn Coordinate System into the direction of the normal
     *-----------------------------------------------------------------------------*/
    /*----From Cartesian to spherical coordinates in order to produce matrixes to turn the vertices--*/
    double phi, theta;
    double normalPolar[3];
    CartesianToPolar(p.GetHitNormal(), normalPolar);
    phi = normalPolar[1];
    theta = normalPolar[0];

    /*----Create matrix matrix and multiply with rotation matrix Y and X----------*/
    set_identity(identity);
    RotateZ(tmpMatrix, identity, phi);
    RotateY(tmpMatrix, identity, theta);

    /*----Multiply vector with matrix to produce new vector-------------------------*/
    double vect[3] = {0.0, 0.0, 0.0};
    double nnn1[3] = {0.0, 0.0, 0.0}; // Nullvektor
    double nnn[3] = {0.0, 0.0, 1.0};  // Einheitsvektor z
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            vect[j] += reflect.GetDirection()[i] * identity[j][i];
            nnn1[j] += nnn[i] * identity[j][i];
        }
    }
    norm(vect);

    reflect.SetParticleCartesian(p.GetHitpoint(), vect, p.GetVelocity());
    //-------------------------------------------------------------------------
}

void Simulation::read_input_file(string filename){
    Json::Value root;
    ifstream in(filename, ios_base::binary);
    if ( !in.good() ){
        cerr << "could not open file: " << filename <<endl;
        exit( EXIT_FAILURE );
    }

    in >> root;

    Simulation::nParticles = root["number_of_particles"].asInt();
    Simulation::job_id = root["job_id"].asInt();
    Simulation::comment = root["comment"].asString();
    Simulation::wall = root["wall_material"].asString();
    Simulation::feedgas = root["thruster_feedgas"].asString();
    Simulation::baffles_flag = root["baffles"].asBool();
    Simulation::source_type = root["source_type"].asString();



}

void Simulation::run() {

    auto RNG = [&]() {
        auto ret = this->dis(this->rng);
        return ret;
    };

    begin_time = std::chrono::high_resolution_clock::now();

    cout << "\n Different rand num? " << RNG() << "\n" << endl;

    Yield y = read_yield("/home/guba/CLionProjects/spuba/input.json", wall, feedgas);

    cout << "Feedgas: " << feedgas << endl;
    cout << "Wall Material: " << wall << endl;

    // TODO someone should find out what this means ;)
    auto tubel = chamber.getGeometry().getGeometry_objects().front().object.c->getPosition()[2];
    Circle wall_sp(OrigRs, 0.0, Rs);
    Segment2D wall_seg(0.0, Rc, tubel, Rc);
    wall_sp.FindInt(wall_seg);
    Vector2D v = wall_sp.GetPoint(1) - wall_sp.GetOrigin();
    double fa = atan2(v.x[1], v.x[0]);
    double df = fa / (NUM_EM_SP_SEGS);
    double dthES = df;

    // TODO make cones an json-input
    if ( baffles_flag ) {
    auto tubeh = chamber.getGeometry().getGeometry_objects()[1].object.t->GetLength();
    setup_cones(chamber.getBaffels(),
                tubeh,
                tubel
    );
    }


    for( auto& thruster : thrusters ){
        thruster.init_all_bins();
    }

    /* Source */
    cout << "Initializing source..." << endl;
    ParticleSource source = read_source_file("/home/guba/CLionProjects/spuba/source.json", source_type);
    cout << "Source with type: \"" << source.getType() << "\" loaded" << endl;

    source.init("0");

    // TODO: make source_position an input
    vector<double> source_position = {0.0, 0.0, 0.0}; //emission source
    for (int i = 0; i < 3; ++i) source_position[i] = thrusters[0].getGeometry().getGeometry_objects().front().object.cyl->getTop()->getPosition()[i];
    cout << "Emission Source at (" << source_position[0] << " , " << source_position[1] << " , " << source_position[3] << " )" << std::endl;
    source.setPosition(source_position); // thruster_01 position
    cout << "Total Flux " << source.GetTotIonFl() << endl;

    cout << std::endl;

    T t;
    /********************************************************************************/

    cout << "-------------- Main cycle started ----------------------" << std::endl;
    for (int i = 0; i < nParticles; i++) {

        for (int j = 0; j < SUBC; j++) {
            if ( (i % 1000000 == 0) && j == 0){
                cout << " " << i << " R x" << NUM_REDEP << " S x" << SUBC << " particles were calculated !!!\n";
            }

            Particle3D p = source.gen_ion3D(rng);
            p.SetENumB(-1);

            fillEmissionBins1D(nBins, thrusters.front().bins.d1.emissionBins, p);

            chamber.intersect(p,baffles_flag,t,dthES);

            if (p.GetType() != DIED) {
                Particle3D reflect(p);
                for (int c = 0; c < NUM_REDEP; c++) {

                    calculate_reflection(p,reflect);

                    double weight = y.interpolate_yield((p.GetEnergy()));
                    if (p.GetENum() >= 0 && p.GetENum() < t.l_v_nReDepCylSeg.size())
                        t.l_v_nReDepCylSeg[p.GetENum()] += weight;
                    if (p.GetENumS() >= 0 && p.GetENumS() < t.l_v_nReDepSpSeg.size())
                        t.l_v_nReDepSpSeg[p.GetENumS()] += weight;

                    if( chamber.intersect_baffels(reflect, baffles_flag) ) {
                        reflect.SetType(DIED);
                    }

                    // Collision of reflected particles with walls
                    if (reflect.GetType() != DIED) {

                        for (int id = 0; id < thrusters.size(); ++id) {
                            auto& thruster = thrusters[id];
                            if( thruster.Intersect(reflect,weight,dthES) ) break;
                        }

                    }
                } // end: for (int c = 0; c < NUM_REDEP; c++)

            } // end: if (p.GetType() != DIED) to generate particle

        } // end: for (unsigned j = 0; j < SUBC; j++)
    } // end: for (unsigned i = 0; i < nParticles; i++)

    cout << "-------------- Main cycle ended ----------------------\n" << std::endl;

}


void Simulation::output() {/***********************************************************************************/

    cout << "Write in files:" << endl;
    cout << "nParticles " << nParticles << endl;
    outputEmission(OUT_FOLDER + "/emission_distribution_"s + to_string(job_id) + ".out", nParticles, nBins, this->thrusters.front().bins.d1.emissionBins);

    // TODO change most of the particle counter to weight counters / output everything in SI units / output depending on geometry
    for (auto& thruster : this->thrusters) {
        Counter& counter = thruster.counter;
        Bins& bin = thruster.bins;
        string thruster_name=thruster.getName();
        outputAbsorbtion(OUT_FOLDER + "/deposition_distribution_"s + thruster_name + "_" + to_string(job_id) + ".out", counter.circleCounter, nParticles, nBins, bin.d1.absorbtionBins);
        outputAbsorbtion2DBot(OUT_FOLDER + "/deposition_distribution_Exit_"s + thruster_name + "_" + to_string(job_id) + ".out", counter.exitCounter, nParticles, nBins2Dphi, nBins2Dr, bin.d2.absorbtion.exit);
        outputAbsorbtion2DChan(OUT_FOLDER + "/deposition_distribution_Chan_"s + thruster_name + "_" + to_string(job_id) + ".out", counter.channelCounter, nParticles, nBins2Dphi, nBins2Dz, bin.d2.absorbtion.chan,
                               thruster.getGeometry().getGeometry_objects().front().object.cyl->GetRadius());
        outputAbsorbtion2DBot(OUT_FOLDER + "/deposition_distribution_Bot_"s + thruster_name + "_" + to_string(job_id) + ".out", counter.bottomCounter, nParticles, nBins2Dphi, nBins2Dr, bin.d2.absorbtion.bot);
        cout << thruster_name << ": ";
        cout << "exit plane was hit by "<< counter.circleCounter << " particles (exitCounter=" << counter.exitCounter << ")";
        cout << " -> thr. channel " << counter.channelCounter << " bottom: " << counter.bottomCounter << " particles" << endl;
    }

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - begin_time);
    cout << "\nTotal run time ~ " << duration.count() << " seconds" << endl;
}

