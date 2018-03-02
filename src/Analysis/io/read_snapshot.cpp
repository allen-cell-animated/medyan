#include "read_snapshot.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "utility.h"

#include "OutputStruct.h"

#include "analysis/io/pdb.h"
#include "Structure/Bead.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Edge.h"
#include "Structure/SurfaceMesh/Membrane.h"
#include "Structure/SurfaceMesh/Vertex.h"

namespace medyan {
namespace analysis {

namespace {
    /// Helper struct to find number of beads (atoms) required for pdb-style output
    struct PdbMaxBead {
        size_t filament = 0;
        size_t linker = 0;
        size_t motor = 0;
        size_t membrane = 0;

        void renew(const OutputStructSnapshot& snapshot) {
            // Beads in filaments
            size_t newFilament = 0;
            for(auto& eachFilament: snapshot.filamentStruct) {
                newFilament += eachFilament.getNumBeads();
            }
            if(newFilament > filament) filament = newFilament;

            // Beads in linkers
            size_t newLinker = snapshot.linkerStruct.size() * 2;
            if(newLinker > linker) linker = newLinker;

            // Beads in motors
            size_t newMotor = snapshot.motorStruct.size() * 2;
            if(newMotor > motor) motor = newMotor;

            // Beads in membranes
            size_t newMembrane = 0;
            for(auto& eachMembrane: snapshot.membraneStruct) {
                /*
                newMembrane += eachMembrane.getNumEdges() * 2;
                */
                newMembrane += eachMembrane.getNumVertices();
            }
            if(newMembrane > membrane) membrane = newMembrane;
        }
    };

}

void SnapshotReader::readAndConvertToVmd(const size_t maxFrames) {
    using namespace std;

    vector<OutputStructSnapshot> snapshots;

    // Read snapshot
    ifstream is(_snapshotFilepath);

    PdbMaxBead maxBead;

    size_t curFrame = 0;

    cout << "Start reading " << _snapshotFilepath << endl;
    string line;
    while(maxFrames == 0 || curFrame < maxFrames) {
        getline(is, line);
        if(!is) break;
        if(line.empty()) continue;

        ++curFrame;
        if (curFrame % 10 == 0)cout << "Frame " << curFrame << endl;

        istringstream iss(line);
        snapshots.emplace_back(0);
        snapshots.back().getFromOutput(is, iss);

        maxBead.renew(snapshots.back());
    }
    cout << "Reading complete. " << curFrame << " frames to be processed." << endl;

    // close input stream
    is.close();

    // Write to pdb
    PdbGenerator pg(_vmdFilepath);

    cout << "Writing to " << _vmdFilepath << endl;
    // Note: in the output, all the coordinates would be 1/10 in size.
    size_t numSnapshots = snapshots.size();
    for(size_t idx = 0; idx < numSnapshots; ++idx) {
        if (idx % 10 == 0) cout << "Generating model " << idx << endl;
        pg.genModel(idx + 1);

        char chain;
        size_t atomSerial;
        size_t atomCount;

        SubSystem s; // Dummy subsystem

        // Filaments
        chain = 'F';
        atomSerial = 0;
        atomCount = 0;
        for(auto& eachFilament: snapshots[idx].filamentStruct) {
            for(auto& eachCoord: eachFilament.getCoords()) {
                ++atomSerial;
                ++atomCount;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++atomSerial; // Separate bonds
        }
        while(atomCount < maxBead.filament) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                0.0, 0.0, 0.0
            );
        }

        // Linkers
        chain = 'L';
        atomSerial = 0;
        atomCount = 0;
        for(auto& eachLinker: snapshots[idx].linkerStruct) {
            for(auto& eachCoord: eachLinker.getCoords()) {
                ++atomSerial;
                ++atomCount;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++atomSerial; // Separate bonds
        }
        while(atomCount < maxBead.linker) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                0.0, 0.0, 0.0
            );
        }

        // Motors
        chain = 'M';
        atomSerial = 0;
        atomCount = 0;
        for(auto& eachMotor: snapshots[idx].motorStruct) {
            for(auto& eachCoord: eachMotor.getCoords()) {
                ++atomSerial;
                ++atomCount;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++atomSerial; // Separate bonds
        }
        while(atomCount < maxBead.motor) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                0.0, 0.0, 0.0
            );
        }

        // Membranes
        chain = 'E';
        atomSerial = 0;
        atomCount = 0;
        for(auto& eachMembrane: snapshots[idx].membraneStruct) {
            /*
            bool buildMembrane = !eachMembrane.getMembrane();
            Membrane* dataMembrane;

            unique_ptr<Membrane> newMembrane;
            if(buildMembrane) {
                newMembrane = make_unique<Membrane>(&s, 0, eachMembrane.getMembraneInfo());
                dataMembrane = newMembrane.get();
            } else {
                dataMembrane = eachMembrane.getMembrane();
            }

            // Output by edges
            for(Edge* e: dataMembrane->getEdgeVector()) {
                for(Vertex* v: e->getVertices()) {
                    ++atomSerial;
                    ++atomCount;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                        v[0], v[1], v[2]
                    );
                }
                ++atomSerial;
            }

            */
            for(auto& vInfo: eachMembrane.getMembraneInfo()) {
                ++atomSerial;
                ++atomCount;
                auto& v = get<0>(vInfo);
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                    v[0], v[1], v[2]
                );
            }
            ++atomSerial;
        }
        while(atomCount < maxBead.membrane) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                0.0, 0.0, 0.0
            );
        }

        // End of model
        pg.genEndmdl();

    } // End of looping through snapshots
    cout << "Writing complete. " << numSnapshots << " models created." << endl;

}

} // namespace analysis
} // namespace medyan
