#include "read_snapshot.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "utility.h"

#include "OutputStruct.h"

#include "analysis/io/pdb.h"
#include "core/log.h"
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
            ++newFilament; // For TER
            if(newFilament > filament) filament = newFilament;

            // Beads in linkers
            size_t newLinker = snapshot.linkerStruct.size() * 2 + 1; // +1 TER
            if(newLinker > linker) linker = newLinker;

            // Beads in motors
            size_t newMotor = snapshot.motorStruct.size() * 2 + 1; // +1 TER
            if(newMotor > motor) motor = newMotor;

            // Beads in membranes
            size_t newMembrane = 0;
            for(auto& eachMembrane: snapshot.membraneStruct) {
                /*
                newMembrane += eachMembrane.getNumEdges() * 2;
                */
                newMembrane += eachMembrane.getNumVertices();
            }
            ++newMembrane; // +1 TER
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

    LOG(STEP) << "Start reading " << _snapshotFilepath;
    string line;
    while(maxFrames == 0 || curFrame < maxFrames) {
        getline(is, line);
        if(!is) break;
        if(line.empty()) continue;

        ++curFrame;
        if (curFrame % 20 == 0) LOG(INFO) << "Frame " << curFrame;

        istringstream iss(line);
        snapshots.emplace_back(0);
        snapshots.back().getFromOutput(is, iss);

        maxBead.renew(snapshots.back());
    }
    LOG(STEP) << "Reading complete. " << curFrame << " frames to be processed.";

    // close input stream
    is.close();

    // Write to pdb
    PdbGenerator pg(_vmdFilepath);

    LOG(STEP) << "Writing to " << _vmdFilepath;
    // Note: in the output, all the coordinates would be 1/10 in size.
    size_t numSnapshots = snapshots.size();
    for(size_t idx = 0; idx < numSnapshots; ++idx) {
        if (idx % 20 == 0) LOG(INFO) << "Generating model " << idx;
        pg.genModel(idx + 1);

        char chain;
        size_t atomSerial;
        size_t atomCount;
        size_t resSeq;

        SubSystem s; // Dummy subsystem

        // Filaments
        chain = 'F';
        atomSerial = 0;
        atomCount = 0;
        resSeq = 0;
        for(auto& eachFilament: snapshots[idx].filamentStruct) {
            for(auto& eachCoord: eachFilament.getCoords()) {
                ++atomSerial;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.filament) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
        }
        pg.genTer(++atomSerial, "ARG", chain, ++resSeq);

        // Linkers
        chain = 'L';
        atomSerial = 0;
        atomCount = 0;
        resSeq = 0;
        for(auto& eachLinker: snapshots[idx].linkerStruct) {
            for(auto& eachCoord: eachLinker.getCoords()) {
                ++atomSerial;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.linker) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
        }
        pg.genTer(++atomSerial, "ARG", chain, ++resSeq);

        // Motors
        chain = 'M';
        atomSerial = 0;
        atomCount = 0;
        resSeq = 0;
        for(auto& eachMotor: snapshots[idx].motorStruct) {
            for(auto& eachCoord: eachMotor.getCoords()) {
                ++atomSerial;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.motor) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
        }
        pg.genTer(++atomSerial, "ARG", chain, ++resSeq);

        // Membranes
        chain = 'E';
        atomSerial = 0;
        atomCount = 0;
        resSeq = 0;
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
                ++resSeq;
                auto& v = get<0>(vInfo);
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    v[0], v[1], v[2]
                );
            }
            ++resSeq;
        }
        while(atomCount < maxBead.membrane) {
            ++atomSerial;
            ++atomCount;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
        }
        pg.genTer(++atomSerial, "ARG", chain, ++resSeq);

        // End of model
        pg.genEndmdl();

    } // End of looping through snapshots
    LOG(STEP) << "Writing complete. " << numSnapshots << " models created.";

}

} // namespace analysis
} // namespace medyan
