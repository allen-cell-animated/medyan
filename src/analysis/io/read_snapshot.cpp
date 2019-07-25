#include "read_snapshot.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "utility.h"

#include "OutputStruct.hpp"

#include "analysis/io/pdb.h"
#include "Core/Globals.hpp"
#include "Util/Io/Log.hpp"
#include "Structure/Bead.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Edge.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"

namespace medyan {
namespace analysis {

namespace {
    /// Helper struct to find number of beads (atoms) required for pdb-style output
    struct PdbMaxBead {
        std::vector<size_t> filament; // Maximum beads in each filament
        size_t linker = 0;
        size_t motor = 0;
        size_t membrane = 0;

        size_t maxBead = 0;

        void renew(const OutputStructSnapshot& snapshot) {
            size_t curMaxBead = 0;

            // Beads in filaments
            for(auto& eachFilament: snapshot.filamentStruct) {
                int curId = eachFilament.getId();
                if(curId >= filament.size()) {
                    filament.resize(curId + 1, 0);
                }
                if(filament[curId] < eachFilament.getNumBeads()) filament[curId] = eachFilament.getNumBeads();
                curMaxBead += filament[curId];
            }

            // Beads in linkers
            size_t newLinker = snapshot.linkerStruct.size() * 2;
            if(newLinker > linker) linker = newLinker;
            curMaxBead += linker;

            // Beads in motors
            size_t newMotor = snapshot.motorStruct.size() * 2;
            if(newMotor > motor) motor = newMotor;
            curMaxBead += motor;

            // Beads in membranes
            size_t newMembrane = 0;
            for(auto& eachMembrane: snapshot.membraneStruct) {
                /*
                newMembrane += eachMembrane.getNumEdges() * 2;
                */
                newMembrane += eachMembrane.getNumVertices();
            }
            if(newMembrane > membrane) membrane = newMembrane;
            curMaxBead += membrane;

            // Total beads
            if(curMaxBead > maxBead) maxBead = curMaxBead;
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
    PdbGenerator pg(_pdbFilepath);
    PsfGenerator psfGen(_psfFilepath); // To generate bond connections for a fixed topology

    LOG(STEP) << "Writing to " << _pdbFilepath << " and " << _psfFilepath;
    psfGen.genHeader();
    psfGen.genNatom(maxBead.maxBead);

    const size_t bondFrame = global().analyzeMembraneBondFrame;

    size_t numSnapshots = snapshots.size();
    for(size_t idx = 0; idx < numSnapshots; ++idx) {
        if (idx % 20 == 0) LOG(INFO) << "Generating model " << idx;
        pg.genModel(idx + 1);

        char chain;
        size_t atomSerial; // Pdb file atom number
        size_t atomId = 0; // Psf file atom number. TER does not increase this variable.
        size_t atomCount; // Temporary counter
        size_t resSeq;

        SubSystem s; // Dummy subsystem

        // Filaments
        chain = 'F';
        atomSerial = 0;
        resSeq = 0;
        auto filamentPtr = snapshots[idx].filamentStruct.begin();
        size_t filamentAlloc = maxBead.filament.size();

        std::vector< std::pair< size_t, size_t >> bondList;

        for(size_t i = 0; i < filamentAlloc; ++i) {
            atomCount = 0;

            if(filamentPtr != snapshots[idx].filamentStruct.end() && filamentPtr->getId() == i) {
                // Filament exists for this id.
                auto& allCoords = filamentPtr->getCoords();
                for(auto& eachCoord: allCoords) {
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                        eachCoord[0], eachCoord[1], eachCoord[2]
                    );
                    if(idx == bondFrame) {
                        psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }
                while(atomCount < maxBead.filament[i]) {
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                        allCoords.back()[0], allCoords.back()[1], allCoords.back()[2]
                    );
                    if(idx == bondFrame) {
                        psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }

                ++filamentPtr;
            } else {
                // Filament does not exist for this id.
                while(atomCount < maxBead.filament[i]) {
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq
                    );
                    if(idx == bondFrame) {
                        psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }
            }

            ++resSeq; // Separate trace.
        }
        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // Linkers
        chain = 'L';
        atomCount = 0;
        resSeq = 0;
        for(auto& eachLinker: snapshots[idx].linkerStruct) {
            for(auto& eachCoord: eachLinker.getCoords()) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
                if(idx == bondFrame) {
                    psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.linker) {
            for(size_t b = 0; b < 2; ++b) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq
                );
                if(idx == bondFrame) {
                    psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq; // Separate bonds
        }
        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // Motors
        chain = 'M';
        atomCount = 0;
        resSeq = 0;
        for(auto& eachMotor: snapshots[idx].motorStruct) {
            for(auto& eachCoord: eachMotor.getCoords()) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
                if(idx == bondFrame) {
                    psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.motor) {
            for(size_t b = 0; b < 2; ++b) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq
                );
                if(idx == bondFrame) {
                    psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq;
        }
        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // Membranes
        chain = 'E';
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
                    ++atomId;
                    ++atomCount;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, atomSerial, ' ',
                        v[0], v[1], v[2]
                    );
                }
                ++atomSerial;
            }

            */

            size_t atomIdSoFar = atomId;

            for(const auto& v: eachMembrane.getMembraneInfo().attributeInitializerInfo.vertexCoordinateList) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, 1/* resSeq */, ' ',
                    v[0], v[1], v[2]
                );
                if(idx == bondFrame) {
                    psfGen.genAtom(atomId, "", 1/* resSeq */, "ARG", "CA", "CA");
                }
            }
            ++resSeq;

            if(idx == bondFrame) {
                const auto& tlist = eachMembrane.getMembraneInfo().triangleVertexIndexList;
                for(const auto& t : tlist) {
                    for(size_t i = 0; i < 3; ++i) {
                        size_t i_next = (i+1) % 3;
                        if(t[i] < t[i_next]) {
                            bondList.emplace_back(
                                atomIdSoFar + t[i] + 1,
                                atomIdSoFar + t[i_next] + 1
                            );
                        }
                    }
                }
            }

        }

        while(atomCount < maxBead.membrane) {
            ++atomSerial;
            ++atomId;
            ++atomCount;
            ++resSeq;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
            if(idx == bondFrame) {
                psfGen.genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
            }
        }

        // Generating bond
        if(idx == bondFrame) {
            size_t numBonds = bondList.size();
            psfGen.genNbond(numBonds);
            psfGen.genBondStart();

            for(const auto& bond : bondList) {
                psfGen.genBond(bond.first, bond.second);
            }

            psfGen.genBondEnd();
            LOG(INFO) << "Bond info generated.";
        }

        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // End of model
        pg.genEndmdl();

    } // End of looping through snapshots
    LOG(STEP) << "Writing complete. " << numSnapshots << " models created.";

}

} // namespace analysis
} // namespace medyan
