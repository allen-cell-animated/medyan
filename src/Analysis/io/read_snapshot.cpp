#include "read_snapshot.hpp"

#include <fstream>
#include <sstream>
#include <vector>

#include "OutputStruct.h"

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
                newMembrane += eachMembrane.getNumEdges() * 2;
            }
            if(newMembrane > membrane) membrane = newMembrane;
        }
    };
}

void SnapshotReader::ReadAndConvertToVmd() {
    std::vector<OutputStructSnapshot> snapshots;
    OutputStructSnapshot dummyMaxHolder(0);

    // Read from output file
    std::ifstream is(_snapshotFilepath);

    PdbMaxBead maxBead;

    std::string line;
    while(true) {
        std::getline(is, line);
        if(!is) break;
        if(line.empty()) continue;

        std::istringstream iss(line);
        snapshots.emplace_back(0);
        snapshots.back().getFromOutput(is, iss);

        maxBead.renew(snapshots.back());
    }

    is.close();

    // Construct the vmd version
    std::ofstream os(_vmdFilepath);

    // TODO: output

    os.close();

}
