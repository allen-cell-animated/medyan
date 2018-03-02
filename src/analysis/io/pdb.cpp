#include "analysis/io/pdb.h"

#include <iomanip>

namespace medyan {
namespace analysis {

using namespace std;

void PdbGenerator::genModel(int serial) {
    if(!_ofs.is_open()) return;

    _ofs
        << "MODEL "             //  1 -  6  Record name
        << setw(4) << serial    // 11 - 14  serial
        << endl;
}
void PdbGenerator::genEndmdl() {
    if(!_ofs.is_open()) return;

    _ofs << "ENDMDL" << endl;
}
void PdbGenerator::genAtom(
    int serial, string name, char altLoc, string resName, char chainID,
    int resSeq, char iCode, double x, double y, double z, double occupancy, double tempFactor,
    string element, string charge
) {
    if(!_ofs.is_open()) return;

    name.resize(4, ' ');
    resName.resize(3, ' ');
    element.resize(2, ' ');
    charge.resize(2, ' ');

    _ofs
        << "ATOM  "                                             //  1 -  6  Record name
        << setw(5) << serial                                    //  7 - 11  serial
        << ' '
        << name                                                 // 13 - 16  name
        << altLoc                                               // 17       altLoc
        << resName                                              // 18 - 20  resName
        << ' '
        << chainID                                              // 22       chainID
        << setw(4) << resSeq                                    // 23 - 26  resSeq
        << iCode                                                // 27       iCode
        << "   "
        << fixed << setw(8) << setprecision(3) << x             // 31 - 38  x
        << fixed << setw(8) << setprecision(3) << y             // 39 - 46  y
        << fixed << setw(8) << setprecision(3) << z             // 47 - 54  z
        << fixed << setw(6) << setprecision(2) << occupancy     // 55 - 60  occupancy
        << fixed << setw(6) << setprecision(2) << tempFactor    // 61 - 66  tempFactor
        << "          "
        << element                                              // 77 - 78  element
        << charge                                               // 79 - 80  charge
        << endl;

}
void PdbGenerator::genTer(
    int serial, string resName, char chainID,
    int resSeq, char iCode
) {
    if(!_ofs.is_open()) return;

    resName.resize(3, ' ');

    _ofs
        << "TER   "                 //  1 -  6  Record name
        << setw(5) << serial        //  7 - 11  serial
        << "      "
        << resName                  // 18 - 20  resName
        << ' '
        << chainID                  // 22       chainID
        << setw(4) << resSeq        // 23 - 26  resSeq
        << iCode                    // 27       iCode
        << endl;
}



} // namespace analysis
} // namespace medyan
