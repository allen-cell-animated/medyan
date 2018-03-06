#ifndef MEDYAN_ANALYSIS_IO_PDB_H
#define MEDYAN_ANALYSIS_IO_PDB_H

#include <iostream>
#include <fstream>
#include <string>

namespace medyan {
namespace analysis {

// Checkout www.wwpdb.org/

class PdbGenerator {
    std::string _pdbFilepath;
    std::ofstream _ofs; ///< Handles the file output

public:
    PdbGenerator(const std::string& filepath): _pdbFilepath(filepath), _ofs(filepath) {}

    void genModel(int serial);
    void genEndmdl();
    void genAtom(int serial, std::string name, char altLoc, std::string resName, char chainID,
        int resSeq, char iCode=' ', double x=0.0, double y=0.0, double z=0.0, double occupancy=1.0, double tempFactor=0.0,
        std::string element="  ", std::string charge="  ");
    void genTer(int serial, std::string resName, char chainID,
        int resSeq, char iCode=' ');

};


} // namespace analysis
} // namespace medyan

#endif