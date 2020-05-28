
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Parser_h
#define MEDYAN_Parser_h

#include <vector>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <ios>

#include "common.h"
#include "SysParams.h"
#include "utility.h"


/// A general parser
/*!
 *  A parser object, when initialized, opens an input file. Upon destruction, it 
 *  closes the file.
 */
class Parser {
protected:
    fstream _inputFile; ///< input file being used
    
public:
    Parser(string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            cout << "There was an error parsing file " << inputFileName
                 << ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Loading file " << inputFileName << endl;
    }
    ~Parser() {_inputFile.close();}
};

/// To parse a system input file, initialized by the Controller.
struct SystemParser {
    
    //@{
    /// Parameter parser. Reads input directly into system parameters
    /// @note - does not check for correctness and consistency here.
    static MechParams    readMechParams(std::istream&);
    static ChemParams    readChemParams(std::istream&);
    static GeoParams     readGeoParams(std::istream&);
    static BoundParams   readBoundParams(std::istream&);
    static DyRateParams  readDyRateParams(std::istream&);
    static SpecialParams readSpecialParams(std::istream&);
    //@}
    
    //@{
    /// Algorithm parser
    static MechParams::MechanicsAlgorithm readMechanicsAlgorithm(std::istream&);
    static ChemParams::ChemistryAlgorithm readChemistryAlgorithm(std::istream&);
    //@}
    
    //@{
    /// Type parser
    static MechParams::MechanicsFFType     readMechanicsFFType(std::istream&);
    static DyRateParams::DynamicRateType   readDynamicRateType(std::istream&);
    static BoundParams::BoundaryType       readBoundaryType(std::istream&);
    static SpecialParams::SpecialSetupType readSpecialSetupType(std::istream&);
    //@}
    
    /// Read Filament information
    static FilamentSetup readFilamentSetup(std::istream&);
    
    /// Read Bubble information
    static BubbleSetup readBubbleSetup(std::istream&);
    
    /// Chemistry information
    static ChemParams::ChemistrySetup readChemistrySetup(std::istream&);
};

/// Used to parse initial Filament information, initialized by the Controller.
struct FilamentParser {
  
    /// Reads filament input file. Returns a vector of tuples containing
    /// filament type and positions (start and end points).
    /// @note - Does not check for coordinate correctness.
    static FilamentData readFilaments(std::istream&);
};

/// Used to parse initial Bubble information, initialized by the Controller.
struct BubbleParser {    
    /// Reads bubble input file. Returns a vector of tuples containing
    /// bubble type and position.
    /// @note - Does not check for coordinate correctness.
    static BubbleData readBubbles(std::istream&);
};


/// Used to parse all chemical information, initialized by the Controller.
struct ChemistryParser {
    /// Reads chemical reactions and species from input file. Returns a
    /// ChemistryData struct containing this data
    /// @note - this does not check for consistency and correctness, the only
    ///         sanity check here is that there are no duplicate species names.
    static ChemistryData readChemistryInput(std::istream&);
};


/// Used to parse pin positions if needed upon restart
class PinRestartParser: public Parser {
    
public:
    PinRestartParser(string inputFileName) : Parser(inputFileName) {}
    ~PinRestartParser() {}
    
    /// Reads pin positions from file, and sets filaments
    void resetPins();
};



#endif
