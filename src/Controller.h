
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Controller_h
#define MEDYAN_Controller_h

#include "common.h"

#include "Output.h"
#include "MController.h"
#include "GController.h"
#include "CController.h"
#include "DRController.h"

//FORWARD DECLARATIONS
class SubSystem;
class Cylinder;
class FilamentBindingManager;

/// Used to initialize, manage, and run an entire simulation.

/*!
 *  The Controller is initialized in the main program, and initializes the SubSystem 
 *  given an initial input directory. After initialization of all member controllers,
 *  the Controller class can run a simulation given the already read input parameters, 
 *  by iteratively switching between mechanical equilibration and stochastic chemical 
 *  steps. The Controller is also responsible for updating positions, reaction rates,
 *  and neighbors lists through its corresponding sub-controllers and the Subsystem.
 */
class Controller {

private:
    double chemistrytime, minimizationtime, neighborlisttime, bindingmanagertime,
        positiontime, rxnratetime;
    string _inputFile; ///< System input file
    
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls

    MController* _mController;   ///< Chemical controller used
    CController* _cController;   ///< Mechanical controller used
    GController* _gController;   ///< Geometry controller used
    DRController* _drController; ///< Dynamic rate controller used
    
    string _inputDirectory;   ///< Input directory being used
    string _outputDirectory;  ///< Output directory being used
    
    vector<Output*> _outputs; ///< Vector of specified outputs
    
    double _runTime;          ///< Total desired runtime for simulation

    double _snapshotTime;     ///< Desired time for each snapshot
    
    double _minimizationTime;  ///< Frequency of mechanical minimization
    double _neighborListTime;  ///< Frequency of neighbor list updates
    
    //@{
    /// Same parameter set as timestep, but in terms of chemical
    /// reaction steps. Useful for small runs and debugging.
    double _runSteps;
    double _snapshotSteps;
    
    double _minimizationSteps;
    double _neighborListSteps;
    ChemistryData _chemData;
    ChemistryAlgorithm _cAlgorithm;
    vector<tuple<short, vector<double>, vector<double>>> fil;
    tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > filaments;
    vector<Compartment*> activatecompartments;
    multimap<int,Compartment*> fCompmap;
    multimap<int,Compartment*> bCompmap;
    //@}
    
    ///INITIALIZATION HELPER FUNCTIONS
    
    /// Set up an initial configuration of a network
    /// For now, only [Bubbles](@ref Bubble) and [Filaments](@ref Filament)
    /// can be initialized before the simulation begins. Any other elements
    /// should be initialized in this function.
    void setupInitialNetwork(SystemParser& p);
    
    /// Setup any special structures needed
    void setupSpecialStructures(SystemParser& p);
    
    ///RUNTIME HELPER FUNCTIONS
    
    /// Move the boundary based on the timestep
    void moveBoundary(double deltaTau);
    ///Activate/deactivate compartments based on the longest filament (along Xaxis).
    void activatedeactivateComp();
    void ControlfrontEndCompobsolete();
    void ControlbackEndCompobsolete();
    void ControlfrontbackEndComp();
    /// Update the positions of all elements in the system
    void updatePositions();
    
#ifdef DYNAMICRATES
    /// Update the reaction rates of all elements in the system
    void updateReactionRates();
#endif
    
    /// Update neighbors lists, called in run
    void updateNeighborLists();

    /// Execute any special protocols needed, for example,
    /// making Linker and Filament species static
    void executeSpecialProtocols();

    /// Reset counters on all elements in the system
    void resetCounters();

    ///Helper function to pin filaments near the boundary
    void pinBoundaryFilaments();
    //Qin
    void pinLowerBoundaryFilaments();
    
public:
    Controller(SubSystem* s);
    ~Controller() {};
    
    ///Initialize the system, given an input and output directory
    void initialize(string inputFile,
                    string inputDirectory,
                    string outputDirectory);
    ///Run the simulation
    void run();
};

#endif
