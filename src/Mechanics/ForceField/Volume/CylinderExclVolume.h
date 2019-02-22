
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CylinderExclVolume_h
#define MEDYAN_CylinderExclVolume_h

#include "common.h"

#include "CylinderVolumeInteractions.h"
#include "NeighborListImpl.h"
#include "HybridNeighborListImpl.h"

#include "SysParams.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

//FORWARD DECLARATIONS
class Cylinder;

/// Represents an excuded volume interaction between two [Cylinders](@ref Cylinder).
template <class CVolumeInteractionType>
class CylinderExclVolume : public CylinderVolumeInteractions {
    
private:
    CVolumeInteractionType _FFType;
    CylinderCylinderNL* _neighborList;  ///< Neighbor list of cylinders
#ifdef HYBRID_NLSTENCILLIST
    HybridCylinderCylinderNL* _HneighborList;
    short _HnlID;
#endif
    ///Array describing the constants in calculation
    int *beadSet;
    double *krep;
    int nint = 0;
#ifdef CUDAACCL
    int * gpu_beadSet = NULL;
    double * gpu_krep = NULL;
    int * gpu_params = NULL;
    CUDAvars cvars;
    double *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    ///Array describing indexed set of interactions
    ///For volume, this is a 4-bead potential
    const static int n = 4;
    static int numInteractions;
    
    ///Constructor
    CylinderExclVolume() {
#ifndef HYBRID_NLSTENCILLIST
        _neighborList = new CylinderCylinderNL(SysParams::Mechanics().VolumeCutoff);
#endif
#ifdef CUDAACCL_NL
        _neighborList->cudacpyforces = true;
#endif
    }
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(double *coord, double *f);

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}

    virtual const string getName() {return "Cylinder Excluded Volume";}

#ifdef HYBRID_NLSTENCILLIST

    virtual void setHNeighborList(HybridCylinderCylinderNL* Hnl) {
        _HneighborList = Hnl;
        _HnlID = Hnl->setneighborsearchparameters(0,0,true,false,SysParams::Mechanics()
                .VolumeCutoff,0.0);
    };

    virtual HybridCylinderCylinderNL* getHNeighborList(){return _HneighborList;};
#endif
};

#endif
