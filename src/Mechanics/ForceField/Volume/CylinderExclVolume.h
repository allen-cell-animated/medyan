
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
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    HybridCylinderCylinderNL* _HneighborList;
    short _HnlID;
#endif
    ///Array describing the constants in calculation
    int *beadSet;
    floatingpoint *krep;
    int nint = 0;
#ifdef CUDAACCL
    int * gpu_beadSet = NULL;
    floatingpoint * gpu_krep = NULL;
    int * gpu_params = NULL;
    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    ///Array describing indexed set of interactions
    ///For volume, this is a 4-bead potential
    const static int n = 4;
    static int numInteractions;

    ///Constructor
    CylinderExclVolume() {
        //If Hybrid NeighborList is not preferred, neighborList is created using Original
        // framework.
#if !defined(HYBRID_NLSTENCILLIST) || !defined(SIMDBINDINGSEARCH)
        _neighborList = new CylinderCylinderNL(SysParams::Mechanics().VolumeCutoff);
#endif
#ifdef CUDAACCL_NL
        _neighborList->cudacpyforces = true;
#endif
    }

    virtual void vectorize();
    virtual void deallocate();

    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}

    virtual const string getName() {return "Cylinder Excluded Volume";}

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)

    virtual void setHNeighborList(HybridCylinderCylinderNL* Hnl) {
        _HneighborList = Hnl;
        _HnlID = Hnl->setneighborsearchparameters(0,0,true,false,SysParams::Mechanics()
                .VolumeCutoff,0.0);
    };

    virtual HybridCylinderCylinderNL* getHNeighborList(){return _HneighborList;};
#endif
};

#endif
