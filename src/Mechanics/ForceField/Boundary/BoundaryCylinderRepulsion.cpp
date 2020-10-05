
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

#include "BoundaryCylinderRepulsion.h"

#include <algorithm> // max

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryElement.h"
#include "BoundaryElementImpl.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"
#include "cross_check.h"
#include "CUDAcommon.h"
#include "BoundaryCylinderRepulsionCUDA.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

using namespace mathfunc;

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {

    //count interactions
    nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
    {

        for(auto &c : _neighborList->getNeighbors(be))
        {
            if(c->isMinusEnd()) nint++;
            nint++;
        }
    }
    CUDAcommon::tmin.numinteractions[9] += nint;
    beadSet = new int[n * nint];
    krep = new floatingpoint[nint];
    slen = new floatingpoint[nint];
    auto beList = BoundaryElement::getBoundaryElements();

    int nbe = BoundaryElement::getBoundaryElements().size();
    int i = 0;
    int ni = 0;
    int bindex = 0;

    nneighbors = new int[nbe];//stores number of interactions per boundary element.
    floatingpoint *beListplane;
    int *nintvec;
    beListplane = new floatingpoint[4 * nbe];
    nintvec = new int[nbe];//stores cumulative number of nneighbors.

    int cumnn=0;
    for (i = 0; i < nbe; i++) {

        auto be = BoundaryElement::getBoundaryElements()[i];//beList[i];
        auto nn = _neighborList->getNeighbors(be).size();

        nneighbors[i] = 0;
        auto idx = 0;

        for (ni = 0; ni < nn; ni++) {
//            auto check=false;
//            auto neighbor = _neighborList->getNeighbors(be)[ni];
            /*std::cout<<"Boundary with cindex "<<neighbor->_dcIndex<<" and ID "
                    ""<<neighbor->getID()<<" with bindices "<<neighbor->getFirstBead()
                    ->getIndex()<<" "<<neighbor->getSecondBead()->getIndex()<<endl;*/
            if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
            {
                bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->getIndex() * 3 + si.bead;
                beadSet[cumnn+idx] = bindex;
                krep[cumnn+idx] = be->getRepulsionConst();
                slen[cumnn+idx] = be->getScreeningLength();
                idx++;
            }
            bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->getIndex() * 3 + si.bead;
            beadSet[cumnn+idx] = bindex;
            krep[cumnn+idx] = be->getRepulsionConst();
            slen[cumnn+idx] = be->getScreeningLength();
            idx++;


//            if (_neighborList->getNeighbors(be)[ni]->isPlusEnd())
//            {bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->getIndex();check=true;}
//            else if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
//            {bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->getIndex();check=true;}
//                if(check){
//                    beadSet[cumnn+idx] = bindex;
//                    krep[cumnn+idx] = be->getRepulsionConst();
//                    slen[cumnn+idx] = be->getScreeningLength();
//                    idx++;
//                }
        }
        nneighbors[i]=idx;
        cumnn+=idx;
        nintvec[i] = cumnn;
#ifdef CUDAACCL
        if(dynamic_cast<PlaneBoundaryElement*>(beList[i])) {
            floatingpoint *x = new floatingpoint[4];
            beList[i]->elementeqn(x);
            beListplane[4 * i] = x[0];
            beListplane[4 * i +1] = x[1];
            beListplane[4 * i +2] = x[2];
            beListplane[4 * i +3] = x[3];
            delete [] x;
        }
        else{
            cout<<"CUDA cannot handle non-plane type boundaries. Exiting..."<<endl;
            exit(EXIT_FAILURE);
        }
#endif
    }
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    _FFType.optimalblocksnthreads(nint, stream);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * nint * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * nint * sizeof(int),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_krep, nint * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_krep, krep, nint * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_slen, nint * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_slen, slen, nint * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nintperbe, nbe * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nintperbe, nintvec, nbe * sizeof(int),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beListplane, 4 * nbe * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beListplane, beListplane, 4 * nbe * sizeof
                                                                                        (floatingpoint), cudaMemcpyHostToDevice));

    //PINNED memory to accelerate data transfer speed
//    CUDAcommon::handleerror(cudaHostAlloc((void**)&U_i, sizeof(floatingpoint), cudaHostAllocDefault), "cudaHOstAlloc",
//                            "BoundaryCylinderRepulsion.cu");
//    CUDAcommon::handleerror(cudaMalloc((void **) &gU, sizeof(floatingpoint)), "cudaMalloc", "BoundaryCylinderRepulsion.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(nint);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += nint;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)),"cuda data"
                                    " transfer",
                            "BoundaryCylinderRepulsion.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "BoundaryCylinderRepulsion.cu");
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"BoundaryCylinderRepulsion.cu",
//                            "vectorizeFF");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif
    delete [] beListplane;
    delete [] nintvec;
}

template<class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] krep;
    delete [] slen;
    delete [] nneighbors;

#ifdef CUDAACCL
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    _FFType.deallocate();
//    if(nint>0) {
//        CUDAcommon::handleerror(cudaStreamDestroy(stream));
        CUDAcommon::handleerror(cudaFree(gpu_beadSet));
        CUDAcommon::handleerror(cudaFree(gpu_krep));
        CUDAcommon::handleerror(cudaFree(gpu_slen));
        CUDAcommon::handleerror(cudaFree(gpu_beListplane));
        CUDAcommon::handleerror(cudaFree(gpu_nintperbe));
        CUDAcommon::handleerror(cudaFree(gpu_params));
//        CUDAcommon::handleerror(cudaFreeHost(U_i));
//        CUDAcommon::handleerror(cudaFree(gU));
//    }
#endif
}

template <class BRepulsionInteractionType>
floatingpoint BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint *coord) {
    floatingpoint U_ii=0.0;

//    std::cout<<"Total boundary nint "<<nint<<endl;
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL
    floatingpoint* gU_i;
//    std::cout<<"Boundary nint "<<nint<<endl;
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_slen, gpu_nintperbe, gpu_beListplane,
//                            gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_slen, gpu_nintperbe, gpu_beListplane,
                            gpu_d, gpu_params);
//    }
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"BoundaryCylinderRepulsion.cu",
//                            "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
#endif
#endif
#ifdef SERIAL
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    U_ii = _FFType.energy(coord, beadSet, krep, slen, nneighbors);

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif
#endif
    return U_ii;
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;

    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_slen, gpu_nintperbe, gpu_beListplane,
                       gpu_params);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_slen, gpu_nintperbe, gpu_beListplane,
                       gpu_params);
    }

    //TODO remove this later need not copy forces back to CPU.
//    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(floatingpoint),
//                                       cudaMemcpyDeviceToHost),"cuda data transfer", "BoundaryCylinderRepulsion.cu");
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef SERIAL
    _FFType.forces(coord, f, beadSet, krep, slen, nneighbors);
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeF += elapsed_runs.count();
#endif
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif
}


namespace {


template< typename InteractionType >
void boundaryCylinderRepulsionLoadForce(
    const InteractionType& interaction, floatingpoint kRep, floatingpoint screenLen,
    const Bead& bo, Bead& bd, typename BoundaryCylinderRepulsion< InteractionType >::LoadForceEnd end,
    BoundaryElement* be
) {
    using LoadForceEnd = typename BoundaryCylinderRepulsion< InteractionType >::LoadForceEnd;

    auto& loadForces = (end == LoadForceEnd::Plus ? bd.loadForcesP : bd.loadForcesM);
    auto& lfi        = (end == LoadForceEnd::Plus ? bd.lfip        : bd.lfim       );

    // Direction of polymerization
    const auto dir = normalizedVector(bd.coordinate() - bo.coordinate());

    // Array of coordinate values to update
    const auto monSize = SysParams::Geometry().monomerSize   [bd.getType()];
    const auto cylSize = SysParams::Geometry().cylinderNumMon[bd.getType()];

    for (int i = 0; i < cylSize; i++) {

        const auto newCoord = bd.coordinate() + (i * monSize) * dir;

        // Projection magnitude ratio on the direction of the cylinder
        // (Effective monomer size) = (monomer size) * proj
        const auto proj = std::max< floatingpoint >(-dot(vector2Vec< 3 >(be->normal(vec2Vector(newCoord))), dir), 0.0);
        const auto loadForce = interaction.loadForces(be->distance(vec2Vector(newCoord)), kRep, screenLen);

        // The load force stored in bead also considers effective monomer size.
        loadForces[i] += proj * loadForce;
    }

    //reset lfi
    lfi = 0;

} // void boundaryRepulsionLoadForce(...)

} // namespace (anonymous)

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeLoadForces() {
//    std::cout<<"BOUNDARY REPULSION LOAD FORCES DOES NOT USE VECTORIZED FORCES/COORDINATES"<<endl;
    for (auto be: BoundaryElement::getBoundaryElements()) {

        for(auto &c : _neighborList->getNeighbors(be)) {


            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();


            //potential acts on second cylinder bead unless this is a minus end
            if(c->isPlusEnd()) {
                boundaryCylinderRepulsionLoadForce(
                    _FFType, kRep, screenLength,
                    *c->getFirstBead(), *c->getSecondBead(), LoadForceEnd::Plus,
                    be
                );

            }

            if(c->isMinusEnd()) {
                boundaryCylinderRepulsionLoadForce(
                    _FFType, kRep, screenLength,
                    *c->getSecondBead(), *c->getFirstBead(), LoadForceEnd::Minus,
                    be
                );

            }

        }

    }
}
template< typename InteractionType >
void BoundaryCylinderRepulsion< InteractionType >::computeLoadForce(Cylinder* c, LoadForceEnd end) const {
    for (auto be : BoundaryElement::getBoundaryElements()) {

        for(auto cyl : _neighborList->getNeighbors(be)) if(c == cyl) {

            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();

            boundaryCylinderRepulsionLoadForce(
                _FFType, kRep, screenLength,
                (end == LoadForceEnd::Plus ? *c->getFirstBead() : *c->getSecondBead()),
                (end == LoadForceEnd::Plus ? *c->getSecondBead() : *c->getFirstBead()),
                end,
                be
            );


            break;

        }

    }
}

// Explicit template instantiations
template class BoundaryCylinderRepulsion< BoundaryCylinderRepulsionExp >;
