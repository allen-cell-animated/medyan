
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
#ifdef HYBRID_NLSTENCILLIST
#include "HybridNeighborListImpl.h"

#include "Bead.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bubble.h"
#include "BoundaryElement.h"

#include "GController.h"
#include "MathFunctions.h"
#include "CUDAcommon.h"

using namespace mathfunc;

short HybridCylinderCylinderNL::totalhybridNL = 0;

void HybridCylinderCylinderNL::updateallcylinderstobin() {
    for(auto cyl:Cylinder::getCylinders())
        updatebin(cyl);
}

void HybridCylinderCylinderNL::assignallcylinderstobin() {
    for(auto cyl:Cylinder::getCylinders())
        assignbin(cyl);
/*    std::cout<<"H Total number of bins "<< _binGrid->getBins().size()<<endl;
    for(auto bin:_binGrid->getBins()){
        std::cout<<bin->getCylinders().size()<<" ";
    }
    std::cout<<endl;*/
}

void HybridCylinderCylinderNL::assignbin(Cylinder* cyl){
    Bin* _bin;
    try {_bin = getBin(cyl->coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    _bin->addCylinder(cyl);
    cyl->_hbinvec.push_back(_bin);
}

void HybridCylinderCylinderNL::unassignbin(Cylinder* cyl, Bin* bin){
    bin->removeCylinder(cyl);
}

void HybridCylinderCylinderNL::updatebin(Cylinder *cyl){
    Bin* _bin;
//    std::cout<<coordinate[0]<<" "<<coordinate[1]<<" "<<coordinate[2]<<endl;
    try {_bin = getBin(cyl->coordinate);}
    catch (exception& e) {
        cout << e.what();
        cyl->printSelf();
        exit(EXIT_FAILURE);
    }
    if(_bin != cyl->_hbinvec.at(_ID)) {
#ifdef CHEMISTRY
        auto oldBin = cyl->_hbinvec.at(_ID);
        auto newBin = _bin;
#endif
        //remove from old compartment, add to new
        oldBin->removeCylinder(cyl);
        cyl->_hbinvec.at(_ID) = newBin;
        _bin->addCylinder(cyl);
    }
}

void HybridCylinderCylinderNL::generateConnections() {
    for(size_t i=0U; i<_grid[0]; ++i) {

        for(size_t j=0U; j<_grid[1]; ++j) {

            for(size_t k=0U; k<_grid[2]; ++k) {
                vector<size_t> indices{i,j,k};
                Bin *target = getBin(indices);//defined in this file.

                vector<double> coordinates =
                        {indices[0] * _binSize[0] + _binSize[0] / 2,
                         indices[1] * _binSize[1] + _binSize[1] / 2,
                         indices[2] * _binSize[2] + _binSize[2] / 2};
                target->setCoordinates(coordinates);
                int stencilcount = 0;

                //Go through all neighbors to get the neighbors list
                short nneighbors = 1;
                for(int ii = -nneighbors; ii <= nneighbors; ii++){
                    for(int jj = -nneighbors; jj <= nneighbors; jj++){
                        for(int kk = -nneighbors; kk <= nneighbors; kk++){
                            //Consider the target bin itself as a neighbor.

                            stencilcount++;
                            int iprime = i+ii;
                            int jprime = j+jj;
                            int kprime = k+kk;

/*                            cout<<ii<<" "<<jj<<" "<<kk<<" "<<iprime<<" "<<jprime<<" "
                                <<kprime<<" "<<_grid[0]<<" "<<_grid[1]<<" "<<_grid[2]<<endl;*/

                            if(iprime<0 or iprime>=int(_grid[0]) or jprime<0 or
                               jprime>=int(_grid[1]) or kprime<0 or
                               kprime>=int(_grid[2]))
                                continue;
                            vector<size_t> currentIndices{size_t(iprime), size_t
                                    (jprime), size_t(kprime)};
                            Bin *neighbor = getBin(currentIndices);
                            target->addNeighbour(neighbor);
                            target->stencilID.push_back(stencilcount-1);//All 125
                            //Only 63 neighbors will be added
                            if(j>=1 || (j==0 && k>0) || (j==0 && k == 0 && i <=0)){
                                target->adduniquepermutationNeighbour(neighbor);
                            }
                        }
                    }
                }
            }
        }
    }


    /*for(size_t i=0U; i<_grid[0]; ++i) {

        for (size_t j = 0U; j < _grid[1]; ++j) {

            for (size_t k = 0U; k < _grid[2]; ++k) {
                vector<size_t> indices{i, j, k};
                Bin *target = getBin(indices);
                std::cout << "Target " << target->coordinates()[0] << " " <<
                          target->coordinates()[1] << " " <<
                          target->coordinates()[2] << " " << endl;
                std::cout<<"Bin size "<<_binSize[0]<<endl;
                for (int ii: {-1, 0, 1}) {
                    for (int jj: {-1, 0, 1}) {
                        for (int kk: {-1, 0, 1}) {
                            int iprime = i + ii;
                            int jprime = j + jj;
                            int kprime = k + kk;
                            if (iprime < 0 or iprime == int(_grid[0]) or jprime < 0 or
                                jprime == int(_grid[1]) or kprime < 0 or
                                kprime == int(_grid[2]))
                                continue;
                            vector<size_t> currentIndices{size_t(iprime), size_t
                                    (jprime), size_t(kprime)};
                            Bin *neighbor = getBin(currentIndices);
                            std::cout << "Neighbor " << neighbor->coordinates()[0]
                                      << " " <<
                                      neighbor->coordinates()[1] << " " <<
                                      neighbor->coordinates()[2] << " " << endl;
                        }
                    }
                }
            }
        }
    }*/
}

void HybridCylinderCylinderNL::initializeBinGrid() {

//    //Initial parameters of system
    auto _nDim = SysParams::Geometry().nDim;
    double searchdist = 1.125 * (sqrt(_largestrMaxsq));
    std::cout<<"H searchdist "<<searchdist<<" rMax "<<sqrt(_largestrMaxsq)<<endl;
    _binSize = {searchdist, searchdist, searchdist};
    if(_nDim >=1) {
        _size.push_back(int(SysParams::Geometry().NX * SysParams::Geometry()
                .compartmentSizeX));
        if( (_size[0]) % int(_binSize[0]) ==0)
            _grid.push_back(_size[0]/_binSize[0]);
        else
            _grid.push_back(_size[0]/_binSize[0] + 1);
        cout<<_grid[0]<<" "<<_size[0]<<" "<<_binSize[0]<<endl;
    }
    if (_nDim >= 2) {
        _size.push_back(int(SysParams::Geometry().NY * SysParams::Geometry()
                .compartmentSizeY));
        if( (_size[1]) % int(_binSize[1]) ==0)
            _grid.push_back(_size[1]/_binSize[1]);
        else
            _grid.push_back(_size[1]/_binSize[1] + 1);
        cout<<_grid[1]<<" "<<_size[1]<<" "<<_binSize[1]<<endl;
    }
    if (_nDim == 3) {
        _size.push_back(int(SysParams::Geometry().NZ * SysParams::Geometry()
                .compartmentSizeZ));
        if( (_size[2]) % int(_binSize[2]) ==0)
            _grid.push_back(_size[2]/_binSize[2]);
        else
            _grid.push_back(_size[2]/_binSize[2] + 1);
        cout<<_grid[2]<<" "<<_size[2]<<" "<<_binSize[2]<<endl;
    }

    //Check that grid and compartmentSize match nDim
    if((_nDim == 3 &&
        _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 &&
        _binSize[0] != 0 &&
        _binSize[1] != 0 &&
        _binSize[2] != 0)){
    }
    else {
        cout << "Bin parameters for CylinderCylinderNeighborLists are invalid. Exiting." <<
             endl;
        exit(EXIT_FAILURE);
    }
    int size = 1;
    for(auto x: _grid) {
        if(x != 0) size*=x;
    }
    cout<<_grid[0]<<" "<<_grid[1]<<" "<<_grid[2]<<endl;
    cout<<size<<endl;
    //Set the instance of this grid with given parameters
    _binGrid = new BinGrid(size, _ID, _binSize);
    //Create connections based on dimensionality
    generateConnections();
}

//You need a vector of all grids so you can loop through and update respective coordinates.
Bin* HybridCylinderCylinderNL::getBin(const vector<double> &coords) {
    //Check if out of bounds
    size_t index = 0;
    size_t i = 0;
    for(auto x: coords)
    {
        //Flatten the coordinates to 1D, get integer index
        if(i == 0) {
            if(x < 0 || x >= (_binSize[0] * _grid[0])) {
                cout<<"get Bin coords x"<<endl;
                throw OutOfBoundsException();
            }

            index += int(x / _binSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_binSize[1] * _grid[1])) {
                cout<<"get Bin coords y"<<endl;
                throw OutOfBoundsException();
            }

            index += int(x / _binSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_binSize[2] * _grid[2])) {
                cout<<"get Bin coords z"<<endl;
                throw OutOfBoundsException();
            }

            index += int(x / _binSize[2]) * _grid[0] * _grid[1];
        }
        i++;
    }

    try {
        return _binGrid->getBin(index);
    }
    catch (exception& e){
        cout << "Bad bin access at..." << endl;
        cout << "Bin index = " << index << endl;
        cout << "Coords = " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
        throw NaNCoordinateException();
    }
}

Bin* HybridCylinderCylinderNL::getBin(const vector<size_t> &indices) {
    size_t index = 0;
    size_t i = 0;
    for(auto x: indices)
    {
        //Flatten the indices to 1D
        if(i == 0) {
            if(x >= _grid[0]) {
                cout<<"get Bin x"<<endl;
                throw OutOfBoundsException();
            }

            index += x;
        }
        else if(i == 1) {
            if(x >= _grid[1]) {
                cout<<"get Bin y"<<endl;
                throw OutOfBoundsException();
            }

            index += x * _grid[0];
        }
        else {
            if(x >= _grid[2]) {
                cout << "get Bin z" << endl;
                throw OutOfBoundsException();
            }

            index += x * _grid[0] * _grid[1];
        }

        i++;
    }
    try {
        return _binGrid->getBin(index);
    }
    catch (exception& e){
        cout << "Bad Bin access at..." << endl;
        cout << "Bin index = " << index << endl;
        cout << "Indices = " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        throw NaNCoordinateException();
    }
}

void HybridCylinderCylinderNL::updateNeighborsbin(Cylinder* currcylinder, bool runtime){
    //clear existing neighbors of currcylinder from all neighborlists
    for(int idx = 0; idx < totaluniquefIDpairs; idx++) {
        int countbounds = _rMaxsqvec[idx].size();
        for (int idx2 = 0; idx2 < countbounds; idx2++) {
            auto HNLID = HNLIDvec[idx][idx2];
            _list4mbinvec[HNLID][currcylinder].clear();
        }
    }
    //get necessary variables
    auto binvec = currcylinder->_hbinvec;//The different hybrid bins that this cylinder
    // belongs to.
    //Check if the cylinder has been assigned a bin. If not, assign.
    if(binvec.size()<=_ID)
        assignbin(currcylinder);
    binvec = currcylinder->_hbinvec;
    //get parent bin corresponding to this hybrid neighbor list.
    auto parentbin =  binvec.at(_ID);
    //get neighboring bins
    vector<Bin*> _neighboringBins = binvec.at(_ID)//Get the bin that belongs to the
                    // current binGrid of interest for this NL.
                                                    ->getNeighbours();
    double *coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    auto cylinderpointervec = CUDAcommon::getSERLvars().cylinderpointervec;
    int cindex = currcylinder->_dcIndex;
    cylinder c = cylindervec[cindex];

    //
    int ncyls2 = 0;
    int tcyl2 = 0;
    int nbincount = 0;
    auto nbinstencil = parentbin->stencilID;// A standard templated numbering of
    // neighboring bins is implemented i.e. based on position w.r.t. bin of interest,
    // neighboring bins are given a particular ID.nbinstencil stores the set of such
    // neighbors that is close to bin of interest. Bins close to the boundary will have
    // < 27 elements in the stencilID vector.
    short ftype1 = c.type; //cylinder type and filament type is one and the
    // same.
    float _largestrMax = sqrt(_largestrMaxsq);
    for (auto &bin : _neighboringBins) {
            bool isbinneeded = _binGrid->iswithincutoff(c.coord,
                                                        parentbin->coordinates(),
                                                        nbinstencil.at(nbincount),
                                                        _largestrMax);
            nbincount++;
            if (isbinneeded) {
                auto cindicesvec = bin->getcindices();
                int numneighbors = cindicesvec.size();
                for (int iter = 0; iter < numneighbors; iter++) {
                    int ncindex = cindicesvec[iter];
                    cylinder ncylinder = cylindervec[ncindex];
                    short ftype2 = ncylinder.type;
//                    //Don't add the same cylinder
//                    if (c.ID == ncylinder.ID) continue;
                    // Testing if a half neighborlist will be stable
                    if(c.ID <= ncylinder.ID) continue;
                    //Don't add if belonging to same parent
                    if (c.filamentID == ncylinder.filamentID) {
                        auto distsep = fabs(c.filamentposition - ncylinder.filamentposition);
                        //if not cross filament, check if not neighboring
                        if (distsep <= 2) continue;
                    }

                    //Loop through all the distance bounds and add to neighborlist
                    for (int idx = 0; idx < totaluniquefIDpairs; idx++) {
                        int countbounds = _rMaxsqvec[idx].size();
                        auto fpairs = _filamentIDvec[idx].data();
                        //Check for cylinder filament types
                        if (ftype1 < ftype2) {
                            if (ftype1 != fpairs[0] || ftype2 != fpairs[1])continue;
                        }
                        else if (ftype1 != fpairs[1] || ftype2 != fpairs[0]) continue;
                        double dist = twoPointDistancesquared(c.coord, ncylinder.coord);
                        if (dist < _smallestrMinsq || dist > _largestrMaxsq) continue;
                        for (int idx2 = 0; idx2 < countbounds; idx2++) {
                            //Dont add if ID is more than cylinder for half-list
                            //if (!_fullstatusvec[idx][idx2] && c.ID <= ncylinder.ID) continue;
                            //Dont add if not within range
                            if (dist > _rMaxsqvec[idx][idx2] ||
                                dist < _rMinsqvec[idx][idx2])
                                continue;
                            short HNLID = HNLIDvec[idx][idx2];
                            //If we got through all of this, add it!
                            Cylinder *Ncylinder = cylinderpointervec[ncindex];
                            _list4mbinvec[HNLID][currcylinder].push_back(Ncylinder);

                            //if runtime, add to other list as well if full
                            /* if (runtime && _fullstatusvec[idx][idx2]) {
                                _list4mbinvec[HNLID][Ncylinder].push_back(currcylinder);
                            }*/
                        }
                    }
                }
            }
    }
}

vector<Cylinder*> HybridCylinderCylinderNL::getNeighborsstencil(short HNLID, Cylinder*
                                                                cylinder) {

    return _list4mbinvec[HNLID][cylinder];
}

void HybridCylinderCylinderNL::addNeighbor(Neighbor* n) {

    //return if not a cylinder!
    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;

    //update neighbors
    updateNeighborsbin(cylinder, true);
}

void HybridCylinderCylinderNL::removeNeighbor(Neighbor* n) {

    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;
    for(int idx = 0; idx < totaluniquefIDpairs; idx++) {
        int countbounds = _rMaxsqvec[idx].size();
        for (int idx2 = 0; idx2 < countbounds; idx2++) {
            auto HNLID = HNLIDvec[idx][idx2];
/*            std::cout << "Removing neighbors of cylinder with cindex " <<
                          cylinder->_dcIndex<<" and ID "<<cylinder->getID() << " from NL " << HNLID << endl;*/
            //Remove from NeighborList
            _list4mbinvec[HNLID].erase(cylinder);
            //Remove from bin
            Bin *bin = cylinder->_hbinvec.at(_ID);
            unassignbin(cylinder, bin);
            //remove from other lists
//            std::cout << "Removed from cylinders ";
            for (auto it = _list4mbinvec[HNLID].begin();
                 it != _list4mbinvec[HNLID].end(); it++) {
                auto cit = find(it->second.begin(), it->second.end(), cylinder);
                {
                    if (cit != it->second.end()) {
                        it->second.erase(cit);
//                        std::cout << it->first->getID() << " ";
                    }
                }
            }
//            std::cout<<endl;
        }
    }

}

void HybridCylinderCylinderNL::reset() {

    //loop through all neighbor keys
    for(int idx = 0; idx < totalhybridNL; idx++) {
        _list4mbinvec[idx].clear();
//        std::cout<<"Hybrid rmin rmax "<<_rMinsqvec[idx]<<" "<<_rMaxsqvec[idx]<<endl;
    }

    /*chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();*/
    //check and reassign cylinders to different bins if needed.
    updateallcylinderstobin();
    _binGrid->updatecindices();
    for(auto cylinder: Cylinder::getCylinders()) {
        updateNeighborsbin(cylinder);
//        for (int idx = 0; idx < totalhybridNL; idx++) {
//            tot[idx] += _list4mbinvec[idx][cylinder].size();
//        }
    }
//    std::cout<<endl;
//    for(int idx = 0; idx < totalhybridNL; idx++)
//        std::cout<<"reset HybridNLSTENCILLIST size "<<" "<<tot[idx]<<endl;
/*    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_sten(mine - mins);
    std::cout<<"Hybrid NLSTEN reset time "<<elapsed_sten.count()<<endl;*/

    //Check if HNLID = 1 is symmetric
/*    short HNLID = 1;
    int idx = 0; int idx2 = _rMaxsqvec[idx].size()-1;
    auto _list4mbin = _list4mbinvec[HNLID];
    std::cout<<"map size = " << _list4mbin.size()<<endl;
    std::cout << "max_size = " << _list4mbin.max_size() <<endl;
    for(auto cylinder: Cylinder::getCylinders()) {
        auto neighbors = _list4mbin[cylinder];
        auto cylinderbin = cylinder->_hbinvec[0];
        auto cylbincoord = cylinderbin->coordinates();
        for(auto ncylinder:neighbors){
            auto ncylinderbin = ncylinder->_hbinvec[0];
            auto ncylbincoord = ncylinderbin->coordinates();
            //look for cylinder in the neighbor list of ncylinder
            auto ncylinderneighbors = _list4mbin[ncylinder];
            std::cout<<"neighborvec size "<<ncylinderneighbors.size()<<" capacity "
                     <<ncylinderneighbors.capacity()<<" max_size "<<ncylinderneighbors
                    .max_size()<<endl;
            if(find(ncylinderneighbors.begin(),ncylinderneighbors.end(),cylinder) ==
                    ncylinderneighbors.end()){
                std::cout<<" cylinder "<<cylinder->getID()<<" from bin "
                        ""<<cylinderbin<<" "
                        "coordinates "<<cylbincoord[0]<<" "<<cylbincoord[1]<<" "
                        ""<<cylbincoord[2]<<" has neighbor cylinder "<<ncylinder->getID()
                         <<" from bin "<<ncylinderbin<<" coordinates "
                        ""<<ncylbincoord[0]<<" "<<ncylbincoord[1]<<" "
                        ""<<ncylbincoord[2]<<endl;
                std::cout<<"But neighbor cylinder does not have cylinder in it's "
                        "neighbors list. Check algorithm."<<endl;
            }
        }
    }*/
}

void HybridCylinderCylinderNL::updateSIMDbindingsites(){

    //check and reassign cylinders to different bins if needed.
    chrono::high_resolution_clock::time_point minscreate, minecreate;
    minscreate = chrono::high_resolution_clock::now();
    updateallcylinderstobin();
    _binGrid->updatecindices();
    _binGrid->createSIMDcoordinates();
    minecreate = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_create(minecreate - minscreate);
    cout<<"SIMD NL create time "<<elapsed_create.count()<<endl;
    short idvecL[2] = {0,0};
    short idvecM[2] = {0,1};
    minscreate = chrono::high_resolution_clock::now();
    for(auto bin:_binGrid->getBins()){
    calculatebspairsLMself<1,true, true>(bin, bspairslinkerself, idvecL);
    calculatebspairsLMenclosed<1,false, true>(bin, bspairslinker,idvecL);
    }
    minecreate = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_calculate(minecreate - minscreate);
    cout<<"SIMD NL calculate time L "<<elapsed_calculate.count()<<endl;
    minscreate = chrono::high_resolution_clock::now();
    for(auto bin:_binGrid->getBins()){
            calculatebspairsLMself<1,true, false>(bin, bspairsmotorself, idvecM);
            calculatebspairsLMenclosed<1,false, false>(bin, bspairsmotor,idvecM);
    }
    minecreate = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_calculateM(minecreate - minscreate);
    cout<<"SIMD NL calculate time M "<<elapsed_calculateM.count()<<endl;
}

template <uint D, bool SELF, bool LinkerorMotor>
void HybridCylinderCylinderNL::calculatebspairsLMself(Bin* bin, dist::dOut<D, SELF>&
bspairsoutSself, short idvec[2]){

    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    auto cylcmp1 = bin->Cyldcindexvec;

    minsfind = chrono::high_resolution_clock::now();
    if(bin->getSIMDcoords<LinkerorMotor>().size()) {
        bspairsoutSself.reset_counters();
        dist::find_distances(bspairsoutSself, bin->getSIMDcoords<LinkerorMotor>(),
                             t_avx_par);
    }
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runfind(minefind - minsfind);
    findtimeV2 += elapsed_runfind.count();
/*    cout<<bin->getSIMDcoords<LinkerorMotor>().size()<<" "<<bspairsoutSself
    .counter[D-1]<<endl;*/
    //@{
    /*if(false) {

        uint N = bspairsoutSself.counter[D-1];
        uint prev_size = getfilID_fpospairs<LinkerorMotor>().size();
        getpairsLinkerorMotor<LinkerorMotor>().resize
                (getpairsLinkerorMotor<LinkerorMotor>().size() + 2*N);
        getfilID_fpospairs<LinkerorMotor>().resize(getfilID_fpospairs<LinkerorMotor>()
                                                           .size() + N);

        minsfind = chrono::high_resolution_clock::now();
        std::vector<std::thread> threads_avx;
        uint nt = nthreads;
        threads_avx.reserve(nt);
        uint prev = 0;
        uint frac = N / nt;
        uint next = frac + N % nt;
        for (uint i = 0; i < nt; ++i) {
            threads_avx.push_back(std::thread(
                    &HybridBindingSearchManager::gatherCylinderIDfromcIndex<D, SELF,
                            LinkerorMotor>, this, std::ref(bspairsoutSself), prev,
                    next, prev_size, _compartment));
            prev = next;
            next = min(N, prev + frac);
        }

        //Join
        for (auto &t : threads_avx)
            t.join();
        threads_avx.clear();

        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_append(minefind - minsfind);
        appendtime += elapsed_append.count();
    }*/
    //@}
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridCylinderCylinderNL::calculatebspairsLMenclosed (Bin* bin, dist::dOut<D,SELF>&
bspairsoutS, short idvec[2]){
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    auto cylcmp1 = bin->Cyldcindexvec;

    for(auto nbin: bin->getuniquepermutationNeighbours()){

        minsfind = chrono::high_resolution_clock::now();

        if(bin->getSIMDcoords<LinkerorMotor>().size() > 0 &&
        nbin->getSIMDcoords<LinkerorMotor>().size() > 0) {
            bspairsoutS.reset_counters();
            dist::find_distances(bspairsoutS, bin->getSIMDcoords<LinkerorMotor>(),
                                 nbin->getSIMDcoords<LinkerorMotor>(), t_avx_par);
        }
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runfind(minefind - minsfind);
        findtimeV2 += elapsed_runfind.count();
        /*cout<<bin->getSIMDcoords<LinkerorMotor>().size()<<" "
            <<nbin->getSIMDcoords<LinkerorMotor>().size()<<" "
            <<bspairsoutS.counter[D-1]<<endl;*/
        //MERGE INTO single vector
        //@{
        /*if(false) {
            minsfind = chrono::high_resolution_clock::now();
            short dim = 0;
            uint N = bspairsoutS.counter[dim];
            uint prev_size = getfilID_fpospairs<LinkerorMotor>().size();
            getpairsLinkerorMotor<LinkerorMotor>().resize
                    (getpairsLinkerorMotor<LinkerorMotor>().size() + 2*N);
            getfilID_fpospairs<LinkerorMotor>().resize(getfilID_fpospairs<LinkerorMotor>()
                                                               .size() + N);

            std::vector<std::thread> threads_avx;
            uint nt = nthreads;
            threads_avx.reserve(nt);
            uint prev = 0;
            uint frac = N / nt;
            uint next = frac + N % nt;
            for (uint i = 0; i < nt; ++i) {
                threads_avx.push_back(std::thread
                                              (&HybridBindingSearchManager::gatherCylinderIDfromcIndex<D, SELF,
                                                       LinkerorMotor>, this,
                                               std::ref(bspairsoutS), prev,
                                               next, prev_size, ncmp));
                prev = next;
                next = min(N, prev + frac);
            }

            //Join
            for (auto &t : threads_avx)
                t.join();
            threads_avx.clear();

            minefind = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_append(minefind - minsfind);
            appendtime += elapsed_append.count();
        }*/
        //@}
    }
}
double HybridCylinderCylinderNL::SIMDtime = 0.0;
double HybridCylinderCylinderNL::HYBDtime = 0.0;
double HybridCylinderCylinderNL::findtime = 0.0;
double HybridCylinderCylinderNL::appendtime = 0.0;
double HybridCylinderCylinderNL::findtimeV2 = 0.0;
double HybridCylinderCylinderNL::SIMDparse1 = 0.0;
double HybridCylinderCylinderNL::SIMDparse2 = 0.0;
double HybridCylinderCylinderNL::SIMDparse3 = 0.0;
double HybridCylinderCylinderNL::SIMDcountbs = 0.0;

dist::dOut<1U,false> HybridCylinderCylinderNL::bspairslinker;
dist::dOut<1U,true> HybridCylinderCylinderNL::bspairslinkerself;
dist::dOut<1U,false> HybridCylinderCylinderNL::bspairsmotor;
dist::dOut<1U,true> HybridCylinderCylinderNL::bspairsmotorself;

#endif