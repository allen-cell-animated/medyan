
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

                vector<floatingpoint> coordinates =
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

                            if(iprime<0 || iprime>=int(_grid[0]) || jprime<0 ||
                               jprime>=int(_grid[1]) || kprime<0 ||
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
    //Creates bins based on the largest rMax in the system.
    auto _nDim = SysParams::Geometry().nDim;
    floatingpoint searchdist = 1.125 * (sqrt(_largestrMaxsq));
    std::cout<<"H searchdist "<<searchdist<<" rMax "<<sqrt(_largestrMaxsq)<<endl;
    _binSize = {searchdist, searchdist, searchdist};

    if(_nDim >=1) {
        _size.push_back(int(SysParams::Geometry().NX * SysParams::Geometry()
                .compartmentSizeX));
        if( (_size[0]) % int(_binSize[0]) ==0)
            _grid.push_back(_size[0]/_binSize[0]);
        else
            _grid.push_back(_size[0]/_binSize[0] + 1);
    }

    if (_nDim >= 2) {
        _size.push_back(int(SysParams::Geometry().NY * SysParams::Geometry()
                .compartmentSizeY));
        if( (_size[1]) % int(_binSize[1]) ==0)
            _grid.push_back(_size[1]/_binSize[1]);
        else
            _grid.push_back(_size[1]/_binSize[1] + 1);
    }

    if (_nDim == 3) {
        _size.push_back(int(SysParams::Geometry().NZ * SysParams::Geometry()
                .compartmentSizeZ));
        if( (_size[2]) % int(_binSize[2]) ==0)
            _grid.push_back(_size[2]/_binSize[2]);
        else
            _grid.push_back(_size[2]/_binSize[2] + 1);
    }

    //Check that grid and compartmentSize match nDim
    if((_nDim == 3 && _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 && _binSize[0] != 0 &&
        _binSize[1] != 0 && _binSize[2] != 0)){
        //Do nothing
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
    cout<<"Grid size "<<_grid[0]<<" "<<_grid[1]<<" "<<_grid[2]<<endl;
    cout<<"Total number of bins "<<size<<endl;
    //Set the instance of this grid with given parameters
    _binGrid = new BinGrid(size, _ID, _binSize);
    //Create connections based on dimensionality
    generateConnections();
    cout<<"connections generated"<<endl;
}

//You need a vector of all grids so you can loop through and update respective coordinates.
Bin* HybridCylinderCylinderNL::getBin(const vector<floatingpoint> &coords) {
    //Check if out of bounds
    size_t index = 0;
    size_t i = 0;
    for(auto x: coords)
    {
        //Flatten the coordinates to 1D, get integer index
        if(i == 0) {
            if(x < 0 || x >= (_binSize[0] * _grid[0])) {

                throw OutOfBoundsException();
            }

            index += int(x / _binSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_binSize[1] * _grid[1])) {

                throw OutOfBoundsException();
            }

            index += int(x / _binSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_binSize[2] * _grid[2])) {

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
    auto cindex = currcylinder->getStableIndex();
    const auto& c = Cylinder::getDbData().value[cindex];

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
                    const auto& ncylinder = Cylinder::getDbData().value[ncindex];
                    short ftype2 = ncylinder.type;
//                    //Don't add the same cylinder
//                    if (c.ID == ncylinder.ID) continue;
                    // Testing if a half neighborlist will be stable
                    if(c.id <= ncylinder.id) continue;
                    //Don't add if belonging to same parent
                    if (c.filamentId == ncylinder.filamentId) {
                        auto distsep = fabs(c.positionOnFilament - ncylinder.positionOnFilament);
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
                        floatingpoint dist = distance2(c.coord, ncylinder.coord);
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
                            Cylinder *Ncylinder = Cylinder::getStableElement(ncindex);
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
}

#endif