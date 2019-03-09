
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

#include "Compartment.h"
#include "MathFunctions.h"
#include "Visitor.h"
#include "Parser.h"
//REMOVE LATER
#include "ChemNRMImpl.h"
#include "Filament.h"
#include "Cylinder.h"
#include <stdint.h>
#include "GController.h"

#ifdef SIMDBINDINGSEARCH

void Compartment::SIMDcoordinates(){
    //bscoords = new dist::Coords;
//    std::cout<<"Address of Coords "<<&bscoords<<endl;
    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    int N = _cylinders.size() * SysParams::Chemistry().maxbindingsitespercylinder;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(
                N), bindsitecoordinatesZ(N);
        cindex_bs.resize(N);
        //cID_bs.resize(N);
        Cyldcindexvec.resize(_cylinders.size());
        CylcIDvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        for (auto cyl:_cylinders) {
            int cindex = cyl->_dcIndex;
            auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
            if (checkftype)
                short _filamentType = cylinderstruct.type;
            auto x1 = cyl->getFirstBead()->coordinate;
            auto x2 = cyl->getSecondBead()->coordinate;
            uint16_t shiftedindex = (i << 4);
            //uint32_t shiftedCID = cyl->getID()<<4;
            Cyldcindexvec[i] = cyl->_dcIndex;
            CylcIDvec[i] = cyl->getID();
            i++;
            uint16_t j = 0;
            for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {

                auto mp = (float) *it / SysParams::Geometry().cylinderNumMon[_filamentType];
                auto coord = midPointCoordinate(x1, x2, mp);
                bindsitecoordinatesX[k] = coord[0];
                bindsitecoordinatesY[k] = coord[1];
                bindsitecoordinatesZ[k] = coord[2];
                //last 4 bits are binding site while first 12 bits are cylinder index.
                cindex_bs[k] = shiftedindex | j;
                //cID_bs[k] = shiftedCID |j;
                k++;
                j++;

            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
//        cout<<"#bsites "<<bindsitecoordinatesX.size()<<endl;
        bscoords.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                             bindsitecoordinatesZ, cindex_bs);
    }
}

void Compartment::SIMDcoordinates4linkersearch(bool isvectorizedgather){
    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 1;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    int N = _cylinders.size() * maxnbs;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(N),
                       bindsitecoordinatesZ(N);
        vector<uint16_t> cindex_bs(N);
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        if(isvectorizedgather) {
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (boundstate[bstatepos][maxnbs * cindex + j]) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
        else{
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if(cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN
                            () == 1.0) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
        bscoordslinker.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                             bindsitecoordinatesZ,
                             cindex_bs);
    }
}

void Compartment::SIMDcoordinates4motorsearch(bool isvectorizedgather){
    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 2;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    int N = _cylinders.size() * maxnbs;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(N),
                bindsitecoordinatesZ(N);
        vector<uint16_t> cindex_bs(N);
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        if(isvectorizedgather) {
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (boundstate[bstatepos][maxnbs * cindex + j]) {
                        int A = boundstate[bstatepos][maxnbs * cindex + j] -
                                cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN
                                ();
                        if(abs(A) != 0){
                            cout<<"OOPS Problematic! Occupies sites are passed to "
                                  "SIMD"<<endl;
                            cout<<"Cylinder "<<cyl->getID()<<" "<<*it<<endl;
                        }
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
        else{
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN
                            () == 1.0) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
        bscoordsmotor.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                             bindsitecoordinatesZ,
                             cindex_bs);
    }
}

void Compartment::SIMDcoordinates_section(){

    for(short i =0; i < 27; i++) {
        partitionedcoordx[i].clear();
        partitionedcoordy[i].clear();
        partitionedcoordz[i].clear();
        cindex_bs_section[i].clear();
    }
    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    int N = _cylinders.size() * SysParams::Chemistry().maxbindingsitespercylinder;
    if(N) {
        Cyldcindexvec.resize(_cylinders.size());
        CylcIDvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint32_t k = 0;

        for (auto cyl:_cylinders) {
            int cindex = cyl->_dcIndex;
            auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
            if (checkftype)
                short _filamentType = cylinderstruct.type;
            auto x1 = cyl->getFirstBead()->coordinate;
            auto x2 = cyl->getSecondBead()->coordinate;
//            uint32_t shiftedindex = (i << 4);
            uint32_t shiftedindex = (cyl->_dcIndex << 4);

            Cyldcindexvec[i] = cyl->_dcIndex;
            CylcIDvec[i] = cyl->getID();
            uint32_t j = 0;
            for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {

                auto mp = (float) *it / SysParams::Geometry().cylinderNumMon[_filamentType];
                auto coord = midPointCoordinate(x1, x2, mp);
                //last 4 bits are binding site while first 12 bits are cylinder index.
                uint32_t index = shiftedindex | j;
                j++;
                int pindices[3];
                getpartition3Dindex(pindices, coord);
                addcoordtopartitons(pindices, coord, index);
                j++;
            }
        }

//    assert(k<65536);
        //Create input vector for SIMD calculations
//        cout<<bscoords.size()<<" "<<partitionedcoordx[0].size()<<endl;
        for(short i =0; i < 27; i++) {
//            cout<<partitionedcoordx[i].size()<<" ";
            bscoords_section[i].init_coords(partitionedcoordx[i], partitionedcoordy[i],
                                 partitionedcoordz[i], cindex_bs_section[i]);
        }
//        cout<<endl;
    }
}

void Compartment::SIMDcoordinates4linkersearch_section(bool isvectorizedgather){

    for(short i =0; i < 27; i++) {
        partitionedcoordx[i].clear();
        partitionedcoordy[i].clear();
        partitionedcoordz[i].clear();
        cindex_bs_section[i].clear();
    }

    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 1;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    double _rMax = 40;
    double searchdist = SysParams::Geometry().largestCylinderSize/2 + _rMax;
    bool rMaxvsCmpSize = (SysParams::Geometry().largestCylinderSize/2 + _rMax) <
                            SysParams::Geometry().largestCompartmentSide/2;

    double coord_bounds[6] = {
            _coords[0] - SysParams::Geometry().compartmentSizeX/2 + searchdist,
            _coords[1] - SysParams::Geometry().compartmentSizeY/2 + searchdist,
            _coords[2] - SysParams::Geometry().compartmentSizeZ/2 + searchdist,
            _coords[0] + SysParams::Geometry().compartmentSizeX/2 - searchdist,
            _coords[1] + SysParams::Geometry().compartmentSizeY/2 - searchdist,
            _coords[2] + SysParams::Geometry().compartmentSizeZ/2 - searchdist};

    int N = _cylinders.size() * maxnbs;
    if(N) {
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
//        cout<<"Cmp coord "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;
        for (auto cyl:_cylinders) {
            uint32_t cindex = cyl->_dcIndex;
            auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
            if (checkftype)
                short _filamentType = cylinderstruct.type;
            auto x1 = cyl->getFirstBead()->coordinate;
            auto x2 = cyl->getSecondBead()->coordinate;
//            uint16_t shiftedindex = (i << 4);
            uint32_t shiftedindex = (cindex << 4);
            Cyldcindexvec[i] = cindex;
            i++;
/*			cout<<cyl->getID()<<" ";
	        auto guessCmp = GController::getCompartment(cyl->coordinate);
	        if(guessCmp != this){
	        	cout<<endl;
		        cout<<"OOPS! Cylinder does not belong to this compartment, Linker "
				"slice"<<endl;
		        cout<<"cylinder compartment "<<guessCmp->coordinates()[0]<<" "
		                                                                   ""<<guessCmp->coordinates()[2]<<endl;
		        cout<<"current compartment "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;
		        exit(EXIT_FAILURE);
	        }*/

            uint32_t j = 0;
            for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                bool state = false;
                if(isvectorizedgather)
                    state = checkoccupancy(boundstate, bstatepos, maxnbs * cindex + j);
                else
                    state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                if (state) {
                    auto mp = (float) *it /
                            SysParams::Geometry().cylinderNumMon[_filamentType];
                    auto coord = midPointCoordinate(x1, x2, mp);
                    //last 4 bits are binding site while first 12 bits are cylinder index.
                    uint32_t index = shiftedindex | j;
                    int pindices[3];
                    if(rMaxvsCmpSize){
                        getpartitionindex<true>(pindices, coord, coord_bounds);
                        addcoordtorMaxbasedpartitons<true>(pindices, coord, index);
                    }
                    else{
                        getpartitionindex<false>(pindices, coord, coord_bounds);
                        addcoordtorMaxbasedpartitons<false>(pindices, coord, index);
                    }
/*                    getpartition3Dindex(pindices, coord, coord_bounds);
                    addcoordtopartitons_smallrmax(pindices, coord, index);
                    getpartition3Dindex(pindices, coord);
                    addcoordtopartitons(pindices, coord, index);*/
                }
                j++;
            }
        }
//        cout<<endl;
//    assert(k<65536);
        //Create input vector for SIMD calculations
//        cout<<"Linker coord size ";
        for(short i =0; i < 27; i++) {
//            cout<<partitionedcoordx[i].size()<<" ";
            bscoords_section_linker[i].init_coords(partitionedcoordx[i],
                    partitionedcoordy[i], partitionedcoordz[i], cindex_bs_section[i]);
        }
//        cout<<endl;
/*        cout<<"Cmp Linker coord "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;
        for(short partition = 0; partition < 27; partition ++) {
            cout<<"Partition "<<partition<<endl;
            for (int i = 0; i < partitionedcoordx[partition].size(); i++) {
                cout << cindex_bs_section[partition][i]<<" "
                                                      ""<<partitionedcoordx[partition][i] << " " <<
                partitionedcoordy[partition][i] <<
                " " << partitionedcoordz[partition][i] <<" ";
            }
            cout<<endl;
            cout << "---------------------" << endl;
        }*/
    }
    else{
	    for(short i =0; i < 27; i++)
		    bscoords_section_linker[i].resize(0);
    }

}

void Compartment::SIMDcoordinates4motorsearch_section(bool isvectorizedgather){

    for(short i =0; i < 27; i++) {
        partitionedcoordx[i].clear();
        partitionedcoordy[i].clear();
        partitionedcoordz[i].clear();
        cindex_bs_section[i].clear();
    }

    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 2;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    double _rMax = 225;
    double searchdist = SysParams::Geometry().largestCylinderSize/2 + _rMax;
    bool rMaxvsCmpSize = (searchdist) <
                         SysParams::Geometry().largestCompartmentSide/2;

//    cout<<"rMaxvsCmpSize "<<rMaxvsCmpSize<<endl;
    double coord_bounds[6] = {
            _coords[0] - SysParams::Geometry().compartmentSizeX/2 + searchdist,
            _coords[1] - SysParams::Geometry().compartmentSizeY/2 + searchdist,
            _coords[2] - SysParams::Geometry().compartmentSizeZ/2 + searchdist,
            _coords[0] + SysParams::Geometry().compartmentSizeX/2 - searchdist,
            _coords[1] + SysParams::Geometry().compartmentSizeY/2 - searchdist,
            _coords[2] + SysParams::Geometry().compartmentSizeZ/2 - searchdist};

/*    cout<<"Cmp corner coords ";
    for(uint i = 0; i < 6; i ++)
        cout<<coord_bounds[i]<<" ";
    cout<<endl;*/

    int N = _cylinders.size() * maxnbs;
    if(N) {
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint32_t k = 0;
            for (auto cyl:_cylinders) {
                uint32_t cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                /*if(cylinderstruct.ID != cyl->getID()){
                	cout<<"OOPS! cylinder is struct vec is inaccurate"<<endl;
                	exit(EXIT_FAILURE);
                }
                auto guessCmp = GController::getCompartment(cyl->coordinate);
                if(guessCmp != this){
                	cout<<"OOPS! Cylinder does not belong to this compartment, Motor slice"
					   ""<<endl;
                	cout<<"cylinder compartment "<<guessCmp->coordinates()[0]<<" "
																			   ""<<guessCmp->coordinates()[2]<<endl;
                	cout<<"current compartment "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;
	                exit(EXIT_FAILURE);
                }*/
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
//                uint16_t shiftedindex = (i << 4);
                uint32_t shiftedindex = (cindex << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint32_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    bool state = false;
                    if(isvectorizedgather)
                        state = checkoccupancy(boundstate, bstatepos,maxnbs * cindex + j);
                    else
                        state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                    if (state) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        uint32_t index = shiftedindex | j;
                        //split and crosscheck
//						cout<<index<<" ";
                        int pindices[3];
                        if(rMaxvsCmpSize){
                            getpartitionindex<true>(pindices, coord, coord_bounds);
                            addcoordtorMaxbasedpartitons<true>(pindices, coord, index);
                        }
                        else{
                            getpartitionindex<false>(pindices, coord, coord_bounds);
                            addcoordtorMaxbasedpartitons<false>(pindices, coord, index);
                        }
/*                        getpartition3Dindex(pindices, coord, coord_bounds);
                        addcoordtopartitons_smallrmax(pindices, coord, index);
                        getpartition3Dindex(pindices, coord);
                        addcoordtopartitons(pindices, coord, index);*/
                    }
                    j++;
                }
            }
//            cout<<endl;
//    assert(k<65536);
        //Create input vector for SIMD calculations
//        cout<<"Motor coords size ";
        for(short i =0; i < 27; i++) {
//            cout<<partitionedcoordx[i].size()<<" ";
            bscoords_section_motor[i].init_coords(partitionedcoordx[i], partitionedcoordy[i],
                    partitionedcoordz[i], cindex_bs_section[i]);
        }

/*        cout<<"indices ";
        for(auto x:cindex_bs_section[0])
        	cout<<x<<" ";
        cout<<endl;*/

//        cout<<endl;
        /*cout<<"Cmp Motor coord "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;
        for(short partition = 0; partition < 27; partition ++) {
            cout<<"Partition "<<partition<<endl;
            for (int i = 0; i < partitionedcoordx[partition].size(); i++) {
                cout << cindex_bs_section[partition][i]<<" "
                                                         ""<<partitionedcoordx[partition][i] << " " <<
                     partitionedcoordy[partition][i] <<
                     " " << partitionedcoordz[partition][i] <<" ";
            }
            cout<<endl;
            cout << "---------------------" << endl;
        }*/
    }
    else{
	    for(short i =0; i < 27; i++)
		    bscoords_section_motor[i].resize(0);
    }
}

void Compartment::getpartition3Dindex(int (&indices)[3], vector<double> coord){
    short i = 0;

    for(auto x:coord){
        indices[i] =(x > _coords[i]);
        i++;
    }
/*    cout<<indices[0]<<" "<<indices[1]<<" "<<indices[2]<<" "<<coord[0]<<" "<<coord[1]<<" "
        <<coord[2]<<" "<<_coords[0]<<" "<<_coords[1]<<" "<<_coords[2]<<endl;*/
}

void Compartment::getpartition3Dindex(int (&indices)[3], vector<double> coord,
                                      double (&cmpcornercoords)[6]){
    short i = 0;

    for(auto x:coord){
        if(x <= cmpcornercoords[i])
            indices[i] = 0;
        else if(x>=cmpcornercoords[i+3])
            indices[i] = 2;
        else
            indices[i] = 1;

        i++;
    }
}

//if rMax+Cylsize/2+delta is less than CmpSize/2
template<>
void Compartment::getpartitionindex<true>(int (&indices)[3], vector<double> coord,
                       double (&cmpcornercoords)[6]){
    short i = 0;

    for(auto x:coord){
        if(x <= cmpcornercoords[i])
            indices[i] = 0;
        else if(x>=cmpcornercoords[i+3])
            indices[i] = 2;
        else
            indices[i] = 1;

        i++;
    }
}

//if rMax+Cylsize/2+delta is greater than CmpSize/2
template<>
void Compartment::getpartitionindex<false>(int (&indices)[3], vector<double> coord,
                                          double (&cmpcornercoords)[6]){
    short i = 0;

    for(auto x:coord){
        if(x >= cmpcornercoords[i])
            indices[i] = 2;
        else if(x<=cmpcornercoords[i+3])
            indices[i] = 0;
        else
            indices[i] = 1;

        i++;
    }
}

void Compartment::addcoord(vector<double> coord, uint32_t index, short i){
    partitionedcoordx[i].push_back(coord[0]);
    partitionedcoordy[i].push_back(coord[1]);
    partitionedcoordz[i].push_back(coord[2]);
    cindex_bs_section[i].push_back(index);
}

/*void Compartment::addcoord(vector<double> coord, uint16_t index, short i){
    partitionedcoordx[i].push_back(coord[0]);
    partitionedcoordy[i].push_back(coord[1]);
    partitionedcoordz[i].push_back(coord[2]);
    cindex_bs_section[i].push_back(index);
}*/

bool Compartment::checkoccupancy(vector<vector<bool>>& boundstate, short bstatepos, int pos){
    return boundstate[bstatepos][pos];

}
bool Compartment::checkoccupancy(Cylinder* cyl, short it, short _filamentType,
                                 short bstatepos){
    bool status;
    if(bstatepos == 1)
        status = cyl->getCCylinder()->getCMonomer(it)->speciesBound(
            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN() == 1.0;
    else if(bstatepos == 2)
        status = cyl->getCCylinder()->getCMonomer(it)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN() == 1.0;
    return status;
}

void Compartment::addcoordtopartitons(int (&pindices)[3], vector<double> coord, uint32_t
                                    index){
    addcoord(coord, index, 0);

    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 26);
    }
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 12);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 24);
    }
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 19);
    }
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 13);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 21);
    }
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 8);
        addcoord(coord, index, 14);
        //Vertex
        addcoord(coord, index, 22);
    }
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 15);
        //Vertex
        addcoord(coord, index, 20);

    }
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 11);
        addcoord(coord, index, 8);
        //Vertex
        addcoord(coord, index, 23);
    }
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 15);
        addcoord(coord, index, 13);
        //Vertex
        addcoord(coord, index, 25);
    }
}

void Compartment::addcoordtopartitons_smallrmax(int (&pindices)[3], vector<double> coord,
                                uint16_t index){
    addcoord(coord, index, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 26);
    }
    //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 18);
        //Vertex
    }
    //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 12);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 24);
    }
    //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 16);
        //Vertex
    }
    //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        //Edge
        //Vertex
    }
    //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 7);
        //Vertex
    }
    //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 19);
    }
    //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 9);
        //Vertex
    }
    //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 13);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 21);
    }
    //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 14);
        //Vertex
    }
    //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 4);
        //Edge
        //Vertex
    }
    //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 12);
        //Vertex
    }
    //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        //Edge
        //Vertex
    }
    //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        //Edge
        //Vertex
    }
    //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 11);
        //Vertex
    }
    //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 3);
        //Edge
        //Vertex
    }
    //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 13);
        //Vertex
    }
    //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 8);
        addcoord(coord, index, 14);
        //Vertex
        addcoord(coord, index, 22);
    }
    //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 10);
        //Vertex
    }
    //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 15);
        //Vertex
        addcoord(coord, index, 20);

    }
    //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 8);
        //Vertex

    }
    //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        //Edge
        //Vertex
    }
    //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 15);
        //Vertex
    }
    //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 11);
        addcoord(coord, index, 8);
        //Vertex
        addcoord(coord, index, 23);
    }
    //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 17);
        //Vertex
    }
    //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 15);
        addcoord(coord, index, 13);
        //Vertex
        addcoord(coord, index, 25);
    }
}

//there are 27 partitions made based on rMax. These partitions are composites of
// underlying 27 unique voxels that are non-verlapping.

//if rMax+Cylsize/2+delta is lesser than CmpSize/2
template<>
void Compartment::addcoordtorMaxbasedpartitons<true>(int (&pindices)[3], vector<double>
        coord, uint32_t index){
    addcoord(coord, index, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 26);
    }
        //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 18);
        //Vertex
    }
        //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 12);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 24);
    }
        //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 16);
        //Vertex
    }
        //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 5);
        //Edge
        //Vertex
    }
        //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 7);
        //Vertex
    }
        //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 19);
    }
        //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 9);
        //Vertex
    }
        //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 13);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 21);
    }
        //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 14);
        //Vertex
    }
        //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 4);
        //Edge
        //Vertex
    }
        //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 12);
        //Vertex
    }
        //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        //Edge
        //Vertex
    }
        //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        //Edge
        //Vertex
    }
        //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 11);
        //Vertex
    }
        //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 3);
        //Edge
        //Vertex
    }
        //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 13);
        //Vertex
    }
        //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 8);
        addcoord(coord, index, 14);
        //Vertex
        addcoord(coord, index, 22);
    }
        //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 10);
        //Vertex
    }
        //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 15);
        //Vertex
        addcoord(coord, index, 20);

    }
        //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 8);
        //Vertex

    }
        //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        //Edge
        //Vertex
    }
        //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 15);
        //Vertex
    }
        //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 11);
        addcoord(coord, index, 8);
        //Vertex
        addcoord(coord, index, 23);
    }
        //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 17);
        //Vertex
    }
        //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 15);
        addcoord(coord, index, 13);
        //Vertex
        addcoord(coord, index, 25);
    }
}

//if rMax+Cylsize/2+delta is greater than CmpSize/2
template<>
void Compartment::addcoordtorMaxbasedpartitons<false>(int (&pindices)[3], vector<double>
        coord, uint32_t index){
    addcoord(coord, index, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        for(int part = 1; part < 27; part++){
            addcoord(coord, index, part);
        }
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 26);
    }
        //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 12);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 24);
        addcoord(coord, index, 26);
    }
        //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 18);
        addcoord(coord, index, 12);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 24);
    }
        //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 26);
    }
        //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 12);
        addcoord(coord, index, 13);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 21);
        addcoord(coord, index, 24);
        addcoord(coord, index, 26);
    }
        //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 9);
        addcoord(coord, index, 12);
        addcoord(coord, index, 13);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 24);
        addcoord(coord, index, 21);
    }
        //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 19);
    }
        //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 13);
        addcoord(coord, index, 16);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 21);
    }
        //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 5);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 9);
        addcoord(coord, index, 13);
        addcoord(coord, index, 7);
        //Vertex
        addcoord(coord, index, 21);
    }
        //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 10);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 22);
        addcoord(coord, index, 26);
    }
        //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 8);
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 14);
        addcoord(coord, index, 15);
        addcoord(coord, index, 16);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 22);
        addcoord(coord, index, 24);
        addcoord(coord, index, 26);
    }
        //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 15);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 24);
    }
        //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 9);
        addcoord(coord, index, 10);
        addcoord(coord, index, 11);
        addcoord(coord, index, 14);
        addcoord(coord, index, 16);
        addcoord(coord, index, 17);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 22);
        addcoord(coord, index, 23);
        addcoord(coord, index, 26);
    }
        //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 9);
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 13);
        addcoord(coord, index, 15);
        addcoord(coord, index, 17);
        addcoord(coord, index, 18);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 21);
        addcoord(coord, index, 24);
        addcoord(coord, index, 25);
    }
        //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 16);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 23);
    }
        //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 8);
        addcoord(coord, index, 9);
        addcoord(coord, index, 11);
        addcoord(coord, index, 13);
        addcoord(coord, index, 15);
        addcoord(coord, index, 16);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 19);
        addcoord(coord, index, 21);
        addcoord(coord, index, 23);
        addcoord(coord, index, 25);
    }
        //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        addcoord(coord, index, 5);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 7);
        addcoord(coord, index, 9);
        addcoord(coord, index, 13);
        addcoord(coord, index, 15);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 21);
        addcoord(coord, index, 25);
    }
        //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 8);
        addcoord(coord, index, 14);
        //Vertex
        addcoord(coord, index, 22);
    }
        //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 14);
        addcoord(coord, index, 15);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 22);
    }
        //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 4);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 15);
        //Vertex
        addcoord(coord, index, 20);

    }
        //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 10);
        addcoord(coord, index, 11);
        addcoord(coord, index, 14);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 22);
        addcoord(coord, index, 23);
    }
        //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 10);
        addcoord(coord, index, 11);
        addcoord(coord, index, 12);
        addcoord(coord, index, 13);
        addcoord(coord, index, 14);
        addcoord(coord, index, 15);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 22);
        addcoord(coord, index, 23);
        addcoord(coord, index, 25);
    }
        //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 3);
        addcoord(coord, index, 4);
        addcoord(coord, index, 6);
        //Edge
        addcoord(coord, index, 10);
        addcoord(coord, index, 12);
        addcoord(coord, index, 13);
        addcoord(coord, index, 15);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 20);
        addcoord(coord, index, 25);
    }
        //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 2);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 11);
        addcoord(coord, index, 8);
        //Vertex
        addcoord(coord, index, 23);
    }
        //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, 1);
        addcoord(coord, index, 2);
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        //Edge
        addcoord(coord, index, 8);
        addcoord(coord, index, 11);
        addcoord(coord, index, 13);
        addcoord(coord, index, 15);
        addcoord(coord, index, 17);
        //Vertex
        addcoord(coord, index, 23);
        addcoord(coord, index, 25);
    }
        //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, 6);
        addcoord(coord, index, 3);
        addcoord(coord, index, 1);
        //Edge
        addcoord(coord, index, 17);
        addcoord(coord, index, 15);
        addcoord(coord, index, 13);
        //Vertex
        addcoord(coord, index, 25);
    }
}

#endif

Compartment& Compartment::operator=(const Compartment &other) {
    
    _species.clear();
    _internal_reactions.clear();
    _diffusion_reactions.clear();
    other.cloneSpecies(this);
    other.cloneReactions(this);
    _diffusion_rates = other._diffusion_rates;
    
    return *this;
}
    
bool Compartment::apply_impl(SpeciesVisitor &v) {
    for(auto &s : _species.species()) {
        v.visit(s.get());
    }
    return true;
}

bool Compartment::apply_impl(ReactionVisitor &v) {
    for(auto &r : _internal_reactions.reactions()) {
        v.visit(r.get());
    }
    return true;
}

vector<ReactionBase*> Compartment::generateDiffusionReactions(Compartment* C)
{
    vector<ReactionBase*> rxns;
    
    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;
        
        if(C->isActivated()) {
            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
            //Diffusion reaction from "this" compartment to C.
            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour},diff_rate);
            this->addDiffusionReaction(R);
            rxns.push_back(R);
        }
    }

    
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

//Qin
vector<ReactionBase*> Compartment::generateScaleDiffusionReactions(Compartment* C)
{
    vector<ReactionBase*> rxns;

    cout << "neighbor: x = " << C->_coords[0] << ", y = " << C->_coords[1] <<endl;
    auto factor = generateScaleFactor(C);
    cout << "factor = " << factor << endl;
    
    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;
        
        if(C->isActivated()) {
            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
            
            auto diff_rate_s = diff_rate * factor;
       
            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour},diff_rate_s);
            this->addDiffusionReaction(R);
            rxns.push_back(R);
        }
        

    }

    
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

//Qin, generate a scaling factor for diffusion constant. For cylinder with 1 compartment in Z direction only
float Compartment::generateScaleFactor(Compartment* C)
{
    vector<ReactionBase*> rxns;
    
    auto lx = SysParams::Geometry().compartmentSizeX;
    auto ly = SysParams::Geometry().compartmentSizeY;
    auto r = SysParams::Boundaries().diameter / 2; //radius
    //float c1;
    
    if((_coords[0] - lx/2) < r && (_coords[0] + lx/2) > r) {
        cout << "Diffusion Scaling failed" << endl;
        return 1;
    }
    
    if((_coords[1] - ly/2) < r && (_coords[1] + ly/2) > r) {
        cout << "Diffusion Scaling failed" << endl;
        return 1;
    }
    
    auto x = _coords[0];
    auto y = _coords[1];
    auto nx = C->_coords[0];
    auto ny = C->_coords[1];
    float c1;
    float c2;
    
    //scale diffusion rate based on compartment area
    //1. find the location of the neighbor compartment
    if(ny == y) {
        
        //2. calculate the interection point
        //if at lower part
        if(y < r) {
            //if at left
            if(nx < x) c1 = x - lx/2;
            //if at right
            else c1 = x + lx/2;
  
            c2 = r - sqrt(r * r - (c1 - r) * (c1 - r));
            
            //3. calculate scaling factor
            //check if interaction is within compartment
            if(c2 < (y + ly/2) && c2 > (y - ly/2)) {
                float factor = (y + ly/2 - c2) / ly;
                return factor;
            }
            else return 1;

        }
        //if at upper part
        else {
            //at left
            if(nx < x) c1 = x - lx/2;

            else c1 = x + lx/2; //right
            
            c2 = r + sqrt(r * r - (c1 - r) * (c1 - r));

            
            //3. calculate scaling factor
            if(c2 < (y + ly/2) && c2 > (y - ly/2)) {
                float factor = (c2 - y + ly/2) / ly;
                return factor;
            }
            else return 1;

        }
    }
    
    else if(nx == x){
        //2. calculate the interection point
        //if at left part
        if(x < r) {
            //if at lower
            if(ny < y) c1 = y - ly/2;

            //if at upper
            else c1 = y + ly/2;
            
            c2 = r - sqrt(r * r - (c1 - r) * (c1 - r));
            
            //3. calculate scaling factor
            //check if interaction is within compartment
            if(c2 < (x + lx/2) && c2 > (x - lx/2)) {
                float factor = (_coords[0] + lx/2 - c2) / lx;
                return factor;
            }
            else return 1;
            
        }
        //if at right part
        else {
            //at lower
            if(ny < y) c1 = y - ly/2;

            else c1 = y + ly/2; //right
            
            c2 = r + sqrt(r * r - (c1 - r) * (c1 - r));
            
            //3. calculate scaling factor
            if(c2 < (x + lx/2) && c2 > (x - lx/2)) {
                float factor = (c2 - x + lx/2) / lx;
                return factor;
            }
            else return 1;
            
        }

    }
    
}


vector<ReactionBase*> Compartment::generateAllDiffusionReactions() {
    
    vector<ReactionBase*> rxns;

    if(_activated) {
        for (auto &C: _neighbours) {
            
            auto newRxns = generateDiffusionReactions(C);
            rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

vector<ReactionBase*> Compartment::generateAllpairsDiffusionReactions() {
    
    vector<ReactionBase*> rxns;
    
    if(_activated) {
        for (auto &C: _neighbours) {
            if(C->isActivated()){
                //generates diffusion reactions both from and to the chosen compartment
            auto newRxns = generateDiffusionReactions(C);
                 rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
            newRxns = C->generateDiffusionReactions(this);
                 rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
               
            }
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

void Compartment::removeDiffusionReactions(ChemSim* chem, Compartment* C)
{
    //look for neighbor's diffusion reactions
    vector<ReactionBase*> to_remove;

    for(auto &r : C->_diffusion_reactions.reactions()) {
        
        auto rs = r.get()->rspecies()[1];
        if(rs->getSpecies().getParent() == this) {

            r->passivateReaction();
            
            chem->removeReaction(r.get());
            
            to_remove.push_back(r.get());
        }

    }
    
    //remove them
    for(auto &r : to_remove)
        C->_diffusion_reactions.removeReaction(r);
    
}

void Compartment::removeAllDiffusionReactions(ChemSim* chem) {
    
    //remove all diffusion reactions that this has ownership of
    for(auto &r : _diffusion_reactions.reactions()) {
        r->passivateReaction();
        chem->removeReaction(r.get());
    }
    
    _diffusion_reactions.clear();
    
    //remove neighboring diffusing reactions with this compartment
    for (auto &C: _neighbours)
        removeDiffusionReactions(chem, C);
}


void Compartment::transferSpecies(int i) {
    //i axis
    //0 X
    //1 Y
    //2 Z
    //3 all directions
    //get active neighbors
    vector<Compartment*> activeNeighbors;
    
    for(auto &neighbor : _neighbours){
        auto ncoord=neighbor->coordinates();

        if(neighbor->isActivated()){
            if(i==3)
                activeNeighbors.push_back(neighbor);
            else if(mathfunc::twoPointDistance(ncoord,_coords)==(abs(_coords[i]-ncoord[i])))
                activeNeighbors.push_back(neighbor);
        }}
    
    assert(activeNeighbors.size() != 0
           && "Cannot transfer species to another compartment... no neighbors are active");
    if(i<3 && activeNeighbors.size()>1){
        cout<<"Error transferring species along an axis. More than 1 neighbor. Exiting. "<< endl;
        exit(EXIT_FAILURE);
    }
    
    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;
    
    for(auto &sp : _species.species()) {
        
        int copyNumber = sp->getN();
        auto nit = activeNeighbors.begin();
        
        if(sp->getFullName().find("Bound") == string::npos){
            while(copyNumber > 0) {
                sp->down();
                
                //choose a random active neighbor
                auto neighbor = *nit;
                
                sp_neighbor = neighbor->findSpeciesByName(sp->getName());
                
                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);
                
                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);
                
                //increase copy number
                
                sp_neighbor->up();
                
                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                copyNumber--;
                
            }
        }
        
        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : _species.species())
            sp->updateReactantPropensities();
    }
}

void Compartment::shareSpecies(int i) {
    //i axis
    //0 X
    //1 Y
    //2 Z
    //3 all directions
    //get active neighbors
    vector<Compartment*> activeNeighbors;
    
    for(auto &neighbor : _neighbours){
        auto ncoord=neighbor->coordinates();
    if(neighbor->isActivated()){
        if(i==3)
            activeNeighbors.push_back(neighbor);
        else if(mathfunc::twoPointDistance(ncoord,_coords)==(abs(_coords[i]-ncoord[i])))
        activeNeighbors.push_back(neighbor);
    }}
    
    assert(activeNeighbors.size() != 0
           && "Cannot share species to another compartment... no neighbors are active");
    if(i<3 && activeNeighbors.size()>1){
        cout<<"Error sharing species along an axis. More than 1 neighbor. Exiting."<< endl;
        exit(EXIT_FAILURE);
    }
    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;
    
    for(auto &sp : _species.species()) {
        auto nit = activeNeighbors.begin();
        auto neighbor = *nit;
        sp_neighbor = neighbor->findSpeciesByName(sp->getName());
        int copyNumber = sp_neighbor->getN();
        int lowerlimit = (int) sp_neighbor->getN()/2;
        if(sp->getFullName().find("Bound") == string::npos){
            while(copyNumber > lowerlimit) {
                sp_neighbor->down();
                
                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);
                
                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);
                
                //increase copy number
                sp->up();
                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                neighbor = *nit;
                sp_neighbor = neighbor->findSpeciesByName(sp->getName());
                copyNumber--;
                
            }
        }

        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : _species.species())
            sp->updateReactantPropensities();
        
    }
}

void Compartment::activate(ChemSim* chem) {
    
    assert(!_activated && "Compartment is already activated.");
    
    //set marker
    _activated = true;
    //add all diffusion reactions
    auto rxns = generateAllpairsDiffusionReactions();
    for(auto &r : rxns) chem->addReaction(r);
    shareSpecies(SysParams::Boundaries().transfershareaxis);
    
    for (auto &C: _neighbours){
        if(C->isActivated()){
            for(auto &r : C->_diffusion_reactions.reactions()) {
                auto rs = r.get()->rspecies()[1];//product
                if(rs->getSpecies().getParent() == this) {
                    auto rs1 = r.get()->rspecies()[0];
                    if(rs1->getN()>0 && r->isPassivated()){
                        r->activateReaction();
                    }
        }
    }
        }
    }
    
}

void Compartment::deactivate(ChemSim* chem) {
    
    //assert no cylinders in this compartment
    assert((_cylinders.size() == 0)
           && "Compartment cannot be deactivated when containing active cylinders.");
    
    assert(_activated && "Compartment is already deactivated.");
    
    //set marker
    _activated = false;
    
    transferSpecies(SysParams::Boundaries().transfershareaxis);
    removeAllDiffusionReactions(chem);
}

bool operator==(const Compartment& a, const Compartment& b) {
    if(a.numberOfSpecies()!=b.numberOfSpecies() or
       a.numberOfInternalReactions()!=b.numberOfInternalReactions())
        return false;
    
    if(typeid(a)!=typeid(b))
        return false;
    
    bool spec_bool = false;
    auto sit_pair = mismatch(a._species.species().begin(),
                             a._species.species().end(),
                             b._species.species().begin(),
            [](const unique_ptr<Species> &A, const unique_ptr<Species> &B)
            {return (*A)==(*B); });
    if(sit_pair.first==a._species.species().end())
        spec_bool=true;
    
    
    bool reac_bool = false;
    auto rit_pair = mismatch(a._internal_reactions.reactions().begin(),
                             a._internal_reactions.reactions().end(),
                             b._internal_reactions.reactions().begin(),
            [](const unique_ptr<ReactionBase> &A, const unique_ptr<ReactionBase> &B)
            {return (*A)==(*B);});
    if(rit_pair.first==a._internal_reactions.reactions().end())
        reac_bool=true;
    
    return spec_bool && reac_bool;
}