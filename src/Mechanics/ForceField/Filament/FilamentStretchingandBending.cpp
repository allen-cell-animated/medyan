
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
#include "FilamentStretchingandBending.h"

#include "FilamentStretchingHarmonicandBendingHarmonic.h"
#include "FilamentStretchingHarmonicandBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "CGMethod.h"
#include "SysParams.h"
#include "SubSystem.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "cross_check.h"
using namespace mathfunc;

template <class FBendingInteractionType>
void FilamentStretchingandBending<FBendingInteractionType>::precomputevars(
		floatingpoint *coord, floatingpoint *cyllengthset,
		floatingpoint *cylbenddotproduct){

    int cylcount = 0;
    int hingecount = 0;
    floatingpoint *coord1, *coord2, *coord3;
    bool fetchcoordsstatus = true;
    for(auto b = 0; b < Bead::getBeads().size()-1;  b++){
        if(fetchcoordsstatus) {
            coord1 = &coord[3 * beadSetall[b]];
            coord2 = &coord[3 * beadSetall[b + 1]];
        }
        else{
            coord1 = coord2;
            coord2 = coord3;
        }
        if (beadpaircyllengthstatus[b]) {
            cyllengthset[cylcount] = twoPointDistance(coord1, coord2);
            cylcount++;
        }
        fetchcoordsstatus = true;
        if(beadtriplethingestatus[b]) {
            coord3 = &coord[3 * beadSetall[b + 2]];
            cylbenddotproduct[hingecount] = scalarProduct(coord1, coord2, coord2, coord3);
            hingecount++;
            //If there are connections between b+1, b+2, and b+3.
            if (beadtriplethingestatus[b + 1]) {
                fetchcoordsstatus = false;
            }
        }
    }

    /*floatingpoint *coord1, *coord2, *coord3;
    int cylcount = 0;
    int dotproductcount = 0;
    for (auto f: Filament::getFilaments()) {
        auto cyl = *f->getCylinderVector().begin();
        coord1 = &coord[3 * cyl->getFirstBead()->getStableIndex()];
        coord2 = &coord[3 * cyl->getSecondBead()->getStableIndex()];
	    cyllengthset[cylcount] = twoPointDistance(coord1,coord2);
        cylcount++;

        for (auto it = f->getCylinderVector().begin()+1;
             it != f->getCylinderVector().end(); it++){
        	coord3 = &coord[3 * (*it)->getSecondBead()->getStableIndex()];
	        cyllengthset[cylcount] = twoPointDistance(coord2,coord3);
	        cylcount++;
	        cylbenddotproduct[dotproductcount] = scalarProduct(coord1,coord2,coord2,coord3);
	        dotproductcount++;
        	coord1 = coord2;
        	coord2 = coord3;
        }
    }*/
}

template <class FBendingInteractionType>
void FilamentStretchingandBending<FBendingInteractionType>::vectorize() {

	// Count number of interactions
	_numInteractions = 0;
	for(auto f : Filament::getFilaments())
		if(f->getCylinderVector().size() > 1) _numInteractions += f->getCylinderVector().size() - 1;

	totalenergy = new floatingpoint[3*SubSystem::tp->numThreads()];
	totalenergy[0] = (floatingpoint)0.0;
	totalenergy[1] = (floatingpoint)0.0;
	totalenergy[2] = (floatingpoint)0.0;
	// The first cylinder in each filament is not considered in the hybrid stretching
	// bending paradigm. Need a function to calculate energies seperately for those
	// cylinders.
	_strnumInteractions = Filament::getFilaments().size();
    //These two numbers should match
/*	cout<<"Str NumInt Method  1 "<<Cylinder::getCylinders().size()<<" Method 2 "
		<<_strnumInteractions+_numInteractions<<endl;*/
	beadSetcylsansbending = new int[nstr * _strnumInteractions];
	kstrsansbending = new floatingpoint[_strnumInteractions];
	eqlsansbending = new floatingpoint[_strnumInteractions];

	kstr = new floatingpoint[_numInteractions];
	eql = new floatingpoint[_numInteractions];
	beadSet = new int[n * _numInteractions];
	kbend = new floatingpoint[_numInteractions];
	eqt = new floatingpoint[_numInteractions];

	//Testing Cylset stores IDs
	cylSet = new int[2*_numInteractions];
	cylSetcylsansbending = new int[_strnumInteractions];

	//Datasets to help calculate precomputed vars.
    auto Totalcyl = Cylinder::getCylinders().size();
    auto Totalbeads = Bead::getBeads().size();
    beadSetall = new int[Totalbeads];
    beadpaircyllengthstatus = new bool[Totalbeads];
    beadtriplethingestatus  = new bool[Totalbeads];
	cyllengthset = new floatingpoint[Totalcyl];
    cylbenddotproduct = new floatingpoint[_numInteractions];
    int beadcount = 0;
    int hingecount = 0;

	int i = 0;

	int istr = 0;

	int cylcount = 0;

	for (auto f: Filament::getFilaments()) {

		auto cyl = *f->getCylinderVector().begin();
		beadSetcylsansbending[nstr * istr] = cyl->getFirstBead()->getStableIndex();
		beadSetcylsansbending[nstr * istr + 1] = cyl->getSecondBead()->getStableIndex();
		kstrsansbending[istr] = cyl->getMCylinder()->getStretchingConst();
		eqlsansbending[istr] = cyl->getMCylinder()->getEqLength();
		cylSetcylsansbending[istr] = cylcount;
		beadSetall[beadcount] = beadSetcylsansbending[nstr * istr];
//        cout<<beadcount<<endl;
		beadcount++;
        beadpaircyllengthstatus[beadcount-1] = 1;
//        cout<<beadcount-1<<endl;
		beadSetall[beadcount] = beadSetcylsansbending[nstr * istr + 1];
		beadcount++;
		cylcount++;
		istr++;

		if (f->getCylinderVector().size() > 1){

			for (auto it = f->getCylinderVector().begin()+1;
			     it != f->getCylinderVector().end(); it++){

				auto it2 = it - 1;
				beadSet[n * i] = (*it2)->getFirstBead()->getStableIndex();
				beadSet[n * i + 1] = (*it)->getFirstBead()->getStableIndex();
				beadSet[n * i + 2] = (*it)->getSecondBead()->getStableIndex();

                cylSet[ncylperint * i] = cylcount - 1;
                cylSet[ncylperint * i + 1] = cylcount;
				cylcount++;

				beadSetall[beadcount] = beadSet[n * i + 2];
				beadpaircyllengthstatus[beadcount-1] = 1;
//                cout<<beadcount-1<<endl;
				beadcount++;

                beadtriplethingestatus[hingecount] = 1;
                hingecount++;

				kbend[i] = (*it)->getMCylinder()->getBendingConst();
				eqt[i]  = (*it)->getMCylinder()->getEqTheta();

				kstr[i] = (*it)->getMCylinder()->getStretchingConst();
				eql[i]  = (*it)->getMCylinder()->getEqLength();

				i++;
			}
		}
		beadpaircyllengthstatus[beadcount-1] = 0;
//        cout<<beadcount-1<<" end"<< endl;
        beadtriplethingestatus[hingecount] = 0;
        beadtriplethingestatus[hingecount+1] = 0;
        hingecount = hingecount + 2;
	}
	/*for(auto f:Filament::getFilaments()){
	    cout<<f->getCylinderVector().size()<<" ";
	}
	cout<<endl;
	for(int i = 0;i < Totalbeads-1; i++){
	    cout<<beadpaircyllengthstatus[i]<<" ";
	}
	cout<<endl;
    for(int i = 0;i < Totalbeads-1; i++){
        cout<<beadtriplethingestatus[i]<<" ";
    }
    cout<<endl;
	cout<<"Ncyl "<<cylcount<<" "<<Cylinder::getCylinders().size()<<" "
                                                                ""<<2*_numInteractions<<endl;
	cout<<endl;*/
}

template<class FBendingInteractionType>
void FilamentStretchingandBending<FBendingInteractionType>::deallocate() {

	delete [] beadSet;
	delete [] kbend;
	delete [] eqt;

	delete [] beadSetcylsansbending;
	delete [] kstr;
	delete [] eql;
	delete [] kstrsansbending;
	delete [] eqlsansbending;
	delete [] totalenergy;

	delete [] cylSetcylsansbending;
	delete [] cylSet;
	delete [] cyllengthset;
	delete [] cylbenddotproduct;
	delete [] beadSetall;
	delete [] beadpaircyllengthstatus;
	delete [] beadtriplethingestatus;
}

//Needs to have a return value.
template <class FStretchingandBendingInteractionType>
floatingpoint FilamentStretchingandBending<FStretchingandBendingInteractionType>::
        computeEnergy(floatingpoint *coord){

	floatingpoint U_ii=0.0;

#ifdef SERIAL

	const int startID = 0;
	int threadID = 0;
	precomputevars(coord, cyllengthset, cylbenddotproduct);
    _FFType.energy(coord, _numInteractions, cylSet, cyllengthset, cylbenddotproduct,
            kstr, kbend, eql, eqt, totalenergy, startID, _numInteractions, threadID);
    _FFType.energy(coord, cylSetcylsansbending, cyllengthset, kstrsansbending,
           eqlsansbending, totalenergy, startID, _strnumInteractions, threadID);

	/*_FFType.energy(coord, _numInteractions, beadSet, kstr, kbend, eql, eqt, totalenergy,
	               startID, _numInteractions, threadID);

    _FFType.energy(coord, beadSetcylsansbending, kstrsansbending, eqlsansbending, totalenergy,
            startID, _strnumInteractions, threadID);*/

#endif
	floatingpoint U = (floatingpoint) 0.0;
	for(int t = 0 ; t < SubSystem::tp->numThreads(); t++){
		for(int j = 0; j <3;j++) {
			if(totalenergy[3*t+j] >= (floatingpoint) 0.0)
				U +=  totalenergy[3*t+j];
			else
				return (floatingpoint) -1.0;
		}
	}

	return U;

}

 template <class FStretchingandBendingInteractionType>
void FilamentStretchingandBending<FStretchingandBendingInteractionType>::computeForces
(floatingpoint *coord, floatingpoint *f) {
#ifdef SERIAL

	_FFType.forces(coord, f, _numInteractions, beadSet, kstr, kbend, eql, eqt);
     const int startID = 0;
     int threadID = 0;
	_FFType.forces(coord, f, beadSetcylsansbending, kstrsansbending, eqlsansbending, startID,
	        _strnumInteractions, threadID);
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

///Template specializations
template floatingpoint FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::computeEnergy(floatingpoint *coord);
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::vectorize();
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::deallocate();
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::precomputevars(
		floatingpoint *coord, floatingpoint *cyllengthset,
		floatingpoint *cylbenddotproduct);


template floatingpoint FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>::computeEnergy(floatingpoint *coord);
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>::computeForces(floatingpoint *coord, floatingpoint *f);
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>::vectorize();
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>::deallocate();
template void FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>::precomputevars(
		floatingpoint *coord, floatingpoint *cyllengthset,
		floatingpoint *cylbenddotproduct);
