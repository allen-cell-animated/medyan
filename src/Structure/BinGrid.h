
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

#ifndef MEDYAN_BinGrid_h
#define MEDYAN_BinGrid_h
#include "common.h"
#include "Composite.h"
#include "Bin.h"
/*!
* BinGrid is a vector of numBins Bins where numBins is provided during instantiation.
* Each Bin in BinGrid is indexed in NeighborListImpl.
* BinGrid pointer belongs to the respective BinGrid.
 * PrototypeBin exists so we can copy properties to all other bins in bingrid.
*/
class BinGrid : public Composite {
private:
//    Bin _prototype_bin; ///< Prototype bin, to be configured
    ///< before initialization
    short _bgType = -1;
public:
    /// Constructor, creates a number of Compartment instances
    BinGrid(int numBins, short bgType): _bgType(bgType) {

        //add children
        for(size_t i=0; i<numBins; ++i)
            addChild(unique_ptr<Component>(new Bin(bgType)));
    }

    /// Get bins that this grid holds
    vector<Bin*> getBins() {

        vector<Bin*> bins;

        for (auto &c : children())
            bins.push_back((Bin*)c.get());

        return bins;
    }

//    /// Get a bin at certain index
    Bin* getBin(int index) {

        return (Bin*)(children()[index].get());
    }

    /// Get name of this bin grid
    virtual string getFullName() const {return string("BinGrid");};

    /// Get the protobin from this grid, in order to configure and then initialize
//    Bin& getProtoBin() {return _prototype_bin;}
//    const Bin& getProtoBin() const {return _prototype_bin;}

    /// Print properties of this grid
    virtual void printSelf() {
        cout << getFullName() << endl;
        cout << "Number of Bin objects: " << numberOfChildren() << endl;
        cout << "Type: " << _bgType <<endl;
        for(auto &c : children())
            c->printSelf();
    }
    ///GetType implementation just returns zero (no CompartmentGrid types yet)
    virtual int getType() {return _bgType;}
};
#endif //MEDYAN_BINGRID_H
