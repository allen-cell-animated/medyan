//
//  example_composite.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 7/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <boost/heap/pairing_heap.hpp>

#include "System.h"
#include "ChemNRMImpl.h"
#include "ChemSim.h"
#include "Signaling.h"

#include "Composite.h"

#include "CompartmentContainer.h"

#include "Visitor.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    
    CompartmentSpatial<3> CC1;
    vector<float> sides{100.0,100.0,100.0};
    CC1.setSides(sides.begin());
    vector<float> coords{12.3,1.2,22.1};
    CC1.setCoords(coords.begin());
    CC1.printSelf();
    cout << endl << endl;
    
    const int NDIM =3;
    const int NGRID = 4;
    CompartmentGrid<NDIM> ccv{NGRID, NGRID, NGRID};
    CompartmentSpatial<NDIM> &Cproto = ccv.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    Cproto.setSides(sides.begin());
    
    ReactionBase *RM1M2 = Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    ReactionBase *RM2M1 = Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    ccv.initialize();
    //ccv.printSelf();
    cout << "Num of Species: " << ccv.countSpecies() << endl;
    cout << "Num of Reactions: " << ccv.countReactions() << endl;
    
    int counter;

    VisitorLambda vl1;
    vl1.setLambda([&counter](Component *c){
        //c->printSelf();
        ++counter;
        return true;});
    
    ccv.apply(vl1);
    cout << "VisitorLambda: visited " << counter << " nodes" << endl;
    
    counter=0;
    VisitorLambda cvl1;
    cvl1.setLambda([&counter](Component *c)
    {
        counter+=c->numberOfSpecies();
        return true;
    });
    ccv.apply(cvl1);
    cout << "VisitorLambda: counted " << counter << " Species" << endl;

    
    counter=0;
    SpeciesVisitorLambda svl1;
    svl1.setLambda( [&counter](Species *s){
        if(counter==23)
            cout << "Visiting Species: i=" << counter << ", " << s << ", Parent ptr=" << s->getParent() << endl;
        ++counter;
        return true;
    });
    ccv.apply(svl1);
    
    cout << endl;

    counter=0;
    svl1.setLambdaPred([](Species *s){return s->getN()==6;});
    svl1.setLambda( [&counter](Species *s){
        if(counter%10==0){
            cout << "Visiting Species: i=" << counter << ", " << s << ", Parent ptr=" << s->getParent() << ", " << *s << endl;
        }
        ++counter;
        return true;
    });
    ccv.apply(svl1);

    counter=0;
    ReactionVisitorLambda rvl1;
    rvl1.setLambdaPred([](ReactionBase *r){return r->getRate()<50.0;});
    rvl1.setLambda( [&counter](ReactionBase *r){
        if(counter%10==0){
            cout << "Visiting Reaction: i=" << counter << ", " << r << ", Parent ptr=" << r->getParent() << ",\n" << *r;
        }
        ++counter;
        return true;
    });
    ccv.apply(rvl1);
    
    cout << endl;

    counter=0;
    TypeVisitorLambda<Compartment> uvl1;
    uvl1.setLambda( [&counter](Component *c){
        if(counter%10==0){
            Compartment *C = static_cast<Compartment*>(c);
            cout << "Visiting Compartment: i=" << counter << ", has " << C->countReactions() << " reactions." << endl;
        }
        ++counter;
        return true;
    });
    ccv.apply(uvl1);
    
    cout << "Main exited..." << endl;
    return 0;
}

