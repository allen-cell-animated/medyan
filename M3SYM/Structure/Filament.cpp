
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include <time.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/range/numeric.hpp>
#include <fstream>

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SubSystem.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<Filament*> Filament::_filaments;


Filament::Filament(SubSystem* s, vector<double>& position,
                                 vector<double>& direction,
                                 bool nucleation, bool branch)

    : Trackable(), _subSystem(s), _ID(_filaments.getID()) {
 
    //create beads
    Bead* b1 = _subSystem->addTrackable<Bead>(position, this, 0);
    
    //choose length
    double length = 0.0;
    
    if(branch)          length = SysParams::Geometry().monomerSize;
    else if(nucleation) length = SysParams::Geometry().minCylinderSize;
    
    auto pos2 = nextPointProjection(position, length, direction);
        
    Bead* b2 = _subSystem->addTrackable<Bead>(pos2, this, 1);
    
    //create cylinder
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b1, b2, 0);
    
    c0->setPlusEnd(true);
    c0->setMinusEnd(true);
    _cylinderVector.push_back(c0);
}


Filament::Filament(SubSystem* s, vector<vector<double> >& position,
                   int numBeads, string projectionType)

    : Trackable(), _subSystem(s), _ID(_filaments.getID()) {
    
    //create a projection of beads
    vector<vector<double>> tmpBeadsCoord;
    
    //straight projection
    if(projectionType == "STRAIGHT")
        tmpBeadsCoord = straightFilamentProjection(position, numBeads);
    //zigzag projection
    else if(projectionType == "ZIGZAG")
        tmpBeadsCoord = zigZagFilamentProjection(position, numBeads);
    //arc projection
    else if(projectionType == "ARC")
        tmpBeadsCoord = arcFilamentProjection(position, numBeads);
   
    //create beads
    auto direction = twoPointDirection(tmpBeadsCoord[0], tmpBeadsCoord[1]);
        
    Bead* b1 = _subSystem->addTrackable<Bead>(tmpBeadsCoord[0], this, 0);
    Bead* b2 = _subSystem->addTrackable<Bead>(tmpBeadsCoord[1], this, 1);
    
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b1, b2, 0,
                                                      false, false, true);
        
    c0->setPlusEnd(true);
    c0->setMinusEnd(true);
    _cylinderVector.push_back(c0);
    
    for (int i = 2; i<numBeads; i++)
        extendFront(tmpBeadsCoord[i]);
}

Filament::~Filament() {
    
    //remove cylinders, beads from system
    for(auto &c : _cylinderVector) {
        //delete bead
        _subSystem->removeTrackable<Bead>(c->getFirstBead());
        
        //remove second bead if last
        if(c->isPlusEnd())
            _subSystem->removeTrackable<Bead>(c->getSecondBead());
        
        //delete cylinder
        _subSystem->removeTrackable<Cylinder>(c);
    }
}


//Extend front for initialization
void Filament::extendFront(vector<double>& coordinates) {
    
    Cylinder* cBack = _cylinderVector.back();
    cBack->setPlusEnd(false);
    
    int lastPositionFilament = cBack->getPositionFilament();
    
    Bead* b2 = cBack->getSecondBead();
    
    //create a new bead
    auto direction = twoPointDirection(b2->coordinate, coordinates);
    
    auto newBeadCoords = nextPointProjection(b2->coordinate,
                         SysParams::Geometry().cylinderSize, direction);
    
    //create bead
    Bead* bNew = _subSystem->addTrackable<Bead>
                 (newBeadCoords, this, b2->getPositionFilament() + 1);
    
    //create cylinder
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>
                   (this, b2, bNew, lastPositionFilament + 1,
                    false, false, true);
    
    c0->setPlusEnd(true);
    _cylinderVector.push_back(c0);
    
}

//Extend back for initialization
void Filament::extendBack(vector<double>& coordinates) {

    Cylinder* cFront = _cylinderVector.front();
    cFront->setMinusEnd(false);
    
    int lastPositionFilament = cFront->getPositionFilament();
    
    Bead* b2 = cFront->getFirstBead();
    
    //create a new bead
    auto direction = twoPointDirection(b2->coordinate, coordinates);
    
    auto newBeadCoords = nextPointProjection(b2->coordinate,
                         SysParams::Geometry().cylinderSize, direction);
    
    //create bead
    Bead* bNew = _subSystem->addTrackable<Bead>
                 (newBeadCoords, this, b2->getPositionFilament() - 1);
    
    //create cylinder
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>
                   (this, bNew, b2, lastPositionFilament - 1,
                    false, false, true);
    
    c0->setMinusEnd(true);
    _cylinderVector.push_front(c0);

}

//extend front at runtime
void Filament::extendFront(short plusEnd) {

    Cylinder* cBack = _cylinderVector.back();
    
    int lastPositionFilament = cBack->getPositionFilament();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    //move last bead of last cylinder forward
    auto direction1 = twoPointDirection(b1->coordinate, b2->coordinate);
    
    auto npp = nextPointProjection(b2->coordinate,
               SysParams::Geometry().monomerSize, direction1);
    
    //create a new bead in same place as b2
    Bead* bNew = _subSystem->addTrackable<Bead>
                 (npp, this, b2->getPositionFilament() + 1);
    
#ifdef MECHANICS
    //transfer the same load force to new bead
    //(approximation until next minimization)
    bNew->loadForce = b2->loadForce;
#endif
    
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>
                   (this, b2, bNew, lastPositionFilament + 1, true);
    
    _cylinderVector.back()->setPlusEnd(false);
    _cylinderVector.push_back(c0);
    _cylinderVector.back()->setPlusEnd(true);
    
#ifdef CHEMISTRY
    //get last cylinder, mark species
    CMonomer* m = _cylinderVector.back()->getCCylinder()->getCMonomer(0);
    m->speciesPlusEnd(plusEnd)->up();
#endif
    
#ifdef DYNAMICRATES
    //update reaction rates
    _cylinderVector.back()->updateReactionRates();
#endif
    
    _deltaPlusEnd++;
}

//extend back at runtime
void Filament::extendBack(short minusEnd) {
    
    Cylinder* cFront = _cylinderVector.front();
    int lastPositionFilament = cFront->getPositionFilament();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();
    
    //move last bead of last cylinder forward
    auto direction1 = twoPointDirection(b1->coordinate, b2->coordinate);
    
    auto npp = nextPointProjection(b2->coordinate,
               SysParams::Geometry().monomerSize, direction1);
    
    //create a new bead in same place as b2
    Bead* bNew = _subSystem->addTrackable<Bead>
                 (npp, this, b2->getPositionFilament() - 1);

#ifdef MECHANICS
    //transfer the same load force to new bead
    //(approximation until next minimization)
    bNew->loadForce = b2->loadForce;
#endif
    
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>
                   (this, bNew, b2, lastPositionFilament - 1, false, true);
    
    _cylinderVector.front()->setMinusEnd(false);
    _cylinderVector.push_front(c0);
    _cylinderVector.front()->setMinusEnd(true);
    
#ifdef CHEMISTRY
    //get first cylinder, mark species
    auto newCCylinder = getCylinderVector().front()->getCCylinder();
    CMonomer* m = newCCylinder->getCMonomer(newCCylinder->getSize() - 1);
    
    m->speciesMinusEnd(minusEnd)->up();
#endif
    
#ifdef DYNAMICRATES
    //update reaction rates
    _cylinderVector.front()->updateReactionRates();
#endif
    
    _deltaMinusEnd++;
}

//Depolymerize front at runtime
void Filament::retractFront() {
    
    Cylinder* retCylinder = _cylinderVector.back();
    _cylinderVector.pop_back();
    
    _subSystem->removeTrackable<Bead>(retCylinder->getSecondBead());
    _subSystem->removeTrackable<Cylinder>(retCylinder);
    
    _cylinderVector.back()->setPlusEnd(true);

#ifdef DYNAMICRATES
    //update rates of new front
    _cylinderVector.back()->updateReactionRates();
#endif
    
    _deltaPlusEnd--;
}

void Filament::retractBack() {
    
    Cylinder* retCylinder = _cylinderVector.front();
    _cylinderVector.pop_front();
    
    _subSystem->removeTrackable<Bead>(retCylinder->getFirstBead());
    _subSystem->removeTrackable<Cylinder>(retCylinder);
    
    _cylinderVector.front()->setMinusEnd(true);
    
#ifdef DYNAMICRATES
    //update rates of new back
    _cylinderVector.front()->updateReactionRates();
#endif
    
    _deltaMinusEnd--;
}

void Filament::polymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    auto direction = twoPointDirection(b1->coordinate, b2->coordinate);
    
    b2->coordinate = nextPointProjection(b2->coordinate,
                     SysParams::Geometry().monomerSize, direction);
    
#ifdef MECHANICS
    //increase eq length, update
    cBack->getMCylinder()->setEqLength(
        cBack->getMCylinder()->getEqLength() +
        SysParams::Geometry().monomerSize);
#endif
}

void Filament::polymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b1 = cFront->getFirstBead();
    Bead* b2 = cFront->getSecondBead();

    auto direction = twoPointDirection(b2->coordinate, b1->coordinate);
    
    b1->coordinate = nextPointProjection(b1->coordinate,
                     SysParams::Geometry().monomerSize, direction);

#ifdef MECHANICS
    //increase eq length, update
    cFront->getMCylinder()->setEqLength(
        cFront->getMCylinder()->getEqLength() +
        SysParams::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();

    auto direction = twoPointDirection(b2->coordinate, b1->coordinate);
    
    b2->coordinate = nextPointProjection(b2->coordinate,
                     SysParams::Geometry().monomerSize, direction);
    
#ifdef MECHANICS
    //decrease eq length, update
    cBack->getMCylinder()->setEqLength(
        cBack->getMCylinder()->getEqLength() -
        SysParams::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b1 = cFront->getFirstBead();
    Bead* b2 = cFront->getSecondBead();
    
    auto direction = twoPointDirection(b1->coordinate, b2->coordinate);
    
    b1->coordinate = nextPointProjection(b1->coordinate,
                     SysParams::Geometry().monomerSize, direction);
    
#ifdef MECHANICS
    //decrease eq length, update
    cFront->getMCylinder()->setEqLength(
        cFront->getMCylinder()->getEqLength() -
        SysParams::Geometry().monomerSize);
#endif
}


void Filament::nucleate(short plusEnd, short filament, short minusEnd) {
    
#ifdef CHEMISTRY
    //chemically initialize species
    CCylinder* cc = _cylinderVector[0]->getCCylinder();
    int monomerPosition = SysParams::Geometry().cylinderIntSize / 2 + 1;
    
    CMonomer* m1 = cc->getCMonomer(monomerPosition - 1);
    CMonomer* m2 = cc->getCMonomer(monomerPosition);
    CMonomer* m3 = cc->getCMonomer(monomerPosition + 1);
    
    //minus end
    m1->speciesMinusEnd(minusEnd)->up();
    
    //filament
    m2->speciesFilament(filament)->up();
    m2->speciesBound(brancherBindingSite)->up();
    m2->speciesBound(linkerBindingSite)->up();
    m2->speciesBound(motorBindingSite)->up();
    
    //plus end
    m3->speciesPlusEnd(plusEnd)->up();
#endif
}


Filament* Filament::sever(int cylinderPosition) {
    
    int vectorPosition = 0;
    
    //loop through cylinder vector, find position
    for(auto &c : _cylinderVector) {
        
        if(c->getPositionFilament() == cylinderPosition) break;
        else vectorPosition++;
    }
    
    //if vector position is zero, we can't sever. return null
    if(vectorPosition == 0) return nullptr;
    
#ifdef CHEMISTRY
    //if any of the cylinders are only one monomer long, we can't sever
    CCylinder* ccBack  = _cylinderVector[vectorPosition - 1]->getCCylinder();
    CMonomer* cm = ccBack->getCMonomer(ccBack->getSize() - 1);
    
    if(cm->activeSpeciesMinusEnd() != -1) return nullptr;
#endif
    
    //create a new filament
    Filament* newFilament = _subSystem->addTrackable<Filament>(_subSystem);
    
    //Split the cylinder vector at position, transfer cylinders to new filament
    for(int i = vectorPosition; i > 0; i--) {
        
        Cylinder* c = _cylinderVector.front();
        _cylinderVector.pop_front();
        
        c->setFilament(newFilament);
        newFilament->_cylinderVector.push_back(c);
    }
    //new front of new filament, back of old
    auto c1 = newFilament->_cylinderVector.back();
    auto c2 = _cylinderVector.front();
    
    ///copy bead at severing point, attach to new filament
    Bead* newB = _subSystem->addTrackable<Bead>
                 (*(c2->getFirstBead()));
    Bead* oldB = c2->getFirstBead();
    
    //offset these beads by a little for safety
    auto msize = SysParams::Geometry().monomerSize;
    
    vector<double> offsetCoord = {(randomInteger(0,1) ? -1 : +1) *
                                   randomDouble(msize, 2 * msize),
                                  (randomInteger(0,1) ? -1 : +1) *
                                   randomDouble(msize, 2 * msize),
                                  (randomInteger(0,1) ? -1 : +1) *
                                   randomDouble(msize, 2 * msize)};
    
    oldB->coordinate[0] += offsetCoord[0];
    oldB->coordinate[1] += offsetCoord[1];
    oldB->coordinate[2] += offsetCoord[2];
    
    newB->coordinate[0] += -offsetCoord[0];
    newB->coordinate[1] += -offsetCoord[1];
    newB->coordinate[2] += -offsetCoord[2];
    
    c1->setSecondBead(newB);
    
    //set plus and minus ends
    c1->setPlusEnd(true);
    c2->setMinusEnd(true);
    
#ifdef CHEMISTRY
    //mark the plus and minus ends of the new and old filament
    CCylinder* cc1 = c1->getCCylinder();
    CCylinder* cc2 = c2->getCCylinder();
    
    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    
    short filamentInt1 = m1->activeSpeciesFilament();
    short filamentInt2 = m2->activeSpeciesFilament();
    
    //plus end
    m1->speciesFilament(filamentInt1)->down();
    m1->speciesPlusEnd(filamentInt1)->up();
    m1->speciesBound(brancherBindingSite)->down();
    m1->speciesBound(linkerBindingSite)->down();
    m1->speciesBound(motorBindingSite)->down();
    
    //minus end
    m2->speciesFilament(filamentInt2)->down();
    m2->speciesMinusEnd(filamentInt2)->up();
    m2->speciesBound(brancherBindingSite)->down();
    m2->speciesBound(linkerBindingSite)->down();
    m2->speciesBound(motorBindingSite)->down();
    
    //remove any cross-cylinder rxns between these two cylinders
    cc1->removeCrossCylinderReactions(cc2);
    cc2->removeCrossCylinderReactions(cc1);
#endif
    
    return newFilament;
}

vector<vector<double>> Filament::straightFilamentProjection(
                       vector<vector<double>>& v, int numBeads) {
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/twoPointDistance(v[1], v[0]);
    tau[0] = invD * ( v[1][0] - v[0][0] );
    tau[1] = invD * ( v[1][1] - v[0][1] );
    tau[2] = invD * ( v[1][2] - v[0][2] );
    
    for (int i = 0; i<numBeads; i++) {
        
        tmpVec[0] = v[0][0] + SysParams::Geometry().cylinderSize * i * tau[0];
        tmpVec[1] = v[0][1] + SysParams::Geometry().cylinderSize * i * tau[1];
        tmpVec[2] = v[0][2] + SysParams::Geometry().cylinderSize * i * tau[2];
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}

vector<vector<double>> Filament::zigZagFilamentProjection(
                       vector<vector<double>>& v, int numBeads){
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/twoPointDistance(v[1], v[0]);
    tau[0] = invD * ( v[1][0] - v[0][0] );
    tau[1] = invD * ( v[1][1] - v[0][1] );
    tau[2] = invD * ( v[1][2] - v[0][2] );
    
    vector<double> perptau = {-tau[1], tau[0], tau[2]};
    
    
    for (int i = 0; i<numBeads; i++) {
        
        if(i%2 == 0) {
            tmpVec[0] = v[0][0] +
            SysParams::Geometry().cylinderSize * i * tau[0];
            tmpVec[1] = v[0][1] +
            SysParams::Geometry().cylinderSize * i * tau[1];
            tmpVec[2] = v[0][2] +
            SysParams::Geometry().cylinderSize * i * tau[2];
        }
        else {
            tmpVec[0] = v[0][0] +
            SysParams::Geometry().cylinderSize * i * perptau[0];
            tmpVec[1] = v[0][1] +
            SysParams::Geometry().cylinderSize * i * perptau[1];
            tmpVec[2] = v[0][2] +
            SysParams::Geometry().cylinderSize * i * perptau[2];
        }
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}

/// Generate a random number
/// @note - created by Aravind 12/2014
double random_g() {
    // HERE
    //return ran3(&time(0));
    return 0.0;
}

/// Create a projection
/// @note - created by Aravind 12/2014
void marsagila(vector<double>&v)
{
    double d1,d2,d3;
    double *x=new double[3];
    d1=2*random_g()-1;
    d2=2*random_g()-1;
    d3=pow(d1,2)+pow(d2,2);
    while(d3>=1)
    {
        d1=2*random_g()-1;
        d2=2*random_g()-1;
        d3=pow(d1,2)+pow(d2,2);
    }
    x[0]=2.0*d1*pow((1.0-d3),0.5);
    x[1]=2.0*d2*pow(1.0-d3,0.5);
    x[2]=1.0-2.0*d3;
    v[0]=2.0*d1*pow((1.0-d3),0.5);
    v[1]=2.0*d2*pow(1.0-d3,0.5);
    v[2]=1.0-2.0*d3;
}

/// Matrix multiply
/// @note - created by Aravind 12/2014
void matrix_mul(boost::numeric::ublas::matrix<double>&X,
                boost::numeric::ublas::matrix<double>&Y,
                boost::numeric::ublas::matrix<double>&Z,
                vector<double>&x,vector<double>&y,
                vector<double>&z,int nbeads,
                vector<vector<double>> &coordinate)
{
    int t,i;
    double dt,length,cyl_length,sum;
    vector<int> id;
    vector<double> dx,dy,dz,dx2,dy2,dz2,length2,dxdy2,dummyy(3);
    using namespace boost::numeric::ublas;
    matrix<double> B(4,4),B2(1,4),temp1(4,1),dummy(1,1),temp2(4,1),temp3(4,1);
    // B
    B(0,0)=1; B(0,1)=0; B(0,2)=0; B(0,3)=0;
    B(1,0)=-3; B(1,1)=3; B(1,2)=0; B(1,3)=0;
    B(2,0)=3; B(2,1)=-6; B(2,2)=3; B(2,3)=0;
    B(3,0)=-1; B(3,1)=3; B(3,2)=-3; B(3,3)=1;
    //
    axpy_prod(B,X,temp1);
    axpy_prod(B,Y,temp2);
    axpy_prod(B,Z,temp3);
    B2(0,0)=1;
    for(t=0;t<=4000;t++)
    {
        dt=0.00025*t;
        B2(0,1)=dt;
        B2(0,2)=dt*dt;
        B2(0,3)=dt*dt*dt;
        axpy_prod(B2,temp1,dummy);
        x.push_back(dummy(0,0));
        axpy_prod(B2,temp2,dummy);
        y.push_back(dummy(0,0));
        axpy_prod(B2,temp3,dummy);
        z.push_back(dummy(0,0));
    }
    adjacent_difference(x.begin(),x.end(),back_inserter(dx));//dx
    adjacent_difference(y.begin(),y.end(),back_inserter(dy));//dy
    adjacent_difference(z.begin(),z.end(),back_inserter(dz));//dz
    
    transform(dx.begin(), dx.end(),dx.begin(),back_inserter(dx2), multiplies<double>());
    transform(dy.begin(), dy.end(),dy.begin(),back_inserter(dy2), multiplies<double>());
    transform(dz.begin(), dz.end(),dz.begin(),back_inserter(dz2), multiplies<double>());
    
    //array of sum(dx^2+dy^2)
    transform(dx2.begin(),dx2.end(),dy2.begin(),
              back_inserter(dxdy2),plus<double>());
    //array of sum(dx^2+dy^2+dz^2)
    transform(dxdy2.begin(),dxdy2.end(),dz2.begin(),
              back_inserter(length2),plus<double>());
    
    std::vector<double> tempLength;
    for(auto x: length2) tempLength.push_back(sqrt(x));
    length2 = tempLength;
    
    length2[0]=0.0;
    
    length = boost::accumulate(length2, 0.0);//arc length.
    //making equal divisions.
    i=0;sum=0.0;id.push_back(0.0);
    cyl_length=length/(nbeads-1);
    while(i<=4000)
    {
        sum+=length2[i];
        if(sum>=cyl_length||i==4000)
        {
            id.push_back(i);
            sum=0.0;
        }
        i++;
    }
    for(i=0;i<id.size();i++)
    {
        dummyy[0]=x[id[i]];
        dummyy[1]=y[id[i]];
        dummyy[2]=z[id[i]];
        coordinate.push_back(dummyy);
    }
}

/// @note - Created by Aravind 12/2014
vector<vector<double>> Filament::arcFilamentProjection(
                       vector<vector<double>>& v, int numBeads) {
    
    using namespace boost::numeric::ublas;
    
    double temp4,mod;
    double temp[3];
    std::vector<double> X3,x3(3),x4(3),X4,x,y,z;
    matrix<double> C(3,3),B(4,4),X(4,1),Y(4,1),Z(4,1);
    std::vector< std::vector<double> > coordinates;
    temp4=0;

    marsagila(x3);
    marsagila(x4);
    
    std::transform(v[0].begin(), v[0].end(), v[1].begin(), temp, std::minus<double>());
    //vector difference.
    std::transform(temp, temp+3, temp, temp, std::multiplies<double>());
    //squaring elements of a vector.
    mod=std::accumulate(temp,temp+3,temp4,std::plus<double>()); //modulus^2.
    mod=pow(mod,0.5); //MODULUS.
    mod=mod*0.5;
    
    std::transform(x3.begin(), x3.end(),x3.begin(),
                   std::bind1st(std::multiplies<double>(), mod));
    std::transform(x4.begin(), x4.end(),x4.begin(),
                   std::bind1st(std::multiplies<double>(), mod));
    //transform point on sphere to the one centred at x1.
    std::transform(x3.begin(), x3.end(), v[0].begin(),
                   std::back_inserter(X3) , std::plus<double>());
    //transform sphere centre to x2.
    std::transform(x4.begin(), x4.end(), v[1].begin(),
                   std::back_inserter(X4) , std::plus<double>());
    
    // The four points that the beizer curve passes through.
    X(0,0)=v[0][0]; X(1,0)=X3[0]; X(2,0)=X4[0]; X(3,0)=v[1][0];
    Y(0,0)=v[0][1]; Y(1,0)=X3[1]; Y(2,0)=X4[1]; Y(3,0)=v[1][1];
    Z(0,0)=v[0][2]; Z(1,0)=X3[2]; Z(2,0)=X4[2]; Z(3,0)=v[1][2];
    //
    matrix_mul(X,Y,Z,x,y,z,numBeads,coordinates);
    return coordinates;
}

void Filament::printInfo() {
    
    cout << endl;
    
    cout << "Filament: ptr = " << this << endl;
    cout << "Filament ID = " << _ID << endl;
    
    cout << endl;
    cout << "Cylinder information..." << endl;
    
    for(auto c : _cylinderVector)
        c->printInfo();
    
    cout << endl;
    
}

bool Filament::isConsistent() {
    
#ifdef CHEMISTRY
    //check consistency of each individual cylinder
    for(auto &c : _cylinderVector) {
        if(!c->getCCylinder()->isConsistent()) {
         
            cout << "Cylinder at position " << c->getPositionFilament()
                 << " is chemically inconsistent" << endl;
            return false;
        }
    }

    //check that it only has one plus end and one minus end
    int numPlusEnd = 0;
    int numMinusEnd = 0;
        
    for(auto &c : _cylinderVector) {
        
        for(int i = 0; i < c->getCCylinder()->getSize(); i++) {
            auto m = c->getCCylinder()->getCMonomer(i);
            
            if(m->activeSpeciesPlusEnd() != -1) numPlusEnd++;
            if(m->activeSpeciesMinusEnd() != -1) numMinusEnd++;
        }
    }
    if(numPlusEnd != 1) {
        cout << "This filament has more than one plus end species." << endl;
        return false;
    }
    if(numMinusEnd != 1) {
        cout << "This filament has more than one minus end species." << endl;
        return false;
    }
     
#endif
    return true;
}


species_copy_t Filament::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto f : _filaments.getElements()) {
        
        //loop through the filament
        for(auto c : f->_cylinderVector) {
            
            for(int i = 0; i < c->getCCylinder()->getSize(); i++) {
                auto m = c->getCCylinder()->getCMonomer(i);
                
                //filament species
                int activeIndex = m->activeSpeciesFilament();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesFilament(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
                //plus end species
                activeIndex = m->activeSpeciesPlusEnd();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesPlusEnd(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
                //minus end species
                activeIndex = m->activeSpeciesMinusEnd();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesMinusEnd(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
            }
        }
    }
    return copyNum;
}



