
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

#ifdef TESTING

#  define DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST
#  ifdef DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST

#    include "gtest/gtest.h"

#    include "common.h"
#    include "test_public.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "GController.h"
#    include "SubSystem.h"
#    include "Component.h"

#    include "Triangle.h"
#    include "MTriangle.h"
#    include "Cylinder.h"

#    include "TriangleCylinderExclVolume.h"
#    include "TriangleCylinderBeadExclVolume.h"

namespace {

    class TriangleCylinderVolumeFFTest: public ::testing::Test {
    protected:
        double radius;
        double height;
        test_public::ComponentDummy dummyParent;
        SubSystem s;
        
        Triangle *t;
        array<Vertex*, 3> tv;
        Cylinder *c;
        array<Bead*, 2> cb;

		TriangleCylinderVolumeFFTest():
            dummyParent(0),
            radius(100), height(20) {
            
            test_public::quickSetupPlayground(&s);

            SysParams::GParams.cylinderNumMon.resize(1, 3);

            SysParams::MParams.MemElasticK.resize(1, 400);
            SysParams::MParams.MemBendingK.resize(1, 100);
            SysParams::MParams.MemEqCurv.resize(1, 0);

            double trans = radius + height; // To ensure that all coordinates of all the beads are greater than 0.

            // Add a triangle
            tv[0] = new Vertex({trans + radius, trans, trans}, &dummyParent, 0);
            tv[1] = new Vertex({trans - radius/2, trans + sqrt(3)*radius/2, trans}, &dummyParent, 0);
            tv[2] = new Vertex({trans - radius/2, trans - sqrt(3)*radius/2, trans}, &dummyParent, 0);
            for(Vertex* eachV: tv) eachV->addToSubSystem();
            t = new Triangle(&dummyParent, tv[0], tv[1], tv[2]);
            t->addToSubSystem();

            // Add a cylinder
            cb[0] = new Bead({trans + radius/20, trans + radius/20, trans + height}, &dummyParent, 0);
            cb[1] = new Bead({trans - radius/20, trans - radius/20, trans + height}, &dummyParent, 0);
            for(Bead* eachB: cb) eachB->addToSubSystem();
            c = new Cylinder(&dummyParent, cb[0], cb[1], 0, 0);
            c->addToSubSystem();
        }
        ~TriangleCylinderVolumeFFTest() {
            SysParams::GParams.cylinderNumMon.resize(0);

            SysParams::MParams.MemElasticK.resize(0);
            SysParams::MParams.MemBendingK.resize(0);
            SysParams::MParams.MemEqCurv.resize(0);

            // Remove the triangle
            t->removeFromSubSystem();
            delete t;
            for(Vertex* eachV: tv) {
                eachV->removeFromSubSystem();
                delete eachV;
            }

            // Remove the cylinder
            c->removeFromSubSystem();
            delete c;
            for(Bead* eachB: cb) {
                eachB->removeFromSubSystem();
                delete eachB;
            }
        }

        void recordCoordinate() {
            for(Vertex* it: tv) it->coordinateP = it->coordinate;
            for(Bead* it: cb) it->coordinateP = it->coordinate;
        }
        void resetCoordinate() {
            for(Vertex* it: tv) it->coordinate = it->coordinateP;
            for(Bead* it: cb) it->coordinate = it->coordinateP;
        }
        void resetForce() {
            for(Vertex* it: tv) it->force.assign(3, 0);
            for(Bead* it: cb) it->force.assign(3, 0);
        }
        void assignRandomForceAuxP(double sigma) {
		normal_distribution<> nd(0, sigma);

            for(Vertex* it: tv) {
                for(double& eachForce: it->forceAuxP) { // forceAuxP is used here simply because it is an empty container.
                    eachForce = nd(Rand::engFixed);
                }
            }
            for(Bead* it: cb) {
                for(double& eachForce: it->forceAuxP) { // forceAuxP is used here simply because it is an empty container.
                    eachForce = nd(Rand::engFixed);
                }
            }
        }
        void moveAlongForceAuxP(double d) {
            for(Vertex* it: tv) {
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    it->coordinate[coordIdx] += it->forceAuxP[coordIdx] * d;
                }
            }
            for(Bead* it: cb) {
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    it->coordinate[coordIdx] += it->forceAuxP[coordIdx] * d;
                }
            }
        }

    };

}

TEST_F(TriangleCylinderVolumeFFTest, Force) {
    TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion> tcv;

    assignRandomForceAuxP(radius/200);
    recordCoordinate();

    // Compare the results with force predictions
    // U(x+h) - U(x-h) = dotProduct(2h, dU/dx) = -dotProduct(2h, F)

    resetCoordinate();
    resetForce();
    // Update position for neighbor list search
    t->updatePosition();
    c->updatePosition();
    t->updateGeometry(true);
    tcv->computeForces();

    // Don't bother updating neighbor list here
    moveAlongForceAuxP(1.0);
    t->updateGeometry(true);
    double U1 = tcv->computeEnergy(0.0);

    resetCoordinate();
    moveAlongForceAuxP(-1.0);
    t->updateGeometry(true);
    double U2 = tcv->computeEnergy(0.0);

    double exDiff = 0.0;
    for(Vertex* v: tv)
        exDiff -= 2 * dotProduct(v->forceAuxP, v->force);
    for(Bead* b: cb)
        exDiff -= 2 * dotProduct(b->forceAuxP, b->force);

    EXPECT_NEAR(U1 - U2, exDiff, abs(exDiff / 1000))
        << tcv->getName() << " force not working properly.";

}


#  endif // DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST
#endif // TESTING

