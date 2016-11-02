
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

//#define DO_THIS_COMPOSITE_TEST
#ifdef DO_THIS_COMPOSITE_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "Composite.h"

TEST(ComponentTest, SimpleMethods) {
    Component cp;
    EXPECT_EQ(nullptr, cp.getParent());
    EXPECT_EQ(true, cp.isRoot());
    EXPECT_EQ(0U, cp.numberOfChildren());
    EXPECT_EQ(0, cp.isComposite());
    EXPECT_EQ(0, cp.isSpeciesContainer());
    EXPECT_EQ(0, cp.isReactionsContainer() );
    EXPECT_EQ(0U, cp.numberOfSpecies());
    EXPECT_EQ(0U, cp.numberOfReactions());
    EXPECT_EQ(0U, cp.countSpecies());
    EXPECT_EQ(0U, cp.countReactions());
    EXPECT_EQ("Component", cp.getFullName());
}

TEST(CompositeTest, SimpleMethods) {
    Composite cp;
    EXPECT_EQ(nullptr, cp.getParent());
    EXPECT_EQ(true, cp.isRoot());
    EXPECT_EQ(0U, cp.numberOfChildren());
    EXPECT_EQ(1, cp.isComposite());
    EXPECT_EQ(0, cp.isSpeciesContainer());
    EXPECT_EQ(0, cp.isReactionsContainer() );
    EXPECT_EQ(0U, cp.numberOfSpecies());
    EXPECT_EQ(0U, cp.numberOfReactions());
    EXPECT_EQ(0U, cp.countSpecies());
    EXPECT_EQ(0U, cp.countReactions());
    EXPECT_EQ("Composite", cp.getFullName());
    EXPECT_EQ(0U, cp.children().size());
}

TEST(CompositeTest, ChildrenMethods) {
    auto X = make_unique<Composite>();
    Component *Xptr = X.get();
    X->addChild(make_unique<Component>());
    Component *chA = X->children(0);
    X->addChild(make_unique<Component>());
    Component *chB = X->children(1);
    X->addChild(make_unique<Component>());
    Component *chC = X->children(2);
    EXPECT_EQ(3U, X->numberOfChildren());
    X->removeChild(chA);
    EXPECT_EQ(2U, X->numberOfChildren());
    EXPECT_EQ(chB, X->children(0));
    EXPECT_EQ(chC, X->children(1));
    
    EXPECT_EQ(true, X->isRoot());
    EXPECT_EQ(X.get(), chB->getParent());
    EXPECT_EQ(X.get(), chC->getParent());
    
    auto Y = make_unique<Composite>();
    Y->addChild(move(X));
    Y->addChild(make_unique<Component>());
    Y->addChild(make_unique<Component>());
    EXPECT_EQ(3U, Y->numberOfChildren());
    EXPECT_EQ(5U, Y->countDescendents());
    
    EXPECT_EQ(Y.get(), chC->getRoot());
    EXPECT_EQ(Y.get(), Xptr->getRoot());
    EXPECT_EQ(Xptr, chB->getParent());
    
    Y->removeChild(Xptr);
    EXPECT_EQ(2U, Y->numberOfChildren());
    EXPECT_EQ(2U, Y->countDescendents());
}

#endif //DO_THIS_COMPOSITE_TEST
#endif //TESTING

