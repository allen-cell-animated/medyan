
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

#ifdef TESTING

#define DO_THIS_REACTION_TEST
#ifdef DO_THIS_REACTION_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "Species.h"
#include "Reaction.h"
#include "Compartment.h"

void rspecies_callback (RSpecies *r, int delta){
    r->getSpecies().setN(33);
}

void reaction_callback (ReactionBase *r){
    r->setRate(0.05);
}

struct ReactionCallback {
    void operator() (ReactionBase *r){
        ++_count;
        r->setRate(1.05);
    }
    int _count;
};

TEST(RSpeciesTest, Main) {
    // Light testing RSpecies by itself, without Reactions
    SpeciesBulk A("A",  10);
    RSpecies& RA(A.getRSpecies());
    EXPECT_EQ(10, RA.getN());
    EXPECT_EQ(&A, &RA.getSpecies());
#ifdef RSPECIES_SIGNALING
    bool is_sig = RA.isSignaling();
    EXPECT_EQ(0, is_sig);
#endif
    EXPECT_EQ(0U, RA.reactantReactions().size());
    EXPECT_EQ(0U, RA.productReactions().size());

    SpeciesBulk B("B",  10);
    RSpecies& RB(B.getRSpecies());
    
    // Using a reaction A->B to see if copy number n is correctly driven
    // (implicitly tests up() and down() methods
    Reaction<1,1> rxn1 = {{&A,&B}, 10.0 };
    rxn1.makeStep();
    rxn1.makeStep();
    EXPECT_EQ(8, RA.getN());
    EXPECT_EQ(12, RB.getN());
    
    // Introducing the back reaction, B->A, to bring the RSpecies back to
    // their original state
    Reaction<1,1> rxn2 = {{&B,&A}, 10.0 };
    rxn2.makeStep();
    rxn2.makeStep();
    EXPECT_EQ(10, RA.getN());
    EXPECT_EQ(10, RB.getN());    
}

TEST(ReactionTest, CTors) {
    SpeciesBulk A("A",  10);
    SpeciesBulk B("B",  10);
    SpeciesBulk C("C",  10);
    SpeciesBulk D("D",  10);
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 }; // A -> B
    Reaction<2,1> rxn2 = { {&A,&B,&C}, 10.0 }; // A+B -> C
    Reaction<1,2> rxn3 = { {&A,&B,&C}, 10.0 }; // A -> B + C
    Reaction<2,2> rxn4 = { {&A,&B,&C,&D}, 10.0 }; // A + B -> C + D

    EXPECT_EQ(1, rxn1.getM()); EXPECT_EQ(1, rxn1.getN());
    EXPECT_EQ(2, rxn2.getM()); EXPECT_EQ(1, rxn2.getN());
    EXPECT_EQ(1, rxn3.getM()); EXPECT_EQ(2, rxn3.getN());
    EXPECT_EQ(2, rxn4.getM()); EXPECT_EQ(2, rxn4.getN());

    // rxn1
    auto it1 = rxn1.rspecies();
    EXPECT_EQ("A{Bulk}",(*it1)->getFullName());
    // rxn2
    auto it2 = rxn2.rspecies();
    EXPECT_EQ("A{Bulk}",(*it2)->getFullName());
    ++it2;
    EXPECT_EQ("B{Bulk}",(*it2)->getFullName());
    EXPECT_EQ("C{Bulk}",(*(rxn2.rspecies()+rxn2.getM()))->getFullName());
    // rxn3
    auto it3 = rxn3.rspecies()+rxn2.getM();
    EXPECT_EQ("A{Bulk}",(*rxn3.rspecies())->getFullName());
    EXPECT_EQ("B{Bulk}",(*it2)->getFullName());
    ++it3;
    EXPECT_EQ("C{Bulk}",(*(rxn2.rspecies()+rxn2.getM()))->getFullName());
    // rxn4
    auto r_it4 = rxn4.rspecies();
    auto p_it4 = rxn4.rspecies()+rxn4.getM();
    EXPECT_EQ("A{Bulk}",(*r_it4)->getFullName());
    ++r_it4;
    EXPECT_EQ("B{Bulk}",(*r_it4)->getFullName());
    EXPECT_EQ("C{Bulk}",(*p_it4)->getFullName());
    ++p_it4;
    EXPECT_EQ("D{Bulk}",(*p_it4)->getFullName());

} 

#ifdef RSPECIES_SIGNALING
TEST(ReactionTest, RSpeciesSignaling) {
    SpeciesBulk A("A",  8);
    A.startSignaling();
    auto c = A.connect(rspecies_callback);
    
    A.getRSpecies().emitSignal(0);
    EXPECT_EQ(33,A.getN());
    
    c.disconnect();
    A.setN(9);
    A.getRSpecies().emitSignal(0);
    EXPECT_EQ(9,A.getN());
}
#endif // of RSPECIES_SIGNALING

#ifdef TRACK_DEPENDENTS
#ifdef TRACK_ZERO_COPY_N
TEST(ReactionTest, Dependents1) {
    SpeciesBulk A("A",  10);
    SpeciesBulk B("B",  10);
    SpeciesBulk C("C",  10);
    RSpecies& RA(A.getRSpecies());
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 };
    Reaction<1,1> rxn3 = { {&A,&C}, 10.0 };
    Reaction<1,1> rxn2 = { {&B,&A}, 10.0 };
    
    rxn1.activateReaction();
    rxn2.activateReaction();
    rxn3.activateReaction();
    
    // note: we have three reactions, (1) A->B, (2) B->A, (3) A->C
    // (1) affects (2) and (3)
    // (2) affects (1) and (3) 
    // (3) affects (1)
    vector<ReactionBase*> ar1 = rxn1.getAffectedReactions();
    vector<ReactionBase*> ar2 = rxn2.getAffectedReactions();
    vector<ReactionBase*> ar3 = rxn3.getAffectedReactions();
    EXPECT_EQ(2U, ar1.size());
    EXPECT_EQ(2U, ar2.size());
    EXPECT_EQ(1U, ar3.size());
    EXPECT_EQ(&rxn1, ar3[0]);// (3) affects (1)
    
    // Testing passivateAssocReacts() and activateAssocReactions()
    // We run A->C 10 times, so A's copy number drops to 0.
    // Hence, the reaction (1) & (3) should be passivated.
    // In the context of this test, where we do not run a Gillespie-like algorithm,
    // this should lead to reaction (2) not having (1) or (3) as dependents - the
    // dependent count of (2) should become zero
    
    for (int i=0; i<10; ++i){
        rxn3.makeStep();
    }
    EXPECT_EQ(0U, A.getN());

    EXPECT_EQ(0U, rxn2.dependents().size());
    
    // But not that the getAffectedReactions() still returns the original
    // depedencices, ignoring the fact that copy number of A is zero.
    ar1 = rxn1.getAffectedReactions();
    ar2 = rxn2.getAffectedReactions();
    ar3 = rxn3.getAffectedReactions();
    EXPECT_EQ(2U, ar1.size());
    EXPECT_EQ(2U, ar2.size());
    EXPECT_EQ(1U, ar3.size());
    EXPECT_EQ(&rxn1, ar3[0]);// (3) affects (1)
    
    // Now let's activate (1) and (3) by moving the copy number of A from 0 to 1
    rxn2.makeStep();
    EXPECT_EQ(1U, RA.getN());
    EXPECT_EQ(2U, rxn2.dependents().size());
}
#endif // of TRACK_ZERO_COPY_N
#endif // of TRACK_UPPER_COPY_N

TEST(ReactionTest, Dependents2) {
    SpeciesBulk A("A",  10);
    SpeciesBulk B("B",  10);
    SpeciesBulk C("C",  10);
    SpeciesBulk D("D",  10);
    SpeciesBulk E("E",  10);
    SpeciesBulk F("F",  10);
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 };
    Reaction<1,1> rxn3 = { {&C,&D}, 10.0 };
    Reaction<1,1> rxn2 = { {&E,&F}, 10.0 };

    rxn1.activateReaction();
    rxn2.activateReaction();
    rxn3.activateReaction();
    
    // Let's add artificial dependencies by hand
    rxn1.registerNewDependent(&rxn2);
    EXPECT_TRUE(rxn1.dependents().find(&rxn2) != rxn1.dependents().end());
    rxn1.registerNewDependent(&rxn3);
    EXPECT_TRUE(rxn1.dependents().find(&rxn3) != rxn1.dependents().end());
    EXPECT_EQ(2U, rxn1.dependents().size());
    rxn1.unregisterDependent(&rxn3);
    EXPECT_EQ(1U, rxn1.dependents().size());
}

TEST(ReactionTest, Propensities) {
    SpeciesBulk A("A",  8);
    SpeciesBulk B("B",  12);
    SpeciesBulk C("C",  14);
    
    Reaction<2,1> rxn = { {&A,&B,&C}, 3.14 }; // A+B -> C
    rxn.activateReaction();
    EXPECT_FLOAT_EQ(8*12*3.14, rxn.computePropensity());
    EXPECT_EQ(8*12, rxn.getProductOfReactants());
    
    Reaction<1,1> rxn2 = { {&A,&B}, 3.14 }; // A->B
    rxn2.activateReaction();
    EXPECT_FLOAT_EQ(8*3.14, rxn2.computePropensity());
    EXPECT_EQ(8, rxn2.getProductOfReactants());

}

#ifdef REACTION_SIGNALING
TEST(ReactionTest, ReactionSignaling) {
    SpeciesBulk A("A",  8);
    SpeciesBulk B("B",  12);
    Reaction<1,1> rxn = { {&A,&B}, 3.14 }; // A->B
    rxn.startSignaling();
    
    // There are at least ways to set up callbacks: 1) functions; 2) functors; 3) lambda
    // functions. In terms of connection management, the return of sm.connect(...) can
    // be captured and used later to temporarily block the callback or permanently
    // disconnect it. See the boost documentation for signals2.
    
    // function<void (Reaction *)> rcb(reaction_callback);

    ConnectionBlock rcb(rxn.connect(reaction_callback), false);
    rxn.emitSignal();
    EXPECT_FLOAT_EQ(0.05, rxn.getRate());
    rcb.block();
    
    ConnectionBlock RCB(rxn.connect(ReactionCallback()), false);
    rxn.emitSignal();
    EXPECT_FLOAT_EQ(1.05, rxn.getRate());
    RCB.block();
    
    boost::signals2::connection clambda =
    rxn.connect([](ReactionBase *r){r->setRate(2.05);});
    ConnectionBlock Rlambda (clambda, false); 
    rxn.emitSignal();
    EXPECT_FLOAT_EQ(2.05, rxn.getRate());
    clambda.disconnect(); // permanently disconnecting this signal

    rcb.unblock();
    rxn.emitSignal();
    EXPECT_FLOAT_EQ(0.05, rxn.getRate());
    rcb.block();
    
    auto c = rxn.connect([](ReactionBase *r){r->setRate(3.05);});
    rxn.emitSignal();
    EXPECT_FLOAT_EQ(3.05, rxn.getRate());
    c.disconnect(); // this permanently disconnects the signal vs blocking
    
}
#endif

TEST(ReactionTest, ReactionCloning) {
    
    Compartment* C1 = new Compartment;
    Compartment* C2 = new Compartment;
    
    Species* ADiff1 = C1->addSpeciesDiffusing("ADiff", 10);
    Species* ADiff2 = C2->addSpeciesDiffusing("ADiff", 10);
    
    Species* BDiff1 = C1->addSpeciesDiffusing("BDiff", 10);
    Species* BDiff2 = C2->addSpeciesDiffusing("BDiff", 10);
    
    ReactionBase* r1 = C1->addInternal<Reaction,1,1>({ADiff1,BDiff1}, 100.0);
#ifdef REACTION_SIGNALING
    auto c = r1->connect([](ReactionBase *r){r->setRate(9.0);});
#endif
    
    ///Clone, check if valid
    ReactionBase* r2 = r1->clone(C2->getSpeciesContainer());
    
    EXPECT_EQ(1, r2->getM());
    EXPECT_EQ(1, r2->getN());
    
    EXPECT_TRUE(r2->containsSpecies(ADiff2));
    EXPECT_TRUE(r2->containsSpecies(BDiff2));
    
    ///Check signal cloning
#ifdef REACTION_SIGNALING
    r2->emitSignal();
    EXPECT_EQ(9.0, r2->getRate());
#endif
    
    ///Clone a reaction where not all species are in compartment
    Compartment* C3 = new Compartment;
    Species* CDiff3 = C3->addSpeciesDiffusing("CDiff", 10);
    Species* ADiff3 = C3->addSpeciesDiffusing("ADiff", 10);
    
    ReactionBase* r3 = C3->addInternal<Reaction,1,1>({ADiff3, CDiff3}, 100.0);
    ReactionBase* r4 = r3->clone(C1->getSpeciesContainer());
    
    ///Should keep cdiff3 in reactants
    EXPECT_TRUE(r4->containsSpecies(CDiff3));
    EXPECT_TRUE(r4->containsSpecies(ADiff1));
    
    ///Check reaction equality
    EXPECT_TRUE(r3->is_equal(*r4));
}


#endif //DO_THIS_REACTION_TEST
#endif //TESTING

