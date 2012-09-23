//
//  test_species_container.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 9/22/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#define DO_THIS_TEST

#ifdef DO_THIS_TEST

#include <iostream>

#include <iostream>
#include "gtest/gtest.h"

#include "SpeciesContainer.h"

using namespace std;
using namespace chem;

TEST(SpeciesPtrContainerVectorTest, Main) {
    SpeciesPtrContainerVector cont;
    Species *arp23 = cont.addSpecies(new SpeciesDiffusing("Arp2/3",55));
    Species *actin = cont.addSpecies<SpeciesDiffusing>("Actin",44);
    Species *profilin = cont.addSpeciesUnique(make_unique<SpeciesDiffusing>("Profilin",33));
    
    EXPECT_EQ(3U, cont.size());
    EXPECT_EQ(3U, cont.species().size());
    
    EXPECT_EQ(arp23,cont.findSpeciesByIndex(0));
    EXPECT_EQ(actin,cont.findSpeciesByIndex(1));
    EXPECT_EQ(profilin,cont.findSpeciesByIndex(2));
    
    EXPECT_EQ(arp23,cont.findSpeciesByName("Arp2/3"));
    EXPECT_EQ(actin,cont.findSpeciesByName("Actin"));
    EXPECT_EQ(profilin,cont.findSpeciesByName("Profilin"));
    
    SpeciesDiffusing another_profilin("Profilin",33);
    EXPECT_EQ(profilin,cont.findSimilarSpecies(another_profilin));
    
    auto molecule = another_profilin.getMolecule();
    EXPECT_EQ(profilin,cont.findSpeciesByMolecule(molecule));
    
    EXPECT_EQ(1,cont.areAllSpeciesUnique());
    
    cont.addSpecies(new SpeciesDiffusing("Arp2/3",66));
    
    EXPECT_EQ(2, cont.removeSpecies("Arp2/3"));
    EXPECT_EQ(2U, cont.size());

    EXPECT_EQ(1, cont.removeSpecies(actin));
    EXPECT_EQ(1U, cont.size());

}


TEST(SpeciesContainerVectorTest, Main) {
    SpeciesContainerVector<SpeciesDiffusing> cont;
    size_t arp23 = cont.addSpecies("Arp2/3",55);
    size_t actin = cont.addSpecies("Actin",44);
    size_t profilin = cont.addSpecies("Profilin",33);
    
    EXPECT_EQ(3U, cont.size());
    EXPECT_EQ(3U, cont.species().size());

    EXPECT_EQ("Arp2/3",cont.findSpecies(arp23).getName());
    EXPECT_EQ("Actin",cont.findSpecies(actin).getName());
    EXPECT_EQ("Profilin",cont.findSpecies(profilin).getName());
    
    EXPECT_EQ("Arp2/3",cont.findSpecies("Arp2/3").getName());
    EXPECT_EQ("Actin",cont.findSpecies("Actin").getName());
    EXPECT_EQ("Profilin",cont.findSpecies("Profilin").getName());

    EXPECT_EQ(profilin,cont.findSpeciesIndex("Profilin"));

    
    SpeciesDiffusing another_profilin("Profilin",33);
    EXPECT_EQ("Profilin",cont.findSimilarSpecies(another_profilin).getName());
    
    auto molecule = another_profilin.getMolecule();
    EXPECT_EQ("Profilin",cont.findSpeciesByMolecule(molecule).getName());
    
    EXPECT_EQ(1,cont.areAllSpeciesUnique());

    cont.addSpecies("Arp2/3",66);

    EXPECT_EQ(2, cont.removeSpecies("Arp2/3"));
    EXPECT_EQ(2U, cont.size());
    
    EXPECT_EQ(1, cont.removeSpecies(actin));
    EXPECT_EQ(1U, cont.size());
    
}


#endif // DO_THIS_TEST