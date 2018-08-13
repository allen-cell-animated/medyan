cmake_minimum_required(VERSION 2.8)
project(Medyan3_2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# see https://cmake.org/cmake/help/latest/module/FindBoost.html
find_package(Boost REQUIRED system thread timer chrono)
include_directories(${Boost_INCLUDE_DIRS})
# include_directories(/usr/local/include/boost)

include_directories(/usr/local/include/)

# include_directories(/usr/local/include/boost)
# include_directories(/usr/local/include/boost/math)
# include_directories(/usr/local/include/boost/pool)
set(CMAKE_CXX_STANDARD 11)
add_definitions(-DCHEMISTRY)
add_definitions(-DMECHANICS)
add_definitions(-DDYNAMICRATES)
add_definitions(-DBOOST_MEM_POOL)
add_definitions(-DBOOL_POOL_NSIZE=65536)
add_definitions(-DTESTING)
add_definitions(-DTRACK_DEPENDENTS)
add_definitions(-DTRACK_ZERO_COPY_N)
add_definitions(-DTRACK_UPPER_COPY_N)
add_definitions(-DREACTION_SIGNALING)
add_definitions(-DRSPECIES_SIGNALING)
#add_definitions(-DDEBUGCONSTANTSEED)
set(SRC_DIR /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src)
set(SOURCE_FILES
        ${SRC_DIR}/Chemistry/BindingManager.cpp
        ${SRC_DIR}/Chemistry/BindingManager.h
        ${SRC_DIR}/Chemistry/ChemCallbacks.h
        ${SRC_DIR}/Chemistry/ChemGillespieImpl.cpp
        ${SRC_DIR}/Chemistry/ChemGillespieImpl.h
        ${SRC_DIR}/Chemistry/ChemManager.cpp
        ${SRC_DIR}/Chemistry/ChemManager.h
        ${SRC_DIR}/Chemistry/ChemNRMImpl.cpp
        ${SRC_DIR}/Chemistry/ChemNRMImpl.h
        ${SRC_DIR}/Chemistry/ChemRNode.h
        ${SRC_DIR}/Chemistry/ChemSim.cpp
        ${SRC_DIR}/Chemistry/ChemSim.h
        ${SRC_DIR}/Chemistry/ChemSimImpl.h
        ${SRC_DIR}/Chemistry/ChemSimpleGillespieImpl.cpp
        ${SRC_DIR}/Chemistry/ChemSimpleGillespieImpl.h
        ${SRC_DIR}/Chemistry/Reaction.cpp
        ${SRC_DIR}/Chemistry/Reaction.h
        ${SRC_DIR}/Chemistry/ReactionBase.cpp
        ${SRC_DIR}/Chemistry/ReactionBase.h
        ${SRC_DIR}/Chemistry/ReactionContainer.h
        ${SRC_DIR}/Chemistry/ReactionTemplate.cpp
        ${SRC_DIR}/Chemistry/ReactionTemplate.h
        ${SRC_DIR}/Chemistry/RSpecies.cpp
        ${SRC_DIR}/Chemistry/RSpecies.h
        ${SRC_DIR}/Chemistry/Species.cpp
        ${SRC_DIR}/Chemistry/Species.h
        ${SRC_DIR}/Chemistry/SpeciesContainer.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryBubbleRepulsion.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryBubbleRepulsion.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryBubbleRepulsionExp.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryBubbleRepulsionExp.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderAttachment.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderAttachment.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderAttachmentHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderAttachmentHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsion.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsion.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionIn.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionIn.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExpIn.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExpIn.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryFF.h
        ${SRC_DIR}/Mechanics/ForceField/Boundary/BoundaryInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingBending.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingBending.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingBendingCosine.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingBendingCosine.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingDihedral.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingDihedral.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingDihedralCosine.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingDihedralCosine.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingFF.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingPosition.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingPosition.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingPositionCosine.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingPositionCosine.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingStretching.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingStretching.h
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingStretchingHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Branching/BranchingStretchingHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleBubbleRepulsionExp.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleBubbleRepulsionExp.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleCylinderRepulsionExp.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleCylinderRepulsionExp.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleFF.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/BubbleInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/MTOCAttachment.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/MTOCAttachment.h
        ${SRC_DIR}/Mechanics/ForceField/Bubble/MTOCAttachmentHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Bubble/MTOCAttachmentHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBending.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBending.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBendingCosine.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBendingCosine.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBendingHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentBendingHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentFF.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentStretching.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentStretching.h
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentStretchingHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Filament/FilamentStretchingHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerFF.h
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerStretching.cpp
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerStretching.h
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerStretchingHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/Linker/LinkerStretchingHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostFF.h
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostStretching.cpp
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostStretching.h
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostStretchingHarmonic.cpp
        ${SRC_DIR}/Mechanics/ForceField/MotorGhost/MotorGhostStretchingHarmonic.h
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderExclVolRepulsion.cpp
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderExclVolume.cpp
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderExclVolume.h
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderVolumeFF.cpp
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderVolumeFF.h
        ${SRC_DIR}/Mechanics/ForceField/Volume/CylinderVolumeInteractions.h
        ${SRC_DIR}/Mechanics/ForceField/ForceField.h
        ${SRC_DIR}/Mechanics/ForceField/ForceFieldManager.cpp
        ${SRC_DIR}/Mechanics/ForceField/ForceFieldManager.h
        ${SRC_DIR}/Mechanics/Minimizer/CGFletcherRievesMethod.cpp
        ${SRC_DIR}/Mechanics/Minimizer/CGFletcherRievesMethod.h
        ${SRC_DIR}/Mechanics/Minimizer/CGMethod.cpp
        ${SRC_DIR}/Mechanics/Minimizer/CGMethod.h
        ${SRC_DIR}/Mechanics/Minimizer/CGPolakRibiereMethod.cpp
        ${SRC_DIR}/Mechanics/Minimizer/CGPolakRibiereMethod.h
        ${SRC_DIR}/Mechanics/Minimizer/CGSteepestDescent.cpp
        ${SRC_DIR}/Mechanics/Minimizer/CGSteepestDescent.h
        ${SRC_DIR}/Mechanics/Minimizer/ConjugateGradient.h
        ${SRC_DIR}/Mechanics/Minimizer/Minimizer.h
        ${SRC_DIR}/Structure/Special/MTOC.cpp
        ${SRC_DIR}/Structure/Special/MTOC.h
        ${SRC_DIR}/Structure/Bead.cpp
        ${SRC_DIR}/Structure/Bead.h
        ${SRC_DIR}/Structure/Boundary.h
        ${SRC_DIR}/Structure/BoundaryElement.cpp
        ${SRC_DIR}/Structure/BoundaryElement.h
        ${SRC_DIR}/Structure/BoundaryElementImpl.cpp
        ${SRC_DIR}/Structure/BoundaryElementImpl.h
        ${SRC_DIR}/Structure/BoundaryImpl.cpp
        ${SRC_DIR}/Structure/BoundaryImpl.h
        ${SRC_DIR}/Structure/BoundarySurface.h
        ${SRC_DIR}/Structure/BoundarySurfaceImpl.cpp
        ${SRC_DIR}/Structure/BoundarySurfaceImpl.h
        ${SRC_DIR}/Structure/BranchingPoint.cpp
        ${SRC_DIR}/Structure/BranchingPoint.h
        ${SRC_DIR}/Structure/Bubble.cpp
        ${SRC_DIR}/Structure/Bubble.h
        ${SRC_DIR}/Structure/CBound.h
        ${SRC_DIR}/Structure/CBranchingPoint.cpp
        ${SRC_DIR}/Structure/CBranchingPoint.h
        ${SRC_DIR}/Structure/CCylinder.cpp
        ${SRC_DIR}/Structure/CCylinder.h
        ${SRC_DIR}/Structure/CLinker.cpp
        ${SRC_DIR}/Structure/CLinker.h
        ${SRC_DIR}/Structure/CMonomer.cpp
        ${SRC_DIR}/Structure/CMonomer.h
        ${SRC_DIR}/Structure/CMotorGhost.cpp
        ${SRC_DIR}/Structure/CMotorGhost.h
        ${SRC_DIR}/Structure/Compartment.cpp
        ${SRC_DIR}/Structure/Compartment.h
        ${SRC_DIR}/Structure/CompartmentGrid.cpp
        ${SRC_DIR}/Structure/CompartmentGrid.h
        ${SRC_DIR}/Structure/Cylinder.cpp
        ${SRC_DIR}/Structure/Cylinder.h
        ${SRC_DIR}/Structure/Database.h
        ${SRC_DIR}/Structure/DynamicNeighbor.h
        ${SRC_DIR}/Structure/Filament.cpp
        ${SRC_DIR}/Structure/Filament.h
        ${SRC_DIR}/Structure/Linker.cpp
        ${SRC_DIR}/Structure/Linker.h
        ${SRC_DIR}/Structure/MBranchingPoint.cpp
        ${SRC_DIR}/Structure/MBranchingPoint.h
        ${SRC_DIR}/Structure/MCylinder.cpp
        ${SRC_DIR}/Structure/MCylinder.h
        ${SRC_DIR}/Structure/MLinker.cpp
        ${SRC_DIR}/Structure/MLinker.h
        ${SRC_DIR}/Structure/MMotorGhost.cpp
        ${SRC_DIR}/Structure/MMotorGhost.h
        ${SRC_DIR}/Structure/MotorGhost.cpp
        ${SRC_DIR}/Structure/MotorGhost.h
        ${SRC_DIR}/Structure/Movable.h
        ${SRC_DIR}/Structure/Neighbor.h
        ${SRC_DIR}/Structure/NeighborList.h
        ${SRC_DIR}/Structure/NeighborListImpl.cpp
        ${SRC_DIR}/Structure/NeighborListImpl.h
        ${SRC_DIR}/Structure/Reactable.h
        ${SRC_DIR}/Structure/SubSystem.cpp
        ${SRC_DIR}/Structure/SubSystem.h
        ${SRC_DIR}/Structure/Trackable.cpp
        ${SRC_DIR}/Structure/Trackable.h
        #        ${SRC_DIR}/TESTS/test_boundaries.cpp
        #        ${SRC_DIR}/TESTS/test_compartments.cpp
        #        ${SRC_DIR}/TESTS/test_composite.cpp
        #        ${SRC_DIR}/TESTS/test_geometry.cpp
        #        ${SRC_DIR}/TESTS/test_gillespie.cpp
        #        ${SRC_DIR}/TESTS/test_gillespie_cross.cpp
        #        ${SRC_DIR}/TESTS/test_nrm.cpp
        #        ${SRC_DIR}/TESTS/test_reactions.cpp
        #        ${SRC_DIR}/TESTS/test_simple_gillespie.cpp
        #        ${SRC_DIR}/TESTS/test_species.cpp
        ${SRC_DIR}/BubbleInitializer.cpp
        ${SRC_DIR}/BubbleInitializer.h
        ${SRC_DIR}/CController.cpp
        ${SRC_DIR}/CController.h
        ${SRC_DIR}/common.cpp
        ${SRC_DIR}/common.h
        ${SRC_DIR}/Component.cpp
        ${SRC_DIR}/Component.h
        ${SRC_DIR}/Composite.cpp
        ${SRC_DIR}/Composite.h
        ${SRC_DIR}/Controller.cpp
        ${SRC_DIR}/Controller.h
        ${SRC_DIR}/DRController.cpp
        ${SRC_DIR}/DRController.h
        ${SRC_DIR}/FilamentInitializer.cpp
        ${SRC_DIR}/FilamentInitializer.h
        ${SRC_DIR}/GController.cpp
        ${SRC_DIR}/GController.h
        ${SRC_DIR}/Histogram.h
        ${SRC_DIR}/MathFunctions.cpp
        ${SRC_DIR}/MathFunctions.h
        ${SRC_DIR}/MController.cpp
        ${SRC_DIR}/MController.h
        ${SRC_DIR}/MEDYAN.cpp
        ${SRC_DIR}/Output.cpp
        ${SRC_DIR}/Output.h
        ${SRC_DIR}/Parser.cpp
        ${SRC_DIR}/Parser.h
        ${SRC_DIR}/Rand.cpp
        ${SRC_DIR}/Rand.h
        ${SRC_DIR}/RateChanger.h
        ${SRC_DIR}/RateChangerImpl.cpp
        ${SRC_DIR}/RateChangerImpl.h
        ${SRC_DIR}/Restart.h
        ${SRC_DIR}/SysParams.cpp
        ${SRC_DIR}/SysParams.h
        ${SRC_DIR}/utility.cpp
        ${SRC_DIR}/utility.h
        ${SRC_DIR}/Visitor.h)

include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Chemistry )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Boundary )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Branching )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Bubble )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Filament )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Linker )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/MotorGhost )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/ForceField/Volume )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Mechanics/Minimizer )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Structure )
include_directories( /home/aravind/Research/CompartmentFrontandBack+turnover+resrart_MEDYAN3_2_CONSTRANDOM/medyan/src/Structure/Special )

#CUDA_ADD_EXECUTABLE(Medyan ${SOURCE_FILES})
add_executable(Medyan3_2 ${SOURCE_FILES})
# Note: a target should be already defined using 'add_executable' or 'add_library'
target_link_libraries(Medyan3_2 ${Boost_LIBRARIES})