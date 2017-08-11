cmake_minimum_required(VERSION 3.8)
project(MEDYAN)

add_definitions(-DCAMKII)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Boost REQUIRED system thread timer chrono)
# see https://cmake.org/cmake/help/latest/module/FindBoost.html
find_package(Boost REQUIRED system thread timer chrono)
include_directories(${Boost_INCLUDE_DIRS})
# include_directories(/usr/local/include/boost)

include_directories(/usr/local/include/)

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

set(SOURCE_FILES
        Chemistry/BindingManager.cpp
        Chemistry/BindingManager.h
        Chemistry/ChemCallbacks.h
        Chemistry/ChemGillespieImpl.cpp
        Chemistry/ChemGillespieImpl.h
        Chemistry/ChemManager.cpp
        Chemistry/ChemManager.h
        Chemistry/ChemNRMImpl.cpp
        Chemistry/ChemNRMImpl.h
        Chemistry/ChemRNode.h
        Chemistry/ChemSim.cpp
        Chemistry/ChemSim.h
        Chemistry/ChemSimImpl.h
        Chemistry/ChemSimpleGillespieImpl.cpp
        Chemistry/ChemSimpleGillespieImpl.h
        Chemistry/Reaction.cpp
        Chemistry/Reaction.h
        Chemistry/ReactionBase.cpp
        Chemistry/ReactionBase.h
        Chemistry/ReactionContainer.h
        Chemistry/ReactionTemplate.cpp
        Chemistry/ReactionTemplate.h
        Chemistry/RSpecies.cpp
        Chemistry/RSpecies.h
        Chemistry/Species.cpp
        Chemistry/Species.h
        Chemistry/SpeciesContainer.h
        Mechanics/ForceField/Boundary/BoundaryBubbleRepulsion.cpp
        Mechanics/ForceField/Boundary/BoundaryBubbleRepulsion.h
        Mechanics/ForceField/Boundary/BoundaryBubbleRepulsionExp.cpp
        Mechanics/ForceField/Boundary/BoundaryBubbleRepulsionExp.h
        Mechanics/ForceField/Boundary/BoundaryCylinderAttachment.cpp
        Mechanics/ForceField/Boundary/BoundaryCylinderAttachment.h
        Mechanics/ForceField/Boundary/BoundaryCylinderAttachmentHarmonic.cpp
        Mechanics/ForceField/Boundary/BoundaryCylinderAttachmentHarmonic.h
        Mechanics/ForceField/Boundary/BoundaryCylinderRepulsion.cpp
        Mechanics/ForceField/Boundary/BoundaryCylinderRepulsion.h
        Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.cpp
        Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.h
        Mechanics/ForceField/Boundary/BoundaryFF.cpp
        Mechanics/ForceField/Boundary/BoundaryFF.h
        Mechanics/ForceField/Boundary/BoundaryInteractions.h
        Mechanics/ForceField/Branching/BranchingBending.cpp
        Mechanics/ForceField/Branching/BranchingBending.h
        Mechanics/ForceField/Branching/BranchingBendingCosine.cpp
        Mechanics/ForceField/Branching/BranchingBendingCosine.h
        Mechanics/ForceField/Branching/BranchingDihedral.cpp
        Mechanics/ForceField/Branching/BranchingDihedral.h
        Mechanics/ForceField/Branching/BranchingDihedralCosine.cpp
        Mechanics/ForceField/Branching/BranchingDihedralCosine.h
        Mechanics/ForceField/Branching/BranchingFF.cpp
        Mechanics/ForceField/Branching/BranchingFF.h
        Mechanics/ForceField/Branching/BranchingInteractions.h
        Mechanics/ForceField/Branching/BranchingPosition.cpp
        Mechanics/ForceField/Branching/BranchingPosition.h
        Mechanics/ForceField/Branching/BranchingPositionCosine.cpp
        Mechanics/ForceField/Branching/BranchingPositionCosine.h
        Mechanics/ForceField/Branching/BranchingStretching.cpp
        Mechanics/ForceField/Branching/BranchingStretching.h
        Mechanics/ForceField/Branching/BranchingStretchingHarmonic.cpp
        Mechanics/ForceField/Branching/BranchingStretchingHarmonic.h
        Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.cpp
        Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.h
        Mechanics/ForceField/Bubble/BubbleBubbleRepulsionExp.cpp
        Mechanics/ForceField/Bubble/BubbleBubbleRepulsionExp.h
        Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.cpp
        Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.h
        Mechanics/ForceField/Bubble/BubbleCylinderRepulsionExp.cpp
        Mechanics/ForceField/Bubble/BubbleCylinderRepulsionExp.h
        Mechanics/ForceField/Bubble/BubbleFF.cpp
        Mechanics/ForceField/Bubble/BubbleFF.h
        Mechanics/ForceField/Bubble/BubbleInteractions.h
        Mechanics/ForceField/Bubble/MTOCAttachment.cpp
        Mechanics/ForceField/Bubble/MTOCAttachment.h
        Mechanics/ForceField/Bubble/MTOCAttachmentHarmonic.cpp
        Mechanics/ForceField/Bubble/MTOCAttachmentHarmonic.h
        Mechanics/ForceField/Camkii/CamkiiBendingCosine.cpp
        Mechanics/ForceField/Camkii/CamkiiBendingCosine.h
        Mechanics/ForceField/Camkii/CamkiiBending.cpp
        Mechanics/ForceField/Camkii/CamkiiBending.h
        Mechanics/ForceField/Camkii/CamkiiFF.cpp
        Mechanics/ForceField/Camkii/CamkiiFF.h
        Mechanics/ForceField/Camkii/CamkiiInteractions.h
        Mechanics/ForceField/Camkii/CamkiiStretching.cpp
        Mechanics/ForceField/Camkii/CamkiiStretching.h
        Mechanics/ForceField/Camkii/CamkiiStretchingHarmonic.cpp
        Mechanics/ForceField/Camkii/CamkiiStretchingHarmonic.h
        Mechanics/ForceField/Filament/FilamentBending.cpp
        Mechanics/ForceField/Filament/FilamentBending.h
        Mechanics/ForceField/Filament/FilamentBendingCosine.cpp
        Mechanics/ForceField/Filament/FilamentBendingCosine.h
        Mechanics/ForceField/Filament/FilamentBendingHarmonic.cpp
        Mechanics/ForceField/Filament/FilamentBendingHarmonic.h
        Mechanics/ForceField/Filament/FilamentFF.cpp
        Mechanics/ForceField/Filament/FilamentFF.h
        Mechanics/ForceField/Filament/FilamentInteractions.h
        Mechanics/ForceField/Filament/FilamentStretching.cpp
        Mechanics/ForceField/Filament/FilamentStretching.h
        Mechanics/ForceField/Filament/FilamentStretchingHarmonic.cpp
        Mechanics/ForceField/Filament/FilamentStretchingHarmonic.h
        Mechanics/ForceField/Linker/LinkerFF.cpp
        Mechanics/ForceField/Linker/LinkerFF.h
        Mechanics/ForceField/Linker/LinkerInteractions.h
        Mechanics/ForceField/Linker/LinkerStretching.cpp
        Mechanics/ForceField/Linker/LinkerStretching.h
        Mechanics/ForceField/Linker/LinkerStretchingHarmonic.cpp
        Mechanics/ForceField/Linker/LinkerStretchingHarmonic.h
        Mechanics/ForceField/MotorGhost/MotorGhostFF.cpp
        Mechanics/ForceField/MotorGhost/MotorGhostFF.h
        Mechanics/ForceField/MotorGhost/MotorGhostInteractions.h
        Mechanics/ForceField/MotorGhost/MotorGhostStretching.cpp
        Mechanics/ForceField/MotorGhost/MotorGhostStretching.h
        Mechanics/ForceField/MotorGhost/MotorGhostStretchingHarmonic.cpp
        Mechanics/ForceField/MotorGhost/MotorGhostStretchingHarmonic.h
        Mechanics/ForceField/Volume/CylinderExclVolRepulsion.cpp
        Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h
        Mechanics/ForceField/Volume/CylinderExclVolume.cpp
        Mechanics/ForceField/Volume/CylinderExclVolume.h
        Mechanics/ForceField/Volume/CylinderVolumeFF.cpp
        Mechanics/ForceField/Volume/CylinderVolumeFF.h
        Mechanics/ForceField/Volume/CylinderVolumeInteractions.h
        Mechanics/ForceField/ForceField.h
        Mechanics/ForceField/ForceFieldManager.cpp
        Mechanics/ForceField/ForceFieldManager.h
        Mechanics/Minimizer/CGFletcherRievesMethod.cpp
        Mechanics/Minimizer/CGFletcherRievesMethod.h
        Mechanics/Minimizer/CGMethod.cpp
        Mechanics/Minimizer/CGMethod.h
        Mechanics/Minimizer/CGPolakRibiereMethod.cpp
        Mechanics/Minimizer/CGPolakRibiereMethod.h
        Mechanics/Minimizer/CGSteepestDescent.cpp
        Mechanics/Minimizer/CGSteepestDescent.h
        Mechanics/Minimizer/ConjugateGradient.h
        Mechanics/Minimizer/Minimizer.h
        Structure/Special/MTOC.cpp
        Structure/Special/MTOC.h
        Structure/Bead.cpp
        Structure/Bead.h
        Structure/Boundary.h
        Structure/BoundaryElement.cpp
        Structure/BoundaryElement.h
        Structure/BoundaryElementImpl.cpp
        Structure/BoundaryElementImpl.h
        Structure/BoundaryImpl.cpp
        Structure/BoundaryImpl.h
        Structure/BoundarySurface.h
        Structure/BoundarySurfaceImpl.cpp
        Structure/BoundarySurfaceImpl.h
        Structure/BranchingPoint.cpp
        Structure/BranchingPoint.h
        Structure/Bubble.cpp
        Structure/Bubble.h
        Structure/Camkii.cpp
        Structure/Camkii.h
        Structure/CBound.h
        Structure/CBranchingPoint.cpp
        Structure/CBranchingPoint.h
        Structure/CCylinder.cpp
        Structure/CCylinder.h
        Structure/CLinker.cpp
        Structure/CLinker.h
        Structure/CMonomer.cpp
        Structure/CMonomer.h
        Structure/CMotorGhost.cpp
        Structure/CMotorGhost.h
        Structure/Compartment.cpp
        Structure/Compartment.h
        Structure/CompartmentGrid.cpp
        Structure/CompartmentGrid.h
        Structure/Cylinder.cpp
        Structure/Cylinder.h
        Structure/Database.h
        Structure/DynamicNeighbor.h
        Structure/Filament.cpp
        Structure/Filament.h
        Structure/Linker.cpp
        Structure/Linker.h
        Structure/MBranchingPoint.cpp
        Structure/MBranchingPoint.h
        Structure/MCylinder.cpp
        Structure/MCylinder.h
        Structure/MLinker.cpp
        Structure/MLinker.h
        Structure/MMotorGhost.cpp
        Structure/MMotorGhost.h
        Structure/MotorGhost.cpp
        Structure/MotorGhost.h
        Structure/Movable.h
        Structure/Neighbor.h
        Structure/NeighborList.h
        Structure/NeighborListImpl.cpp
        Structure/NeighborListImpl.h
        Structure/Reactable.h
        Structure/SubSystem.cpp
        Structure/SubSystem.h
        Structure/Trackable.cpp
        Structure/Trackable.h
        BubbleInitializer.cpp
        BubbleInitializer.h
        CController.cpp
        CController.h
        common.cpp
        common.h
        Component.cpp
        Component.h
        Composite.cpp
        Composite.h
        Controller.cpp
        Controller.h
        TestController.cpp
        TestController.h
        DRController.cpp
        DRController.h
        FilamentInitializer.cpp
        FilamentInitializer.h
        GController.cpp
        GController.h
        Histogram.h
        MathFunctions.cpp
        MathFunctions.h
        MController.cpp
        MController.h
        MEDYAN.cpp
        Output.cpp
        Output.h
        Parser.cpp
        Parser.h
        Rand.cpp
        Rand.h
        RateChanger.h
        RateChangerImpl.cpp
        RateChangerImpl.h
        Restart.h
        SysParams.cpp
        SysParams.h
        utility.cpp
        utility.h
        Visitor.h)

add_executable(MEDYAN ${SOURCE_FILES})
include_directories( ./ )
include_directories( ./Chemistry )
include_directories( ./Mechanics )
include_directories( ./Mechanics/ForceField )
include_directories( ./Mechanics/ForceField/Boundary )
include_directories( ./Mechanics/ForceField/Branching )
include_directories( ./Mechanics/ForceField/Bubble )
include_directories( ./Mechanics/ForceField/Camkii )
include_directories( ./Mechanics/ForceField/Filament )
include_directories( ./Mechanics/ForceField/Linker )
include_directories( ./Mechanics/ForceField/MotorGhost )
include_directories( ./Mechanics/ForceField/Volume )
include_directories( ./Mechanics/Minimizer )
include_directories( ./Structure )
include_directories( ./Structure/Special )

# check if we need CUDA TBD
# FIND_PACKAGE(CUDA)

# CUDA_ADD_EXECUTABLE(MEDYAN ${SOURCE_FILES})
# add_executable(src ${SOURCE_FILES})
# Note: a target should be already defined using 'add_executable' or 'add_library'
target_link_libraries(MEDYAN ${Boost_LIBRARIES})