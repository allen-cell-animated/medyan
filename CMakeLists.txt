cmake_minimum_required(VERSION 3.12)

project(medyan)

#######################################
# Command line variables
#######################################

#######################################
# Source files
#######################################

set(src_dir          "src")
set(src_dir_external "external")

# List sources
file(GLOB_RECURSE src_list LIST_DIRECTORIES true CONFIGURE_DEPENDS "${src_dir}/*.cpp")

#######################################
# Compiling configs
#######################################

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -ffast-math -fno-finite-math-only -mtune=native -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -funroll-loops -flto")
endif()

# Add the executable
add_executable(medyan ${src_list})

# Include paths
target_include_directories(medyan PUBLIC
    ${src_dir}
    ${src_dir}/Chemistry
    ${src_dir}/Mechanics
    ${src_dir}/Mechanics/ForceField
    ${src_dir}/Mechanics/ForceField/Boundary
    ${src_dir}/Mechanics/ForceField/Branching
    ${src_dir}/Mechanics/ForceField/Bubble
    ${src_dir}/Mechanics/ForceField/Filament
    ${src_dir}/Mechanics/ForceField/Linker
    ${src_dir}/Mechanics/ForceField/MotorGhost
    ${src_dir}/Mechanics/ForceField/Volume
    ${src_dir}/Mechanics/Minimizer
    ${src_dir}/Structure
    ${src_dir}/Structure/Special
    ${src_dir_external}
)

# Preprocessor macros
target_compile_definitions(medyan PRIVATE
    # Basic system
    CHEMISTRY
    MECHANICS
    DYNAMICRATES

    # Chemistry algorithm
    TRACK_DEPENDENTS # Needed for Gillespie and NRM algorithm
    # Track zero copy and max copy number.
    # For passivating and activating reactions accordingly
    TRACK_ZERO_COPY_N
    TRACK_UPPER_COPY_N
    # Reaction signaling. Needed for all Filament, Linker,
    # Motor, and BranchingPoint reactions in system
    REACTION_SIGNALING
    RSPECIES_SIGNALING

    # Energy minimization algorithm
    MOVEBEADSLINESEARCH
    SERIAL

    # Neighbors list calculations (choose 1 of 3)
    #NLORIGINAL
    #HYBRID_NLSTENCILLIST
    SIMDBINDINGSEARCH
    NPROCS=1 # deprecated

    # Boost memory
    BOOST_MEM_POOL
    BOOL_POOL_NSIZE=65536

    # Debug
    #PLOSFEEDBACK
    CHECKFORCES_INF_NAN
)

# Additional link directories
if(${MEDYAN_ADDITIONAL_LINK_DIRS})
    target_link_libraries(medyan PRIVATE ${MEDYAN_ADDITIONAL_LINK_DIRS})
endif()
# RPath specification
if(${MEDYAN_RPATH})
    set_target_properties(medyan PROPERTIES
        BUILD_RPATH ${MEDYAN_RPATH}
    )
endif()


#######################################
# External dependencies
#######################################

# boost
if("${MEDYAN_BOOST_INSTALL_MODE}" STREQUAL "")
    # Default boost install using vcpkg

elseif("${MEDYAN_BOOST_INSTALL_MODE}" STREQUAL "find")
    # Automatically try to find boost
    find_package(Boost 1.67.0 REQUIRED)
    target_include_directories(medyan PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(medyan PRIVATE ${Boost_LIBRARIES})

elseif("${MEDYAN_BOOST_INSTALL_MODE}" STREQUAL "manual")
    # Manually specify boost location
    target_include_directories(medyan PUBLIC "${MEDYAN_BOOST_INCLUDE_DIR}")
    target_link_directories(medyan PRIVATE "${MEDYAN_BOOST_LIBRARY_DIR}")
    target_link_libraries(medyan PRIVATE "${MEDYAN_BOOST_LIBRARIES}")

endif()

# catch2
find_package(Catch2 CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE Catch2::Catch2)

# eigen3
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE Eigen3::Eigen)

