# 4.1 (Released 2019)

## Breaking changes

## New features
- Added a thread pool implementation, which may facilitate multi-threading computations.

## Refactoring and optimizations

## Bug fixes

# 4.0 (Released 2019-07-05)

## New features
- Dropped support for pre C++14 compilers.
- Added support for MSVC compilers and added Visual Studio solution file.
- Added cylindrical boundary type.
- Added support for runtime modification of filament polymerization/depolymerization rate.
- Added support for energy dissipation tracking.

## Optimizations
- Used flattened storage for force field computations.
- Improved performance for neighbor list and binding site search. Binding site search now gains SIMD support. (See [file](./docs/Medyan4.0.pdf) for detailed performance benchmark)

## Bug fixes
- Fixed a bug in boundary repulsion potential.
- Fixed numerical errors that might occur in cylinder volume exclusion potential.

# 3.2.1 (Released 2018-08-23)

## Bug fixes
- Rewrote `Makefile` and fixed the incorrect dependency generation. Also added `PLOSFEEDBACK` to the macro list.
- Fixed a typo in the source.
- Fixed some bugs related to mechanochemical feedback when `PLOSFEEDBACK` is defined.

# 3.2 (Released 2018-08-14)

## Bug fixes
- Filament
    - Fixed a bug on filament severing by cofilin.

        Colfilin reactions required certain callback modifications to effectively simulate them without crashing the program. Those changes are included in this version.

    - Load force computation now considers directionality of filaments.

        Previous models of Brownian ratchet feedback in Medyan considered just the distance from boundary. In this version, both the distance and angle are considered. This leads to a more accurate model for brownian ratchet.

    - Fixed a bug on load force clearing.

        In the previous implementations of MEDYAN, brownian ratchet forces were cleared improperly before being overwritten. As a result, the filament polymerization/depolymerization rates were inaccurate. This is especially true in high χ parameter cases.

- Motors and Linkers
    - Fixed a bug on stretch force computation for motors and linkers.
- Restart protocol can now effectively restart networks with shorter filaments compared to the earlier version.
- Other
    - Improved arrangement of main controller procedure.
    - Bug fix on retrieving compartment grid center coordinate. (used rarely in MEDYAN)
    - Reviewed example files to ensure accuracy of parameters mentioned.
    - `MOTORWALKINGFORWARD` string was used instead of `MOTORWALKINGBACKWARD`. This bug only affects simulations with bi-directional motors.
    - Fixed diffusion bug in `moveBoundary` function.

        Previous version of MEDYAN had a bug which resulted in partial activation of diffusion reactions when a new comaprtment is activated. This has been addressed in this version.

## New features
- Simulations can now control chemical reaction volume based on network span if necessary. Please refer to Usage guide for more details.
- Added output files for trajectory of concentration, plus end type / coordinate and chemical reaction counter on filament.
- Added feature to visualize plus ends in trajectories using VMD. Please refer to UsageGuide for more information.
- Added plotting function to the analysis script in `AnalyzeTrajectory.py`.
- Added example input files for different systems.
- Added a nucleation example input files.
- Added treadmilling tracking.

    Simulations can now track treadmilling in each filament by tracking the total number of polymerization, and depolymerization reactions explicitly in plus and minus ends in addition to nucleation reactions during each chemical evolution cycle.
 
- Added macro `DEBUGCONSTANTSEED`.

    It helps simulate trajectories without stochasticity thereby helping in faster debugging/cross comparison of codes. Simulations might be slower in this mode and so it is better to use smaller systems while debugging.

- Added macro `PLOSFEEDBACK`.

    As many mechanochemical models have been used in publications from Papoian lab, this macro ensures that the feedback model is the same as mentioned in PLOS computatinal biology paper.

    Please refer to Model Guide on Mechanochemical feedback for more information.

    Note: macro changes require edition of `Makefile` and so code should be recompiled to effect the updates.
