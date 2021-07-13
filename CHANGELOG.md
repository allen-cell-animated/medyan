# 5.0.0 (Released TBD)

## New features
- Added membrane simulation.


# 4.3.0 (Released 2021-07-12)

## New features
- Dropped support for pre C++17 compilers.
- Added simulation configuration generator, and enhanced input file parser to allow structured input using S-Expressions (4d71a27).
- Added `ReactionDy` type that supports reactions with dynamic number of reactants and products (e37feb3).

## Refactoring
- Refactored the data structure in the energy minimization, such that it allows for non-coordinate degrees of freedom (9a02e0f6d).
- Improved reference semantics for `Vec`-like structures as function arguments (27b296ad4).


# 4.2.0 (Released 2021-07-11)

## Enhancements
- Force components that reach >1e15 will trigger an alarm in addition to NaN and Inf values (39ce574).
- Removed thread pool usage (c78fb5e).
- Optimized filament and branching force field calculations (df76443, 4d523e7).
- Added source groups corresponding to the file tree structure in `CMakeLists.txt` (8fc9ead).
- `OPTIMOUT` macro offers more detailed timings of various parts of the code so the user has better understanding of the rate limiting steps in the simulation (873c814).
- Optimized species searching when cloning reactions (52f4c47).
- Optimized addition of dependent reactions (7b86dcd).
- Hessian analysis
    - Eigenvalue computation can be optionally disabled with `EIGENTRACKING: OFF` while computing Hessian matrix (0c455b5).
    - Add output to visualize all eigenmodes in Hessian tracking (16155a6).

## Bug fixes

- Restart
    - During restart phase, motor walking rate was not turned off under a couple of conditions. Fixed that (39ce574).
    - During restart phase, reactions for motor/linker is now turned off (39ce574).
    - During restart phase, `Cylinder::_filID` was not set. Fixed that (106c1e2).
    - Fixed restarted trajectories not having diffusion of chemical species (2903f54).

- Force fields
    - Fixed branching interactions when the offspring filament is attached to the plus end of parent filament (c87d3d9, 840de95).
    - Fixed initial force computation for combined-filament-stretching-and-bending force field, which is no longer enabled by default (99c9bd6).
    - Fixed incorrect forces for branching position cosine (4d523e7).
    - Fixed wrong constant used in branching bending energy computation (a68e859).
    - Fixed incorrect volume exclusion force expressions (49dd9d1).
    - Volume exclusion potential now considers cylinder equilibrium lengths (5d12d11).

- Others
    - After chemistry, when the diffusing species copy numbers are set, now use `up()` and `down()` functions instead of only `setN()` to trigger crosschecks and signals (4b8a5e5).
    - Replaced `std::bind2nd` with lambda function (805b221).
    - Fixed reaction rate scaling by volume in newly created compartments (f6c7783).
    - Fixed the coordinate calculation for binding sites using the `Cylinder::adjustedrelativeposition` function. However, to ensure backward compatibility, it is turned off by default (dbf7527, dbecef5).
    - The cloned reaction is fixed to be non-proto-compartment (e0e6fac).
    - Fixed incorrect filament nucleation when bubbles exist (78a8f8c).
    - Fixed incorrect creation of `FilamentCreationCallback` (7b77123).
    - Fixed build and test workflows on GitHub Actions (bf63f6f).


# 4.1.2 (Released 2020-06-01)

## Enhancements
- New building procedure using `CMake` is supported and recommended (a1b28e254).
- Added a force field that combines filament stretching and bending (772133585).
- Added quadratic line search in gradient descent method (772133585).

## Bug fixes
- The copy number of chemical species are now assigned before dynamic rate initialization (7175987be).
- Fixed branching dihedral force field cosine form (772133585).
- Myosin duty ratio is now calculated at runtime (772133585).

# 4.1.1 (Released 2020-01-30)

## New features
- The system input file no longer requires species number of types specification (b887ba2).
- Rearranged examples and added new ones (48993f8).

## Refactoring and optimizations
- Refactored the logger (a2c69f9).
- Integrated the tests (8de0a0f).
- Increased readability of on-screen outputs.
- `Makefile` default optimization changed from `-Os` to `-O2` (b8468c3).

## Bug fixes
- Fixed distribution generation in chemistry algorithm (5897eed).
- Fixed bugs related to SIMD macros (13fe728).
- Fixed volume exclusion not including short filaments (d0961d4).
- Fixed Myosin motor duty ratio calculation. Now it is based on binding/unbinding rates. (6f57531)

# 4.1.0 (Released 2019-10-29)

## New features
- Added a thread pool implementation, which may facilitate multi-threading computations (fd69e22).
- Added new MTOC functions (e286a3d).
- Added AFM pulling simulation (b581d1a).
- Added mechanical Hessian analysis (b581d1a).
- Added branching dihedral quadratic force field (37b6173).
- Added a cell list data structure for elements contained in the compartments (586f8f5).
- Used a double linked list data structure instead of `std::unordered_set` for `Reactable`s and `Movable`s (0f32b73).
- The dissipation tracking can now have higher resolution, reflecting the mechanical energy change for every force field (5e7eed8).

## Refactoring and optimizations
- Refactored the `Database` class (bc0222c, 3c71316).
- Removed stretched energy computations (ee6220f).
- Various improvements on the restarting procedure (0f32b73).

## Bug fixes
- Fixed motor stretch force not reset before minimization.
- Distance between the boundary and the filaments can be negative now.
- Cylinder neighbor list for volume exclusion will generate multiple neighbor list if there is more than one filament type.
- Filament creation by nucleation and branching is not allowed in partially activated compartment with volume fraction < 0.5. It potentially prevents the infinite while loop during filament creation.
- Fixed `BoundaryCylinderRepulsionIn` force field.
- Fixed `MTOCAttachment` potential.
- Fixed mechanochemical feedback on motor initialization under `PLOSFEEDBACK` version (1fcd9cf).
- Fixed cases where 2 ends a motor binds on the same binding site.
- Fixed compartment transfer/share axis (befaf74).
- Fixed the way of counting number of interactions for filament bending (75bc733).
- Fixed the newly created filament not having mechanochemical feedback (9dc5d27).
- Fixed the problem that the type `floatingpoint` cannot be aliased to `double` (87fad88).
- Fixed `dist_avx.h` not compatible with AVX mode (ce314fa).
- Fixed the issue with changing number of binding sites per cylinder in SIMD binding search (9a39874).
- Fixed the boundary pinning force field (49e254b).
- Fixed minor portability issues on different compilers (eb10ca1).
- Fixed ending simulation message not displayed when simulation finishes (4a2f610).
- `PLOSFEEDBACK` is now turned off by default (3ff0837).
- Fixed the documents for MEDYAN 4.0 speed up comparison (3ff0837).

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
