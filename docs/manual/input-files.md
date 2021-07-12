# Input files

All physical quantities in the input files will use the following base units, unless otherwise specified.

| name      | unit |
|-----------|------|
| time      | s    |
| length    | nm   |
| mass      | g    |

Some common derived units are displayed below.

| name      | unit |
|-----------|------|
| force     | pN   |
| energy    | pN⋅nm |

## System input files

TODO others

### Chemistry

The following chemical parameters can be set. It should be noted that the number of parameters listed for each chemical species type that resides on a filament must match the number of filament types, specified in the system input file. This must be consistent for all filament types. To set multiple parameters corresponding to multiple filaments, list the parameters with space in between after the parameter qualifier.

All chemical parameters must be set unless otherwise noted in the description. For the motor parameters, the number of parameters must match the number of motor species in the system. For more information on chemical algorithms, see [Popov et al (2016)](https://doi.org/10.1371/journal.pcbi.1004877).

An alternate set of parameters can be specified in replacement of `RUNTIME` for smaller systems in which simulation time is based on explicit reaction steps; if `RUNTIME` is not initialized or set to zero, the parameter `RUNSTEPS` and its associated chemical step-based parameter set will be used if provided.

| item | type | description |
|------|------|-------------|
| `CHEMISTRYFILE` | string | Input chemistry file. Should be in the input directory. |
| `CALGORITHM` | `{GILLESPIE, NRM}` | Chemistry algorithm used. |
| `RUNSTEPS` | int | Number of total chemical steps. If `RUNTIME` is set, will not be used. |
| `RUNTIME` | double | Total runtime of simulation. |
| `SNAPSHOTSTEPS` | int | Number of steps per snapshot. If `SNAPSHOTTIME` is set, will not be used. |
| `SNAPSHOTTIME` | double | Time of each snapshot. |
| `MINIMIZATIONSTEPS` | int | Number of chemical steps per mechanical equilibration. If `MINIMIZATIONTIME` is set, will not be used. |
| `MINIMIZATIONTIME` | double | Time between each mechanical equilibration. |
| `NEIGHBORLISTSTEPS` | int | Number of chemical steps per neighbor list update. This includes updating chemical reactions as well as force fields which rely on neighbor lists. If `NEIGHBORLISTTIME` is set, will not be used. |
| `NEIGHBORLISTTIME` | int | Time between each neighbor list update. |
| `NUMDIFFUSINGSPECIES` | int | Diffusing species in system. |
| `NUMBULKSPECIES` | int | Bulk species in system. |
| `NUMFILAMENTTYPES` | int | Number of different filament types. |
| `NUMFILAMENTSPECIES` | int | Filament species in system for each filament type defined. |
| `NUMPLUSENDSPECIES` | int | Plus end species in system for each filament type defined. |
| `NUMMINUSENDSPECIES` | int | Minus end species in system for each filament type defined. |
| `NUMBOUNDSPECIES` | int | Bound species in system for each filament type defined. |
| `NUMLINKERSPECIES` | int | Cross-linker species in system for each filament type defined.  |
| `NUMMOTORSPECIES` | int | Motor species in system for each filament type defined. |
| `NUMBRANCHERSPECIES` | int | Brancher species in system for each filament type defined. |
| `NUMBINDINGSITES` | int | Number of binding sites per cylinder for each filament type defined. This will set binding sites for cross-linkers, motors, and other binding molecules. |
| `NUMMOTORHEADSMIN` | int | Minimum number of motor heads per motor species defined. |
| `NUMMOTORHEADSMAX` | int | Maximum number of motor heads per motor species defined. |
| `MOTORSTEPSIZE` | double | Single motor head step size. |
|  `DISSIPATIONTRACKING` | `{OFF, ON}` | Whether to switch on the dissipation tracking feature. |
|  `LINKERBINDINGSKIP` | int | Switches on the different binding tracking feature to allow motors to have more binding spots per cylinder than linkers. The specified integer is the number of binding sites that the cross-linkers will skip before accepting a possible binding site. |
|  `EVENTTRACKING` | `{OFF, ON}` | Whether to switch on the event tracking feature. |
       
### Dynamic rates

The following dynamic rate forms and parameters can be set. These parameters are characteristic lengths and amplitudes of the rate changing equations outlined in [Popov et al (2016)](https://doi.org/10.1371/journal.pcbi.1004877). These can be tuned to mimic the stall and unbinding mechanochemical coupling of cross-linkers and myosin II motors. Note that if dynamic rates are enabled, the number of dynamic rate forms for each type of reaction must match the number of species of that type specified in the system input file, i.e. the number of forms for cross-linker unbinding must match the number of cross-linker species, etc.

The number of parameters specified for each type of dynamic rate form must match the number of parameters required for those forms. See below for details, and see [Popov et al (2016)](https://doi.org/10.1371/journal.pcbi.1004877) for more information on the explicit forms. Parameters must be listed in order of the form that they correspond to, also corresponding to the species that they represent.

| item | type | description |
|------|------|-------------|
| `DFPOLYMERIZATIONTYPE` | `{BROWRATCHET}` | Filament polymerization dynamic rate form. |
| `DFPOLYMERIZATIONLEN` | double | Characteristic length for filament polymerization dynamic rate form. |
| `DLUNBINDINGTYPE` | `{CATCHSLIP, SLIP}` | Cross-linker unbinding dynamic rate form. If `CATCHSLIP`, two parameters for `DLUNBINDINGLEN` and `DLUNBINDINGAMP` are needed to define the functional form. If `SLIP`, one `DLUNBIDINGLEN` is needed to define the functional form. |
| `DLUNBINDINGLEN` | double | Characteristic length of cross-linker unbinding dynamic rate form. |
| `DLUNBINDINGAMP` | double | Amplitude of cross-linker unbinding dynamic rate form. |
| `DMUNBINDINGTYPE` | `{LOWDUTYCATCHSLIP, LOWDUTYSLIP}` | Myosin II unbinding dynamic rate form. If `LOWDUTYCATCHSLIP`, two parameters for `DMUNBINDINGFORCE` are needed to define the functional form. If `LOWDUTYSLIP`, one `DMUNBIDINGFORCE` is needed to define the functional form. |
| `DMUNBINDINGFORCE` | double | Characteristic force of myosin II unbinding dynamic rate form. |
| `DMWALKINGTYPE` | `{LOWDUTYSTALL}` | Myosin II walking dynamic rate form. | 

### Starting filament configuration

These parameters define the initial configuration and length of filaments in the system. It is noted that at least one filament, plus end, and minus end chemical species must be initialized in the chemistry input file, or a startup error will result.

| item | type | description |
|------|------|-------------|
| `FILAMENTFILE` | string | Name of filament initialization file. This is not required. |
| `NUMFILAMENTS` | int | Number of random filaments to initialize. These filaments will be randomly distributed in the system volume. |
| `FILAMENTLENGTH` | int | Number of cylinders per filament to initialize, defining the initial length of the filaments. |
| `FILAMENTTYTPE` | int | Filament type to initialize. |
| `PROJECTIONTYPE` | `{STRAIGHT, ZIGZAG, ARC, PREDEFINED}` | Specifies how the beads are sampled between two ends of a filament. |

Projection type meaning for a filament specified by a set of 3D Cartesian coordinates `[v0, v1, v2 ...]`, and the number of beads.
- `STRAIGHT` - Creates a filament with minus end at `v0`, and extends number-of-bead full-size cylinders in `v0 -> v1` direction. `v2, v3 ...` are ignored.
- `ZIGZAG` - Creates a filament with minus end at `v0`, and extends number-of-bead full-size cylinders in (1) `v0 -> v1` direction (2) another different direction alternatively. `v2, v3 ...` are ignored.
- `ARC` - Not sure. (TODO: NEED DOCS)
- `PREDEFINED` - Creates a filament with bead coordinates `[v0, v1, v2 ...]`. WARNING: Each cylinder in the filament is treated as a full-sized cylinder, no matter what the initial bead coordinates are. That means each cylinder in the filament can be compressed or stretched initially.

### TODO others
