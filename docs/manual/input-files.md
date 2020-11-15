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

### Starting filament configuration

These parameters define the initial configuration and length of filaments in the system. It is noted that at least one filament, plus end, and minus end chemical species must be initialized in the chemistry input file, or a startup error will result.

| item | type | description |
|------|------|-------------|
| `FILAMENTFILE` | string | Name of filament initialization file. This is not required. |
| `NUMFILAMENTS` | int | Number of random filaments to initialize. These filaments will be randomly distributed in the system volume. |
| `FILAMENTLENGTH` | int | Number of cylinders per filament to initialize, defining the initial length of the filaments. |
| `FILAMENTTYTPE` | int | Filament type to initialize. |
| `PROJECTIONTYPE` | string | `{STRAIGHT, ZIGZAG, ARC, PREDEFINED}` Specifies how the beads are sampled between two ends of a filament. |

Projection type meaning for a filament specified by a set of 3D Cartesian coordinates `[v0, v1, v2 ...]`, and the number of beads.
- `STRAIGHT` - Creates a filament with minus end at `v0`, and extends number-of-bead full-size cylinders in `v0 -> v1` direction. `v2, v3 ...` are ignored.
- `ZIGZAG` - Creates a filament with minus end at `v0`, and extends number-of-bead full-size cylinders in (1) `v0 -> v1` direction (2) another different direction alternatively. `v2, v3 ...` are ignored.
- `ARC` - Not sure. (TODO: NEED DOCS)
- `PREDEFINED` - Creates a filament with bead coordinates `[v0, v1, v2 ...]`. WARNING: Each cylinder in the filament is treated as a full-sized cylinder, no matter what the initial bead coordinates are. That means each cylinder in the filament can be compressed or stretched initially.
