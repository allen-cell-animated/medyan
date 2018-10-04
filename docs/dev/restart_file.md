## Introduction

During MEDYAN simulation, the file containing all the sufficient information to rebuild the whole system will be generated and overwritten. MEDYAN is able to initiate the system from a restart file so that the system is recovered to the point that the file is last updated.

The file is in binary format, and information will be stored in a certain structure which can then be read and parsed by MEDYAN. The details of the structure will be illustrated below.

The information in the input files need not to be included in the restart file as they will still be read prior to reading the restart file. However, some information in the input files will be overridden by the restart file (and thus useless). For example, the number of filaments might be different from what's specified in the system input.

## File structure

The information is stored hierarchically. Each layer in the hierarchy possesses a header and a body, and is possibly the body of higher level layers. The headers normally contain useful information (metadata) of the bodies, and the length information might also be included.

- Global header

### Filament

#### Header

```
     0                   1                   2                   3
     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    |                    Number of beads (uint32)                   |
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    |        Filament Type          |           Checksum            |
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
```

#### Content

The data contains information of bead coordinates, which should be in the format of triplets of 64-bit floats. So the total number of bytes used must be `24 * N`. Including header, the total number of bytes is `24 * N + 8`.
