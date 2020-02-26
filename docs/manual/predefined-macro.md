# Predefined macros

A list of all the predefined macros and their usages

| Macro | Description |
|-------|-------------|
| `BOOST_MEM_POOL` | Enable boost memory pool optimizations. |
| `BOOL_POOL_NSIZE` | Set boost memory pool size. |
| `CHECKFORCES_INF_NAN` | Enable checks for `inf` or `NaN` forces. |
| `CHEMISTRY` | Enable system chemistry. |
| `DYNAMICRATES` | Enable dynamic rate changing. This macro can only be defined if both `CHEMISTRY` and `MECHANICS` are defined. |
| `HYBRID_NLSTENCILLIST` | An optimized neighbor list implementation. Conflicts with `NLORIGINAL` and `SIMDBINDINGSEARCH`. |
| `MECHANICS` | Enable system mechanics. |
| `NLORIGINAL` | The neighbor list implementation similar to MEDYAN v3.2. Conflicts with `HYBRID_NLSTENCILLIST` and `SIMDBINDINGSEARCH`. |
| `REACTION_SIGNALING` | Enable reaction callback signaling. |
| `RSPECIES_SIGNALING` | Enable species callback signaling. |
| `SERIAL`    | Enable serial implementation of energy minimizer. |
| `SIMDBINDINGSEARCH` | Enable SIMD-based neighbor list and binding site pair search protocol. Conflicts with `NLORIGINAL` and `HYBRID_NLSTENCILLIST`. |
| `TRACK_DEPENDENTS` | Track reaction dependents in system. |
| `TRACK_ZERO_COPY_N` | Passivate reactions with zero copy number. |
| `TRACK_UPPER_COPY_N` | Passivate reactions with big copy number. |
