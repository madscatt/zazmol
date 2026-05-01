# PDB Fixed-Column Reference

This is a reference note for the standalone C++ PDB parser. The official wwPDB
format is useful for column ranges, but Python `zazmol` remains the behavior
oracle for tolerated legacy input, field names, and defaults.

Primary references:

- wwPDB PDB Format v3.3 Coordinate Section:
  <https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html>
- wwPDB PDB Format v3.3 Introduction:
  <https://www.wwpdb.org/documentation/file-format-content/format33/sect1.html>
- wwPDB PDB Format v2.3 Coordinate Section:
  <https://www.wwpdb.org/documentation/file-format-content/format23/sect9.html>
- wwPDB PDB Format v2.3 CONECT Section:
  <https://www.wwpdb.org/documentation/file-format-content/format23/sect10.html>

## ATOM/HETATM Columns

The current C++ parser should slice the same fixed columns used by Python
`zazmol`:

| Python/SasMol field | PDB columns | Python slice | Notes |
| --- | ---: | --- | --- |
| `record` / `atom` | 1-6 | `line[0:6]` | `ATOM` or `HETATM` after trim. |
| `original_index` | 7-11 | `line[6:11]` | Original atom serial. |
| `name` | 13-16 | `line[12:16]` | Trimmed atom name. |
| `loc` | 17 | `line[16]` | Preserved only in `pdbscan` mode; otherwise blank. |
| `resname` | 18-21 in Python | `line[17:21]` | Python intentionally reads four columns. |
| `chain` | 22 | `line[21]` | Single-character chain identifier. |
| `resid` | 23-26 | `line[22:26]` | Parsed as integer, original text also retained. |
| `rescode` | 27 | `line[26]` | Insertion code. |
| `x` | 31-38 | `line[30:38]` | `coord_type` storage. |
| `y` | 39-46 | `line[38:46]` | `coord_type` storage. |
| `z` | 47-54 | `line[46:54]` | `coord_type` storage. |
| `occupancy` | 55-60 | `line[54:60]` | Python default differs by path; preserve behavior. |
| `beta` | 61-66 | `line[60:66]` | Python default differs by path; preserve behavior. |
| `segname` | 73-76 | `line[72:76]` | Legacy/SASMOL segment identifier. |
| `element` | 77-78 | `line[76:78]` | Resolved later when blank. |
| `charge` | 79-80 | `line[78:80]` | Optional. |

## MODEL/ENDMDL/END

The wwPDB spec requires paired `MODEL`/`ENDMDL` records for multi-model
entries, and an `END` record terminates the file. Python `zazmol` is more
tolerant for useful legacy cases:

- single-frame files with no final `END` are accepted
- END-separated pseudo-trajectories are accepted when atom counts match
- MODEL/ENDMDL trajectories require consistent frame structure

## CONECT

The legacy CONECT record uses 5-column atom serial fields after the record name.
Python `zazmol` stores connectivity against original atom serials and remaps to
current indices on output. The C++ code should preserve that behavior rather
than inventing a new connectivity schema.
