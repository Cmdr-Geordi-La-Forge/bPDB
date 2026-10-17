# Binary PDB (bPDB) Format Specification

A SIMD-optimized, heavily compressed binary format for macromolecular coordinate data.

## 1. File Architecture Overview

The bPDB format utilizes a strictly 8-byte aligned, linear TLV (Type-Length-Value) architecture. It eliminates random-access pointer chasing in favor of sequential stream processing, optimizing hardware prefetching. 

1. **Fixed Global Header (8 Bytes)**
2. **Extensible Metadata Block (TLV)**
3. **Coordinate Stream (Chunked)**

---

## 2. Global Header & Metadata

### 2.1 Fixed Global Header (8 Bytes)
The file strictly begins with this minimal 8-byte block.

| Offset | Type      | Description                     |
|:-------|:----------|:--------------------------------|
| `+0`   | `char[4]` | Magic String: `"bPDB"`          |
| `+4`   | `uint8`   | Format Version (e.g., `0x01`)   |
| `+5`   | `uint8`   | Deposition Day (1-31)           |
| `+6`   | `uint8`   | Deposition Month (1-12)         |
| `+7`   | `uint8`   | Deposition Year (`YYYY - 1973`) |

### 2.2 Extensible Metadata Block (TLV)
Allows forward-compatible metadata addition. Older parsers can skip unknown 8-byte codes by advancing the stream by `section_len` + padding.

| Type       | Description                                                          |
|:-----------|:---------------------------------------------------------------------|
| `uint16`   | Number of Metadata Sections (`num_meta`)                             |
| `uint16[]` | Array of `num_meta` 2-byte Codes (e.g., `0x0001` ID, `0x0002` TITLE) |
| padding    | Padded with `0x00` to the nearest 8-byte boundary                    |

For each of the `num_meta` sections declared above, the following payload structure strictly applies in order:

| Type      | Description                                            |
|:----------|:-------------------------------------------------------|
| `uint64`  | Unpadded length of the section payload (`section_len`) |
| `uint8[]` | The raw payload data                                   |
| padding   | Padded with `0x00` to the nearest 8-byte boundary      |

---

## 3. Coordinate Stream (Chunked)

The coordinate stream runs from the end of the Metadata Block until EOF. It operates purely sequentially. The stream is divided into "Chunks" of residues that share the same File and Chain origin. 

### 3.1 Chunk Header (16 Bytes)
Each Chunk begins with a 16-byte header that explicitly identifies its origin.

| Offset | Type       | Description                                               |
|:-------|:-----------|:----------------------------------------------------------|
| `+0`   | `uint16`   | Number of Residues contained in this chunk                |
| `+2`   | `uint8`    | Event Flags (`0x04` = Has TER, `0x08` = Has TER + OXT)    |
| `+3`   | `uint8`    | Fragment ID (e.g., `1` for AlphaFold `-F1`. `0` for none) |
| `+4`   | `char[2]`  | Chain ID (e.g., `" A"`, `"AA"`)                           |
| `+6`   | `char[10]` | File ID (e.g., `"1V0Q"`, `"A0A017SE85"`)                  |

### 3.2 Residue & Atom Blocks
Every residue in a chunk strictly consists of 1 Residue Block followed by a static number of Atom Blocks (dictated by the Superset Dictionary).

#### Residue Block (8 Bytes)
Establishes the base absolute coordinates and sequence metadata.

| Byte Offset | Type    | Description                                                          |
|:------------|:--------|:---------------------------------------------------------------------|
| `+0`        | `uint8` | Residue Type ID (0-36 from dictionary).                              |
| `+1`        | `uint8` | **B-factor MSB**: Top 8 bits of absolute 21-bit B-factor for Atom 1. |
| `+2`        | `uint8` | **resSeq LSB**: Bottom 8 bits of 16-bit Residue Sequence Number.     |
| `+3`        | `uint8` | **X-base MSB**: Top 8 bits of absolute 21-bit X coordinate.          |
| `+4`        | `uint8` | **resSeq MSB**: Top 8 bits of 16-bit Residue Sequence Number.        |
| `+5`        | `uint8` | **Y-base MSB**: Top 8 bits of absolute 21-bit Y coordinate.          |
| `+6`        | `char`  | **iCode**: Insertion code character (space or null if none).         |
| `+7`        | `uint8` | **Z-base MSB**: Top 8 bits of absolute 21-bit Z coordinate.          |

#### Atom Block (8 Bytes)
Interpreted as four 16-bit lanes mapped directly to SIMD pipelines.

| 16-bit Lane | Component | Bits 3-15 (13-bit MSB Math)                      | Bits 0-2 (3-bit LSB Metadata)             |
|:------------|:----------|:-------------------------------------------------|:------------------------------------------|
| `0`         | B-factor  | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bits 0-2                        |
| `1`         | X Coord   | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bits 3-5                        |
| `2`         | Y Coord   | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bit 6 \| AltLoc \| Charge Bit 0 |
| `3`         | Z Coord   | Delta from Atom N-1 (or absolute LSB for Atom 1) | Charge Bits 1-3                           |
