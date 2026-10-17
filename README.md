# Binary PDB (bPDB) Format Specification

A SIMD-optimized, heavily compressed binary format for macromolecular coordinate data.

## 1. File Architecture Overview

The bPDB format utilizes a strictly 8-byte aligned, sequentially ordered architecture. Mandatory blocks occur in a guaranteed order without type codes, allowing the decoder to stream linearly through the file while optimizing hardware prefetching. Future extensions are handled by a dedicated Metadata TLV block immediately preceding the coordinate stream.

1. **Fixed Global Header (32 Bytes)**
2. **Mandatory Title Block**
3. **Mandatory File ID Array**
4. **Mandatory Chain Name Array**
5. **Extensible Metadata Block (TLV)**
6. **Coordinate Stream**

---

## 2. Sequential Blocks

### 2.1 Fixed Global Header (32 Bytes)
The file strictly begins with this 32-byte block.

| Offset | Type       | Description                                   |
|:-------|:-----------|:----------------------------------------------|
| `+0`   | `char[4]`  | Magic String: `"bPDB"`                        |
| `+4`   | `uint8`    | Format Version (e.g., `0x01`)                 |
| `+5`   | `uint8`    | Deposition Day (1-31)                         |
| `+6`   | `uint8`    | Deposition Month (1-12)                       |
| `+7`   | `uint8`    | Deposition Year (`YYYY - 1973`)               |
| `+8`   | `char[16]` | `pdbID` (PDB ID or AlphaFold ID, null-padded) |
| `+24`  | `uint32`   | Number of Files (`num_files`)                 |
| `+28`  | `uint32`   | Number of Chains (`num_chains`)               |

### 2.2 Title Block
| Type     | Description                                       |
|:---------|:--------------------------------------------------|
| `uint32` | Unpadded length of the Title string (`title_len`) |
| `char[]` | Title string                                      |
| padding  | Padded with `0x00` to the nearest 8-byte boundary |

### 2.3 File ID Array
| Type         | Description                                                              |
|:-------------|:-------------------------------------------------------------------------|
| `uint32`     | Number of concatenated files (`f_count`)                                 |
| `char[][16]` | Array of `f_count` 16-byte strings identifying split/concatenated files. |
| padding      | Padded with `0x00` to the nearest 8-byte boundary                        |

### 2.4 Chain Name Array
| Type | Description |
| :--- | :--- |
| `uint32` | Number of unique chains (`c_count`) |
| `char[][4]` | Array of `c_count` 4-byte strings identifying unique chains. |
| padding | Padded with `0x00` to the nearest 8-byte boundary |

### 2.5 Extensible Metadata Block (TLV)
Allows forward-compatible metadata addition. Older parsers can safely skip unknown codes by advancing the stream by `section_len` + padding.

| Type | Description |
| :--- | :--- |
| `uint16` | Number of Metadata Sections (`num_meta`) |
| `uint16[]`| Array of `num_meta` 2-byte Codes (e.g., `0x0001` Resolution, `0x0002` CONECT) |
| padding | Padded with `0x00` to the nearest 8-byte boundary |

For each of the `num_meta` sections declared above, the following payload structure strictly applies in order:

| Type | Description |
| :--- | :--- |
| `uint64` | Unpadded length of the section payload (`section_len`) |
| `uint8[]` | The raw payload data |
| padding | Padded with `0x00` to the nearest 8-byte boundary |

---

## 3. Coordinate Stream (Chunked)

The coordinate stream runs from the end of the Metadata Block until EOF. It operates purely sequentially without dictionary offset lookups. 

The stream is divided into "Chunks" of residues. Each Chunk begins with an 8-byte header that defines the properties of the upcoming sequence.

### 3.1 Chunk Header (8 Bytes)
| Byte  | Type     | Description                                |
|:------|:---------|:-------------------------------------------|
| `0-3` | `uint32` | Number of Residues contained in this chunk |
| `4-7` | `uint32` | 1-Bit Event Flags (See below)              |

### 3.2 Chunk Flags
The flags dictate state changes before and after the chunk is processed.

| Bit | Hex    | Description                                                       |
|:----|:-------|:------------------------------------------------------------------|
| 0   | `0x01` | Next Chain (Advances internal `chain_names` index *before* chunk) |
| 1   | `0x02` | Next File (Advances internal `file_ids` index *before* chunk)     |
| 2   | `0x04` | Has TER (Outputs a TER line *after* parsing residues)             |
| 3   | `0x08` | Has OXT (Parses an 8-byte OXT terminal block *after* residues)    |

### 3.3 Residue & Atom Blocks

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

| 16-bit Lane | Component | Bits 3-15 (13-bit MSB Math) | Bits 0-2 (3-bit LSB Metadata) |
| :--- | :--- | :--- | :--- |
| `0` | B-factor | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bits 0-2 |
| `1` | X Coord | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bits 3-5 |
| `2` | Y Coord | Delta from Atom N-1 (or absolute LSB for Atom 1) | Occupancy Bit 6 \| AltLoc \| Charge Bit 0 |
| `3` | Z Coord | Delta from Atom N-1 (or absolute LSB for Atom 1) | Charge Bits 1-3 |

*   **Math Precision:** The 13-bit deltas provide a signed range of ±4096 (±4.096 Å or Å²).
*   **Metadata Unpacking:** Occupancy is 7 bits unsigned (0-100). AltLoc is 1 bit (Blank or 'A'). Charge is 4 bits signed (-8 to +7). 
*   **Missing Atoms:** Any atom missing from the experimental data is injected as a "Dummy Atom" (Occupancy = 0) with a coordinate delta of 0 to maintain static block sizing without branching.

#### OXT Terminal Block (8 Bytes)
If the `Has OXT` flag is set in the Chunk Header, this block appears immediately after the residue's final Atom Block. It is encoded identically to standard Atom Blocks, but its deltas are calculated relative to the absolute coordinates of Atom 0 of the preceding residue to prevent overflow.
