# Binary PDB (bPDB) Format Specification (V0)

The bPDB V0 architecture utilizes a Data-Oriented (SoA), completely unmasked coordinate atom stream. It achieves extreme compression and zero-shift SIMD parsing speeds by employing native 13-bit integer deltas and completely symmetric 16-byte absolute Anchor blocks. The geometric running sum is perfectly unbroken across the entire chain. Models are isolated into independent files to support concurrent streaming.

Every structure and section in the format is strictly **8-byte aligned**.

## 1. File Architecture Overview

A valid bPDB file consists of four sequential sections:
1. **Fixed Global Header** (64 Bytes)
2. **Metadata Bitmask & Payload** (Extensible, 8-byte aligned)
3. **Chain Clusters** (Repeated: Header -> Bytecode Sequence -> Coordinate Stream)
4. **CONECT Stream** (Flat 32-bit integer array to EOF)

## 2. Global Header & Metadata

### 2.1 Fixed Global Header (64 Bytes)
| Offset | Type       | Description                                                   |
|:-------|:-----------|:--------------------------------------------------------------|
| `+0`   | `char[4]`  | Magic String: `"bPDB"`                                        |
| `+4`   | `uint8`    | Format Version (`0x00`)                                       |
| `+5`   | `uint8`    | Deposition Day (1-31)                                         |
| `+6`   | `uint8`    | Deposition Month (1-12)                                       |
| `+7`   | `uint8`    | Deposition Year (`YYYY - 1973`)                               |
| `+8`   | `char[12]` | `pdbID` (null-padded, stripped of AF DB prefixes/suffixes)    |
| `+20`  | `uint16`   | `total_models_minus_1` (0 indicates a single-model structure) |
| `+22`  | `uint16`   | `current_model_minus_1` (0 indicates the first model)         |
| `+24`  | `char[40]` | `classification` (HEADER classification string, null-padded)  |

### 2.2 Global 64-bit Metadata Bitmask
Immediately following the 64-byte fixed header is a single 8-byte (`uint64_t`) metadata header.
*   **Bits 0-7:** File version duplicate (`0x00`).
*   **Bits 16-63:** The 48-bit Metadata Presence Mask.

The decoder reads this mask to determine exactly which metadata payloads follow. Payloads are written in bit-sequential order, heavily padded to 8-byte boundaries:

*   **Bit 0 (`META_AF_MODE`):** If `1`, indicates an AlphaFold/Predicted model. Occupancy is statically assumed to be `1.00` and Charge `0`.
*   **Bit 1 (`RESOLUTION`):** `float` (4 bytes) + `0x00` pad (4 bytes).
*   **Bit 2 (`TITLE`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 3 (`EXPDTA`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 4 (`PROGRAM`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 5 (`HET_DICT`):** `uint64_t` dict_length + `uint8_t` custom_count + inline dictionary definitions + padding.
*   **Bits 6-10:** Binary arrays for `HELIX`, `SHEET`, `SSBOND`, `CISPEP`, and `SITE`. Each starts with an 8-byte `length` followed by raw `struct` dumps, aligned to 8 bytes.
*   **Bit 11 (`HYDROGENS`):** If `1`, indicates the coordinate stream natively interleaves explicit hydrogen atoms.

## 3. The Chain Clusters

Coordinates are organized into sequential "Chains." Each chain declares its topology and length up front, followed by an unbroken stream of atom blocks.

### 3.1 Chain Header (8 Bytes)
| Offset | Type     | Description                                                                                                                     |
|:-------|:---------|:--------------------------------------------------------------------------------------------------------------------------------|
| `+0`   | `uint32` | **num_residues** (14-bit MSB) & **num_atoms** (18-bit LSB). `num_residues` = `0x3FFF` signals EOF, `0x3FFE` signals a fragment. |
| `+4`   | `uint8`  | **chain_id 1 & N-Term flag**: Bits `0-6` encode the 7-bit ASCII chain char. **Bit 7 (`0x80`) is the `is_n_term` flag**.         |
| `+5`   | `uint8`  | **chain_id 2 & C-Term flag**: Bits `0-6` encode the second ASCII chain char. **Bit 7 (`0x80`) is the `is_c_term` flag**.        |
| `+6`   | `int16`  | **start_resSeq & TER flag**: Signed starting Residue Sequence Number shifted left (`resSeq << 1`). **Bit 0 is the `TER` flag**. |

### 3.2 Sequence Track (Bytecode Machine)
A 1-dimensional array of `uint8_t` immediately following the Chain Header.
*   **`0x00` (GAP):** Indicates a missing residue. Increments `resSeq` by 1. Parses no atoms.
*   **Bits 0-6 (Type ID):** Encodes Valid Residues mapped directly to the internal 256-element dictionary array:
    *   **`0x01-0x1F`: Standard Nucleic Acids** (e.g., A, C, G, U, DA)
    *   **`0x20-0x5F`: Standard Amino Acids** (e.g., ALA, ARG, SEC, PYL)
    *   **`0x60-0x7F`: Solvents** (e.g., HOH)
    *   **`0x80+`:** Custom `HET_DICT` molecules.
*   **Bit 7 (`0x80` - iCode Flag):** If set, the *immediately following byte* in the array is an ASCII Insertion Code character (e.g., `'A'`).

*(The Sequence Track is padded with `0x00` to the nearest 8-byte boundary).*

## 4. Coordinate Atom Track (Zero-Shift SIMD Array)

Math operations execute natively on perfectly aligned blocks. Properties like Occupancy, Charge, and Alternate Locations are evaluated natively within the 8-byte atom lanes, dropping state-tracking requirements.

### 4.1 Unbroken Difference Accumulation (`Control == 0`)
A standard 8-byte block handles native intra-chain coordinate deltas up to **± 4.095 Å** (scaled by 1000). The accumulator is a pure unbroken running sum, ignoring residue boundaries entirely.

| Bytes | Type     | Description                                                                                                    |
|:------|:---------|:---------------------------------------------------------------------------------------------------------------|
| `0-1` | `uint16` | Bits 3-15: $\Delta B$ (13-bit). Bit 2: Occupancy bit 6. Bit 1: AltLoc Presence Flag. **Bit 0: Control (`0`)**. |
| `2-3` | `uint16` | Bits 3-15: $\Delta X$ (13-bit). Bits 0-2: Occupancy bits 3-5.                                                  |
| `4-5` | `uint16` | Bits 3-15: $\Delta Y$ (13-bit). Bits 0-2: Occupancy bits 0-2.                                                  |
| `6-7` | `uint16` | Bits 3-15: $\Delta Z$ (13-bit). Bits 0-2: Charge (3-bit unsigned, evaluate as `val - 3`).                      |

### 4.2 Absolute New Anchors (`Control == 1`)
If a spatial jump exceeds ± 4.095 Å (chain gaps, isolated HETATMs, waters), the stream emits a 16-byte absolute New Anchor. The Anchor is **100% symmetric** with the standard atom block lane order (B, X, Y, Z), forming a precise 28-bit coordinate space.

All explicit metadata (Occupancy, Charge, AltLoc) is offloaded entirely to the remaining 12 bits of the LSB block.

**1. The MSB Block (Bytes 0-7, `Control == 1`)**
Contains the 15-bit Most Significant Bits of the absolute values.
*   `+0-1`: 15-bit B-factor MSB `<< 1` | `1` (Control Flag).
*   `+2-3`: 15-bit Absolute X MSB `<< 1` | `0`.
*   `+4-5`: 15-bit Absolute Y MSB `<< 1` | `0`.
*   `+6-7`: 15-bit Absolute Z MSB `<< 1` | `0`.

**2. The LSB Block (Bytes 8-15, `Control == 0`)**
Immediately follows the MSB block and exactly mimics the standard $\Delta$ block layout, carrying the 13-bit LSBs alongside the native metadata (Occupancy and Charge).

### 4.3 Implicit Alternate Location Engine
Explicit alternate location characters (`0x7F DUP` flags) are completely removed from the format. The decoder relies strictly on topological grouping natively stored in the stream. The `AltLoc Presence Flag` temporarily pauses the decoder's sequence dictionary iterator while natively driving an implicit alphabetical state machine (`A -> B -> C`) for infinite overlapping conformers without overhead.

### 4.4 Split-Topology Dictionary & Hydrogen Stream
Heavy atoms and hydrogen atoms are strictly partitioned into parallel data tracks (`h_atoms`, `hy_atoms`) with dedicated topological bounds (`start_norm`, `end_norm`, `end_cterm`). The decoder effortlessly snaps its local dictionary iteration bounds to include or exclude terminal variants based on the header flags.

## 5. End of File & Connectivity

### 5.1 Fragment Boundary & End of Model Markers
*   **Fragment Boundary (`0x3FFE`):** Indicates a seamless transition between split input files (e.g., AF-DB fragments).
*   **End of Model (`0x3FFF`):** Indicates the absolute end of the spatial coordinate stream for the current model, triggering the decoder to break the coordinate loop.

### 5.2 CONECT Stream
Following the `0x3FFF` EOF marker is a flat array of 32-bit integers representing atomic linkages.
1. `int32_t N` (Number of target atoms)
2. `int32_t source_atom_id`
3. $N$ sequential `int32_t target_atom_ids`

*(Padded with `0x00` to an 8-byte boundary at the physical End of File).*
