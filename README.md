# Binary PDB (bPDB) Format Specification (V0)

The bPDB V0 architecture utilizes a Data-Oriented (SoA), completely unmasked coordinate atom stream. It achieves extreme compression and zero-shift SIMD parsing speeds by employing native 13-bit (or 15-bit for AFDB) integer deltas, 1-cycle integer string matching, and completely symmetric 16-byte absolute Anchor blocks. The geometric running sum is perfectly unbroken across the entire chain. Models are isolated into independent files to support concurrent streaming.

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
Immediately following the 64-byte fixed header is a single 8-byte (`uint64_t`) metadata header representing a flat 64-bit Presence Mask. The decoder reads this mask to determine exactly which metadata payloads follow. Payloads are written in bit-sequential order, perfectly mirroring the official PDB sections, heavily padded to 8-byte boundaries:

**Engine & Format Flags (Bits 0-15):**
*   **Bit 0 (`AF_MODE`):** Indicates an AlphaFold/Predicted model (`--af`). Occupancy is statically assumed to be `1.00` and Charge `0`.
*   **Bit 1 (`AFDB`):** Indicates an official AlphaFold Database file (`--afdb`). Drops DeepMind `REMARK 1` copyright disclaimers and `CRYST1` tags during encoding, seamlessly injecting them at decode.
*   **Bit 2 (`HYDROGENS`):** Indicates the coordinate stream natively interleaves explicit hydrogen atoms.

**Section 2: Title Section (Bits 16-31):**
*   **Bit 16 (`TITLE`):** `uint64_t` length + payload string + padding.
*   **Bit 17 (`COMPND`):** TLV-encoded COMPND blocks. EC numbers are aggressively packed into a single 32-bit integer.
*   **Bit 18 (`SOURCE`):** TLV-encoded SOURCE blocks. Integrates zero-copy NCBI taxonomy IDs to dynamically reconstruct scientific/common names.
*   **Bit 19 (`EXPDTA`):** `uint64_t` length + payload string + padding.
*   **Bit 20 (`RESOLUTION`):** `float` (4 bytes) + `0x00` pad (4 bytes).
*   **Bit 21 (`PROGRAM`):** `uint64_t` length + payload string + padding.

**Section 3: Primary Structure Section (Bits 32-35):**
*   **Bit 32 (`DBREF`):** Binary array for Database References.
*   **Bit 33 (`SEQRES`):** Contains the exact sequence chain blocks. Suppressed completely for `--afdb` files.

**Section 4: Heterogen Section (Bits 36-39):**
*   **Bit 36 (`HET_DICT`):** `uint64_t` dict_length + `uint8_t` custom_count + inline dictionary definitions + padding.

**Section 5: Secondary Structure Section (Bits 40-43):**
*   **Bit 40 (`HELIX`):** Binary array starting with an 8-byte `length` followed by raw `struct` dumps.
*   **Bit 41 (`SHEET`):** Binary array starting with an 8-byte `length` followed by raw `struct` dumps.

**Section 6: Connectivity Annotation Section (Bits 44-47):**
*   **Bit 44 (`SSBOND`):** Binary array starting with an 8-byte `length` followed by raw `struct` dumps.
*   **Bit 45 (`CISPEP`):** Binary array starting with an 8-byte `length` followed by raw `struct` dumps.

**Section 7: Miscellaneous Features Section (Bits 48-51):**
*   **Bit 48 (`SITE`):** Binary array starting with an 8-byte `length` followed by raw `struct` dumps.

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
*   **Bits 0-6 (Type ID):** Encodes Valid Residues mapped directly to the internal 256-element dictionary array.
*   **Bit 7 (`0x80` - iCode Flag):** If set, the *immediately following byte* in the array is an ASCII Insertion Code character.

## 4. Coordinate Atom Track (Dynamic Vector SIMD Array)

Math operations execute natively on perfectly aligned blocks using architecture-specific dynamic shifts (`_mm_sra_epi16`, `vshl_s16`).

### 4.1 Unbroken Difference Accumulation (`Control == 0`)
A standard 8-byte block handles native intra-chain coordinate deltas up to **± 4.095 Å** (13-bit standard) or **± 16.383 Å** (15-bit AFDB). The accumulator is a pure unbroken running sum.

| Bytes | Type     | Description                                                                                                    |
|:------|:---------|:---------------------------------------------------------------------------------------------------------------|
| `0-1` | `uint16` | Bits 3-15: $\Delta B$ (13-bit). Bit 2: Occupancy bit 6. Bit 1: AltLoc Presence Flag. **Bit 0: Control (`0`)**. |
| `2-3` | `uint16` | Bits 3-15: $\Delta X$ (13-bit). Bits 0-2: Occupancy bits 3-5.                                                  |
| `4-5` | `uint16` | Bits 3-15: $\Delta Y$ (13-bit). Bits 0-2: Occupancy bits 0-2.                                                  |
| `6-7` | `uint16` | Bits 3-15: $\Delta Z$ (13-bit). Bits 0-2: Charge (3-bit unsigned, evaluate as `val - 3`).                      |

*(For AFDB mode, the lowest 3 bits of all vectors are natively merged into a 15-bit delta stream, overriding the explicit occupancy/charge packing).*

### 4.2 Absolute New Anchors (`Control == 1`)
If a spatial jump exceeds the delta limit, the stream emits a 16-byte absolute New Anchor. The Anchor is **100% symmetric** with the standard atom block lane order (B, X, Y, Z), forming a precise 28-bit coordinate space.

All explicit metadata (Occupancy, Charge, AltLoc) is offloaded entirely to the remaining 12 bits of the LSB block.

## 5. End of File & Connectivity

### 5.1 Fragment Boundary & End of Model Markers
*   **Fragment Boundary (`0x3FFE`):** Indicates a seamless transition between split input files (e.g., AF-DB fragments).
*   **End of Model (`0x3FFF`):** Indicates the absolute end of the spatial coordinate stream for the current model.

### 5.2 CONECT Stream
Following the `0x3FFF` EOF marker is a flat array of 32-bit integers representing atomic linkages.
1. `int32_t N` (Number of target atoms)
2. `int32_t source_atom_id`
3. $N$ sequential `int32_t target_atom_ids`

## 6. Command Line Interface

*   **Encode:** `./bpdb encode [--af] [--afdb] [-h|--hydrogens] <in.pdb> <out.bpdb>`
*   **Decode:** `./bpdb decode <in.bpdb> <out.pdb>`
*   **Benchmark:** `./bpdb bench <in.bpdb>`
*   **Statistics:** `./bpdb stat <in.bpdb>`
