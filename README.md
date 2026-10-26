# Binary PDB (bPDB) Format Specification (V5)

The bPDB V5 architecture utilizes a Data-Oriented (SoA), completely unmasked coordinate atom stream. It achieves extreme compression and zero-shift SIMD parsing speeds by employing native 16-bit signed integer deltas, perfectly interleaved 16-byte New Anchor blocks, and an unbroken geometric running sum across the entire chain.

Every structure and section in the format is strictly **8-byte aligned**.

---

## 1. File Architecture Overview

A valid bPDB file consists of four sequential sections:
1. **Fixed Global Header** (64 Bytes)
2. **Metadata Bitmask & Payload** (Extensible, 8-byte aligned)
3. **Chain Clusters** (Repeated: Header -> Bytecode Sequence -> Coordinate Stream)
4. **CONECT Stream** (Flat 32-bit integer array to EOF)

---

## 2. Global Header & Metadata

### 2.1 Fixed Global Header (64 Bytes)
| Offset | Type       | Description                                                  |
|:-------|:-----------|:-------------------------------------------------------------|
| `+0`   | `char[4]`  | Magic String: `"bPDB"`                                       |
| `+4`   | `uint8`    | Format Version (`0x05`)                                      |
| `+5`   | `uint8`    | Deposition Day (1-31)                                        |
| `+6`   | `uint8`    | Deposition Month (1-12)                                      |
| `+7`   | `uint8`    | Deposition Year (`YYYY - 1973`)                              |
| `+8`   | `char[16]` | `pdbID` (null-padded)                                        |
| `+24`  | `char[40]` | `classification` (HEADER classification string, null-padded) |

### 2.2 Global 64-bit Metadata Bitmask
Immediately following the 64-byte fixed header is a single 8-byte (`uint64_t`) metadata header.
*   **Bits 0-7:** File version duplicate (`0x05`).
*   **Bits 16-63:** The 48-bit Metadata Presence Mask.

The decoder reads this mask to determine exactly which metadata payloads follow. Payloads are written in bit-sequential order, heavily padded to 8-byte boundaries:

*   **Bit 0 (`META_AF_MODE`):** If `1`, indicates an AlphaFold/Predicted model. Occupancy is statically assumed to be `1.00` and Charge `0`, bypassing all related metadata evaluation. (No payload bytes).
*   **Bit 1 (`RESOLUTION`):** `float` (4 bytes) + `0x00` pad (4 bytes).
*   **Bit 2 (`TITLE`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 3 (`EXPDTA`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 4 (`PROGRAM`):** `uint64_t` length + payload string + 8-byte alignment padding.
*   **Bit 5 (`HET_DICT`):** `uint64_t` dict_length + `uint8_t` custom_count + inline dictionary definitions + padding.
*   **Bits 6-10:** Binary arrays for `HELIX`, `SHEET`, `SSBOND`, `CISPEP`, and `SITE`. Each starts with an 8-byte `length` followed by raw `struct` dumps, aligned to 8 bytes.

---

## 3. The Chain Clusters

Coordinates are organized into sequential "Chains." Each chain declares its topology and length up front, followed by an unbroken stream of atom blocks.

### 3.1 Chain Header (8 Bytes)
| Offset | Type     | Description |
|:-------|:---------|:------------|
| `+0`   | `uint16` | **num_residues**: Number of 1-byte links in the Sequence Track. (`0xFFFF` signals EOF). |
| `+2`   | `uint16` | **num_atoms**: Total number of spatial atoms in the stream. |
| `+4`   | `uint8`  | **chain_id 1 & TER flag**: Bits `0-6` encode the 7-bit ASCII chain char. **Bit 7 (`0x80`) is the `TER` flag** (`1` = Print `TER` at chain end). |
| `+5`   | `char`   | **chain_id 2**: Second chain character (if any, else space). |
| `+6`   | `uint16` | **start_resSeq**: The starting Residue Sequence Number. |

### 3.2 Sequence Track (Bytecode Machine)
A 1-dimensional array of `uint8_t` immediately following the Chain Header. The dictionary operates on a **1-based index**.

*   **`0x00` (GAP):** Indicates a missing residue. Increments `resSeq` by 1. Parses no atoms.
*   **Bits 0-6 (Type ID):** Encodes Valid Residues (1 to 126):
    *   `1-37`: Standard Amino/Nucleic Acids.
    *   `38+`: Special termini (`NH2`) and Custom `HET_DICT` molecules.
*   **Bit 7 (`0x80` - iCode Flag):** If set, the *immediately following byte* in the array is an ASCII Insertion Code character (e.g., `'A'`).

*(The Sequence Track is padded with `0x00` to the nearest 8-byte boundary).*

---

## 4. Coordinate Atom Track (Zero-Shift SIMD Array)

To achieve native SIMD load speeds, standard atom coordinates no longer use per-lane bit-packing. Math operations execute natively on perfectly aligned signed 16-bit integers. Properties like Occupancy, Charge, and Alternate Locations are evaluated natively, and explicitly written into the stream via an interleaved 16-byte New Anchor Block only when necessary.

### 4.1 Unbroken Difference Accumulation (`Control == 0`)
A standard 8-byte block handles native intra-chain coordinate deltas up to **$\pm 32.767$ Å** (scaled by 1000). The accumulator is a pure unbroken running sum, ignoring residue boundaries entirely. Standard blocks natively imply a sparse charge of `0`.

| Bytes | Type          | Description |
|:------|:--------------|:------------|
| `0-1` | `uint16_t`    | Bits 2-15: Absolute B-factor (0.00 to 163.83). Bit 1: AltLoc Presence Flag. **Bit 0: Control Flag (`0`)**. |
| `2-3` | `int16_t`     | $\Delta X$ relative to previous atom. |
| `4-5` | `int16_t`     | $\Delta Y$ relative to previous atom. |
| `6-7` | `int16_t`     | $\Delta Z$ relative to previous atom. |

### 4.2 Absolute New Anchors & Sparse Charge (`Control == 1`)
If a spatial jump exceeds $\pm 32.767$ Å (chain gaps, isolated HETATMs, waters), OR if the implicit atom states (Occupancy or physical Charge) change, the stream emits a 16-byte absolute New Anchor. 

Because ~99% of biological atoms are uncharged, charge is not a persistent state variable. If an atom has a non-zero charge, it is forced to use a New Anchor. The immediate next standard block intrinsically drops the charge back to 0 without requiring another jump.

This leverages 31-bit native geometry via direct memory concatenation: The Base block MSBs combine perfectly with the 16-bit LSBs from the immediate next Atom block to stably construct absolute coordinates up to **$\pm 8,388.6$ Å** using native pointer loading. A 3-bit charge (-3 to +4) is compressed into the lowest bit of each coordinate MSB byte.

**1. The Base Jump Block (8 Bytes, `Control == 1`)**
Appears first. Contains perfectly interleaved 15-bit MSBs and sparse charge encoding.
*   `+0`: Absolute Occupancy (0-100) `<< 1` | **Control Flag (`1`)**.
*   `+1`: Absolute B-factor MSB.
*   `+2-3`: 15-bit Absolute X-coordinate MSB `<< 1` | Charge Bit 0.
*   `+4-5`: 15-bit Absolute Y-coordinate MSB `<< 1` | Charge Bit 1.
*   `+6-7`: 15-bit Absolute Z-coordinate MSB `<< 1` | Charge Bit 2.

**2. The Concatenated Atom Block (8 Bytes)**
Immediately follows the Base Block. Contains 16-bit LSBs.
*   `+0-1`: Absolute B-factor LSBs (14 bits) `<< 2` | AltLoc Presence Flag `<< 1` | Control Flag (`0`).
*   `+2-3`: X-coordinate LSBs (16 bits).
*   `+4-5`: Y-coordinate LSBs (16 bits).
*   `+6-7`: Z-coordinate LSBs (16 bits).

### 4.3 Implicit Alternate Location Engine
Explicit alternate location characters (`0x7F DUP` flags) are completely removed from the format. The decoder relies strictly on topological grouping natively stored in the stream. The `AltLoc Presence Flag` (`bit 1` of the Atom Block) temporarily pauses the decoder's sequence dictionary iterator while natively driving an implicit alphabetical state machine (`A -> B -> C`) for infinite overlapping conformers without overhead.

### 4.4 Terminal OXT Inference
There is no `OXT` instruction in the sequence bytecode. Upon sequence track exhaustion, the decoder processes the final remaining atoms in the chain dynamically as the terminal `OXT` atoms.

---

## 5. End of File & Connectivity

### 5.1 EOF Marker
The coordinate stream ends when a Chain Header presents a `num_residues` value of `0xFFFF`.

### 5.2 CONECT Stream
Following the EOF marker is a flat array of 32-bit integers representing atomic linkages.
1. `int32_t N` (Number of target atoms)
2. `int32_t source_atom_id`
3. $N$ sequential `int32_t target_atom_ids`

*(Padded with `0x00` to an 8-byte boundary at the physical End of File).*
