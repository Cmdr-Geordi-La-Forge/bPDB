# Binary PDB (bPDB) Format Specification (V3)

The bPDB V3 architecture utilizes a completely self-regulating, Data-Oriented (SoA) atom stream. It achieves extreme compression and branchless SIMD parsing speeds by employing coordinate difference accumulation, cross-residue geometry anchoring, and in-stream control bits, eliminating the need for rigid structural chunks.

Every structure and section in the format is strictly **8-byte aligned**.

---

## 1. File Architecture Overview

A valid bPDB file consists of four sequential sections:
1. **Fixed Global Header** (64 Bytes)
2. **Metadata Block (TLV)** (Extensible, 8-byte aligned)
3. **Chain Clusters** (Repeated: Header -> Bytecode Sequence -> Coordinate Stream)
4. **CONECT Stream** (Flat 32-bit integer array to EOF)

---

## 2. Global Header & Metadata

### 2.1 Fixed Global Header (64 Bytes)
| Offset | Type       | Description                                                  |
|:-------|:-----------|:-------------------------------------------------------------|
| `+0`   | `char[4]`  | Magic String: `"bPDB"`                                       |
| `+4`   | `uint8`    | Format Version (`0x03`)                                      |
| `+5`   | `uint8`    | Deposition Day (1-31)                                        |
| `+6`   | `uint8`    | Deposition Month (1-12)                                      |
| `+7`   | `uint8`    | Deposition Year (`YYYY - 1973`)                              |
| `+8`   | `char[16]` | `pdbID` (null-padded)                                        |
| `+24`  | `char[40]` | `classification` (HEADER classification string, null-padded) |

### 2.2 Extensible Metadata Block (TLV)
* `uint16_t num_meta`: Total number of metadata sections present.
* `uint16_t codes[num_meta]`: Array of section IDs (`0x0001` TITLE, `0x0002` HET_DICT, etc).
* *[Padded with `0x00` to an 8-byte boundary]*
* For each section defined above, sequentially: 
  * `uint64_t section_len` (Unpadded length of the payload).
  * `uint8_t[] payload` (The raw metadata bytes).
  * *[Padded with `0x00` to an 8-byte boundary]*

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
*   **`0x7F` (DUP):** Flags an Alternate Conformer. The *next* sequence byte belongs to the same `resSeq` and advances the `AltLoc` char (e.g., `'A'` to `'B'`).
*   **Bits 0-6 (Type ID):** Encodes Valid Residues (1 to 126):
    *   `1-37`: Standard Amino/Nucleic Acids.
    *   `38+`: Special termini (`NH2`) and Custom `HET_DICT` molecules.
*   **Bit 7 (`0x80` - iCode Flag):** If set, the *immediately following byte* in the array is an ASCII Insertion Code character (e.g., `'A'`).

*(The Sequence Track is padded with `0x00` to the nearest 8-byte boundary).*

---

## 4. Coordinate Atom Track (64-bit SIMD Array)

Immediately following the sequence padding is the coordinate stream. An atom entry is an 8-byte (`uint64_t`) block. The CPU evaluates it as four 16-bit lanes mapped directly to SIMD pipelines.

*   **Missing Atoms:** Fully unresolved atoms are encoded as exactly `0x0000000000000000`. They execute as a 1-cycle bypass, contributing a mathematically safe $0$ $\Delta$ to the accumulator.

### 4.1 Metadata LSB Packing
The 3 Least Significant Bits (LSBs) of each 16-bit lane (12 bits total) are reserved for inline metadata.

| Lane (16-bit) | Bits 3-15 (13-bit MSB Math Payload) | Bits 0-2 (3-bit LSB Metadata) |
|:--------------|:------------------------------------|:------------------------------|
| `0` (B-fac)   | Spatial Math payload                | Occupancy Bits `[0:2]`        |
| `1` (X Coord) | Spatial Math payload                | Occupancy Bits `[3:5]`        |
| `2` (Y Coord) | Spatial Math payload                | Occupancy Bit `[6]` \| Charge Bit `[0]` \| **Control Bit** |
| `3` (Z Coord) | Spatial Math payload                | Charge Bits `[1:2]` \| AltLoc Present Bit |

*   **Occupancy (7 bits):** Absolute value `0` to `100`.
*   **Charge (3 bits):** Absolute mapping (`Actual_Charge + 4`), storing values `0` to `7`.
*   **AltLoc (1 bit):** Indicates presence of an alternate conformer.
*   **Control (1 bit):** Triggers absolute spatial anchoring.

### 4.2 Relative Difference Accumulation (`Control == 0`)
The standard 8-byte atom block contains a signed 13-bit $\Delta$ relative to a cached anchor. Because 13 bits natively support a range of $\pm 4.096$ Å, the stream traces physical covalent bonds rather than generic array indices:Intra-Residue (Backbone): $\Delta$ is evaluated relative to the previously parsed atom (e.g., `N` $\rightarrow$ `CA` $\rightarrow$ `C` $\rightarrow$ `O`).Intra-Residue (Sidechain): The $\Delta$ for the `CB` atom is evaluated strictly against a cached `CA` accumulator. Subsequent sidechain atoms accumulate linearly off `CB`.Inter-Residue (Polymer Bond): The first atom of a new residue (`N` or `P`) is anchored by a cached accumulator from the previous residue's active linkage (`C` or `O3'`).

### 4.2 Relative Difference Accumulation (`Control == 0`)
The standard 8-byte atom block contains a signed 13-bit $\Delta$ (scaled by $1000$ for coords, $100$ for B-factor).
Because 13 bits provide a signed range of $\pm 4096$, the stream natively bridges spatial gaps up to **$\pm 4.096$ Å**.
*   **Intra-Residue:** $\Delta$ from the previously parsed atom.
*   **Inter-Residue:** When crossing a residue boundary, the $\Delta$ is evaluated against a cached anchor to bridge the physical polymer bond:
    *   *Amino Acids:* Anchored to the 3rd atom (`C`).
    *   *Nucleic Acids:* Anchored to the 9th atom (`O3'`).

### 4.3 In-Stream Absolute Anchor (`Control == 1`)
If a distance exceeds $\pm 4.096$ Å (e.g., chain starts or isolated `HOH` molecules), the encoder explicitly sets the Control Bit to `1` on the atom. The CPU reconstructs absolute coordinates using exactly **16 bytes**:

**1. The High/Metadata Block (8 Bytes)**
This block carries the atom's valid LSB metadata. However, its 13-bit Math Payload MSBs act as the **High Bits** of the absolute coordinates.

**2. The Low/Base Block (8 Bytes)**
The immediate next 8 bytes contain the low bits of the absolute coordinates.
*   `Byte 0`: Low 8 bits of the B-factor.
*   `Byte 1`: Reserved / Empty.
*   `Bytes 2-3`: Low 16 bits of the X coordinate.
*   `Bytes 4-5`: Low 16 bits of the Y coordinate.
*   `Bytes 6-7`: Low 16 bits of the Z coordinate.

**Execution:** The absolute spatial anchor is seamlessly reconstructed via bitwise concatenation:
*   `Abs_X = (High_Block_X_13bit << 16) | Low_Block_X_16bit` (Result: 29-bit coordinate)
*   `Abs_B = (High_Block_B_13bit << 8) | Low_Block_B_8bit` (Result: 21-bit B-factor)

*(Note: The first atom of any `HETATM` implicitly forces `Control = 1` to cleanly severe it from polymer backbone geometry).*

### 4.4 Terminal OXT Inference
There is no `OXT` instruction in the sequence bytecode. Instead, the decoder dynamically tallies the total atoms expected from the sequence dictionary. 
If the Chain Header's `num_atoms` is exactly **`expected_atoms + 1`**, the decoder automatically processes the final spatial atom block as the terminal `OXT` atom, anchoring it to the 3rd atom (`C`) of the final residue.

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
