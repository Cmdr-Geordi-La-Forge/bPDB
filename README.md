# Binary PDB (bPDB) Format Specification (V3)

The V3 architecture completely drops the rigid "Chunk" boundaries in favor of an entirely self-regulating, Data-Oriented atom stream. It utilizes coordinate difference accumulation, cross-residue physical anchoring, and inline bit-flags to achieve maximum SIMD parsing speed and file compression. 

Every structure in the format is strictly **8-byte aligned**.

## 1. File Architecture Overview

The file consists of a linear progression of four sections:
1. **Fixed Global Header** (64 Bytes)
2. **Metadata Block (TLV)** (Extensible, 8-byte aligned)
3. **Chain Clusters** (Sequence bytecode followed by unbroken atom arrays)
4. **CONECT Stream** (Flat 32-bit integer array to EOF)

---

## 2. Chain Clusters

Atoms are clustered in chains. Each cluster begins with a sequence bytecode track, immediately followed by the coordinate stream for that sequence. 

### 2.1 Chain Header
The residue array is prepended by an 8-byte header establishing the chain parameters.

| Offset | Type      | Description                                                                                                                        |
|:-------|:----------|:-----------------------------------------------------------------------------------------------------------------------------------|
| `+0`   | `uint16`  | **num_residues**: The number of 1-byte links in the sequence array. (`0xFFFF` signals End of Coordinate Stream / Start of CONECT). |
| `+2`   | `uint16`  | **num_atoms**: Total atoms in this chain (used for OXT detection).                                                                 |
| `+4`   | `char[2]` | **chain_id**: The chain identifier (e.g., `'A'`, `\0`).                                                                            |
| `+6`   | `uint16`  | **start_resSeq**: The starting Residue Sequence Number.                                                                            |

### 2.2 Sequence Track (The 1-Byte Chain)
Following the header is a contiguous array of 1-byte elements, tightly encoding the chain's topology. The array is padded with `0x00` to the nearest 8-byte boundary.

*   **`0x00` (Missing Residue):** Increments `resSeq` by 1 without requesting any atoms (encodes gaps).
*   **Bits 0-6 (Type ID):** Encodes 127 possible residues:
    *   `1-37`: Standard Amino and Nucleic Acids.
    *   `38`: Amidated termination (`NH2` / `HN1`, `HN2`).
    *   `39`: Water (`HOH`).
    *   `40+`: Custom molecules defined in the `HET_DICT` metadata block.
*   **Bit 7 (Insertion Indicator):** If set to `1`, the *immediately following byte* in the array is the ASCII character for the insertion code (e.g., `'A'`).

---

## 3. Coordinate Atom Track (64-bit SIMD Array)

Immediately following the padded Sequence Track is the coordinate stream. An atom entry is an 8-byte value interpreted as four 16-bit lanes mapped to the XMM/NEON SIMD pipelines.

### 3.1 12-Bit Metadata Packing
The 3 Least Significant Bits (LSBs) of each 16-bit lane are strictly reserved for inline metadata, leaving exactly 13 bits for mathematical accumulation.

| 16-bit Lane | Bits 3-15 (13-bit MSB Math) | Bits 0-2 (3-bit LSB Metadata)                                |
|:------------|:----------------------------|:-------------------------------------------------------------|
| `0` (B-fac) | Math Payload                | Occupancy Bits `[0:2]`                                       |
| `1` (X)     | Math Payload                | Occupancy Bits `[3:5]`                                       |
| `2` (Y)     | Math Payload                | Occupancy Bit `[6]` \| Charge Bit `[0]` \| Control Bit `[0]` |
| `3` (Z)     | Math Payload                | Charge Bits `[1:2]` \| AltLoc Bit `[0]`                      |

*   **Occupancy (7 bits):** Stores absolute values from `0` to `100`.
*   **Charge (3 bits):** Stores absolute values mapped from `-4` to `+3` (by adding `+4` during encoding).
*   **AltLoc (1 bit):** If `1`, indicates an alternate conformer.
*   **Control Bit (1 bit):** If `1`, triggers the In-Stream Absolute Anchor mechanism (see below).

### 3.2 Standard Difference Accumulation (`Control == 0`)
For standard atoms, the 13 MSBs contain a signed $\Delta$ representing the physical distance from an anchor atom. 
Because the 13-bit Math MSBs support a strict signed range of $\pm 4096$, they natively support a maximum coordinate jump of **$\pm 4.096$ Å**. 
*   **Intra-Residue:** The $\Delta$ is evaluated relative to the *previous atom of the same residue*.
*   **Inter-Residue (AA):** The first atom of a new Amino Acid residue is anchored by the 3rd atom (`C`) of the previous residue.
*   **Inter-Residue (NA):** The first atom of a new Nucleic Acid residue is anchored by the 9th atom (`O3'`) of the previous residue.

### 3.3 In-Stream Absolute Anchor (`Control == 1`)
When the physical distance between two atoms exceeds $\pm 4.096$ Å (e.g., the very first atom of a chain, or isolated `HOH` molecules), the stream anchors itself natively without breaking into chunk headers. 

If an Atom block has its **Control Bit set to `1`**, the standard tracking is paused, and the next two 8-byte blocks establish a new absolute coordinate base:

**1. The Base Block (8 Bytes)**
The immediate next 8 bytes contain the low bits of the 29-bit coordinate integers and 21-bit B-factor.
*   `Byte 0`: Lower 1 byte of the 21-bit B-factor.
*   `Byte 1`: Reserved / Pad.
*   `Bytes 2-3`: Lower 2 bytes (16 bits) of the 29-bit X coordinate.
*   `Bytes 4-5`: Lower 2 bytes (16 bits) of the 29-bit Y coordinate.
*   `Bytes 6-7`: Lower 2 bytes (16 bits) of the 29-bit Z coordinate.

**2. The Anchor Atom Block (8 Bytes)**
The following 8-byte block represents the actual atom. Its metadata LSBs are read normally, but its 13-bit "differences" are treated as the high bits. They are concatenated with the Base Block to form the exact absolute integers:
*   `Absolute_X = (Atom_X_13bit << 16) | Base_X_16bit`
*   `Absolute_B = (Atom_B_13bit << 8) | Base_B_8bit`

### 3.4 Terminal OXT Inference
The OXT termination is implied purely by counting. The decoder counts down the atoms requested by the sequence array. If `num_atoms` (declared in the Chain Header) is exactly 1 greater than the expected dictionary sum, the final 8-byte block of the chain is implicitly parsed as the `OXT` atom. No explicit sequence bytecode or EOF terminator is required within the chain.
