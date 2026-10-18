# Binary PDB (bPDB) Format Specification

A SIMD-optimized, heavily compressed binary format for macromolecular coordinate data.

## 1. File Architecture Overview

The bPDB format utilizes a strictly 8-byte aligned, linear TLV architecture. It eliminates random-access pointer chasing in favor of sequential stream processing, optimizing hardware prefetching.

1. **Fixed Global Header (64 Bytes)**
2. **Extensible Metadata Block (TLV)**
3. **Coordinate Stream (Chunked)**

---

## 2. Global Header & Metadata

### 2.1 Fixed Global Header (64 Bytes)
The file strictly begins with this 64-byte block, containing the mandatory identity of the structure.

| Offset | Type       | Description                                                  |
|:-------|:-----------|:-------------------------------------------------------------|
| `+0`   | `char[4]`  | Magic String: `"bPDB"`                                       |
| `+4`   | `uint8`    | Format Version (e.g., `0x01`)                                |
| `+5`   | `uint8`    | Deposition Day (1-31)                                        |
| `+6`   | `uint8`    | Deposition Month (1-12)                                      |
| `+7`   | `uint8`    | Deposition Year (`YYYY - 1973`)                              |
| `+8`   | `char[16]` | `pdbID` (PDB ID or AlphaFold ID, null-padded)                |
| `+24`  | `char[40]` | `classification` (HEADER classification string, null-padded) |

### 2.2 Extensible Metadata Block (TLV)
Allows forward-compatible metadata addition. Older parsers can skip unknown 8-byte codes by advancing the stream by `section_len` + padding.

| Type       | Description                                                                                                                                                                          |
|:-----------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `uint16`   | Number of Metadata Sections (`num_meta`)                                                                                                                                             |
| `uint16[]` | Array of `num_meta` 2-byte Codes (`0x0001` TITLE, `0x0002` HET_DICT, `0x0003` EXPDTA, `0x0004` RESOLUTION, `0x0005` PROGRAM, `0x0006` SSBOND, `0x0007` CISPEP, `0x0008` HELIX, `0x0009` SHEET, `0x000A` SITE) |
| padding    | Padded with `0x00` to the nearest 8-byte boundary                                                                                                                                    |

For each of the `num_meta` sections declared above, the following payload structure strictly applies in order:

| Type      | Description                                            |
|:----------|:-------------------------------------------------------|
| `uint64`  | Unpadded length of the section payload (`section_len`) |
| `uint8[]` | The raw payload data                                   |
| padding   | Padded with `0x00` to the nearest 8-byte boundary      |

### 2.3 Metadata Payload Formats

#### META_HET_DICT (0x0002)
Defines file-specific HETATMs so they can be encoded as compact 8-byte blocks exactly like standard residues. The custom IDs begin dynamically at `37`.
* `uint8`: Number of custom residues defined.
* For each custom residue:
  * `char[4]`: Residue Name (e.g., "NAG ")
  * `uint8`: Number of Atoms
  * For each atom:
    * `char[4]`: Atom Name
    * `char[2]`: Element

#### General String Data (0x0001, 0x0003, 0x0005)
* **META_TITLE (0x0001):** Raw string payload containing the concatenated title.
* **META_EXPDTA (0x0003):** Raw string payload containing the experimental data string.
* **META_PROGRAM (0x0005):** Fixed 56-byte `char` array containing the combined refinement program and version (e.g., `"NULL vNULL"`). Does not require padding calculations as `56 % 8 == 0`.

#### META_RESOLUTION (0x0004)
* `float`: 4-byte resolution value. 
* `uint32`: 4 bytes of `0x00` padding to complete the 8-byte alignment.

### 2.4 Structural Metadata & `PdbResId` Anchor
Secondary structures and linkages (`0x0006` through `0x000A`) utilize a highly compressed 8-byte anchor to represent specific residues without string duplication:

```c
typedef struct {
    char resName[4]; // e.g., "ALA\0"
    char chain;      // e.g., 'A'
    char iCode;      // e.g., ' '
    uint16_t resSeq; // e.g., 145
} PdbResId;
```

* **META_SSBOND (0x0006):** Contains `PdbResId res1`, `PdbResId res2`, `char sym1[8]`, `char sym2[8]`, `float length`.
* **META_CISPEP (0x0007):** Contains `PdbResId res1`, `PdbResId res2`, `int32 modNum`, `float angle`.
* **META_HELIX (0x0008):** Contains `PdbResId start`, `PdbResId end`, `char helixID[4]`, `uint8 helix_class`, `uint8 pad`, `uint16 length`, `char comment[32]`.
* **META_SHEET (0x0009):** Contains `PdbResId start`, `PdbResId end`, `char sheetID[4]`, `int8 num_strands`, `int8 strand_num`, `int8 sense`, `uint8 pad`.
* **META_SITE (0x000A):** A header containing `char siteID[4]` and `uint32 num_res`, followed by `num_res` sequential `PdbResId` structs.

---

## 3. Coordinate Stream (Chunked)

The coordinate stream runs from the end of the Metadata Block until EOF. It operates purely sequentially. The stream is divided into "Chunks".

### 3.1 Chunk Header (32 Bytes)
Each Chunk begins with a 32-byte header that explicitly identifies its origin and establishes a 32-bit Absolute Global Anchor to prevent coordinate overflow in massive macromolecular structures.

| Offset | Type       | Description                                               |
|:-------|:-----------|:----------------------------------------------------------|
| `+0`   | `uint16`   | Number of Residues contained in this chunk                |
| `+2`   | `uint8`    | Event Flags (`0x04` = Has TER, `0x08` = Has TER + OXT)    |
| `+3`   | `uint8`    | Fragment ID (e.g., `1` for AlphaFold `-F1`. `0` for none) |
| `+4`   | `char[2]`  | Chain ID (e.g., `" A"`, `"AA"`)                           |
| `+6`   | `char[14]` | File ID (e.g., `"1V0Q"`, `"A0A017SE85"`)                  |
| `+20`  | `int32`    | **Anchor X**: Absolute X coordinate of the first atom     |
| `+24`  | `int32`    | **Anchor Y**: Absolute Y coordinate of the first atom     |
| `+28`  | `int32`    | **Anchor Z**: Absolute Z coordinate of the first atom     |

### 3.2 Residue & Atom Blocks

Every residue (standard ATOM or file-specific HETATM) strictly consists of 1 Residue Block followed by a static number of Atom Blocks.

#### Residue Block (8 Bytes)
Establishes the base coordinate relative to the 32-bit Chunk Anchor.

| Byte Offset | Type    | Description                                                           |
|:------------|:--------|:----------------------------------------------------------------------|
| `+0`        | `uint8` | Residue Type ID (0-36 standard, 37+ from HET_DICT).                   |
| `+1`        | `uint8` | **B-factor MSB**: Top 8 bits of absolute 22-bit B-factor for Atom 1.  |
| `+2`        | `uint8` | **resSeq LSB**: Bottom 8 bits of 16-bit Residue Sequence Number.      |
| `+3`        | `int8`  | **X-base MSB**: Top 8 bits of signed 22-bit distance to Chunk Anchor. |
| `+4`        | `uint8` | **resSeq MSB**: Top 8 bits of 16-bit Residue Sequence Number.         |
| `+5`        | `int8`  | **Y-base MSB**: Top 8 bits of signed 22-bit distance to Chunk Anchor. |
| `+6`        | `char`  | **iCode**: Insertion code character (space or null if none).          |
| `+7`        | `int8`  | **Z-base MSB**: Top 8 bits of signed 22-bit distance to Chunk Anchor. |

#### Atom Block (8 Bytes)
Interpreted as four 16-bit lanes mapped directly to a uniform 2-bit shift SIMD pipeline.

| 16-bit Lane | Component | Bits 2-15 (14-bit MSB Math)                           | Bits 0-1 (2-bit LSB Metadata)     |
|:------------|:----------|:------------------------------------------------------|:----------------------------------|
| `0`         | B-factor  | Delta from Atom N-1 (or absolute LSB for Atom 1)      | Occupancy Bits 0-1                |
| `1`         | X Coord   | Delta from Atom N-1 (or anchor offset LSB for Atom 1) | Occupancy Bits 2-3                |
| `2`         | Y Coord   | Delta from Atom N-1 (or anchor offset LSB for Atom 1) | Occupancy Bits 4-5                |
| `3`         | Z Coord   | Delta from Atom N-1 (or anchor offset LSB for Atom 1) | Occupancy Bit 6 \| AltLoc (1 bit) |

*   **Math Precision:** The 14-bit deltas provide an extremely robust signed range of ±8192 (±8.192 Å or Å²), absorbing long bounds and disordered loops natively.
*   **Global Bounds:** The 22-bit Residue Base Offset provides a ±2097.152 Å radius from the Chunk Anchor.
*   **Metadata Unpacking:** Occupancy is 7 bits unsigned (0-100). AltLoc is 1 bit (Blank or 'A').

#### OXT Terminal Block (8 Bytes)
If the `Has OXT` flag (`0x08`) is set in the Chunk Header, this block appears immediately after the residue's final Atom Block. It is encoded identically to standard Atom Blocks, but its deltas are calculated relative to the absolute coordinates of Atom 0 of the preceding residue to prevent intra-residue overflow.

---

## 4. Connectivity Section (CONECT)
Immediately following the final coordinate chunk, an EOF Marker Chunk signals the end of the coordinate stream and the beginning of the `CONECT` payload.
* **EOF Marker (32 Bytes):** A standard chunk header where `num_residues` is explicitly set to `0xFFFF`.

To guarantee extremely fast parsing and zero runtime memory allocations, the connectivity data is stored as a flat, sequential 32-bit integer stream directly mirroring the PDB structure.

#### Connectivity Stream 
For every source atom with connections, the sequence is written sequentially:
1. `int32_t num_targets` (The number of target atoms $N$)
2. `int32_t source_atom_id`
3. `int32_t target_1`
4. ... `int32_t target_N`

The stream is padded with `0x00` to the nearest 8-byte boundary at the physical End of File.
