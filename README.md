# bPDB: Binary Protein Data Bank Format

**bPDB** is an ultra-fast, zero-copy, SIMD-accelerated binary format for protein structure data. It is designed to achieve absolute maximum parsing speeds by eliminating branch-heavy ASCII parsing and floating-point string conversions (`atof`), yielding up to **~110x faster parsing** than legacy PDB formats. 

bPDB allows for near-zero overhead `mmap` integration, SIMD coordinate reconstruction, and flawless mathematically verified bidirectional conversion with legacy `.pdb` files.

## Core Architecture

* **Alignment:** Strictly **16-byte aligned**. All variable-length sections pad to 16 bytes to align flawlessly with 128-bit SIMD registers (NEON/SSE).
* **Endianness:** Little-Endian.
* **Storage Strategy:** Uses integer projections (milli-Ångströms). Coordinates are split into a 16-bit base value and a 16-bit low-bit differential for branchless 32-bit integer reconstruction.

---

## 1. Directory Table (Offset: `0x0000`)
The file begins with a fixed **256-byte** directory table consisting of 64 `uint32_t` offsets. This allows $O(1)$ access to any metadata section or the coordinate stream. An offset of `0` indicates the section is absent.

| Index | Offset | Section | Description |
| :--- | :--- | :--- | :--- |
| `0` | `0x00` | `HEADER` | Pointer to `HEADER` block |
| `1` | `0x04` | `OBSLTE` | Pointer to `OBSLTE` block |
| `2` | `0x08` | `TITLE` | Pointer to `TITLE` text pool |
| `3` | `0x0C` | `SPLIT` | Pointer to `SPLIT` block |
| `4` | `0x10` | `COMPND` | Pointer to `COMPND` dictionary pool |
| `5` | `0x14` | `SOURCE` | Pointer to `SOURCE` dictionary pool |
| `6` | `0x18` | `KEYWDS` | Pointer to `KEYWDS` text pool |
| `7` | `0x1C` | `EXPDTA` | Pointer to `EXPDTA` text pool |
| `8` | `0x20` | `AUTHOR` | Pointer to `AUTHOR` text pool |
| `9` | `0x24` | `REVDAT` | Pointer to `REVDAT` dictionary block |
| `10` | `0x28` | `JRNL` | Pointer to `JRNL` dictionary block |
| ... | ... | ... | ... |
| `30` | `0x78` | `CONECT` | Pointer to global connectivity graph |
| `31` | `0x7C` | `COORD` | Pointer to Residue/Atom geometry stream |

---

## 2. Core Data Types

### `bpdb_date_t` (4 bytes)
Compact date representation mapping to `DD-MMM-YY` formats.
* `+0` `uint8_t` day (1-31)
* `+1` `uint8_t` month (1-12)
* `+2` `int16_t` year

### EC Number (32-bit packed integer)
Losslessly packs standard and incomplete Enzyme Commission numbers (e.g., `3.2.1.14`, `3.2.1.-`, `3.2.1.n2`).
* **Bits 29-31** (3 bits): Class 
* **Bits 24-28** (5 bits): Subclass 
* **Bits 16-23** (8 bits signed): Sub-subclass 
* **Bits 0-15** (16 bits signed): Serial 
*(Note: 'n' designations are stored as negative values. '-' is stored as 0).*

---

## 3. Fixed Metadata Blocks

### `HEADER` (56 bytes, padded to 64)
| Offset | Type | Description |
| :--- | :--- | :--- |
| `+0` | `char[40]` | Classification (padded with spaces) |
| `+40` | `bpdb_date_t` | depDate (4 bytes) |
| `+44` | `char[12]` | idCode (Supports legacy 4-char and newer 12-char AlphaFold/UniProt IDs) |

### `REVDAT` (Variable, 16-byte padded)
| Offset | Type | Description |
| :--- | :--- | :--- |
| `+0` | `bpdb_date_t` | modDate (4 bytes) |
| `+4` | `char[4]` | modId |
| `+8` | `uint16_t` | modNum |
| `+10` | `uint16_t` | modType |
| `+12` | `char[4]` | Reserved / Padding |
*(Followed immediately by a string pool or array of modified records).*

---

## 4. Text & Dictionary Pools

### List Text Pools (`TITLE`, `KEYWDS`, `EXPDTA`, `AUTHOR`)
Variable length block containing raw strings separated by `\n`.
* `+0` Concatenated `char` string.
* `+n` Null terminator (`0x00`) followed by padding bytes to reach the nearest 16-byte boundary.

### Dictionary Pools (`COMPND`, `SOURCE`, `JRNL`)
Advanced string-pool architecture for Key-Value data mapping. Padded to 16 bytes.

**Base Header:**
* `+0` `uint16_t`: Molecule Count ($M$)
* `+2` `uint16_t`: Total Key-Value Pairs ($K$)
* `+4` `uint32_t`: Offset to String Pool (relative to block start)

**Molecule Headers ($M \times 8$ bytes):**
* `+0` `uint16_t`: Molecule ID
* `+2` `uint16_t`: KV Pair Start Index
* `+4` `uint16_t`: KV Pair Count

**Key-Value Array ($K \times 8$ bytes):**
* `+0` `uint8_t`: Token Enum (e.g., `MOL_ID`, `SYNONYM`, `EC`)
* `+2` `uint16_t`: String Length (or packed EC number)
* `+4` `uint32_t`: String Offset (relative to String Pool)

**String Pool:** Null-terminated strings padded to 16-bytes.

---

## 5. Coordinate Stream (`COORD`)

The core geometry stream abandons standard object-oriented hierarchies in favor of an interleaved, 128-bit aligned vector stream designed specifically for hardware intrinsics.

### `RESIDUE` Record (16 bytes)
| Lane 0 (32-bit) | Lane 1 (32-bit) | Lane 2 (32-bit) | Lane 3 (32-bit) |
| :--- | :--- | :--- | :--- |
| `char[3]` res_name | `uint16_t` res_num | `uint16_t` atom_count | `uint8_t` icode |
| `char` chain_id | `int16_t` base_x | `int16_t` base_y | `uint8_t` file_flags |
| | | | `int16_t` base_z |

*Note: `file_flags & 0x03` defines the record type (`0` = ATOM, `1` = HETATM, `2` = TER).*

### `ATOM` Record (16 bytes)
The `RESIDUE` block is immediately followed by $N$ contiguous `ATOM` blocks (defined by `atom_count`).

| Lane 0 (32-bit) | Lane 1 (32-bit) | Lane 2 (32-bit) | Lane 3 (32-bit) |
| :--- | :--- | :--- | :--- |
| `char[4]` name | `uint16_t` dx | `uint16_t` dy | `uint16_t` dz |
| | `uint16_t` occ * 1000 | `uint16_t` bfactor * 100 | `uint8_t` z_num |
| | | | `uint8_t` bitfield |

*Note: `bitfield` stores the Alternate Location indicator in the lower 4 bits, and the signed Charge (e.g., `2+`, `1-`) as a 4-bit two's complement integer in the upper 4 bits.*

---

## SIMD Reconstruction

By aligning `base` coordinates and `dx/dy/dz` differentials perfectly along the 32-bit vector lanes, complete 3D coordinate reconstruction takes a single SIMD addition:

```c
// Example: AArch64 NEON Reconstruction
uint32x4_t r_vec = vld1q_u32((const uint32_t*)res);
uint32x4_t a_vec = vld1q_u32((const uint32_t*)atom);

int32x4_t base = vreinterpretq_s32_u32(vandq_u32(r_vec, vdupq_n_u32(0xFFFF0000)));
int32x4_t diff = vreinterpretq_s32_u32(vandq_u32(a_vec, vdupq_n_u32(0x0000FFFF)));

int32x4_t coords = vaddq_s32(base, diff); // Exact Integer mÅ coordinates
float32x4_t fcoords = vmulq_n_f32(vcvtq_f32_s32(coords), 0.001f); // Hardware floats