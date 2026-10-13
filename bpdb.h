#ifndef BPDB_H
#define BPDB_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
#elif defined(__SSE2__) || defined(__x86_64__)
    #include <emmintrin.h>
#endif

#define DIR_HEADER    0
#define DIR_OBSLTE    1
#define DIR_TITLE     2
#define DIR_SPLIT     3
#define DIR_COMPND    4
#define DIR_SOURCE    5
#define DIR_KEYWDS    6
#define DIR_EXPDTA    7
#define DIR_AUTHOR    8
#define DIR_REVDAT    9
#define DIR_JRNL     10
#define DIR_CAVEAT   11
#define DIR_NUMMDL   12
#define DIR_MDLTYP   13
#define DIR_SPRSDE   14
#define DIR_REMARKS  15
#define DIR_REMARK_1 16
#define DIR_REMARK_2 17
#define DIR_REMARK_3 18
#define DIR_REMARK_4 19
#define DIR_REMARK_5 20
#define DIR_CONECT   30
#define DIR_COORD    31

#define BPDB_MOL_FIELD_MOLECULE             0
#define BPDB_MOL_FIELD_CHAIN                1
#define BPDB_MOL_FIELD_FRAGMENT             2
#define BPDB_MOL_FIELD_SYNONYM              3
#define BPDB_MOL_FIELD_EC                   4
#define BPDB_MOL_FIELD_ENGINEERED           5
#define BPDB_MOL_FIELD_MUTATION             6
#define BPDB_MOL_FIELD_COMPND_OTHER         7
#define BPDB_MOL_FIELD_SYNTHETIC            8
#define BPDB_MOL_FIELD_ORG_SCIENTIFIC       9
#define BPDB_MOL_FIELD_ORG_COMMON           10
#define BPDB_MOL_FIELD_ORG_TAXID            11
#define BPDB_MOL_FIELD_STRAIN               12
#define BPDB_MOL_FIELD_VARIANT              13
#define BPDB_MOL_FIELD_CELL_LINE            14
#define BPDB_MOL_FIELD_ATCC                 15
#define BPDB_MOL_FIELD_ORGAN                16
#define BPDB_MOL_FIELD_TISSUE               17
#define BPDB_MOL_FIELD_CELL                 18
#define BPDB_MOL_FIELD_ORGANELLE            19
#define BPDB_MOL_FIELD_SECRETION            20
#define BPDB_MOL_FIELD_CELL_LOC             21
#define BPDB_MOL_FIELD_PLASMID              22
#define BPDB_MOL_FIELD_GENE                 23
#define BPDB_MOL_FIELD_EXP_SYSTEM           24
#define BPDB_MOL_FIELD_EXP_SYS_COMMON       25
#define BPDB_MOL_FIELD_EXP_SYS_TAXID        26
#define BPDB_MOL_FIELD_EXP_SYS_STRAIN       27
#define BPDB_MOL_FIELD_EXP_SYS_VARIANT      28
#define BPDB_MOL_FIELD_EXP_SYS_CELL_LINE    29
#define BPDB_MOL_FIELD_EXP_SYS_ATCC         30
#define BPDB_MOL_FIELD_EXP_SYS_ORGAN        31
#define BPDB_MOL_FIELD_EXP_SYS_TISSUE       32
#define BPDB_MOL_FIELD_EXP_SYS_CELL         33
#define BPDB_MOL_FIELD_EXP_SYS_ORGANELLE    34
#define BPDB_MOL_FIELD_EXP_SYS_CELL_LOC     35
#define BPDB_MOL_FIELD_EXP_SYS_VECTOR_TYPE  36
#define BPDB_MOL_FIELD_EXP_SYS_VECTOR       37
#define BPDB_MOL_FIELD_EXP_SYS_PLASMID      38
#define BPDB_MOL_FIELD_EXP_SYS_GENE         39
#define BPDB_MOL_FIELD_SOURCE_OTHER         40
#define BPDB_MOL_FIELD_MAX                  41

#define BPDB_JRNL_FIELD_AUTH 0
#define BPDB_JRNL_FIELD_TITL 1
#define BPDB_JRNL_FIELD_EDIT 2
#define BPDB_JRNL_FIELD_REF  3
#define BPDB_JRNL_FIELD_PUBL 4
#define BPDB_JRNL_FIELD_REFN 5
#define BPDB_JRNL_FIELD_PMID 6
#define BPDB_JRNL_FIELD_DOI  7
#define BPDB_JRNL_FIELD_MAX  8

#pragma pack(push, 1)

typedef struct {
    uint8_t day;
    uint8_t month;
    int16_t year;
} bpdb_date_t;

typedef struct {
    char classification[40];
    bpdb_date_t depDate;
    char idCode[12];
} bpdb_header_t;

typedef struct {
    uint16_t num_molecules;
    uint16_t num_fields; 
} bpdb_mol_meta_header_t;

typedef struct {
    char name[3];
    char chain_id;
    uint16_t res_num;
    int16_t base_x;
    uint16_t atom_count;
    int16_t base_y;
    uint8_t icode;
    uint8_t file_flags;
    int16_t base_z;
} bpdb_residue_t;

typedef struct {
    char name[4];
    uint16_t dx;
    uint16_t occupancy;
    uint16_t dy;
    uint16_t bfactor;
    uint16_t dz;
    uint8_t z_num;
    uint8_t bitfield;
} bpdb_atom_t;

typedef struct {
    bpdb_date_t date;
    char idCode[12];
    uint16_t count;
} bpdb_id_list_t;

typedef struct {
    char idCode[12];
} bpdb_caveat_t;

typedef struct {
    bpdb_date_t modDate;
    char modId[4];
    uint16_t modNum;
    uint16_t modType;
    uint16_t num_records;
} bpdb_revdat_rec_t;

#pragma pack(pop)

static inline uint32_t pack_ec(int c, int sc, int ssc, int ser) {
    return ((c & 0x07) << 29) | ((sc & 0x1F) << 24) | 
           (((uint32_t)(ssc & 0xFF)) << 16) | ((uint32_t)(ser & 0xFFFF));
}

static inline uint32_t parse_ec_string(const char* str) {
    int p[4] = {0, 0, 0, 0}, idx = 0;
    while (*str && idx < 4) {
        if (*str == '-') { p[idx++] = 0; str++; } 
        else if (*str == 'n') { p[idx++] = -atoi(str + 1); while (*str && *str != '.') str++; } 
        else if (*str >= '0' && *str <= '9') { p[idx++] = atoi(str); while (*str >= '0' && *str <= '9') str++; } 
        else str++;
    }
    return pack_ec(p[0], p[1], p[2], p[3]);
}

static inline void unpack_ec_array(const uint32_t* ec_ptr, char* out_buf) {
    uint32_t ec_count = ec_ptr[0];
    char* p = out_buf;
    for (uint32_t k = 0; k < ec_count; k++) {
        uint32_t ec = ec_ptr[k+1];
        int c = (ec >> 29) & 0x07;
        int sc = (ec >> 24) & 0x1F;
        int ssc = (int8_t)((ec >> 16) & 0xFF);
        int ser = (int16_t)(ec & 0xFFFF);
        
        if (k > 0) {
            *p++ = ','; *p++ = ' ';
        }
        p += sprintf(p, "%d.%d.", c, sc);
        
        if (ssc == 0) p += sprintf(p, "-.");
        else if (ssc < 0) p += sprintf(p, "n%d.", -ssc);
        else p += sprintf(p, "%d.", ssc);
        
        if (ser == 0) p += sprintf(p, "-");
        else if (ser < 0) p += sprintf(p, "n%d", -ser);
        else p += sprintf(p, "%d", ser);
    }
    *p = '\0';
}

static inline void unpack_taxid_array(const uint32_t* tax_ptr, char* out_buf) {
    uint32_t count = tax_ptr[0];
    char* p = out_buf;
    for (uint32_t k = 0; k < count; k++) {
        if (k > 0) {
            *p++ = ','; *p++ = ' ';
        }
        p += sprintf(p, "%u", tax_ptr[k+1]);
    }
    *p = '\0';
}

static inline void decode_coords_simd(const bpdb_residue_t* res, const bpdb_atom_t* atom, 
                                      float* fx, float* fy, float* fz) {
#if defined(__ARM_NEON) || defined(__aarch64__)
    uint32x4_t r_vec = vld1q_u32((const uint32_t*)res);
    uint32x4_t a_vec = vld1q_u32((const uint32_t*)atom);
    uint32x4_t mask_high = vdupq_n_u32(0xFFFF0000);
    uint32x4_t mask_low  = vdupq_n_u32(0x0000FFFF);
    int32x4_t base = vreinterpretq_s32_u32(vandq_u32(r_vec, mask_high));
    int32x4_t diff = vreinterpretq_s32_u32(vandq_u32(a_vec, mask_low));
    int32x4_t coords = vaddq_s32(base, diff);
    float32x4_t fcoords = vmulq_n_f32(vcvtq_f32_s32(coords), 0.001f);
    *fx = vgetq_lane_f32(fcoords, 1);
    *fy = vgetq_lane_f32(fcoords, 2);
    *fz = vgetq_lane_f32(fcoords, 3);
#elif defined(__SSE2__) || defined(__x86_64__)
    __m128i r_vec = _mm_loadu_si128((const __m128i*)res);
    __m128i a_vec = _mm_loadu_si128((const __m128i*)atom);
    __m128i mask_high = _mm_set1_epi32(0xFFFF0000);
    __m128i mask_low  = _mm_set1_epi32(0x0000FFFF);
    __m128i base = _mm_and_si128(r_vec, mask_high);
    __m128i diff = _mm_and_si128(a_vec, mask_low);
    __m128i coords = _mm_add_epi32(base, diff);
    __m128 fcoords = _mm_cvtepi32_ps(coords);
    fcoords = _mm_mul_ps(fcoords, _mm_set1_ps(0.001f));
    float c[4];
    _mm_storeu_ps(c, fcoords);
    *fx = c[1]; *fy = c[2]; *fz = c[3];
#else
    int32_t x = (((int32_t)res->base_x) << 16) | atom->dx;
    int32_t y = (((int32_t)res->base_y) << 16) | atom->dy;
    int32_t z = (((int32_t)res->base_z) << 16) | atom->dz;
    *fx = x * 0.001f; *fy = y * 0.001f; *fz = z * 0.001f;
#endif
}

static inline void decode_coords_int(const bpdb_residue_t* res, const bpdb_atom_t* atom, 
                                     int32_t* ix, int32_t* iy, int32_t* iz) {
    *ix = (((int32_t)res->base_x) << 16) | atom->dx;
    *iy = (((int32_t)res->base_y) << 16) | atom->dy;
    *iz = (((int32_t)res->base_z) << 16) | atom->dz;
}
#endif
