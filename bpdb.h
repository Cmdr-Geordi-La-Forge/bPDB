#ifndef BPDB_H
#define BPDB_H

#include <stdint.h>

#if defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
#elif defined(__SSE2__) || defined(__x86_64__)
    #include <emmintrin.h>
#endif

#define DIR_HEADER 0
#define DIR_OBSLTE 1
#define DIR_TITLE  2
#define DIR_SPLIT  3
#define DIR_COMPND 4
#define DIR_SOURCE 5
#define DIR_KEYWDS 6
#define DIR_EXPDTA 7
#define DIR_AUTHOR 8
#define DIR_REVDAT 9
#define DIR_JRNL   10
#define DIR_CONECT 30
#define DIR_COORD  31

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
    bpdb_date_t modDate;
    char modId[4];
    uint16_t modNum;
    uint16_t modType;
} bpdb_revdat_t;

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

#pragma pack(pop)

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
