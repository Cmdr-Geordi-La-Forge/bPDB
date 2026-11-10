#include "bpdb.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static const char* COMPND_KEYS[] = {"", "MOL_ID", "MOLECULE", "CHAIN", "FRAGMENT", "SYNONYM", "EC", "ENGINEERED", "MUTATION", "OTHER_DETAILS"};
static const char* SOURCE_KEYS[] = {"", "MOL_ID", "SYNTHETIC", "FRAGMENT", "ORGANISM_SCIENTIFIC", "ORGANISM_COMMON", "ORGANISM_TAXID",
                             "STRAIN", "VARIANT", "CELL_LINE", "ATCC", "ORGAN", "TISSUE", "CELL", "ORGANELLE", "SECRETION",
                             "CELLULAR_LOCATION", "PLASMID", "GENE", "EXPRESSION_SYSTEM", "EXPRESSION_SYSTEM_COMMON",
                             "EXPRESSION_SYSTEM_TAXID", "EXPRESSION_SYSTEM_STRAIN", "EXPRESSION_SYSTEM_VARIANT",
                             "EXPRESSION_SYSTEM_CELL_LINE", "EXPRESSION_SYSTEM_ATCC_NUMBER", "EXPRESSION_SYSTEM_ORGAN",
                             "EXPRESSION_SYSTEM_TISSUE", "EXPRESSION_SYSTEM_CELL", "EXPRESSION_SYSTEM_ORGANELLE",
                             "EXPRESSION_SYSTEM_CELLULAR_LOCATION", "EXPRESSION_SYSTEM_VECTOR_TYPE", "EXPRESSION_SYSTEM_VECTOR",
                             "EXPRESSION_SYSTEM_PLASMID", "EXPRESSION_SYSTEM_GENE", "OTHER_DETAILS"};

static const char AFDB_REMARK_BLOCK[] =
  "REMARK   1                                                                      \n"
  "REMARK   1 REFERENCE 1                                                          \n"
  "REMARK   1  AUTH   JOHN JUMPER, RICHARD EVANS, ALEXANDER PRITZEL, TIM GREEN,    \n"
  "REMARK   1  AUTH 2 MICHAEL FIGURNOV, OLAF RONNEBERGER, KATHRYN TUNYASUVUNAKOOL, \n"
  "REMARK   1  AUTH 3 RUSS BATES, AUGUSTIN ZIDEK, ANNA POTAPENKO, ALEX BRIDGLAND,  \n"
  "REMARK   1  AUTH 4 CLEMENS MEYER, SIMON A A KOHL, ANDREW J BALLARD,             \n"
  "REMARK   1  AUTH 5 ANDREW COWIE, BERNARDINO ROMERA-PAREDES, STANISLAV NIKOLOV,  \n"
  "REMARK   1  AUTH 6 RISHUB JAIN, JONAS ADLER, TREVOR BACK, STIG PETERSEN,        \n"
  "REMARK   1  AUTH 7 DAVID REIMAN, ELLEN CLANCY, MICHAL ZIELINSKI,                \n"
  "REMARK   1  AUTH 8 MARTIN STEINEGGER, MICHALINA PACHOLSKA, TAMAS BERGHAMMER,    \n"
  "REMARK   1  AUTH 9 DAVID SILVER, ORIOL VINYALS, ANDREW W SENIOR,                \n"
  "REMARK   1  AUTH10 KORAY KAVUKCUOGLU, PUSHMEET KOHLI, DEMIS HASSABIS            \n"
  "REMARK   1  TITL   HIGHLY ACCURATE PROTEIN STRUCTURE PREDICTION WITH ALPHAFOLD  \n"
  "REMARK   1  REF    NATURE                        V. 596   583 2021              \n"
  "REMARK   1  REFN                   ISSN 0028-0836                               \n"
  "REMARK   1  PMID   34265844                                                     \n"
  "REMARK   1  DOI    10.1038/s41586-021-03819-2                                   \n"
  "REMARK   1                                                                      \n"
  "REMARK   1 DISCLAIMERS                                                          \n"
  "REMARK   1 ALPHAFOLD DATA, COPYRIGHT (2021) DEEPMIND TECHNOLOGIES LIMITED. THE  \n"
  "REMARK   1 INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY AND CAUTION SHOULD\n"
  "REMARK   1 BE EXERCISED IN ITS USE. IT IS PROVIDED \"AS-IS\" WITHOUT ANY WARRANTY \n"
  "REMARK   1 OF ANY KIND, WHETHER EXPRESSED OR IMPLIED. NO WARRANTY IS GIVEN THAT \n"
  "REMARK   1 USE OF THE INFORMATION SHALL NOT INFRINGE THE RIGHTS OF ANY THIRD    \n"
  "REMARK   1 PARTY. THE INFORMATION IS NOT INTENDED TO BE A SUBSTITUTE FOR        \n"
  "REMARK   1 PROFESSIONAL MEDICAL ADVICE, DIAGNOSIS, OR TREATMENT, AND DOES NOT   \n"
  "REMARK   1 CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. IT IS AVAILABLE FOR \n"
  "REMARK   1 ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE.           \n";

static const char AFDB_CRYST_BLOCK[] =
  "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n"
  "ORIGX1      1.000000  0.000000  0.000000        0.00000                         \n"
  "ORIGX2      0.000000  1.000000  0.000000        0.00000                         \n"
  "ORIGX3      0.000000  0.000000  1.000000        0.00000                         \n"
  "SCALE1      1.000000  0.000000  0.000000        0.00000                         \n"
  "SCALE2      0.000000  1.000000  0.000000        0.00000                         \n"
  "SCALE3      0.000000  0.000000  1.000000        0.00000                         \n";

#if defined(__AVX2__) || defined(__SSE4_1__)
#include <immintrin.h>
typedef __m128i accum_t;
#define GET_SIMD_B(acc) _mm_extract_epi32((acc), 0)
#define GET_SIMD_X(acc) _mm_extract_epi32((acc), 1)
#define GET_SIMD_Y(acc) _mm_extract_epi32((acc), 2)
#define GET_SIMD_Z(acc) _mm_extract_epi32((acc), 3)
#elif defined(__aarch64__) || defined(__ARM_NEON)
#include <arm_neon.h>
typedef int32x4_t accum_t;
#define GET_SIMD_B(acc) vgetq_lane_s32((acc), 0)
#define GET_SIMD_X(acc) vgetq_lane_s32((acc), 1)
#define GET_SIMD_Y(acc) vgetq_lane_s32((acc), 2)
#define GET_SIMD_Z(acc) vgetq_lane_s32((acc), 3)
#else
typedef struct { int32_t b, x, y, z; } accum_t;
#define GET_SIMD_B(acc) ((acc).b)
#define GET_SIMD_X(acc) ((acc).x)
#define GET_SIMD_Y(acc) ((acc).y)
#define GET_SIMD_Z(acc) ((acc).z)
#endif

// Unifies memory reads and creates a 100% symmetric SIMD phase
#if defined(__AVX2__) || defined(__SSE4_1__)
#define SIMD_DECODE_PHASE1(ptr, out_accum) do {                         \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      __m128i msb = _mm_cvtepi16_epi32(_mm_cvtsi64_si128(ctrl));        \
      out_accum = _mm_slli_epi32(_mm_srai_epi32(msb, 1), 13);           \
      __m128i lsb = _mm_cvtepi16_epi32(_mm_cvtsi64_si128(*(uint64_t*)(ptr + 8))); \
      out_accum = _mm_add_epi32(out_accum, _mm_srai_epi32(lsb, 3));     \
      ptr += 16;                                                        \
    } else {                                                            \
      __m128i loaded = _mm_cvtsi64_si128(ctrl);                         \
      __m128i deltas = _mm_srai_epi16(loaded, 3);                       \
      out_accum = _mm_add_epi32(out_accum, _mm_cvtepi16_epi32(deltas)); \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#elif defined(__aarch64__) || defined(__ARM_NEON)
#define SIMD_DECODE_PHASE1(ptr, out_accum) do {                         \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      int16x4_t msb = vreinterpret_s16_u8(vcreate_u8(ctrl));            \
      int32x4_t msb_32 = vshlq_n_s32(vmovl_s16(vshr_n_s16(msb, 1)), 13);\
      int16x4_t lsb = vreinterpret_s16_u8(vcreate_u8(*(uint64_t*)(ptr + 8))); \
      out_accum = vaddq_s32(msb_32, vmovl_s16(vshr_n_s16(lsb, 3)));     \
      ptr += 16;                                                        \
    } else {                                                            \
      int16x4_t loaded = vreinterpret_s16_u8(vcreate_u8(ctrl));         \
      int16x4_t deltas = vshr_n_s16(loaded, 3);                         \
      out_accum = vaddq_s32(out_accum, vmovl_s16(deltas));              \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#else
#define SIMD_DECODE_PHASE1(ptr, out_accum) do {                         \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      out_accum.b = (((int16_t)(uint16_t)ctrl) >> 1) << 13;             \
      out_accum.x = (((int16_t)(uint16_t)(ctrl >> 16)) >> 1) << 13;     \
      out_accum.y = (((int16_t)(uint16_t)(ctrl >> 32)) >> 1) << 13;     \
      out_accum.z = (((int16_t)(uint16_t)(ctrl >> 48)) >> 1) << 13;     \
      uint64_t lsb_ctrl = *(uint64_t*)(ptr + 8);                        \
      out_accum.b += ((int16_t)(uint16_t)lsb_ctrl) >> 3;                \
      out_accum.x += ((int16_t)(uint16_t)(lsb_ctrl >> 16)) >> 3;        \
      out_accum.y += ((int16_t)(uint16_t)(lsb_ctrl >> 32)) >> 3;        \
      out_accum.z += ((int16_t)(uint16_t)(lsb_ctrl >> 48)) >> 3;        \
      ptr += 16;                                                        \
    } else {                                                            \
      out_accum.b += ((int16_t)(uint16_t)ctrl) >> 3;                    \
      out_accum.x += ((int16_t)(uint16_t)(ctrl >> 16)) >> 3;            \
      out_accum.y += ((int16_t)(uint16_t)(ctrl >> 32)) >> 3;            \
      out_accum.z += ((int16_t)(uint16_t)(ctrl >> 48)) >> 3;            \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#endif

#if defined(__AVX2__) || defined(__SSE4_1__)
#define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, alt_out) do {   \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      __m128i msb = _mm_cvtepi16_epi32(_mm_cvtsi64_si128(ctrl));        \
      out_accum = _mm_slli_epi32(_mm_srai_epi32(msb, 1), 13);           \
      uint64_t lsb_ctrl = *(uint64_t*)(ptr + 8);                        \
      uint16_t b_lsb = (uint16_t)lsb_ctrl;                              \
      uint16_t x_lsb = (uint16_t)(lsb_ctrl >> 16);                      \
      uint16_t y_lsb = (uint16_t)(lsb_ctrl >> 32);                      \
      uint16_t z_lsb = (uint16_t)(lsb_ctrl >> 48);                      \
      alt_out = (b_lsb >> 1) & 1;                                       \
      c_occ = ((b_lsb & 4) << 4) | ((x_lsb & 7) << 3) | (y_lsb & 7);    \
      c_chg = (z_lsb & 7) - 3;                                          \
      __m128i lsb = _mm_cvtepi16_epi32(_mm_cvtsi64_si128(lsb_ctrl));    \
      out_accum = _mm_add_epi32(out_accum, _mm_srai_epi32(lsb, 3));     \
      ptr += 16;                                                        \
    } else {                                                            \
      uint16_t b_ctrl = (uint16_t)ctrl;                                 \
      uint16_t x_ctrl = (uint16_t)(ctrl >> 16);                         \
      uint16_t y_ctrl = (uint16_t)(ctrl >> 32);                         \
      uint16_t z_ctrl = (uint16_t)(ctrl >> 48);                         \
      alt_out = (b_ctrl >> 1) & 1;                                      \
      c_occ = ((b_ctrl & 4) << 4) | ((x_ctrl & 7) << 3) | (y_ctrl & 7); \
      c_chg = (z_ctrl & 7) - 3;                                         \
      __m128i loaded = _mm_cvtsi64_si128(ctrl);                         \
      __m128i deltas = _mm_srai_epi16(loaded, 3);                       \
      out_accum = _mm_add_epi32(out_accum, _mm_cvtepi16_epi32(deltas)); \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#elif defined(__aarch64__) || defined(__ARM_NEON)
#define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, alt_out) do {   \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      int16x4_t msb = vreinterpret_s16_u8(vcreate_u8(ctrl));            \
      int32x4_t msb_32 = vshlq_n_s32(vmovl_s16(vshr_n_s16(msb, 1)), 13);\
      uint64_t lsb_ctrl = *(uint64_t*)(ptr + 8);                        \
      uint16_t b_lsb = (uint16_t)lsb_ctrl;                              \
      uint16_t x_lsb = (uint16_t)(lsb_ctrl >> 16);                      \
      uint16_t y_lsb = (uint16_t)(lsb_ctrl >> 32);                      \
      uint16_t z_lsb = (uint16_t)(lsb_ctrl >> 48);                      \
      alt_out = (b_lsb >> 1) & 1;                                       \
      c_occ = ((b_lsb & 4) << 4) | ((x_lsb & 7) << 3) | (y_lsb & 7);    \
      c_chg = (z_lsb & 7) - 3;                                          \
      int16x4_t lsb = vreinterpret_s16_u8(vcreate_u8(lsb_ctrl));        \
      out_accum = vaddq_s32(msb_32, vmovl_s16(vshr_n_s16(lsb, 3)));     \
      ptr += 16;                                                        \
    } else {                                                            \
      uint16_t b_ctrl = (uint16_t)ctrl;                                 \
      uint16_t x_ctrl = (uint16_t)(ctrl >> 16);                         \
      uint16_t y_ctrl = (uint16_t)(ctrl >> 32);                         \
      uint16_t z_ctrl = (uint16_t)(ctrl >> 48);                         \
      alt_out = (b_ctrl >> 1) & 1;                                      \
      c_occ = ((b_ctrl & 4) << 4) | ((x_ctrl & 7) << 3) | (y_ctrl & 7); \
      c_chg = (z_ctrl & 7) - 3;                                         \
      int16x4_t loaded = vreinterpret_s16_u8(vcreate_u8(ctrl));         \
      int16x4_t deltas = vshr_n_s16(loaded, 3);                         \
      out_accum = vaddq_s32(out_accum, vmovl_s16(deltas));              \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#else
#define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, alt_out) do {   \
    uint64_t ctrl = *(uint64_t*)(ptr);                                  \
    if ((uint32_t)ctrl & 0x01) {                                        \
      out_accum.b = (((int16_t)(uint16_t)ctrl) >> 1) << 13;             \
      out_accum.x = (((int16_t)(uint16_t)(ctrl >> 16)) >> 1) << 13;     \
      out_accum.y = (((int16_t)(uint16_t)(ctrl >> 32)) >> 1) << 13;     \
      out_accum.z = (((int16_t)(uint16_t)(ctrl >> 48)) >> 1) << 13;     \
      uint64_t lsb_ctrl = *(uint64_t*)(ptr + 8);                        \
      uint16_t b_lsb = (uint16_t)lsb_ctrl;                              \
      uint16_t x_lsb = (uint16_t)(lsb_ctrl >> 16);                      \
      uint16_t y_lsb = (uint16_t)(lsb_ctrl >> 32);                      \
      uint16_t z_lsb = (uint16_t)(lsb_ctrl >> 48);                      \
      alt_out = (b_lsb >> 1) & 1;                                       \
      c_occ = ((b_lsb & 4) << 4) | ((x_lsb & 7) << 3) | (y_lsb & 7);    \
      c_chg = (z_lsb & 7) - 3;                                          \
      out_accum.b += ((int16_t)(uint16_t)lsb_ctrl) >> 3;                \
      out_accum.x += ((int16_t)(uint16_t)(lsb_ctrl >> 16)) >> 3;        \
      out_accum.y += ((int16_t)(uint16_t)(lsb_ctrl >> 32)) >> 3;        \
      out_accum.z += ((int16_t)(uint16_t)(lsb_ctrl >> 48)) >> 3;        \
      ptr += 16;                                                        \
    } else {                                                            \
      uint16_t b_ctrl = (uint16_t)ctrl;                                 \
      uint16_t x_ctrl = (uint16_t)(ctrl >> 16);                         \
      uint16_t y_ctrl = (uint16_t)(ctrl >> 32);                         \
      uint16_t z_ctrl = (uint16_t)(ctrl >> 48);                         \
      alt_out = (b_ctrl >> 1) & 1;                                      \
      c_occ = ((b_ctrl & 4) << 4) | ((x_ctrl & 7) << 3) | (y_ctrl & 7); \
      c_chg = (z_ctrl & 7) - 3;                                         \
      out_accum.b += ((int16_t)b_ctrl) >> 3;                            \
      out_accum.x += ((int16_t)x_ctrl) >> 3;                            \
      out_accum.y += ((int16_t)y_ctrl) >> 3;                            \
      out_accum.z += ((int16_t)z_ctrl) >> 3;                            \
      ptr += 8;                                                         \
    }                                                                   \
  } while(0)
#endif

volatile int32_t benchmark_sink = 0;

static inline void unpack_ec_number(uint32_t packed, char* out_str) {
  int class_val = (packed >> 29) & 0x07;
  int sub = (packed >> 24) & 0x1F;
  int8_t subsub = (int8_t)((packed >> 16) & 0xFF);
  int16_t serial = (int16_t)(packed & 0xFFFF);

  char s3[16], s4[16];
  if (subsub == 0) strcpy(s3, "-"); else if (subsub < 0) snprintf(s3, sizeof(s3), "n%d", -subsub); else snprintf(s3, sizeof(s3), "%d", subsub);
  if (serial == 0) strcpy(s4, "-"); else if (serial < 0) snprintf(s4, sizeof(s4), "n%d", -serial); else snprintf(s4, sizeof(s4), "%d", serial);

  snprintf(out_str, 32, "%d.%d.%s.%s", class_val, sub, s3, s4);
}

void write_tlv_to_pdb(char** w_ptr, const uint8_t* tlv_buf, size_t tlv_len, int is_source, TaxDict* taxdict, int af_mode) {
  size_t ptr = 0; char tbuf[128];
  const char* pfx = is_source ? "SOURCE" : "COMPND";
  const char** keys = is_source ? SOURCE_KEYS : COMPND_KEYS;
  int cont_num = 1;

  while (ptr < tlv_len) {
    uint16_t mlen = *(uint16_t*)(tlv_buf + ptr); ptr += 2;
    size_t m_end = ptr + mlen;

    while (ptr < m_end) {
      uint8_t token_id = tlv_buf[ptr++];
      uint8_t vlen = tlv_buf[ptr++];
      int is_last = (ptr + vlen == tlv_len);
      const char* semi = is_last ? "" : ";";

      if (!is_source && token_id == COMPND_EC) {
        uint32_t packed = *(uint32_t*)(tlv_buf + ptr);
        char ec_str[32]; unpack_ec_number(packed, ec_str);

        char cont_str[8]; if (cont_num == 1) strcpy(cont_str, "    "); else snprintf(cont_str, 8, "%4d ", cont_num);
        int len = snprintf(tbuf, 128, "%-6s%s%s: %s%s", pfx, cont_str, keys[token_id], ec_str, semi);
        memcpy(*w_ptr, tbuf, len); *w_ptr += 81;
        ptr += 4;
        cont_num++;
      } else if (is_source && (token_id == SOURCE_ORGANISM_TAXID || token_id == SOURCE_EXPRESSION_SYSTEM_TAXID)) {
        uint32_t taxid = *(uint32_t*)(tlv_buf + ptr); ptr += 4;

        if (af_mode) {
          const char* sci = get_tax_name(taxdict, taxid, 0);
          if (sci && strlen(sci) > 0) {
            char cont_str[8]; if (cont_num == 1) strcpy(cont_str, "    "); else snprintf(cont_str, 8, "%4d ", cont_num);
            int sci_tok = (token_id == SOURCE_ORGANISM_TAXID) ? SOURCE_ORGANISM_SCIENTIFIC : SOURCE_EXPRESSION_SYSTEM;
            int len = snprintf(tbuf, 128, "%-6s%s%s: %s;", pfx, cont_str, keys[sci_tok], sci);
            memcpy(*w_ptr, tbuf, len); *w_ptr += 81;
            cont_num++;
          }
        }

        char cont_str[8]; if (cont_num == 1) strcpy(cont_str, "    "); else snprintf(cont_str, 8, "%4d ", cont_num);
        int len = snprintf(tbuf, 128, "%-6s%s%s: %u%s", pfx, cont_str, keys[token_id], taxid, semi);
        memcpy(*w_ptr, tbuf, len); *w_ptr += 81;
        cont_num++;
      } else {
        char cont_str[8]; if (cont_num == 1) strcpy(cont_str, "    "); else snprintf(cont_str, 8, "%4d ", cont_num);
        int len = snprintf(tbuf, 128, "%-6s%s%s: %.*s%s", pfx, cont_str, keys[token_id], vlen, (char*)(tlv_buf + ptr), semi);
        memcpy(*w_ptr, tbuf, len); *w_ptr += 81;
        ptr += vlen;
        cont_num++;
      }
    }
  }
}

void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out, int is_first_model, int is_last_model) {
  if (size < 64) return;

  size_t offset = 64;
  uint64_t meta_header = *(uint64_t*)(buffer + offset); offset += 8;
  uint64_t mask = meta_header >> 16;
  int af_mode = (mask & (1ULL << META_AF_MODE)) ? 1 : 0;
  int has_hydrogens = (mask & (1ULL << META_HYDROGENS)) ? 1 : 0;

  size_t dict_scan_off = offset;
  if (mask & (1ULL << META_RESOLUTION)) dict_scan_off += 8;
  if (mask & (1ULL << META_TITLE))   { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_EXPDTA))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_PROGRAM)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_HET_DICT)) {
    uint64_t dict_len = *(uint64_t*)(buffer + dict_scan_off);
    dict_scan_off += 8 + dict_len + ((8 - (dict_len % 8)) % 8);
  }
  if (mask & (1ULL << META_HELIX))   { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_SHEET))   { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_SSBOND))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_CISPEP))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_SITE))    { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_COMPND))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_SOURCE))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_DBREF))   { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << META_SEQRES))  { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }

  offset = dict_scan_off;

  if (!out) {
    int32_t local_sink = 0;
    while (offset < size) {
      uint32_t r_a = *(uint32_t*)(buffer + offset);
      uint32_t num_residues = r_a >> 18;
      uint32_t num_atoms = r_a & 0x3FFFF;

      if (num_residues == 0x3FFE) { offset += 8; continue; }
      if (num_residues == 0x3FFF) { offset += 8; break; }

      offset += 8;

      const uint8_t* seq_ptr = buffer + offset;
      int seq_bytes = num_residues;
      const uint8_t* search_ptr = seq_ptr;
      const uint8_t* end_ptr = seq_ptr + seq_bytes;

      while (search_ptr < end_ptr) {
        const uint8_t* found = (const uint8_t*)memchr(search_ptr, 0x7C, end_ptr - search_ptr);
        if (!found) break;
        seq_bytes++; end_ptr++;
        search_ptr = found + 2;
      }

      offset += seq_bytes + ((8 - (seq_bytes & 7)) & 7);
      const uint8_t* atom_ptr = buffer + offset;

      accum_t accum;
#if defined(__AVX2__) || defined(__SSE4_1__)
      accum = _mm_setzero_si128();
#elif defined(__aarch64__) || defined(__ARM_NEON)
      accum = vdupq_n_s32(0);
#else
      accum.x = 0; accum.y = 0; accum.z = 0; accum.b = 0;
#endif

      for (uint32_t a = 0; a < num_atoms; a++) {
        SIMD_DECODE_PHASE1(atom_ptr, accum);
      }
      local_sink += GET_SIMD_X(accum); offset = (size_t)(atom_ptr - buffer);
    }
    benchmark_sink += local_sink;
    return;
  }

  ResidueDict* UDICT = (ResidueDict*)calloc(256, sizeof(ResidueDict));
  for (size_t i = 0; i < NUM_DICT_RES; i++) {
    UDICT[DICT[i].type_id] = DICT[i];
  }

  int day = buffer[5], month = buffer[6], year = 1973 + buffer[7];
  char pdbID[13] = {0}; memcpy(pdbID, buffer + 8, 12); trim_string(pdbID);
  uint16_t total_models_minus_1 = *(uint16_t*)(buffer + 20);
  uint16_t current_model_minus_1 = *(uint16_t*)(buffer + 22);
  char classification[41] = {0}; memcpy(classification, buffer + 24, 40); trim_string(classification);

  offset = 64 + 8;
  float resolution = -1.0f;
  const char* p_title   = "";      uint64_t title_len  = 0;
  const char* p_expdta  = "";      uint64_t expdta_len = 0;
  const char* p_program = "";

  const uint8_t* p_ssbonds = NULL; uint64_t num_ssbonds = 0;
  const uint8_t* p_cispeps = NULL; uint64_t num_cispeps = 0;
  const uint8_t* p_helices = NULL; uint64_t num_helices = 0;
  const uint8_t* p_sheets  = NULL; uint64_t num_sheets  = 0;
  const uint8_t* p_sites   = NULL; uint64_t len_sites   = 0;

  const uint8_t* p_compnd  = NULL; uint64_t compnd_len  = 0;
  const uint8_t* p_source  = NULL; uint64_t source_len  = 0;
  const uint8_t* p_dbrefs  = NULL; uint64_t num_dbrefs  = 0;
  const uint8_t* p_seqres  = NULL; uint64_t seqres_len  = 0;

  if (mask & (1ULL << META_RESOLUTION)) { resolution = *(float*)(buffer + offset); offset += 8; }
  if (mask & (1ULL << META_TITLE)) { title_len = *(uint64_t*)(buffer + offset); offset += 8; p_title = (const char*)(buffer + offset); offset += title_len; offset += (8 - (title_len % 8)) % 8; }
  if (mask & (1ULL << META_EXPDTA)) { expdta_len = *(uint64_t*)(buffer + offset); offset += 8; p_expdta = (const char*)(buffer + offset); offset += expdta_len; offset += (8 - (expdta_len % 8)) % 8; }
  if (mask & (1ULL << META_PROGRAM)) { uint64_t p_len = *(uint64_t*)(buffer + offset); offset += 8; p_program = (const char*)(buffer + offset); offset += p_len; offset += (8 - (p_len % 8)) % 8; }
  if (mask & (1ULL << META_HET_DICT)) {
    uint64_t dict_len = *(uint64_t*)(buffer + offset); offset += 8;
    const uint8_t* p_het_dict = buffer + offset;
    uint8_t custom_count = p_het_dict[0]; size_t ptr = 1;
    for (int c = 0; c < custom_count; c++) {
      int u_idx = 0x81 + c; memcpy(UDICT[u_idx].resName, p_het_dict + ptr, 3); UDICT[u_idx].resName[3] = '\0'; ptr += 4;
      UDICT[u_idx].h_end_norm = p_het_dict[ptr]; UDICT[u_idx].h_end_cterm = p_het_dict[ptr]; ptr++;
      for (int a = 0; a < UDICT[u_idx].h_end_norm; a++) {
        memcpy(UDICT[u_idx].h_atoms[a], p_het_dict + ptr, 4); UDICT[u_idx].h_atoms[a][4] = '\0'; ptr += 4;
        memcpy(UDICT[u_idx].h_elements[a], p_het_dict + ptr, 2); UDICT[u_idx].h_elements[a][2] = '\0'; ptr += 2;
      }
    }
    offset += dict_len; offset += (8 - (dict_len % 8)) % 8;
  }
  if (mask & (1ULL << META_HELIX)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_helices = buffer + offset; num_helices = L / sizeof(MetaHelix); offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_SHEET)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_sheets = buffer + offset; num_sheets = L / sizeof(MetaSheet); offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_SSBOND)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_ssbonds = buffer + offset; num_ssbonds = L / sizeof(MetaSSBond); offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_CISPEP)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_cispeps = buffer + offset; num_cispeps = L / sizeof(MetaCispep); offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_SITE)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_sites = buffer + offset; len_sites = L; offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_COMPND)) { compnd_len = *(uint64_t*)(buffer + offset); offset += 8; p_compnd = buffer + offset; offset += compnd_len; offset += (8 - (compnd_len % 8)) % 8; }
  if (mask & (1ULL << META_SOURCE)) { source_len = *(uint64_t*)(buffer + offset); offset += 8; p_source = buffer + offset; offset += source_len; offset += (8 - (source_len % 8)) % 8; }
  if (mask & (1ULL << META_DBREF)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_dbrefs = buffer + offset; num_dbrefs = L / sizeof(MetaDBRef); offset += L; offset += (8 - (L % 8)) % 8; }
  if (mask & (1ULL << META_SEQRES)) { seqres_len = *(uint64_t*)(buffer + offset); offset += 8; p_seqres = buffer + offset; offset += seqres_len; offset += (8 - (seqres_len % 8)) % 8; }

  size_t out_lines = 0;

  if (is_first_model) {
    out_lines += 1;
    if (compnd_len > 0) {
      size_t p = 0; while(p < compnd_len) { uint16_t mlen = *(uint16_t*)(p_compnd+p); p += 2; size_t mend = p + mlen; while(p < mend) { p += 2 + p_compnd[p+1]; out_lines++; } }
    }
    if (source_len > 0) {
      TaxDict t = load_taxdict();
      size_t p = 0; while(p < source_len) {
        uint16_t mlen = *(uint16_t*)(p_source+p); p += 2; size_t mend = p + mlen;
        while(p < mend) {
          uint8_t tok = p_source[p]; uint8_t vlen = p_source[p+1]; p += 2;
          if (tok == SOURCE_ORGANISM_TAXID || tok == SOURCE_EXPRESSION_SYSTEM_TAXID) {
            uint32_t taxid = *(uint32_t*)(p_source+p); out_lines++; // +1 for the TAXID line itself
            if (af_mode) {
              const char* sci = get_tax_name(&t, taxid, 0);
              if (sci && strlen(sci) > 0) out_lines++; // +1 dynamically generated Scientific Name
            }
            p += 4;
          } else { p += vlen; out_lines++; }
        }
      }
    }
    if (title_len > 0) out_lines += (title_len + 69) / 70;
    if (mask & (1ULL << META_AFDB)) {
      out_lines += 35;
      uint32_t first_r_a = *(uint32_t*)(buffer + offset);
      uint32_t af_res_count = first_r_a >> 18;
      out_lines += (af_res_count + 12) / 13;
    }
    out_lines += num_dbrefs;
    if (seqres_len > 0) {
      size_t scan = 0;
      while (scan < seqres_len) {
        uint16_t num_res = *(uint16_t*)(p_seqres + scan + 1);
        out_lines += (num_res + 12) / 13; // 13 residues per line
        scan += 3 + num_res;
      }
    }
    if (expdta_len > 0) out_lines++;
    if (resolution > 0.0f) out_lines += 2;
    if (strlen(p_program) > 0) out_lines += 3;
    out_lines += num_helices + num_sheets + num_ssbonds + num_cispeps;
    if (len_sites > 0) {
      size_t s_ptr = 0;
      while (s_ptr < len_sites) {
        uint32_t num_res = *(uint32_t*)(p_sites + s_ptr + 4);
        out_lines += (num_res + 3) / 4;
        s_ptr += 8 + num_res * sizeof(PdbResId);
      }
    }
    out_lines += 1;
  }

  if (total_models_minus_1 > 0 || (mask & (1ULL << META_AFDB))) out_lines += 2;

  size_t scan_off = offset, conect_offset = 0;
  while (scan_off < size) {
    uint32_t r_a = *(uint32_t*)(buffer + scan_off);
    uint32_t num_residues = r_a >> 18;
    uint32_t num_atoms = r_a & 0x3FFFF;
    if (num_residues == 0x3FFE) { scan_off += 8; continue; }
    if (num_residues == 0x3FFF) { scan_off += 8; conect_offset = scan_off; break; }

    int16_t seq_raw = *(int16_t*)(buffer + scan_off + 6);
    int has_ter = seq_raw & 1;
    scan_off += 8;

    const uint8_t* scan_seq_ptr = buffer + scan_off;
    int scan_seq_bytes = num_residues;
    const uint8_t* search_ptr = scan_seq_ptr;
    const uint8_t* end_ptr = scan_seq_ptr + scan_seq_bytes;
    while (search_ptr < end_ptr) {
      const uint8_t* found = (const uint8_t*)memchr(search_ptr, 0x7C, end_ptr - search_ptr);
      if (!found) break;
      scan_seq_bytes++; end_ptr++;
      search_ptr = found + 2;
    }
    scan_off += scan_seq_bytes + ((8 - (scan_seq_bytes & 7)) & 7);

    const uint8_t* atom_ptr = buffer + scan_off;
    for (uint32_t a = 0; a < num_atoms; a++) {
      if (*(uint64_t*)(atom_ptr) != 0) out_lines++;
      if ((atom_ptr)[0] & 0x01) atom_ptr += 16; else atom_ptr += 8;
    }
    if (has_ter) out_lines++;
    scan_off = (size_t)(atom_ptr - buffer);
  }

  if (is_last_model && conect_offset > 0) {
    scan_off = conect_offset;
    while (scan_off + 4 <= size) {
      int32_t N = *(int32_t*)(buffer + scan_off);
      if (N <= 0 || N > 4 || scan_off + 8 + (N << 2) > size) break;
      out_lines++; scan_off += 8 + (N << 2);
    }
  }

  if (out_lines == 0) { free(UDICT); return; }
  size_t out_size = out_lines * 81 - (af_mode ? 1 : 0);
  char* out_data = (char*)malloc(out_size); memset(out_data, ' ', out_size);
  for (size_t i = 0; i < out_lines - (af_mode ? 1 : 0); i++) out_data[i * 81 + 80] = '\n';

  char* w_ptr = out_data; char tbuf[128];

  if (is_first_model) {
    memcpy(w_ptr, "HEADER    ", 10);
    size_t c_len = strlen(classification); if (c_len > 0) memcpy(w_ptr + 10, classification, c_len > 40 ? 40 : c_len);
    int dlen = snprintf(tbuf, 128, "%02d-%3s-%02d", day, MONTHS[month], year % 100); memcpy(w_ptr + 50, tbuf, dlen);
    size_t id_len = strlen(pdbID); if (id_len > 0 && !af_mode) memcpy(w_ptr + 62, pdbID, id_len > 18 ? 18 : id_len); w_ptr += 81;

    if (title_len > 0) {
      for (uint64_t i = 0; i < title_len; i += 70) {
        memcpy(w_ptr, "TITLE   ", 8); if (i > 0) { int ln = snprintf(tbuf, 128, "%2d", (int)(i/70)+1); memcpy(w_ptr + 8, tbuf, ln); }
        size_t copy_len = title_len - i; if (copy_len > 70) copy_len = 70;
        memcpy(w_ptr + 10, p_title + i, copy_len); w_ptr += 81;
      }
    }

    TaxDict taxdict = load_taxdict();
    if (compnd_len > 0) write_tlv_to_pdb(&w_ptr, p_compnd, compnd_len, 0, &taxdict, af_mode);
    if (source_len > 0) write_tlv_to_pdb(&w_ptr, p_source, source_len, 1, &taxdict, af_mode);
    if (mask & (1ULL << META_AFDB)) {
      memcpy(w_ptr, AFDB_REMARK_BLOCK, 28 * 81); w_ptr += 28 * 81;
    }

    for (uint32_t i=0; i<num_dbrefs; i++) {
      MetaDBRef* d = (MetaDBRef*)(p_dbrefs + i * sizeof(MetaDBRef));
      int len = snprintf(tbuf, 128, "DBREF  %4s %c %4d%c %4d%c %-6s %-8s %-12s %5d%c %5d%c",
                         (af_mode ? "XXXX" : pdbID), d->chain, d->seq_begin, d->iCode_begin, d->seq_end, d->iCode_end,
                         d->db_name, d->db_accession, d->db_id_code, d->db_seq_begin, d->db_iCode_begin, d->db_seq_end, d->db_iCode_end);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    if (mask & (1ULL << META_AFDB)) {
      uint32_t first_r_a = *(uint32_t*)(buffer + offset);
      uint32_t num_res = first_r_a >> 18;
      const uint8_t* af_seq_ptr = buffer + offset + 8;

      int row = 1;
      for (uint32_t r = 0; r < num_res; r += 13) {
        int printed = snprintf(tbuf, 128, "SEQRES %3d A %4d ", row++, num_res);
        for (int t = 0; t < 13 && r + t < num_res; t++) {
          uint8_t type_id = af_seq_ptr[r + t];
          const char* rname = (type_id > 0) ? UDICT[type_id].resName : "UNK";
          printed += snprintf(tbuf + printed, 128 - printed, " %3s", rname);
        }
        memcpy(w_ptr, tbuf, printed); w_ptr += 81;
      }
    }

    if (seqres_len > 0) {
      size_t scan = 0;
      while (scan < seqres_len) {
        char chain = (char)p_seqres[scan];
        uint16_t num_res = *(uint16_t*)(p_seqres + scan + 1);
        scan += 3;

        int row = 1;
        for (uint16_t r = 0; r < num_res; r += 13) {
          int printed = snprintf(tbuf, 128, "SEQRES %3d %c %4d ", row++, chain, num_res);
          for (int t = 0; t < 13 && r + t < num_res; t++) {
            uint8_t type_id = p_seqres[scan++];
            const char* rname = (type_id > 0) ? UDICT[type_id].resName : "UNK";
            printed += snprintf(tbuf + printed, 128 - printed, " %3s", rname);
          }
          memcpy(w_ptr, tbuf, printed); w_ptr += 81;
        }
      }
    }
    if (expdta_len > 0) { memcpy(w_ptr, "EXPDTA    ", 10); memcpy(w_ptr + 10, p_expdta, expdta_len > 70 ? 70 : expdta_len); w_ptr += 81; }
    if (resolution > 0.0f) { memcpy(w_ptr, "REMARK   2", 10); w_ptr += 81; int len = snprintf(tbuf, 128, "REMARK   2 RESOLUTION. %7.2f ANGSTROMS.", resolution); memcpy(w_ptr, tbuf, len); w_ptr += 81; }
    if (strlen(p_program) > 0) { memcpy(w_ptr, "REMARK   3", 10); w_ptr += 81; memcpy(w_ptr, "REMARK   3 REFINEMENT.", 22); w_ptr += 81; int len = snprintf(tbuf, 128, "REMARK   3   PROGRAM     : %s", p_program); memcpy(w_ptr, tbuf, len); w_ptr += 81; }

    for (uint32_t i=0; i<num_helices; i++) {
      MetaHelix* h = (MetaHelix*)(p_helices + i * sizeof(MetaHelix));
      int len = snprintf(tbuf, 128, "HELIX  %3d %3s %-3s %c %4d%c %-3s %c %4d%c%2d%-30s %5d",
                         i+1, h->helixID[0] ? h->helixID : "   ", h->start.resName, h->start.chain, (int16_t)h->start.resSeq, h->start.iCode == 0 ? ' ' : h->start.iCode, h->end.resName, h->end.chain, (int16_t)h->end.resSeq, h->end.iCode == 0 ? ' ' : h->end.iCode, h->helix_class, h->comment[0] ? h->comment : "", h->length);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    for (uint32_t i=0; i<num_sheets; i++) {
      MetaSheet* s = (MetaSheet*)(p_sheets + i * sizeof(MetaSheet));
      int len = snprintf(tbuf, 128, "SHEET  %3d %3s%2d %-3s %c%4d%c %-3s %c%4d%c%2d",
                         s->strand_num, s->sheetID[0] ? s->sheetID : "   ", s->num_strands, s->start.resName, s->start.chain, (int16_t)s->start.resSeq, s->start.iCode == 0 ? ' ' : s->start.iCode, s->end.resName, s->end.chain, (int16_t)s->end.resSeq, s->end.iCode == 0 ? ' ' : s->end.iCode, s->sense);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    for (uint32_t i=0; i<num_ssbonds; i++) {
      MetaSSBond* b = (MetaSSBond*)(p_ssbonds + i * sizeof(MetaSSBond));
      int len = snprintf(tbuf, 128, "SSBOND %3d %-3s %c %4d%c   %-3s %c %4d%c",
                         i+1, b->res1.resName, b->res1.chain, (int16_t)b->res1.resSeq, b->res1.iCode == 0 ? ' ' : b->res1.iCode, b->res2.resName, b->res2.chain, (int16_t)b->res2.resSeq, b->res2.iCode == 0 ? ' ' : b->res2.iCode);
      if (b->sym1[0] || b->sym2[0] || b->length > 0.0f) { memset(tbuf + len, ' ', 59 - len); len = 59; if (b->length > 0.0f) len += snprintf(tbuf + len, 128 - len, "%6s %6s %5.2f", b->sym1, b->sym2, b->length); else len += snprintf(tbuf + len, 128 - len, "%6s %6s", b->sym1, b->sym2); }
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    for (uint32_t i=0; i<num_cispeps; i++) {
      MetaCispep* c = (MetaCispep*)(p_cispeps + i * sizeof(MetaCispep));
      int len = snprintf(tbuf, 128, "CISPEP %3d %-3s %c %4d%c   %-3s %c %4d%c       %3d       %6.2f",
                         i+1, c->res1.resName, c->res1.chain, (int16_t)c->res1.resSeq, c->res1.iCode == 0 ? ' ' : c->res1.iCode, c->res2.resName, c->res2.chain, (int16_t)c->res2.resSeq, c->res2.iCode == 0 ? ' ' : c->res2.iCode, c->modNum, c->angle);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    if (len_sites > 0) {
      size_t s_ptr = 0; int s_seq = 1; char last_site[4] = {0};
      while (s_ptr < len_sites) {
        MetaSiteHeader* sh = (MetaSiteHeader*)(p_sites + s_ptr); s_ptr += sizeof(MetaSiteHeader);
        if (strncmp(last_site, sh->siteID, 3) != 0) { s_seq = 1; memcpy(last_site, sh->siteID, 3); }
        for (uint32_t b = 0; b < sh->num_res; b += 4) {
          int printed = snprintf(tbuf, 128, "SITE   %3d %-3s %2d ", s_seq++, sh->siteID, sh->num_res);
          for (uint32_t t = 0; t < 4 && b+t < sh->num_res; t++) { PdbResId* r = (PdbResId*)(p_sites + s_ptr + (b+t)*sizeof(PdbResId)); printed += snprintf(tbuf + printed, 128 - printed, "%-3s %c%4d%c ", r->resName, r->chain, (int16_t)r->resSeq, r->iCode == 0 ? ' ' : r->iCode); }
          memcpy(w_ptr, tbuf, printed); w_ptr += 81;
        }
        s_ptr += sh->num_res * sizeof(PdbResId);
      }
    }
  }

  if (mask & (1ULL << META_AFDB)) {
    memcpy(w_ptr, AFDB_CRYST_BLOCK, 7 * 81);
    w_ptr += 7 * 81;
  }

  int atom_serial = 1;
  if (total_models_minus_1 > 0 || (mask & (1ULL << META_AFDB))) {
    int len = snprintf(tbuf, 128, "MODEL     %4d", current_model_minus_1 + 1);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }

  int has_anchor = 0;
  while (offset < size) {
    uint32_t r_a = *(uint32_t*)(buffer + offset);
    uint32_t num_residues = r_a >> 18;
    uint32_t num_atoms = r_a & 0x3FFFF;
    if (num_residues == 0x3FFE) { offset += 8; continue; }
    if (num_residues == 0x3FFF) { break; }

    uint8_t c1_raw = buffer[offset + 4];
    uint8_t c2_raw = buffer[offset + 5];
    int16_t seq_raw = *(int16_t*)(buffer + offset + 6);

    char c1 = c1_raw & 0x7F; if (c1 == 0) c1 = ' ';
    char c2 = c2_raw & 0x7F; if (c2 == 0) c2 = ' ';

    int is_n_term = (c1_raw & 0x80) != 0;
    int is_c_term = (c2_raw & 0x80) != 0;
    int has_ter = seq_raw & 1;
    int16_t start_resSeq = seq_raw >> 1;

    offset += 8;

    const uint8_t* seq_ptr = buffer + offset;
    int seq_bytes = 0; uint32_t expected_atoms = 0;
    int first_seen = 0;

    for (uint32_t r = 0; r < num_residues; r++) {
      uint8_t type = seq_ptr[seq_bytes++];
      if (type == 0x7C) {
        seq_bytes++;
      } else if (type != 0x00 && type != 0x7F) {
        int is_first = (!first_seen) && is_n_term;
        int is_last  = (r == (uint32_t)(num_residues - 1)) && is_c_term;
        first_seen = 1;

        expected_atoms += (is_last ? UDICT[type].h_end_cterm : UDICT[type].h_end_norm) - (is_first ? 0 : UDICT[type].h_start_norm);
        if (has_hydrogens) {
          expected_atoms += (is_last ? UDICT[type].hy_end_cterm : UDICT[type].hy_end_norm) - (is_first ? 0 : UDICT[type].hy_start_norm);
        }
      }
    }
    offset += seq_bytes + ((8 - (seq_bytes & 7)) & 7);
    const uint8_t* atom_ptr = buffer + offset;

    int32_t current_resSeq = start_resSeq; int first_res = 1;
    uint8_t last_valid_type = 0; char iCode = ' ';

    int32_t current_occ = 100;
    int8_t current_charge = 0;
    char active_alt_char = 'A';
    int first_decoded = 0;

    accum_t accum;
    accum_t anchor_accum;
    has_anchor = 0;
#if defined(__AVX2__) || defined(__SSE4_1__)
    accum = _mm_setzero_si128();
    anchor_accum = _mm_setzero_si128();
#elif defined(__aarch64__) || defined(__ARM_NEON)
    accum = vdupq_n_s32(0);
    anchor_accum = vdupq_n_s32(0);
#else
    accum.x = 0; accum.y = 0; accum.z = 0; accum.b = 0;
    anchor_accum.x = 0; anchor_accum.y = 0; anchor_accum.z = 0; anchor_accum.b = 0;
#endif

    int atoms_left = 0;
    int atom_dict_idx = 0;
    uint32_t parsed_base_atoms = 0;
    int s_idx = 0;
    int no_inc = 0;

    int processing_hydrogens = 0;
    int h_start = 0, h_end = 0;
    int hy_start = 0, hy_end = 0;

    for (uint32_t block_idx = 0; block_idx < num_atoms; block_idx++) {
      if (atom_serial > 99999) atom_serial = 0;

      if (atoms_left == 0) {
        if (has_anchor) {
          accum = anchor_accum;
          has_anchor = 0;
        }
        if (parsed_base_atoms < expected_atoms) {
          uint8_t meta;
          iCode = ' ';
          do {
            if (s_idx >= seq_bytes) break;
            meta = seq_ptr[s_idx++];
            if (meta == 0x00) { current_resSeq++; first_res = 0; no_inc = 0;}
            else if (meta == 0x7F) { no_inc = 1; }
            else if (meta == 0x7C) { iCode = seq_ptr[s_idx++]; }
          } while (meta == 0x00 || meta == 0x7F || meta == 0x7C);
          if (s_idx <= seq_bytes) {
            last_valid_type = meta & 0xFF;
            if (!first_res && !no_inc) { current_resSeq++; }
            first_res = 0; no_inc = 0;

            int is_first = (!first_decoded) && is_n_term;
            int is_last  = (s_idx == seq_bytes) && is_c_term;
            first_decoded = 1;

            h_start = is_first ? 0 : UDICT[last_valid_type].h_start_norm;
            h_end   = is_last  ? UDICT[last_valid_type].h_end_cterm : UDICT[last_valid_type].h_end_norm;

            hy_start = 0; hy_end = 0;
            if (has_hydrogens) {
                hy_start = is_first ? 0 : UDICT[last_valid_type].hy_start_norm;
                hy_end   = is_last  ? UDICT[last_valid_type].hy_end_cterm : UDICT[last_valid_type].hy_end_norm;
            }

            atoms_left = (h_end - h_start) + (hy_end - hy_start);
            atom_dict_idx = h_start;
            processing_hydrogens = 0;
          } else {
            atoms_left = 1;
          }
        } else {
          atoms_left = 1;
        }
        active_alt_char = 'A';
      }

      int alt_flag = 0;

      if (*(uint64_t*)atom_ptr != 0) {
        SIMD_DECODE_MACRO(atom_ptr, accum, current_occ, current_charge, alt_flag);
        int32_t v_b = GET_SIMD_B(accum);

        int anchor_idx = (last_valid_type >= 0x20 && last_valid_type < 0x60) ? 2 : 8;
        if (last_valid_type < 0x60 && atom_dict_idx == anchor_idx && !processing_hydrogens) {
          if (active_alt_char == 'A') {
            anchor_accum = accum;
            has_anchor = 1;
          }
        }

        char altLoc = ' ';
        if (alt_flag || active_alt_char > 'A') {
          altLoc = active_alt_char++;
        }

        const char* a_name;
        const char* e_name;

        if (!processing_hydrogens && atom_dict_idx >= h_start && atom_dict_idx < h_end) {
            a_name = UDICT[last_valid_type].h_atoms[atom_dict_idx];
            e_name = UDICT[last_valid_type].h_elements[atom_dict_idx];
        } else if (processing_hydrogens && atom_dict_idx >= hy_start && atom_dict_idx < hy_end) {
            a_name = UDICT[last_valid_type].hy_atoms[atom_dict_idx];
            e_name = UDICT[last_valid_type].hy_elements[atom_dict_idx];
        } else {
            a_name = " UNK";
            e_name = " X";
        }

        const char* prefix = (last_valid_type < 0x60) ? "ATOM  " : "HETATM";

        if (af_mode) {
          int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f  1.00%6.2f          %2s  ",
                             prefix, atom_serial++, a_name, altLoc, UDICT[last_valid_type].resName,
                             c1, c2, current_resSeq, iCode, GET_SIMD_X(accum) / 1000.0f, GET_SIMD_Y(accum) / 1000.0f, GET_SIMD_Z(accum) / 1000.0f, v_b / 100.0f, e_name);
          if (len > 80) len = 80;
          memcpy(w_ptr, tbuf, len); w_ptr += 81;
        } else {
          char chg_str[3] = "  ";
          if (current_charge > 0) { chg_str[0] = '0' + current_charge; chg_str[1] = '+'; }
          else if (current_charge < 0) { chg_str[0] = '0' - current_charge; chg_str[1] = '-'; }

          int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                             prefix, atom_serial++, a_name, altLoc, UDICT[last_valid_type].resName,
                             c1, c2, current_resSeq, iCode, GET_SIMD_X(accum) / 1000.0f, GET_SIMD_Y(accum) / 1000.0f, GET_SIMD_Z(accum) / 1000.0f, current_occ / 100.0f, v_b / 100.0f, e_name, chg_str);
          if (len > 80) len = 80;
          memcpy(w_ptr, tbuf, len); w_ptr += 81;
        }
      } else {
        atom_ptr += 8;
      }

      if (alt_flag == 0) {
        atoms_left--;
        parsed_base_atoms++;
        atom_dict_idx++;
        active_alt_char = 'A';

        if (!processing_hydrogens && atom_dict_idx == h_end) {
            processing_hydrogens = 1;
            atom_dict_idx = hy_start;
        }
      }
    }

    if (has_ter) {
      int len = snprintf(tbuf, 128, "TER   %5d      %3s%c%c%4d%c", atom_serial++, UDICT[last_valid_type].resName, c1, c2, current_resSeq, iCode);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
    offset = (size_t)(atom_ptr - buffer);
  }

  if (total_models_minus_1 > 0 || (mask & (1ULL << META_AFDB))) {
    memcpy(w_ptr, "ENDMDL                                                                          ", 80); w_ptr += 81;
  }

  if (is_last_model && conect_offset > 0) {
    size_t conect_scan = conect_offset;
    while (conect_scan + 4 <= size) {
      int32_t N = *(int32_t*)(buffer + conect_scan);
      if (N <= 0 || N > 4 || conect_scan + 8 + (N << 2) > size) break;
      conect_scan += 4; int32_t u = *(int32_t*)(buffer + conect_scan); conect_scan += 4;
      int printed = snprintf(tbuf, 128, "CONECT%5d", u);
      for (int32_t t = 0; t < N; t++) { printed += snprintf(tbuf + printed, 128 - printed, "%5d", *(int32_t*)(buffer + conect_scan)); conect_scan += 4; }
      memcpy(w_ptr, tbuf, printed); w_ptr += 81;
    }
  }
  memcpy(w_ptr, "END                                                                             ", 80);
  fwrite(out_data, 1, out_size, out); free(out_data); free(UDICT);
}
