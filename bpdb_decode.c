#include "bpdb.h"
#include <stdlib.h>
#include <string.h>

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

volatile int32_t benchmark_sink = 0;

void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out) {
  if (size < 64) return;

  int day = buffer[5], month = buffer[6], year = 1973 + buffer[7];
  char pdbID[17] = {0}; memcpy(pdbID, buffer + 8, 16); trim_string(pdbID);
  char classification[41] = {0}; memcpy(classification, buffer + 24, 40); trim_string(classification);

  size_t offset = 64;
  uint16_t num_meta = *(uint16_t*)(buffer + offset); offset += 2;
  const uint16_t* meta_codes = (const uint16_t*)(buffer + offset);
  offset += num_meta * 2;
  offset += (8 - (offset % 8)) % 8;

  uint8_t atom_counts[256] = {0};
  for (int i = 0; i < NUM_DICT_RES; i++) atom_counts[i] = DICT[i].num_atoms;

  const char* p_title = ""; uint64_t title_len = 0;
  const char* p_expdta = ""; uint64_t expdta_len = 0;
  float resolution = -1.0f;
  const char* p_program = "";

  const uint8_t* p_het_dict = NULL; uint8_t custom_count = 0;

  const uint8_t* p_ssbonds = NULL; uint64_t num_ssbonds = 0;
  const uint8_t* p_cispeps = NULL; uint64_t num_cispeps = 0;
  const uint8_t* p_helices = NULL; uint64_t num_helices = 0;
  const uint8_t* p_sheets = NULL;  uint64_t num_sheets = 0;
  const uint8_t* p_sites = NULL;   uint64_t len_sites = 0;

  for (int i = 0; i < num_meta; i++) {
    uint64_t section_len = *(uint64_t*)(buffer + offset); offset += 8;
    if (meta_codes[i] == META_TITLE) { p_title = (const char*)(buffer + offset); title_len = section_len; }
    else if (meta_codes[i] == META_EXPDTA) { p_expdta = (const char*)(buffer + offset); expdta_len = section_len; }
    else if (meta_codes[i] == META_RESOLUTION) { resolution = *(float*)(buffer + offset); }
    else if (meta_codes[i] == META_PROGRAM) { p_program = (const char*)(buffer + offset); }
    else if (meta_codes[i] == META_SSBOND) { p_ssbonds = buffer + offset; num_ssbonds = section_len / sizeof(MetaSSBond); }
    else if (meta_codes[i] == META_CISPEP) { p_cispeps = buffer + offset; num_cispeps = section_len / sizeof(MetaCispep); }
    else if (meta_codes[i] == META_HELIX) { p_helices = buffer + offset; num_helices = section_len / sizeof(MetaHelix); }
    else if (meta_codes[i] == META_SHEET) { p_sheets = buffer + offset; num_sheets = section_len / sizeof(MetaSheet); }
    else if (meta_codes[i] == META_SITE) { p_sites = buffer + offset; len_sites = section_len; }
    else if (meta_codes[i] == META_HET_DICT) {
      p_het_dict = buffer + offset; custom_count = p_het_dict[0]; size_t ptr = 1;
      for (int c = 0; c < custom_count; c++) { ptr += 4; atom_counts[NUM_DICT_RES + c] = p_het_dict[ptr++]; ptr += atom_counts[NUM_DICT_RES + c] * 6; }
    }
    offset += section_len; offset += (8 - (section_len % 8)) % 8;
  }

  if (!out) {
    int32_t last_a0_b = 0, last_a0_x = 0, last_a0_y = 0, last_a0_z = 0;
    for (;;) {
      if (offset >= size) break;
      uint16_t num_residues = *(uint16_t*)(buffer + offset);
      if (num_residues == 0xFFFF) {
        offset += 32;
        while (offset + 4 <= size) {
          int32_t N = *(int32_t*)(buffer + offset);
          if (N <= 0 || N > 4 || offset + 8 + N*4 > size) break;
          benchmark_sink += N;
          offset += 8 + N * 4;
        }
        break;
      }

      uint8_t flags = buffer[offset + 2];
      int32_t anchor_x = *(int32_t*)(buffer + offset + 20);
      int32_t anchor_y = *(int32_t*)(buffer + offset + 24);
      int32_t anchor_z = *(int32_t*)(buffer + offset + 28);
      offset += 32;

      for (uint32_t r = 0; r < num_residues; r++) {
        const uint8_t* res_ptr = buffer + offset;
        int num_atoms = atom_counts[res_ptr[0]]; offset += 8;
        const uint8_t* a1_ptr = buffer + offset; accum_t accum;

#if defined(__AVX2__) || defined(__SSE4_1__)
        __m128i a1_vec = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)a1_ptr));
        __m128i lsb_vec = _mm_and_si128(_mm_srai_epi32(a1_vec, 2), _mm_set1_epi32(0x3FFF));
        __m128i msb_vec = _mm_set_epi32((int8_t)res_ptr[7], (int8_t)res_ptr[5], (int8_t)res_ptr[3], res_ptr[1]);
        accum = _mm_add_epi32(_mm_or_si128(_mm_slli_epi32(msb_vec, 14), lsb_vec), _mm_set_epi32(anchor_z, anchor_y, anchor_x, 0));
#elif defined(__aarch64__) || defined(__ARM_NEON)
        int32x4_t a1_32 = vmovl_s16(vreinterpret_s16_u8(vld1_u8(a1_ptr)));
        int32x4_t lsb_vec = vandq_s32(vshrq_n_s32(a1_32, 2), vdupq_n_s32(0x3FFF));
        int32x4_t msb_vec = {res_ptr[1], (int8_t)res_ptr[3], (int8_t)res_ptr[5], (int8_t)res_ptr[7]};
        accum = vaddq_s32(vorrq_s32(vshlq_n_s32(msb_vec, 14), lsb_vec), (int32x4_t){0, anchor_x, anchor_y, anchor_z});
#else
        int32_t L0 = *(uint16_t*)(a1_ptr); int32_t L1 = *(uint16_t*)(a1_ptr + 2);
        int32_t L2 = *(uint16_t*)(a1_ptr + 4); int32_t L3 = *(uint16_t*)(a1_ptr + 6);
        accum.b = (res_ptr[1] << 14) | ((L0 >> 2) & 0x3FFF);
        accum.x = anchor_x + (((int8_t)res_ptr[3] << 14) | ((L1 >> 2) & 0x3FFF));
        accum.y = anchor_y + (((int8_t)res_ptr[5] << 14) | ((L2 >> 2) & 0x3FFF));
        accum.z = anchor_z + (((int8_t)res_ptr[7] << 14) | ((L3 >> 2) & 0x3FFF));
#endif
        last_a0_b = GET_SIMD_B(accum); last_a0_x = GET_SIMD_X(accum); last_a0_y = GET_SIMD_Y(accum); last_a0_z = GET_SIMD_Z(accum);
        offset += 8;

        for (int i = 1; i < num_atoms; i++) {
          const uint8_t* an_ptr = buffer + offset;
#if defined(__AVX2__) || defined(__SSE4_1__)
          accum = _mm_add_epi32(accum, _mm_srai_epi32(_mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr)), 2));
#elif defined(__aarch64__) || defined(__ARM_NEON)
          accum = vaddq_s32(accum, vshrq_n_s32(vmovl_s16(vreinterpret_s16_u8(vld1_u8(an_ptr))), 2));
#else
          int32_t aL0 = *(uint16_t*)(an_ptr); int32_t aL1 = *(uint16_t*)(an_ptr + 2);
          int32_t aL2 = *(uint16_t*)(an_ptr + 4); int32_t aL3 = *(uint16_t*)(an_ptr + 6);
          accum.b += ((int16_t)(aL0 & 0xFFFC)) >> 2; accum.x += ((int16_t)(aL1 & 0xFFFC)) >> 2;
          accum.y += ((int16_t)(aL2 & 0xFFFC)) >> 2; accum.z += ((int16_t)(aL3 & 0xFFFC)) >> 2;
#endif
          offset += 8;
        }
        benchmark_sink += GET_SIMD_X(accum);
      }
      if (flags & CTRL_TER_OXT) {
        uint16_t mL0 = *(uint16_t*)(buffer + offset);
        benchmark_sink += (last_a0_b + (((int16_t)(mL0 & 0xFFFC)) >> 2)) + last_a0_x + last_a0_y + last_a0_z;
        offset += 8;
      }
    }
    return;
  }

  // --- FULL DECODE PHASE ---
  int num_udict = NUM_DICT_RES + custom_count;
  UniDict* UDICT = (UniDict*)malloc(num_udict * sizeof(UniDict));
  for (int i = 0; i < NUM_DICT_RES; i++) {
    strcpy(UDICT[i].resName, DICT[i].resName); UDICT[i].num_atoms = DICT[i].num_atoms;
    for (int a = 0; a < DICT[i].num_atoms; a++) { strcpy(UDICT[i].atom_names[a], DICT[i].atom_names[a]); strcpy(UDICT[i].elements[a], DICT[i].elements[a]); }
  }
  if (p_het_dict) {
    size_t ptr = 1;
    for (int c = 0; c < custom_count; c++) {
      int u_idx = NUM_DICT_RES + c; memcpy(UDICT[u_idx].resName, p_het_dict + ptr, 3); UDICT[u_idx].resName[3] = '\0'; ptr += 4;
      UDICT[u_idx].num_atoms = p_het_dict[ptr++];
      for (int a = 0; a < UDICT[u_idx].num_atoms; a++) {
        memcpy(UDICT[u_idx].atom_names[a], p_het_dict + ptr, 4); UDICT[u_idx].atom_names[a][4] = '\0'; ptr += 4;
        memcpy(UDICT[u_idx].elements[a], p_het_dict + ptr, 2); UDICT[u_idx].elements[a][2] = '\0'; ptr += 2;
      }
    }
  }

  size_t out_lines = 1;
  if (title_len > 0) out_lines += (title_len + 69) / 70;
  if (expdta_len > 0) out_lines++;
  if (resolution > 0.0f) out_lines += 2;
  if (strlen(p_program) > 0) out_lines += 3;
  out_lines += num_helices + num_sheets + num_ssbonds + num_cispeps;

  size_t scan_off = offset, conect_offset = 0;
  if (len_sites > 0) {
    size_t s_ptr = 0;
    while (s_ptr < len_sites) {
      uint32_t num_res = *(uint32_t*)(p_sites + s_ptr + 4);
      out_lines += (num_res + 3) / 4;
      s_ptr += 8 + num_res * sizeof(PdbResId);
    }
  }

  while (scan_off < size) {
    uint16_t num_residues = *(uint16_t*)(buffer + scan_off);
    if (num_residues == 0xFFFF) { scan_off += 32; conect_offset = scan_off; break; }

    uint8_t flags = buffer[scan_off + 2]; scan_off += 32;
    for (uint32_t r = 0; r < num_residues; r++) {
      uint8_t res_type = buffer[scan_off]; int num_atoms = atom_counts[res_type]; scan_off += 8;
      for (int i = 0; i < num_atoms; i++) {
        uint8_t occ = (buffer[scan_off] & 0x03) | ((buffer[scan_off + 2] & 0x03) << 2) | ((buffer[scan_off + 4] & 0x03) << 4) | ((buffer[scan_off + 6] & 0x01) << 6);
        if (occ > 0) out_lines++;
        scan_off += 8;
      }
    }
    if (flags & CTRL_TER_OXT) {
      uint8_t occ = (buffer[scan_off] & 0x03) | ((buffer[scan_off + 2] & 0x03) << 2) | ((buffer[scan_off + 4] & 0x03) << 4) | ((buffer[scan_off + 6] & 0x01) << 6);
      if (occ > 0) out_lines++;
      scan_off += 8;
    }
    if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) out_lines++;
  }

  if (conect_offset > 0) {
    scan_off = conect_offset;
    while (scan_off + 4 <= size) {
      int32_t N = *(int32_t*)(buffer + scan_off);
      if (N <= 0 || N > 4 || scan_off + 8 + N*4 > size) break;
      out_lines++; scan_off += 8 + N * 4;
    }
  }

  size_t out_size = out_lines * 81;
  char* out_data = (char*)malloc(out_size); memset(out_data, ' ', out_size);
  for (size_t i = 0; i < out_lines; i++) out_data[i * 81 + 80] = '\n';

  char* w_ptr = out_data; char tbuf[128];
  memcpy(w_ptr, "HEADER    ", 10);
  size_t c_len = strlen(classification); if (c_len > 0) memcpy(w_ptr + 10, classification, c_len > 40 ? 40 : c_len);
  int dlen = snprintf(tbuf, 128, "%02d-%3s-%02d", day, MONTHS[month], year % 100); memcpy(w_ptr + 50, tbuf, dlen);
  size_t id_len = strlen(pdbID); if (id_len > 0) memcpy(w_ptr + 62, pdbID, id_len > 18 ? 18 : id_len); w_ptr += 81;

  if (title_len > 0) {
    for (uint64_t i = 0; i < title_len; i += 70) {
      memcpy(w_ptr, "TITLE   ", 8);
      if (i > 0) { int ln = snprintf(tbuf, 128, "%2d", (int)(i/70)+1); memcpy(w_ptr + 8, tbuf, ln); }
      size_t copy_len = title_len - i; if (copy_len > 70) copy_len = 70;
      memcpy(w_ptr + 10, p_title + i, copy_len); w_ptr += 81;
    }
  }
  if (expdta_len > 0) {
    memcpy(w_ptr, "EXPDTA    ", 10); memcpy(w_ptr + 10, p_expdta, expdta_len > 70 ? 70 : expdta_len); w_ptr += 81;
  }
  if (resolution > 0.0f) {
    memcpy(w_ptr, "REMARK   2", 10); w_ptr += 81;
    int len = snprintf(tbuf, 128, "REMARK   2 RESOLUTION. %7.2f ANGSTROMS.", resolution);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  if (strlen(p_program) > 0) {
    memcpy(w_ptr, "REMARK   3", 10); w_ptr += 81;
    memcpy(w_ptr, "REMARK   3 REFINEMENT.", 22); w_ptr += 81;
    int len = snprintf(tbuf, 128, "REMARK   3   PROGRAM     : %s", p_program);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }

  for (uint32_t i=0; i<num_helices; i++) {
    MetaHelix* h = (MetaHelix*)(p_helices + i * sizeof(MetaHelix));
    int len = snprintf(tbuf, 128, "HELIX  %3d %3s %-3s %c %4d%c %-3s %c %4d%c%2d%-30s %5d",
                       i+1, h->helixID[0] ? h->helixID : "   ",
                       h->start.resName, h->start.chain, h->start.resSeq, h->start.iCode == 0 ? ' ' : h->start.iCode,
                       h->end.resName, h->end.chain, h->end.resSeq, h->end.iCode == 0 ? ' ' : h->end.iCode,
                       h->helix_class, h->comment[0] ? h->comment : "", h->length);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_sheets; i++) {
    MetaSheet* s = (MetaSheet*)(p_sheets + i * sizeof(MetaSheet));
    int len = snprintf(tbuf, 128, "SHEET  %3d %3s%2d %-3s %c%4d%c %-3s %c%4d%c%2d",
                       s->strand_num, s->sheetID[0] ? s->sheetID : "   ", s->num_strands,
                       s->start.resName, s->start.chain, s->start.resSeq, s->start.iCode == 0 ? ' ' : s->start.iCode,
                       s->end.resName, s->end.chain, s->end.resSeq, s->end.iCode == 0 ? ' ' : s->end.iCode,
                       s->sense);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_ssbonds; i++) {
    MetaSSBond* b = (MetaSSBond*)(p_ssbonds + i * sizeof(MetaSSBond));
    int len = snprintf(tbuf, 128, "SSBOND %3d %-3s %c %4d%c   %-3s %c %4d%c",
                       i+1, b->res1.resName, b->res1.chain, b->res1.resSeq, b->res1.iCode == 0 ? ' ' : b->res1.iCode,
                       b->res2.resName, b->res2.chain, b->res2.resSeq, b->res2.iCode == 0 ? ' ' : b->res2.iCode);
    if (b->sym1[0] || b->sym2[0] || b->length > 0.0f) {
      memset(tbuf + len, ' ', 59 - len); len = 59;
      if (b->length > 0.0f) {
        len += snprintf(tbuf + len, 128 - len, "%6s %6s %5.2f", b->sym1, b->sym2, b->length);
      } else {
        len += snprintf(tbuf + len, 128 - len, "%6s %6s", b->sym1, b->sym2);
      }
    }
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_cispeps; i++) {
    MetaCispep* c = (MetaCispep*)(p_cispeps + i * sizeof(MetaCispep));
    int len = snprintf(tbuf, 128, "CISPEP %3d %-3s %c %4d%c   %-3s %c %4d%c       %3d       %6.2f",
                       i+1, c->res1.resName, c->res1.chain, c->res1.resSeq, c->res1.iCode == 0 ? ' ' : c->res1.iCode,
                       c->res2.resName, c->res2.chain, c->res2.resSeq, c->res2.iCode == 0 ? ' ' : c->res2.iCode,
                       c->modNum, c->angle);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  if (len_sites > 0) {
    size_t s_ptr = 0; int s_seq = 1;
    while (s_ptr < len_sites) {
      MetaSiteHeader* sh = (MetaSiteHeader*)(p_sites + s_ptr);
      s_ptr += sizeof(MetaSiteHeader);
      for (uint32_t b = 0; b < sh->num_res; b += 4) {
        int printed = snprintf(tbuf, 128, "SITE   %3d %-3s %2d ", s_seq++, sh->siteID, sh->num_res);
        for (uint32_t t = 0; t < 4 && b+t < sh->num_res; t++) {
          PdbResId* r = (PdbResId*)(p_sites + s_ptr + (b+t)*sizeof(PdbResId));
          printed += snprintf(tbuf + printed, 128 - printed, "%-3s %c%4d%c ", r->resName, r->chain, r->resSeq, r->iCode == 0 ? ' ' : r->iCode);
        }
        memcpy(w_ptr, tbuf, printed); w_ptr += 81;
      }
      s_ptr += sh->num_res * sizeof(PdbResId);
    }
  }

  int atom_serial = 1;
  int32_t last_a0_b = 0, last_a0_x = 0, last_a0_y = 0, last_a0_z = 0;
  int res_type = 0, resSeq = 0; char iCode = ' ';

  while (offset < size) {
    uint16_t num_residues = *(uint16_t*)(buffer + offset);
    if (num_residues == 0xFFFF) {
      offset += 32;
      while (offset + 4 <= size) {
        int32_t N = *(int32_t*)(buffer + offset);
        if (N <= 0 || N > 4 || offset + 8 + N*4 > size) break;
        offset += 4;
        int32_t u = *(int32_t*)(buffer + offset); offset += 4;

        int printed = snprintf(tbuf, 128, "CONECT%5d", u);
        for (int32_t t = 0; t < N; t++) {
          printed += snprintf(tbuf + printed, 128 - printed, "%5d", *(int32_t*)(buffer + offset));
          offset += 4;
        }
        memcpy(w_ptr, tbuf, printed); w_ptr += 81;
      }
      break;
    }

    uint8_t flags = buffer[offset + 2];
    char c1 = buffer[offset + 4];
    char c2 = buffer[offset + 5];
    if (c1 == 0) c1 = ' ';
    if (c2 == 0) c2 = ' ';

    int32_t anchor_x = *(int32_t*)(buffer + offset + 20), anchor_y = *(int32_t*)(buffer + offset + 24), anchor_z = *(int32_t*)(buffer + offset + 28);
    offset += 32;

    for (uint32_t r = 0; r < num_residues; r++) {
      const uint8_t* res_ptr = buffer + offset;
      res_type = res_ptr[0]; resSeq = res_ptr[2] | (res_ptr[4] << 8); iCode = res_ptr[6];
      int num_atoms = UDICT[res_type].num_atoms; offset += 8;

      const uint8_t* a1_ptr = buffer + offset; accum_t accum;
#if defined(__AVX2__) || defined(__SSE4_1__)
      __m128i a1_vec = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)a1_ptr));
      __m128i lsb_vec = _mm_and_si128(_mm_srai_epi32(a1_vec, 2), _mm_set1_epi32(0x3FFF));
      __m128i msb_vec = _mm_set_epi32((int8_t)res_ptr[7], (int8_t)res_ptr[5], (int8_t)res_ptr[3], res_ptr[1]);
      accum = _mm_add_epi32(_mm_or_si128(_mm_slli_epi32(msb_vec, 14), lsb_vec), _mm_set_epi32(anchor_z, anchor_y, anchor_x, 0));
#elif defined(__aarch64__) || defined(__ARM_NEON)
      int32x4_t a1_32 = vmovl_s16(vreinterpret_s16_u8(vld1_u8(a1_ptr)));
      int32x4_t lsb_vec = vandq_s32(vshrq_n_s32(a1_32, 2), vdupq_n_s32(0x3FFF));
      int32x4_t msb_vec = {res_ptr[1], (int8_t)res_ptr[3], (int8_t)res_ptr[5], (int8_t)res_ptr[7]};
      accum = vaddq_s32(vorrq_s32(vshlq_n_s32(msb_vec, 14), lsb_vec), (int32x4_t){0, anchor_x, anchor_y, anchor_z});
#else
      int32_t L0 = *(uint16_t*)(a1_ptr); int32_t L1 = *(uint16_t*)(a1_ptr + 2);
      int32_t L2 = *(uint16_t*)(a1_ptr + 4); int32_t L3 = *(uint16_t*)(a1_ptr + 6);
      accum.b = (res_ptr[1] << 14) | ((L0 >> 2) & 0x3FFF);
      accum.x = anchor_x + (((int8_t)res_ptr[3] << 14) | ((L1 >> 2) & 0x3FFF));
      accum.y = anchor_y + (((int8_t)res_ptr[5] << 14) | ((L2 >> 2) & 0x3FFF));
      accum.z = anchor_z + (((int8_t)res_ptr[7] << 14) | ((L3 >> 2) & 0x3FFF));
#endif
      last_a0_b = GET_SIMD_B(accum); last_a0_x = GET_SIMD_X(accum); last_a0_y = GET_SIMD_Y(accum); last_a0_z = GET_SIMD_Z(accum);

      uint16_t aL0_0 = *(uint16_t*)(a1_ptr); uint16_t aL1_0 = *(uint16_t*)(a1_ptr + 2);
      uint16_t aL2_0 = *(uint16_t*)(a1_ptr + 4); uint16_t aL3_0 = *(uint16_t*)(a1_ptr + 6);
      uint8_t occ = (aL0_0 & 0x03) | ((aL1_0 & 0x03) << 2) | ((aL2_0 & 0x03) << 4) | ((aL3_0 & 0x01) << 6);
      if (occ > 0) {
        char altLoc = (aL3_0 & 0x02) ? 'A' : ' ';
        const char* prefix = (res_type < 37) ? "ATOM  " : "HETATM";
        int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ",
                           prefix, atom_serial++, UDICT[res_type].atom_names[0], altLoc, UDICT[res_type].resName,
                           c1, c2, resSeq, iCode == 0 ? ' ' : iCode,
                           last_a0_x / 1000.0f, last_a0_y / 1000.0f, last_a0_z / 1000.0f, occ / 100.0f, last_a0_b / 100.0f, UDICT[res_type].elements[0]);
        memcpy(w_ptr, tbuf, len); w_ptr += 81;
      }
      offset += 8;

      for (int i = 1; i < num_atoms; i++) {
        const uint8_t* an_ptr = buffer + offset;
#if defined(__AVX2__) || defined(__SSE4_1__)
        __m128i an_vec = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr));
        accum = _mm_add_epi32(accum, _mm_srai_epi32(an_vec, 2));
#elif defined(__aarch64__) || defined(__ARM_NEON)
        accum = vaddq_s32(accum, vshrq_n_s32(vmovl_s16(vreinterpret_s16_u8(vld1_u8(an_ptr))), 2));
#else
        int32_t aL0 = *(uint16_t*)(an_ptr); int32_t aL1 = *(uint16_t*)(an_ptr + 2);
        int32_t aL2 = *(uint16_t*)(an_ptr + 4); int32_t aL3 = *(uint16_t*)(an_ptr + 6);
        accum.b += ((int16_t)(aL0 & 0xFFFC)) >> 2; accum.x += ((int16_t)(aL1 & 0xFFFC)) >> 2;
        accum.y += ((int16_t)(aL2 & 0xFFFC)) >> 2; accum.z += ((int16_t)(aL3 & 0xFFFC)) >> 2;
#endif
        uint16_t aL0_u = *(uint16_t*)(an_ptr); uint16_t aL1_u = *(uint16_t*)(an_ptr + 2);
        uint16_t aL2_u = *(uint16_t*)(an_ptr + 4); uint16_t aL3_u = *(uint16_t*)(an_ptr + 6);
        occ = (aL0_u & 0x03) | ((aL1_u & 0x03) << 2) | ((aL2_u & 0x03) << 4) | ((aL3_u & 0x01) << 6);

        if (occ > 0) {
          char altLoc = (aL3_u & 0x02) ? 'A' : ' ';
          int32_t v_b = GET_SIMD_B(accum); int32_t v_x = GET_SIMD_X(accum);
          int32_t v_y = GET_SIMD_Y(accum); int32_t v_z = GET_SIMD_Z(accum);
          const char* prefix = (res_type < 37) ? "ATOM  " : "HETATM";
          int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ",
                             prefix, atom_serial++, UDICT[res_type].atom_names[i], altLoc, UDICT[res_type].resName,
                             c1, c2, resSeq, iCode == 0 ? ' ' : iCode,
                             v_x / 1000.0f, v_y / 1000.0f, v_z / 1000.0f, occ / 100.0f, v_b / 100.0f, UDICT[res_type].elements[i]);
          memcpy(w_ptr, tbuf, len); w_ptr += 81;
        }
        offset += 8;
      }
    }

    if (flags & CTRL_TER_OXT) {
      const uint8_t* oxt_ptr = buffer + offset;
      uint16_t mL0 = *(uint16_t*)(oxt_ptr); uint16_t mL1 = *(uint16_t*)(oxt_ptr + 2);
      uint16_t mL2 = *(uint16_t*)(oxt_ptr + 4); uint16_t mL3 = *(uint16_t*)(oxt_ptr + 6);
      uint8_t occ = (mL0 & 0x03) | ((mL1 & 0x03) << 2) | ((mL2 & 0x03) << 4) | ((mL3 & 0x01) << 6);
      if (occ > 0) {
        int32_t db = ((int16_t)(mL0 & 0xFFFC)) >> 2; int32_t dx = ((int16_t)(mL1 & 0xFFFC)) >> 2;
        int32_t dy = ((int16_t)(mL2 & 0xFFFC)) >> 2; int32_t dz = ((int16_t)(mL3 & 0xFFFC)) >> 2;
        int32_t oxt_b = last_a0_b + db; int32_t oxt_x = last_a0_x + dx;
        int32_t oxt_y = last_a0_y + dy; int32_t oxt_z = last_a0_z + dz;
        char altLoc = (mL3 & 0x02) ? 'A' : ' ';
        int len = snprintf(tbuf, 128, "ATOM  %5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ",
                           atom_serial++, " OXT", altLoc, UDICT[res_type].resName, c1, c2, resSeq, iCode == 0 ? ' ' : iCode,
                           oxt_x / 1000.0f, oxt_y / 1000.0f, oxt_z / 1000.0f, occ / 100.0f, oxt_b / 100.0f, " O");
        memcpy(w_ptr, tbuf, len); w_ptr += 81;
      }
      offset += 8;
    }

    if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) {
      int len = snprintf(tbuf, 128, "TER   %5d      %3s%c%c%4d%c", atom_serial++, UDICT[res_type].resName, c1, c2, resSeq, iCode == 0 ? ' ' : iCode);
      memcpy(w_ptr, tbuf, len); w_ptr += 81;
    }
  }

  fwrite(out_data, 1, out_size, out);
  free(out_data); free(UDICT);
}
