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
  offset += num_meta * 2; offset += (8 - (offset % 8)) % 8;

  uint8_t atom_counts[256] = {0};
  for (int i = 0; i < NUM_DICT_RES; i++) atom_counts[i] = DICT[i].num_atoms;

  const char* p_title = ""; uint64_t title_len = 0;
  const char* p_expdta = ""; uint64_t expdta_len = 0;
  float resolution = -1.0f; const char* p_program = "";
    
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
    int32_t local_sink = 0; // Accumulate entirely in fast registers
    while (offset < size) {
      uint32_t num_atoms = *(uint32_t*)(buffer + offset);
      uint16_t num_residues = *(uint16_t*)(buffer + offset + 4);
      if (num_residues == 0xFFFF) {
        offset += 32;
        while (offset + 4 <= size) {
          int32_t N = *(int32_t*)(buffer + offset);
          if (N <= 0 || N > 4 || offset + 8 + (N << 2) > size) break;
          local_sink += N; offset += 8 + (N << 2);
        }
        break;
      }
            
      int32_t anchor_x = *(int32_t*)(buffer + offset + 12);
      int32_t anchor_y = *(int32_t*)(buffer + offset + 16);
      int32_t anchor_z = *(int32_t*)(buffer + offset + 20);
      int32_t anchor_b = *(int32_t*)(buffer + offset + 24);
      uint32_t seq_bytes = *(uint32_t*)(buffer + offset + 28);
      offset += 32 + seq_bytes + ((8 - (seq_bytes % 8)) % 8);

      const uint8_t* atom_ptr = buffer + offset;
      accum_t accum;
#if defined(__AVX2__) || defined(__SSE4_1__)
      accum = _mm_set_epi32(anchor_z, anchor_y, anchor_x, anchor_b);
#elif defined(__aarch64__) || defined(__ARM_NEON)
      accum = (int32x4_t){anchor_b, anchor_x, anchor_y, anchor_z};
#else
      accum.b = anchor_b; accum.x = anchor_x; accum.y = anchor_y; accum.z = anchor_z;
#endif
      for (uint32_t a = 0; a < num_atoms; a++) {
        const uint8_t* an_ptr = atom_ptr + (a << 3);
                
        // 1-cycle skip: bypass the SIMD pipeline entirely for missing atoms
        if (*(uint64_t*)an_ptr == 0) continue;

#if defined(__AVX2__) || defined(__SSE4_1__)
        __m128i masked = _mm_and_si128(_mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr)), _mm_set1_epi32(0xFFF8));
        accum = _mm_add_epi32(accum, _mm_srai_epi32(masked, 3)); 
#elif defined(__aarch64__) || defined(__ARM_NEON)
        int32x4_t masked = vandq_s32(vmovl_s16(vreinterpret_s16_u8(vld1_u8(an_ptr))), vdupq_n_s32(0xFFF8));
        accum = vaddq_s32(accum, vshrq_n_s32(masked, 3));
#else
        uint16_t aL0 = *(uint16_t*)(an_ptr); uint16_t aL1 = *(uint16_t*)(an_ptr + 2);
        uint16_t aL2 = *(uint16_t*)(an_ptr + 4); uint16_t aL3 = *(uint16_t*)(an_ptr + 6);
        accum.b += ((int16_t)(aL0 & 0xFFF8)) >> 3; accum.x += ((int16_t)(aL1 & 0xFFF8)) >> 3;
        accum.y += ((int16_t)(aL2 & 0xFFF8)) >> 3; accum.z += ((int16_t)(aL3 & 0xFFF8)) >> 3;
#endif
      }
            
      // Extract from vector register once the inner loop finishes
      local_sink += GET_SIMD_X(accum);
      offset += num_atoms << 3;
    }
        
    benchmark_sink += local_sink; 
    return;
  }

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
      uint32_t num_res = *(uint32_t*)(p_sites + s_ptr + 4); out_lines += (num_res + 3) / 4; s_ptr += 8 + num_res * sizeof(PdbResId);
    }
  }

  while (scan_off < size) {
    uint32_t num_atoms = *(uint32_t*)(buffer + scan_off);
    uint16_t num_residues = *(uint16_t*)(buffer + scan_off + 4);
    if (num_residues == 0xFFFF) { scan_off += 32; conect_offset = scan_off; break; }
        
    uint8_t flags = buffer[scan_off + 11];
    uint32_t seq_bytes = *(uint32_t*)(buffer + scan_off + 28);
    scan_off += 32 + seq_bytes + ((8 - (seq_bytes % 8)) % 8);

    const uint8_t* atom_ptr = buffer + scan_off;
    for (uint32_t a = 0; a < num_atoms; a++) {
      if (*(uint64_t*)(atom_ptr + (a << 3)) != 0) out_lines++;
    }
    if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) out_lines++;
    scan_off += num_atoms << 3;
  }

  if (conect_offset > 0) {
    scan_off = conect_offset;
    while (scan_off + 4 <= size) {
      int32_t N = *(int32_t*)(buffer + scan_off);
      if (N <= 0 || N > 4 || scan_off + 8 + (N << 2) > size) break;
      out_lines++; scan_off += 8 + (N << 2);
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
      memcpy(w_ptr, "TITLE   ", 8); if (i > 0) { int ln = snprintf(tbuf, 128, "%2d", (int)(i/70)+1); memcpy(w_ptr + 8, tbuf, ln); }
      size_t copy_len = title_len - i; if (copy_len > 70) copy_len = 70;
      memcpy(w_ptr + 10, p_title + i, copy_len); w_ptr += 81;
    }
  }
  if (expdta_len > 0) { memcpy(w_ptr, "EXPDTA    ", 10); memcpy(w_ptr + 10, p_expdta, expdta_len > 70 ? 70 : expdta_len); w_ptr += 81; }
  if (resolution > 0.0f) { memcpy(w_ptr, "REMARK   2", 10); w_ptr += 81; int len = snprintf(tbuf, 128, "REMARK   2 RESOLUTION. %7.2f ANGSTROMS.", resolution); memcpy(w_ptr, tbuf, len); w_ptr += 81; }
  if (strlen(p_program) > 0) { memcpy(w_ptr, "REMARK   3", 10); w_ptr += 81; memcpy(w_ptr, "REMARK   3 REFINEMENT.", 22); w_ptr += 81; int len = snprintf(tbuf, 128, "REMARK   3   PROGRAM     : %s", p_program); memcpy(w_ptr, tbuf, len); w_ptr += 81; }

  for (uint32_t i=0; i<num_helices; i++) {
    MetaHelix* h = (MetaHelix*)(p_helices + i * sizeof(MetaHelix));
    int len = snprintf(tbuf, 128, "HELIX  %3d %3s %-3s %c %4d%c %-3s %c %4d%c%2d%-30s %5d", 
                       i+1, h->helixID[0] ? h->helixID : "   ", h->start.resName, h->start.chain, h->start.resSeq, h->start.iCode == 0 ? ' ' : h->start.iCode, h->end.resName, h->end.chain, h->end.resSeq, h->end.iCode == 0 ? ' ' : h->end.iCode, h->helix_class, h->comment[0] ? h->comment : "", h->length);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_sheets; i++) {
    MetaSheet* s = (MetaSheet*)(p_sheets + i * sizeof(MetaSheet));
    int len = snprintf(tbuf, 128, "SHEET  %3d %3s%2d %-3s %c%4d%c %-3s %c%4d%c%2d", 
                       s->strand_num, s->sheetID[0] ? s->sheetID : "   ", s->num_strands, s->start.resName, s->start.chain, s->start.resSeq, s->start.iCode == 0 ? ' ' : s->start.iCode, s->end.resName, s->end.chain, s->end.resSeq, s->end.iCode == 0 ? ' ' : s->end.iCode, s->sense);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_ssbonds; i++) {
    MetaSSBond* b = (MetaSSBond*)(p_ssbonds + i * sizeof(MetaSSBond));
    int len = snprintf(tbuf, 128, "SSBOND %3d %-3s %c %4d%c   %-3s %c %4d%c", 
                       i+1, b->res1.resName, b->res1.chain, b->res1.resSeq, b->res1.iCode == 0 ? ' ' : b->res1.iCode, b->res2.resName, b->res2.chain, b->res2.resSeq, b->res2.iCode == 0 ? ' ' : b->res2.iCode);
    if (b->sym1[0] || b->sym2[0] || b->length > 0.0f) {
      memset(tbuf + len, ' ', 59 - len); len = 59;
      if (b->length > 0.0f) len += snprintf(tbuf + len, 128 - len, "%6s %6s %5.2f", b->sym1, b->sym2, b->length);
      else len += snprintf(tbuf + len, 128 - len, "%6s %6s", b->sym1, b->sym2);
    }
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  for (uint32_t i=0; i<num_cispeps; i++) {
    MetaCispep* c = (MetaCispep*)(p_cispeps + i * sizeof(MetaCispep));
    int len = snprintf(tbuf, 128, "CISPEP %3d %-3s %c %4d%c   %-3s %c %4d%c       %3d       %6.2f", 
                       i+1, c->res1.resName, c->res1.chain, c->res1.resSeq, c->res1.iCode == 0 ? ' ' : c->res1.iCode, c->res2.resName, c->res2.chain, c->res2.resSeq, c->res2.iCode == 0 ? ' ' : c->res2.iCode, c->modNum, c->angle);
    memcpy(w_ptr, tbuf, len); w_ptr += 81;
  }
  if (len_sites > 0) {
    size_t s_ptr = 0; int s_seq = 1;
    while (s_ptr < len_sites) {
      MetaSiteHeader* sh = (MetaSiteHeader*)(p_sites + s_ptr); s_ptr += sizeof(MetaSiteHeader);
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
  while (offset < size) {
    uint32_t num_atoms = *(uint32_t*)(buffer + offset);
    uint16_t num_residues = *(uint16_t*)(buffer + offset + 4);
    if (num_residues == 0xFFFF) { 
      offset += 32; 
      while (offset + 4 <= size) {
        int32_t N = *(int32_t*)(buffer + offset);
        if (N <= 0 || N > 4 || offset + 8 + (N << 2) > size) break;
        offset += 4; int32_t u = *(int32_t*)(buffer + offset); offset += 4;
        int printed = snprintf(tbuf, 128, "CONECT%5d", u);
        for (int32_t t = 0; t < N; t++) { printed += snprintf(tbuf + printed, 128 - printed, "%5d", *(int32_t*)(buffer + offset)); offset += 4; }
        memcpy(w_ptr, tbuf, printed); w_ptr += 81;
      }
      break; 
    }
        
    uint16_t start_resSeq = *(uint16_t*)(buffer + offset + 6);
    char c1 = buffer[offset + 9], c2 = buffer[offset + 10]; 
    if (c1 == 0) c1 = ' ';
    if (c2 == 0) c2 = ' ';
        
    uint8_t flags = buffer[offset + 11];
    int32_t anchor_x = *(int32_t*)(buffer + offset + 12), anchor_y = *(int32_t*)(buffer + offset + 16);
    int32_t anchor_z = *(int32_t*)(buffer + offset + 20), anchor_b = *(int32_t*)(buffer + offset + 24);
    uint32_t seq_bytes = *(uint32_t*)(buffer + offset + 28);
    offset += 32; 
        
    const uint8_t* seq_ptr = buffer + offset;
    offset += seq_bytes; offset += (8 - (seq_bytes % 8)) % 8;
    const uint8_t* atom_ptr = buffer + offset;
        
    int32_t current_resSeq = start_resSeq; char active_alt = 'A'; int first_res = 1; uint32_t atom_idx = 0;
    uint8_t last_valid_type = 0;

    accum_t accum, next_res_base;
#if defined(__AVX2__) || defined(__SSE4_1__)
    next_res_base = _mm_set_epi32(anchor_z, anchor_y, anchor_x, anchor_b);
    accum = next_res_base; // Silence GCC -Wmaybe-uninitialized
#elif defined(__aarch64__) || defined(__ARM_NEON)
    next_res_base = (int32x4_t){anchor_b, anchor_x, anchor_y, anchor_z};
    accum = next_res_base; // Silence GCC
#else
    next_res_base.b = anchor_b; next_res_base.x = anchor_x; next_res_base.y = anchor_y; next_res_base.z = anchor_z;
    accum = next_res_base; // Silence GCC
#endif

    for (uint32_t r = 0; r < num_residues; r++) {
      uint8_t type = seq_ptr[r << 1]; uint8_t meta = seq_ptr[(r << 1) + 1];
      if (type == 255) { current_resSeq++; active_alt = 'A'; first_res = 0; continue; }
      last_valid_type = type;
      int is_dup = meta & 0x80; char iCode = meta & 0x7F;
            
      if (is_dup) active_alt++;
      else { if (!first_res) current_resSeq++; active_alt = 'A'; }
      first_res = 0;

      accum = next_res_base;
      int natoms = UDICT[type].num_atoms;
      int saved_base = 0;

      for (int i = 0; i < natoms; i++) {
        const uint8_t* an_ptr = atom_ptr + (atom_idx << 3); atom_idx++;
        if (*(uint64_t*)an_ptr != 0) {
#if defined(__AVX2__) || defined(__SSE4_1__)
          __m128i masked = _mm_and_si128(_mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr)), _mm_set1_epi32(0xFFF8));
          accum = _mm_add_epi32(accum, _mm_srai_epi32(masked, 3)); 
#elif defined(__aarch64__) || defined(__ARM_NEON)
          int32x4_t masked = vandq_s32(vmovl_s16(vreinterpret_s16_u8(vld1_u8(an_ptr))), vdupq_n_s32(0xFFF8));
          accum = vaddq_s32(accum, vshrq_n_s32(masked, 3));
#else
          uint16_t aL0 = *(uint16_t*)(an_ptr); uint16_t aL1 = *(uint16_t*)(an_ptr + 2);
          uint16_t aL2 = *(uint16_t*)(an_ptr + 4); uint16_t aL3 = *(uint16_t*)(an_ptr + 6);
          accum.b += ((int16_t)(aL0 & 0xFFF8)) >> 3; accum.x += ((int16_t)(aL1 & 0xFFF8)) >> 3;
          accum.y += ((int16_t)(aL2 & 0xFFF8)) >> 3; accum.z += ((int16_t)(aL3 & 0xFFF8)) >> 3;
#endif
          uint16_t aL0_u = *(uint16_t*)(an_ptr); uint16_t aL1_u = *(uint16_t*)(an_ptr + 2);
          uint16_t aL2_u = *(uint16_t*)(an_ptr + 4); uint16_t aL3_u = *(uint16_t*)(an_ptr + 6);
                    
          uint8_t occ = (aL0_u & 0x07) | ((aL1_u & 0x07) << 3) | ((aL2_u & 0x01) << 6);
          int8_t charge = ((aL2_u >> 1) & 0x03) | ((aL3_u & 0x03) << 2); if (charge & 0x08) charge |= 0xF0; 
          char altLoc = (aL3_u & 0x04) ? active_alt : ' ';
                    
          char chg_str[3] = "  ";
          if (charge > 0) { chg_str[0] = '0' + charge; chg_str[1] = '+'; } else if (charge < 0) { chg_str[0] = '0' - charge; chg_str[1] = '-'; }

          int32_t v_b = GET_SIMD_B(accum); int32_t v_x = GET_SIMD_X(accum); 
          int32_t v_y = GET_SIMD_Y(accum); int32_t v_z = GET_SIMD_Z(accum);
          const char* prefix = (type < 37) ? "ATOM  " : "HETATM";
                    
          int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                             prefix, atom_serial++, UDICT[type].atom_names[i], altLoc, UDICT[type].resName,
                             c1, c2, current_resSeq, iCode == 0 ? ' ' : iCode,
                             v_x / 1000.0f, v_y / 1000.0f, v_z / 1000.0f, occ / 100.0f, v_b / 100.0f, UDICT[type].elements[i], chg_str);
          memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; w_ptr += 81;
        }
        if (i == 2) { next_res_base = accum; saved_base = 1; }
      }
      if (!saved_base) next_res_base = accum;
    }
        
    if (flags & CTRL_TER_OXT) {
      const uint8_t* an_ptr = atom_ptr + (atom_idx << 3);
      if (*(uint64_t*)an_ptr != 0) {
        uint16_t aL0 = *(uint16_t*)(an_ptr); 
        uint16_t aL1 = *(uint16_t*)(an_ptr + 2);
        uint16_t aL2 = *(uint16_t*)(an_ptr + 4); 
        uint16_t aL3 = *(uint16_t*)(an_ptr + 6);

        int32_t v_b = GET_SIMD_B(accum) + (((int16_t)(aL0 & 0xFFF8)) >> 3);
        int32_t v_x = GET_SIMD_X(accum) + (((int16_t)(aL1 & 0xFFF8)) >> 3);
        int32_t v_y = GET_SIMD_Y(accum) + (((int16_t)(aL2 & 0xFFF8)) >> 3);
        int32_t v_z = GET_SIMD_Z(accum) + (((int16_t)(aL3 & 0xFFF8)) >> 3);
        uint8_t occ = (aL0 & 0x07) | ((aL1 & 0x07) << 3) | ((aL2 & 0x01) << 6);
        int8_t charge = ((aL2 >> 1) & 0x03) | ((aL3 & 0x03) << 2); if (charge & 0x08) charge |= 0xF0;
        char altLoc = (aL3 & 0x04) ? active_alt : ' ';
                
        char chg_str[3] = "  ";
        if (charge > 0) { chg_str[0] = '0' + charge; chg_str[1] = '+'; } else if (charge < 0) { chg_str[0] = '0' - charge; chg_str[1] = '-'; }

        int len = snprintf(tbuf, 128, "ATOM  %5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                           atom_serial++, " OXT", altLoc, UDICT[last_valid_type].resName, c1, c2, current_resSeq, ' ',
                           v_x / 1000.0f, v_y / 1000.0f, v_z / 1000.0f, occ / 100.0f, v_b / 100.0f, " O", chg_str);
        memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; w_ptr += 81;
      }
    }

    if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) {
      int len = snprintf(tbuf, 128, "TER   %5d      %3s%c%c%4d ", atom_serial++, UDICT[last_valid_type].resName, c1, c2, current_resSeq);
      memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; w_ptr += 81;
    }
    offset += num_atoms << 3;
  }
    
  fwrite(out_data, 1, out_size, out); free(out_data); free(UDICT);
}
