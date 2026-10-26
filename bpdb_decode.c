#include "bpdb.h"
#include <stdlib.h>
#include <string.h>

#if defined(__AVX2__) || defined(__SSE4_1__)
#include <immintrin.h>
    typedef __m128i accum_t;
    #define GET_SIMD_X(acc) _mm_extract_epi32((acc), 1)
    #define GET_SIMD_Y(acc) _mm_extract_epi32((acc), 2)
    #define GET_SIMD_Z(acc) _mm_extract_epi32((acc), 3)
#elif defined(__aarch64__) || defined(__ARM_NEON)
#include <arm_neon.h>
    typedef int32x4_t accum_t;
    #define GET_SIMD_X(acc) vgetq_lane_s32((acc), 1)
    #define GET_SIMD_Y(acc) vgetq_lane_s32((acc), 2)
    #define GET_SIMD_Z(acc) vgetq_lane_s32((acc), 3)
#else
    typedef struct { int32_t unused, x, y, z; } accum_t;
    #define GET_SIMD_X(acc) ((acc).x)
    #define GET_SIMD_Y(acc) ((acc).y)
    #define GET_SIMD_Z(acc) ((acc).z)
#endif

// Zero-Shift Math: Natively casts 15-bit jump MSBs.
#if defined(__AVX2__) || defined(__SSE4_1__)
    #define SIMD_DECODE_PHASE1(ptr, out_accum) do { \
        if ((ptr)[0] & 0x01) { \
            int32_t vx = (((int16_t)*(uint16_t*)((ptr) + 2)) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            int32_t vy = (((int16_t)*(uint16_t*)((ptr) + 4)) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            int32_t vz = (((int16_t)*(uint16_t*)((ptr) + 6)) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            out_accum = _mm_set_epi32(vz, vy, vx, 0); \
            ptr += 16; \
        } else { \
            out_accum = _mm_add_epi32(out_accum, _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)(ptr)))); \
            ptr += 8; \
        } \
    } while(0)
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #define SIMD_DECODE_PHASE1(ptr, out_accum) do { \
        if ((ptr)[0] & 0x01) { \
            int32_t vx = (((int16_t)*(uint16_t*)((ptr) + 2)) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            int32_t vy = (((int16_t)*(uint16_t*)((ptr) + 4)) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            int32_t vz = (((int16_t)*(uint16_t*)((ptr) + 6)) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            out_accum = (int32x4_t){0, vx, vy, vz}; \
            ptr += 16; \
        } else { \
            out_accum = vaddq_s32(out_accum, vmovl_s16(vreinterpret_s16_u8(vld1_u8(ptr)))); \
            ptr += 8; \
        } \
    } while(0)
#else
    #define SIMD_DECODE_PHASE1(ptr, out_accum) do { \
        if ((ptr)[0] & 0x01) { \
            out_accum.x = (((int16_t)*(uint16_t*)((ptr) + 2)) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            out_accum.y = (((int16_t)*(uint16_t*)((ptr) + 4)) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            out_accum.z = (((int16_t)*(uint16_t*)((ptr) + 6)) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            ptr += 16; \
        } else { \
            out_accum.x += *(int16_t*)((ptr)+2); \
            out_accum.y += *(int16_t*)((ptr)+4); \
            out_accum.z += *(int16_t*)((ptr)+6); \
            ptr += 8; \
        } \
    } while(0)
#endif

#if defined(__AVX2__) || defined(__SSE4_1__)
    #define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, b_out, alt_out) do { \
        if ((ptr)[0] & 0x01) { \
            uint16_t xa = *(uint16_t*)((ptr) + 2); \
            uint16_t ya = *(uint16_t*)((ptr) + 4); \
            uint16_t za = *(uint16_t*)((ptr) + 6); \
            c_chg = (xa & 1) | ((ya & 1) << 1) | ((za & 1) << 2); \
            c_chg -= 3; \
            c_occ = (ptr)[0] >> 1; \
            int32_t b_msb = (ptr)[1]; \
            int32_t vx = (((int16_t)xa) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            int32_t vy = (((int16_t)ya) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            int32_t vz = (((int16_t)za) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            out_accum = _mm_set_epi32(vz, vy, vx, 0); \
            uint16_t b_ctrl = *(uint16_t*)((ptr) + 8); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = (b_msb << 14) | (b_ctrl >> 2); \
            ptr += 16; \
        } else { \
            c_chg = 0; \
            uint16_t b_ctrl = *(uint16_t*)(ptr); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = b_ctrl >> 2; \
            out_accum = _mm_add_epi32(out_accum, _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)(ptr)))); \
            ptr += 8; \
        } \
    } while(0)
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, b_out, alt_out) do { \
        if ((ptr)[0] & 0x01) { \
            uint16_t xa = *(uint16_t*)((ptr) + 2); \
            uint16_t ya = *(uint16_t*)((ptr) + 4); \
            uint16_t za = *(uint16_t*)((ptr) + 6); \
            c_chg = (xa & 1) | ((ya & 1) << 1) | ((za & 1) << 2); \
            c_chg -= 3; \
            c_occ = (ptr)[0] >> 1; \
            int32_t b_msb = (ptr)[1]; \
            int32_t vx = (((int16_t)xa) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            int32_t vy = (((int16_t)ya) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            int32_t vz = (((int16_t)za) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            out_accum = (int32x4_t){0, vx, vy, vz}; \
            uint16_t b_ctrl = *(uint16_t*)((ptr) + 8); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = (b_msb << 14) | (b_ctrl >> 2); \
            ptr += 16; \
        } else { \
            c_chg = 0; \
            uint16_t b_ctrl = *(uint16_t*)(ptr); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = b_ctrl >> 2; \
            int32x4_t loaded = vmovl_s16(vreinterpret_s16_u8(vld1_u8(ptr))); \
            loaded = vandq_s32(loaded, (int32x4_t){0, -1, -1, -1}); \
            out_accum = vaddq_s32(out_accum, loaded); \
            ptr += 8; \
        } \
    } while(0)
#else
    #define SIMD_DECODE_MACRO(ptr, out_accum, c_occ, c_chg, b_out, alt_out) do { \
        if ((ptr)[0] & 0x01) { \
            uint16_t xa = *(uint16_t*)((ptr) + 2); \
            uint16_t ya = *(uint16_t*)((ptr) + 4); \
            uint16_t za = *(uint16_t*)((ptr) + 6); \
            c_chg = (xa & 1) | ((ya & 1) << 1) | ((za & 1) << 2); \
            c_chg -= 3; \
            c_occ = (ptr)[0] >> 1; \
            int32_t b_msb = (ptr)[1]; \
            out_accum.x = (((int16_t)xa) >> 1 << 16) | *(uint16_t*)((ptr) + 10); \
            out_accum.y = (((int16_t)ya) >> 1 << 16) | *(uint16_t*)((ptr) + 12); \
            out_accum.z = (((int16_t)za) >> 1 << 16) | *(uint16_t*)((ptr) + 14); \
            uint16_t b_ctrl = *(uint16_t*)((ptr) + 8); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = (b_msb << 14) | (b_ctrl >> 2); \
            ptr += 16; \
        } else { \
            c_chg = 0; \
            uint16_t b_ctrl = *(uint16_t*)(ptr); \
            alt_out = (b_ctrl >> 1) & 0x01; \
            b_out = b_ctrl >> 2; \
            out_accum.x += *(int16_t*)((ptr)+2); \
            out_accum.y += *(int16_t*)((ptr)+4); \
            out_accum.z += *(int16_t*)((ptr)+6); \
            ptr += 8; \
        } \
    } while(0)
#endif

volatile int32_t benchmark_sink = 0;

void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out) {
    if (size < 64) return;
    
    size_t offset = 64;
    uint64_t meta_header = *(uint64_t*)(buffer + offset); offset += 8;
    uint64_t mask = meta_header >> 16;
    int af_mode = (mask & (1ULL << 0)) ? 1 : 0;

    static uint8_t base_natoms[256] = {0};
    static int dict_init = 0;
    if (!dict_init) {
        for (int i = 0; i < NUM_DICT_RES; i++) base_natoms[i+1] = DICT[i].num_atoms;
        dict_init = 1;
    }

    uint8_t r_natoms[256];
    memcpy(r_natoms, base_natoms, 256);

    size_t dict_scan_off = offset;
    if (mask & (1ULL << 1)) dict_scan_off += 8;
    if (mask & (1ULL << 2)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 3)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 4)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 5)) {
        uint64_t dict_len = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8;
        const uint8_t* p_het_dict = buffer + dict_scan_off;
        uint8_t custom_count = p_het_dict[0];
        size_t ptr = 1;
        for (int c = 0; c < custom_count; c++) {
            int u = NUM_DICT_RES + 1 + c; ptr += 4;
            r_natoms[u] = p_het_dict[ptr++]; ptr += r_natoms[u] * 6;
        }
        dict_scan_off += dict_len + ((8 - (dict_len % 8)) % 8);
    }
    if (mask & (1ULL << 6)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 7)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 8)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 9)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 10)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }

    offset = dict_scan_off;

    if (!out) {
        int32_t local_sink = 0;
        while (offset < size) {
            uint16_t num_residues = *(uint16_t*)(buffer + offset);
            if (num_residues == 0xFFFF) break;
            
            uint16_t num_atoms = *(uint16_t*)(buffer + offset + 2);
            offset += 8;
            
            const uint8_t* seq_ptr = buffer + offset;
            int seq_bytes = 0; 
            for (uint32_t r = 0; r < num_residues; r++) if (seq_ptr[seq_bytes++] & 0x80) seq_bytes++;
            
            offset += seq_bytes + ((8 - (seq_bytes & 7)) & 7);
            const uint8_t* atom_ptr = buffer + offset;
            
            accum_t accum;
#if defined(__AVX2__) || defined(__SSE4_1__)
            accum = _mm_setzero_si128();
#elif defined(__aarch64__) || defined(__ARM_NEON)
            accum = vdupq_n_s32(0);
#else
            accum.x = 0; accum.y = 0; accum.z = 0;
#endif
            
            for (uint32_t a = 0; a < num_atoms; a++) {
                if (*(uint64_t*)atom_ptr != 0) {
                    SIMD_DECODE_PHASE1(atom_ptr, accum);
                } else {
                    atom_ptr += 8;
                }
            }
            local_sink += GET_SIMD_X(accum); offset = (size_t)(atom_ptr - buffer);
        }
        benchmark_sink += local_sink; 
        return;
    }

    int day = buffer[5], month = buffer[6], year = 1973 + buffer[7];
    char pdbID[17] = {0}; memcpy(pdbID, buffer + 8, 16); trim_string(pdbID);
    char classification[41] = {0}; memcpy(classification, buffer + 24, 40); trim_string(classification);

    offset = 64 + 8;
    float resolution = -1.0f;
    const char* p_title = ""; uint64_t title_len = 0;
    const char* p_expdta = ""; uint64_t expdta_len = 0;
    const char* p_program = "";
    
    const uint8_t* p_ssbonds = NULL; uint64_t num_ssbonds = 0;
    const uint8_t* p_cispeps = NULL; uint64_t num_cispeps = 0;
    const uint8_t* p_helices = NULL; uint64_t num_helices = 0;
    const uint8_t* p_sheets = NULL;  uint64_t num_sheets = 0;
    const uint8_t* p_sites = NULL;   uint64_t len_sites = 0;

    UniDict* UDICT = (UniDict*)malloc(256 * sizeof(UniDict));
    for (int i = 0; i < NUM_DICT_RES; i++) {
        strcpy(UDICT[i+1].resName, DICT[i].resName); UDICT[i+1].num_atoms = DICT[i].num_atoms;
        for (int a = 0; a < DICT[i].num_atoms; a++) { strcpy(UDICT[i+1].atom_names[a], DICT[i].atom_names[a]); strcpy(UDICT[i+1].elements[a], DICT[i].elements[a]); }
    }

    if (mask & (1ULL << 1)) { resolution = *(float*)(buffer + offset); offset += 8; }
    if (mask & (1ULL << 2)) { title_len = *(uint64_t*)(buffer + offset); offset += 8; p_title = (const char*)(buffer + offset); offset += title_len; offset += (8 - (title_len % 8)) % 8; }
    if (mask & (1ULL << 3)) { expdta_len = *(uint64_t*)(buffer + offset); offset += 8; p_expdta = (const char*)(buffer + offset); offset += expdta_len; offset += (8 - (expdta_len % 8)) % 8; }
    if (mask & (1ULL << 4)) { uint64_t p_len = *(uint64_t*)(buffer + offset); offset += 8; p_program = (const char*)(buffer + offset); offset += p_len; offset += (8 - (p_len % 8)) % 8; }
    if (mask & (1ULL << 5)) {
        uint64_t dict_len = *(uint64_t*)(buffer + offset); offset += 8;
        const uint8_t* p_het_dict = buffer + offset;
        uint8_t custom_count = p_het_dict[0]; size_t ptr = 1;
        for (int c = 0; c < custom_count; c++) {
            int u_idx = NUM_DICT_RES + 1 + c; memcpy(UDICT[u_idx].resName, p_het_dict + ptr, 3); UDICT[u_idx].resName[3] = '\0'; ptr += 4;
            UDICT[u_idx].num_atoms = p_het_dict[ptr]; ptr++;
            for (int a = 0; a < UDICT[u_idx].num_atoms; a++) {
                memcpy(UDICT[u_idx].atom_names[a], p_het_dict + ptr, 4); UDICT[u_idx].atom_names[a][4] = '\0'; ptr += 4;
                memcpy(UDICT[u_idx].elements[a], p_het_dict + ptr, 2); UDICT[u_idx].elements[a][2] = '\0'; ptr += 2;
            }
        }
        offset += dict_len; offset += (8 - (dict_len % 8)) % 8;
    }
    if (mask & (1ULL << 6)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_helices = buffer + offset; num_helices = L / sizeof(MetaHelix); offset += L; offset += (8 - (L % 8)) % 8; }
    if (mask & (1ULL << 7)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_sheets = buffer + offset; num_sheets = L / sizeof(MetaSheet); offset += L; offset += (8 - (L % 8)) % 8; }
    if (mask & (1ULL << 8)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_ssbonds = buffer + offset; num_ssbonds = L / sizeof(MetaSSBond); offset += L; offset += (8 - (L % 8)) % 8; }
    if (mask & (1ULL << 9)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_cispeps = buffer + offset; num_cispeps = L / sizeof(MetaCispep); offset += L; offset += (8 - (L % 8)) % 8; }
    if (mask & (1ULL << 10)) { uint64_t L = *(uint64_t*)(buffer + offset); offset += 8; p_sites = buffer + offset; len_sites = L; offset += L; offset += (8 - (L % 8)) % 8; }

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
        uint16_t num_residues = *(uint16_t*)(buffer + scan_off);
        if (num_residues == 0xFFFF) { scan_off += 8; conect_offset = scan_off; break; }
        uint16_t num_atoms = *(uint16_t*)(buffer + scan_off + 2);
        char c1_raw = buffer[scan_off + 4]; int has_ter = (c1_raw & 0x80) != 0;
        scan_off += 8;
        
        int scan_seq_bytes = 0;
        const uint8_t* scan_seq_ptr = buffer + scan_off;
        for (uint32_t r = 0; r < num_residues; r++) if (scan_seq_ptr[scan_seq_bytes++] & 0x80) scan_seq_bytes++; 
        scan_off += scan_seq_bytes + ((8 - (scan_seq_bytes & 7)) & 7);

        const uint8_t* atom_ptr = buffer + scan_off;
        for (uint32_t a = 0; a < num_atoms; a++) {
            if (*(uint64_t*)(atom_ptr) != 0) out_lines++;
            if ((atom_ptr)[0] & 0x01) atom_ptr += 16; else atom_ptr += 8;
        }
        if (has_ter) out_lines++; 
        scan_off = (size_t)(atom_ptr - buffer);
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

    int atom_serial = 1;
    while (offset < size) {
        uint16_t num_residues = *(uint16_t*)(buffer + offset);
        if (num_residues == 0xFFFF) { break; }
        
        uint16_t num_atoms = *(uint16_t*)(buffer + offset + 2);
        char c1_raw = buffer[offset + 4];
        char c2 = buffer[offset + 5]; 
        char c1 = c1_raw & 0x7F;
        int has_ter = (c1_raw & 0x80) != 0;
        if (c1 == 0) c1 = ' '; 
        if (c2 == 0) c2 = ' ';
        int16_t start_resSeq = *(int16_t*)(buffer + offset + 6);
        offset += 8; 
        
        const uint8_t* seq_ptr = buffer + offset;
        int seq_bytes = 0; uint32_t expected_atoms = 0;
        for (uint32_t r = 0; r < num_residues; r++) { 
            uint8_t type = seq_ptr[seq_bytes++];
            if (type == 0x7E || type == 0x7D) {
                expected_atoms++;
            } else if (type != 0x00 && type != 0x7F) {
                expected_atoms += r_natoms[type & 0x7F];
            }
            if (type & 0x80) seq_bytes++; 
        }
        offset += seq_bytes + ((8 - (seq_bytes & 7)) & 7);
        const uint8_t* atom_ptr = buffer + offset;
        
        int32_t current_resSeq = start_resSeq; int first_res = 1; 
        uint8_t last_valid_type = 0; char iCode = ' ';

        int32_t current_occ = 10000;
        int8_t current_charge = 0;
        char active_alt_char = 'A';

        accum_t accum;
#if defined(__AVX2__) || defined(__SSE4_1__)
        accum = _mm_setzero_si128();  
#elif defined(__aarch64__) || defined(__ARM_NEON)
        accum = vdupq_n_s32(0);  
#else
        accum.x = 0; accum.y = 0; accum.z = 0;  
#endif

        int atoms_left = 0;
        int atom_dict_idx = 0;
        uint32_t parsed_base_atoms = 0;
        int s_idx = 0;
        int no_inc = 0;

        for (uint32_t block_idx = 0; block_idx < num_atoms; block_idx++) {
            if (atoms_left == 0) {
                if (parsed_base_atoms < expected_atoms) {
                    uint8_t meta;
                    int has_5p = 0, has_3p = 0;
                    do {
                        if (s_idx >= seq_bytes) break;
                        meta = seq_ptr[s_idx++];
                        if (meta == 0x00) { current_resSeq++; first_res = 0; no_inc = 0; }
                        else if (meta == 0x7F) { no_inc = 1; }
                        else if (meta == 0x7E) { has_5p = 1; }
                        else if (meta == 0x7D) { has_3p = 1; }
                        else {
                            if (meta & 0x80) iCode = seq_ptr[s_idx++]; else iCode = ' ';
                        }
                    } while (meta == 0x00 || meta == 0x7F || meta == 0x7E || meta == 0x7D);
                    
                    if (s_idx <= seq_bytes) {
                        last_valid_type = meta & 0x7F;
                        if (!first_res && !no_inc) current_resSeq++;
                        first_res = 0;
                        no_inc = 0;
                        atoms_left = r_natoms[last_valid_type] + has_5p + has_3p;
                        atom_dict_idx = has_5p ? -1 : 0;
                    } else {
                        atoms_left = 1; 
                    }
                } else {
                    atoms_left = 1;
                }
                active_alt_char = 'A';
            }

            int i = atom_dict_idx;
            int alt_flag = 0;
            
            if (*(uint64_t*)atom_ptr != 0) {
                int32_t v_b;
                SIMD_DECODE_MACRO(atom_ptr, accum, current_occ, current_charge, v_b, alt_flag);
                
                char altLoc = ' ';
                if (alt_flag || active_alt_char > 'A') {
                    altLoc = active_alt_char++;
                }
                
                const char* a_name;
                const char* e_name;
                
                if (i == -1) {
                    a_name = " OP3"; e_name = " O";
                } else if (i == r_natoms[last_valid_type]) {
                    a_name = " OXT"; e_name = " O";
                } else if (i >= 0 && i < r_natoms[last_valid_type]) {
                    a_name = UDICT[last_valid_type].atom_names[i];
                    e_name = UDICT[last_valid_type].elements[i];
                } else {
                    a_name = " UNK"; e_name = " X"; 
                }

                const char* prefix = (last_valid_type <= 37 || i == -1 || i == r_natoms[last_valid_type]) ? "ATOM  " : "HETATM";
                
                if (af_mode) {
                    int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f  1.00%6.2f          %2s  ",
                            prefix, atom_serial++, a_name, altLoc, UDICT[last_valid_type].resName,
                            c1, c2, current_resSeq, iCode, GET_SIMD_X(accum) / 1000.0f, GET_SIMD_Y(accum) / 1000.0f, GET_SIMD_Z(accum) / 1000.0f, v_b / 100.0f, e_name);
                    memcpy(w_ptr, tbuf, len); w_ptr += 81;
                } else {
                    char chg_str[3] = "  ";
                    if (current_charge > 0) { chg_str[0] = '0' + current_charge; chg_str[1] = '+'; } 
                    else if (current_charge < 0) { chg_str[0] = '0' - current_charge; chg_str[1] = '-'; }

                    int len = snprintf(tbuf, 128, "%s%5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                            prefix, atom_serial++, a_name, altLoc, UDICT[last_valid_type].resName,
                            c1, c2, current_resSeq, iCode, GET_SIMD_X(accum) / 1000.0f, GET_SIMD_Y(accum) / 1000.0f, GET_SIMD_Z(accum) / 1000.0f, current_occ / 100.0f, v_b / 100.0f, e_name, chg_str);
                    memcpy(w_ptr, tbuf, len); w_ptr += 81;
                }
            } else {
                atom_ptr += 8;
            }

            if (alt_flag == 0) {
                atoms_left--;
                atom_dict_idx++;
                parsed_base_atoms++;
                active_alt_char = 'A';
            }
        }

        if (has_ter) {
            int len = snprintf(tbuf, 128, "TER   %5d      %3s%c%c%4d%c", atom_serial++, UDICT[last_valid_type].resName, c1, c2, current_resSeq, iCode);
            memcpy(w_ptr, tbuf, len); w_ptr += 81;
        }
        offset = (size_t)(atom_ptr - buffer);
    }
    
    if (conect_offset > 0) {
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
    
    fwrite(out_data, 1, out_size, out); free(out_data); free(UDICT);
}
