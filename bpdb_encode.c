#include "bpdb.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

static const char* COMPND_KEYS[] = {"", "MOL_ID", "MOLECULE", "CHAIN", "FRAGMENT", "SYNONYM", "EC", "ENGINEERED", "MUTATION", "OTHER_DETAILS"};
static const char* SOURCE_KEYS[] = {"", "MOL_ID", "SYNTHETIC", "FRAGMENT", "ORGANISM_SCIENTIFIC", "ORGANISM_COMMON", "ORGANISM_TAXID",
                                    "STRAIN", "VARIANT", "CELL_LINE", "ATCC", "ORGAN", "TISSUE", "CELL", "ORGANELLE", "SECRETION",
                                    "CELLULAR_LOCATION", "PLASMID", "GENE", "EXPRESSION_SYSTEM", "EXPRESSION_SYSTEM_COMMON",
                                    "EXPRESSION_SYSTEM_TAXID", "EXPRESSION_SYSTEM_STRAIN", "EXPRESSION_SYSTEM_VARIANT",
                                    "EXPRESSION_SYSTEM_CELL_LINE", "EXPRESSION_SYSTEM_ATCC_NUMBER", "EXPRESSION_SYSTEM_ORGAN",
                                    "EXPRESSION_SYSTEM_TISSUE", "EXPRESSION_SYSTEM_CELL", "EXPRESSION_SYSTEM_ORGANELLE",
                                    "EXPRESSION_SYSTEM_CELLULAR_LOCATION", "EXPRESSION_SYSTEM_VECTOR_TYPE", "EXPRESSION_SYSTEM_VECTOR",
                                    "EXPRESSION_SYSTEM_PLASMID", "EXPRESSION_SYSTEM_GENE", "OTHER_DETAILS"};

// --- StreamBuffer: Reusable 64KB dynamic buffers ---
void stream_init(StreamBuffer* buf) {
  buf->cap = 65536; // Start small (64KB), avoids 16MB upfront allocation
  buf->data = (uint8_t*)malloc(buf->cap);
  buf->len = 0;
}
void stream_reset(StreamBuffer* buf) { buf->len = 0; } // Zero-alloc reset
void stream_write(StreamBuffer* buf, const void* src, size_t size) {
  if (buf->len + size > buf->cap) {
    buf->cap = (buf->cap == 0) ? 65536 : (buf->cap * 2);
    if (buf->len + size > buf->cap) buf->cap = buf->len + size;
    buf->data = (uint8_t*)realloc(buf->data, buf->cap);
  }
  memcpy(buf->data + buf->len, src, size);
  buf->len += size;
}
void stream_free(StreamBuffer* buf) { free(buf->data); }

int parse_month(const char* m) { for (int i = 1; i <= 12; i++) if (strncmp(m, MONTHS[i], 3) == 0) return i; return 1; }

void trim_string(char* str) {
  if (!str) return;
  int start = 0;
  while (str[start] && isspace((unsigned char)str[start])) start++;
  int end = strlen(str) - 1;
  while (end >= start && isspace((unsigned char)str[end])) end--;
  int len = end - start + 1;
  if (len > 0) {
    memmove(str, str + start, len);
    str[len] = '\0';
  } else {
    str[0] = '\0';
  }
}

// --- Strict Column Parsing Helpers (Bounds-checked for MMAP) ---
static void ext_str(char* dst, const char* line, int line_len, int col, int width) {
  int start = col - 1;
  if (start >= line_len) { dst[0] = '\0'; return; }
  int cp = width;
  if (start + cp > line_len) cp = line_len - start;
  while(cp > 0 && (line[start+cp-1] == '\n' || line[start+cp-1] == '\r')) cp--;
  if (cp > 0) { memcpy(dst, line+start, cp); dst[cp] = '\0'; } else dst[0] = '\0';
  trim_string(dst);
}

static void ext_str_rtrim(char* dst, const char* line, int line_len, int col, int width) {
  int start = col - 1;
  if (start >= line_len) { dst[0] = '\0'; return; }
  int cp = width;
  if (start + cp > line_len) cp = line_len - start;
  while(cp > 0 && isspace((unsigned char)line[start+cp-1])) cp--;
  if (cp > 0) { memcpy(dst, line+start, cp); dst[cp] = '\0'; } else dst[0] = '\0';
}

static int ext_int(const char* line, int line_len, int col, int width) {
  char buf[32]; ext_str(buf, line, line_len, col, width);
  return strlen(buf) > 0 ? atoi(buf) : 0;
}

static float ext_float(const char* line, int line_len, int col, int width) {
  char buf[32]; ext_str(buf, line, line_len, col, width);
  return strlen(buf) > 0 ? atof(buf) : 0.0f;
}

static char ext_char(const char* line, int line_len, int col) {
  int start = col - 1;
  if (start >= line_len) return ' ';
  return line[start] == '\n' || line[start] == '\r' ? ' ' : line[start];
}

// --- Zero-alloc integer parsing (Replaces atof) ---
static inline int32_t fast_parse_coord(const char* ptr) {
  int32_t val = 0, sign = 1;
  for (int i = 0; i < 8; i++) {
    char c = ptr[i];
    if (c >= '0' && c <= '9') val = val * 10 + (c - '0');
    else if (c == '-') sign = -1;
  }
  return sign * val;
}

static inline int32_t fast_parse_f62(const char* ptr) {
  int32_t val = 0, sign = 1;
  for (int i = 0; i < 6; i++) {
    char c = ptr[i];
    if (c >= '0' && c <= '9') val = val * 10 + (c - '0');
    else if (c == '-') sign = -1;
  }
  return sign * val;
}

// --- Lookup Table Hashing ---
static inline uint32_t res_hash(const char* res) {
  uint32_t c0 = res[0], c1 = res[1], c2 = res[2];
  if (c0 < 32 || c0 > 95) c0 = 32;
  if (c1 < 32 || c1 > 95) c1 = 32;
  if (c2 < 32 || c2 > 95) c2 = 32;
  return ((c0 - 32) << 12) | ((c1 - 32) << 6) | (c2 - 32);
}

// --- Zero-alloc AltArray ---
#define MAX_ALTS 16
typedef struct {
  int32_t b[MAX_ALTS], x[MAX_ALTS], y[MAX_ALTS], z[MAX_ALTS], occ[MAX_ALTS];
  int8_t charge[MAX_ALTS]; char alt[MAX_ALTS];
  int count;
} AltArray;

void alt_init(AltArray* a) { a->count = 0; }
void alt_push(AltArray* a, int32_t b, int32_t x, int32_t y, int32_t z, int32_t occ, int8_t charge, char alt) {
  if (a->count < MAX_ALTS) {
    int i = a->count++;
    a->b[i] = b; a->x[i] = x; a->y[i] = y; a->z[i] = z; a->occ[i] = occ; a->charge[i] = charge; a->alt[i] = alt;
  }
}

typedef struct {
  int resSeq; char iCode; int type;
  AltArray h_atoms[256]; int h_present[256];
  AltArray hy_atoms[256]; int hy_present[256];
} EncoderState;

typedef struct {
  uint8_t* seq_track; uint8_t* atom_track;
  int seq_bytes; int atom_bytes; int num_atoms; int num_residues;
  int32_t current_x, current_y, current_z, current_b;
  int start_resSeq; int last_resSeq; char last_iCode;
  char c1, c2; int has_ter, has_ter_raw;
  int is_n_term, is_c_term;
  int last_type;
} ChainBuilder;

void cb_init(ChainBuilder* cb) {
  cb->seq_track = (uint8_t*)malloc(240000); cb->atom_track = (uint8_t*)malloc(120000 * 8 + 100);
  cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0;
  cb->start_resSeq = 0; cb->has_ter = 0; cb->has_ter_raw = 0;
  cb->is_n_term = 0; cb->is_c_term = 0;
  cb->current_b = 0; cb->last_type = 0;
}
void cb_reset(ChainBuilder* cb) {
  cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0;
  cb->has_ter = 0; cb->has_ter_raw = 0; cb->is_n_term = 0; cb->is_c_term = 0; cb->last_type = 0;
}
void cb_free(ChainBuilder* cb) { free(cb->seq_track); free(cb->atom_track); }

void flush_chain(ChainBuilder* cb, StreamBuffer* buf) {
  if (cb->num_atoms == 0) { cb->has_ter = 0; cb->has_ter_raw = 0; cb->is_n_term = 0; cb->is_c_term = 0; return; }
  uint8_t header[8] = {0};
  *(uint32_t*)(header + 0) = (cb->num_residues << 18) | (cb->num_atoms & 0x3FFFF);
  header[4] = (cb->c1 & 0x7F) | (cb->is_n_term ? 0x80 : 0);
  header[5] = (cb->c2 & 0x7F) | (cb->is_c_term ? 0x80 : 0);
  *(int16_t*)(header + 6) = (cb->start_resSeq << 1) | (cb->has_ter ? 1 : 0);

  stream_write(buf, header, 8);
  stream_write(buf, cb->seq_track, cb->seq_bytes);
  uint32_t pad = (8 - (cb->seq_bytes & 7)) & 7;
  if (pad > 0) { uint64_t zero = 0; stream_write(buf, &zero, pad); }
  stream_write(buf, cb->atom_track, cb->atom_bytes);
  cb_reset(cb);
}

static inline void split_msb_lsb(int32_t val, int16_t* msb, int16_t* lsb) {
  int32_t base = (val + 4096) >> 13; *msb = (int16_t)base; *lsb = (int16_t)(val - (base << 13));
}

static inline void encode_alt_array(ChainBuilder* cb, AltArray* arr, int force_jump, int af_mode) {
  for (int a = 0; a < arr->count; a++) {
    int32_t x = arr->x[a], y = arr->y[a], z = arr->z[a], b_int = arr->b[a];
    
    if (b_int > 4194303) {
      b_int = 4194303;
    }
    if (b_int < 0) {
      b_int = 0;
    }
    
    if (cb->num_atoms == 0) { cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int; }
    int32_t dx = x - cb->current_x, dy = y - cb->current_y, dz = z - cb->current_z, db = b_int - cb->current_b;
    int alt_flag = (a < arr->count - 1) ? 1 : 0;
    int32_t occ_int = arr->occ[a]; if (occ_int > 100) occ_int = 100; else if (occ_int < 0) occ_int = 0;
    uint8_t charge_u = (arr->charge[a] + 3) & 0x07;
    int limit = af_mode ? 16383 : 4095;
    int needs_jump = (cb->num_atoms == 0) || (force_jump && a == 0) ||
      (dx > limit || dx < -(limit+1) || dy > limit || dy < -(limit+1) || dz > limit || dz < -(limit+1) || db > limit || db < -(limit+1));

    uint8_t* ptr = &cb->atom_track[cb->atom_bytes];
    if (needs_jump) {
      int16_t b_msb, b_lsb, x_msb, x_lsb, y_msb, y_lsb, z_msb, z_lsb;
      split_msb_lsb(b_int, &b_msb, &b_lsb); split_msb_lsb(x, &x_msb, &x_lsb); split_msb_lsb(y, &y_msb, &y_lsb); split_msb_lsb(z, &z_msb, &z_lsb);
      *(uint16_t*)(ptr)     = (b_msb << 1) | 1; *(uint16_t*)(ptr + 2) = (x_msb << 1) | 0; *(uint16_t*)(ptr + 4) = (y_msb << 1) | 0; *(uint16_t*)(ptr + 6) = (z_msb << 1) | 0;
      if (af_mode) {
        *(uint16_t*)(ptr + 8) = (b_lsb << 1) | 0; *(uint16_t*)(ptr + 10) = (x_lsb << 1) | 0; *(uint16_t*)(ptr + 12) = (y_lsb << 1) | 0; *(uint16_t*)(ptr + 14) = (z_lsb << 1) | 0;
      } else {
        *(uint16_t*)(ptr + 8)  = (b_lsb << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
        *(uint16_t*)(ptr + 10) = (x_lsb << 3) | (((occ_int >> 4) & 3) << 1) | ((charge_u >> 2) & 1);
        *(uint16_t*)(ptr + 12) = (y_lsb << 3) | (((occ_int >> 2) & 3) << 1) | ((charge_u >> 1) & 1);
        *(uint16_t*)(ptr + 14) = (z_lsb << 3) | ((occ_int & 3) << 1) | (charge_u & 1);
      }
      cb->atom_bytes += 16;
    } else {
      if (af_mode) {
        *(uint16_t*)(ptr) = (db << 1) | 0; *(uint16_t*)(ptr + 2) = (dx << 1) | 0; *(uint16_t*)(ptr + 4) = (dy << 1) | 0; *(uint16_t*)(ptr + 6) = (dz << 1) | 0;
      } else {
        *(uint16_t*)(ptr)     = (db << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
        *(uint16_t*)(ptr + 2) = (dx << 3) | (((occ_int >> 4) & 3) << 1) | ((charge_u >> 2) & 1);
        *(uint16_t*)(ptr + 4) = (dy << 3) | (((occ_int >> 2) & 3) << 1) | ((charge_u >> 1) & 1);
        *(uint16_t*)(ptr + 6) = (dz << 3) | ((occ_int & 3) << 1) | (charge_u & 1);
      }
      cb->atom_bytes += 8;
    }
    cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int; cb->num_atoms++;
  }
}

void flush_residue(ChainBuilder* cb, EncoderState* st, ResidueDict* UDICT, int is_first, int is_last, int include_hydrogens, int af_mode) {
  if (st->type <= 0) return;
  if (cb->num_residues == 0) cb->start_resSeq = st->resSeq;

  if (cb->num_residues > 0) {
    int diff = st->resSeq - cb->last_resSeq;
    if (diff == 0) { cb->seq_track[cb->seq_bytes++] = 0x7F; cb->num_residues++; }
    else if (diff > 1 && diff <= 100) { for (int g = 1; g < diff; g++) { cb->seq_track[cb->seq_bytes++] = 0x00; cb->num_residues++; } }
  }

  if (st->iCode != ' ' && st->iCode != 0) { cb->seq_track[cb->seq_bytes++] = 0x7C; cb->seq_track[cb->seq_bytes++] = st->iCode; cb->num_residues++; }
  cb->seq_track[cb->seq_bytes++] = st->type; cb->last_type = st->type; cb->num_residues++;
  
  int is_hetatm = (st->type >= 0x60);
  int h_start = is_first ? 0 : UDICT[st->type].h_start_norm;
  int h_end   = is_last  ? UDICT[st->type].h_end_cterm : UDICT[st->type].h_end_norm;

  for (int i = h_start; i < h_end; i++) {
    if (!st->h_present[i] || st->h_atoms[i].count == 0) { memset(&cb->atom_track[cb->atom_bytes], 0, 8); cb->atom_bytes += 8; cb->num_atoms++; continue; }
    encode_alt_array(cb, &st->h_atoms[i], (is_hetatm && i == h_start), af_mode);
  }

  if (include_hydrogens) {
    int hy_start = is_first ? 0 : UDICT[st->type].hy_start_norm;
    int hy_end   = is_last  ? UDICT[st->type].hy_end_cterm : UDICT[st->type].hy_end_norm;
    for (int i = hy_start; i < hy_end; i++) {
      if (!st->hy_present[i] || st->hy_atoms[i].count == 0) { memset(&cb->atom_track[cb->atom_bytes], 0, 8); cb->atom_bytes += 8; cb->num_atoms++; continue; }
      encode_alt_array(cb, &st->hy_atoms[i], 0, af_mode);
    }
  }

  if (!is_hetatm) {
    int anchor_idx = (st->type >= 0x20 && st->type < 0x60) ? 2 : 8;
    if (anchor_idx >= h_start && anchor_idx < h_end && st->h_present[anchor_idx]) {
      cb->current_x = st->h_atoms[anchor_idx].x[0]; cb->current_y = st->h_atoms[anchor_idx].y[0];
      cb->current_z = st->h_atoms[anchor_idx].z[0]; cb->current_b = st->h_atoms[anchor_idx].b[0];
    }
  }
  cb->last_resSeq = st->resSeq; cb->last_iCode = st->iCode;
}

// --- Pre-allocated Global Encoder Context ---
typedef struct {
  char title[4096];
  char pdbID[13];
  char classification[128];
  char expdta[128];
  char program[56];
  char raw_compnd[16384];
  char raw_source[16384];

  float resolution;
  int day, month, year;
  int id_from_header;
  int parsed_num_models;
  int counted_models;
  int custom_het_idx;

  char current_seqres_chain;
  uint16_t current_seqres_count;
  int has_seqres_chain[256];

  StreamBuffer seqres_buf, ssbonds, cispeps, helices, sheets, sites, conects, dbrefs;
} EncodeMetadata;

typedef struct {
  StreamBuffer coords;
  EncoderState st;
  ChainBuilder cb;
  char current_c1, current_c2;
  uint16_t current_model_idx;
  uint16_t total_models_minus_1;
} CoordinateContext;

typedef struct {
  ResidueDict UDICT[256];
  uint8_t res_lut[262144];
  struct { uint8_t map[8192]; } h_hash[256], hy_hash[256];

  StreamBuffer meta_buffer;

  EncodeMetadata m;
  CoordinateContext cctx;
  int initialized;
} BpdbEncoderContext;

static BpdbEncoderContext gctx;

static void init_gctx() {
  if (gctx.initialized) return;
  memset(&gctx, 0, sizeof(BpdbEncoderContext));
  for (size_t i = 0; i < NUM_DICT_RES; i++) {
    gctx.UDICT[DICT[i].type_id] = DICT[i];
    uint32_t r_hash = res_hash(DICT[i].resName);
    gctx.res_lut[r_hash] = DICT[i].type_id;

    memset(gctx.h_hash[DICT[i].type_id].map, 0xFF, 8192);
    for (int a = 0; a < DICT[i].h_end_cterm; a++) {
      uint32_t val = *(uint32_t*)(DICT[i].h_atoms[a]);
      uint32_t h = (val ^ (val >> 13) ^ (val >> 21)) & 0x1FFF;
      gctx.h_hash[DICT[i].type_id].map[h] = a;
    }

    memset(gctx.hy_hash[DICT[i].type_id].map, 0xFF, 8192);
    for (int a = 0; a < DICT[i].hy_end_cterm; a++) {
      uint32_t val = *(uint32_t*)(DICT[i].hy_atoms[a]);
      uint32_t h = (val ^ (val >> 13) ^ (val >> 21)) & 0x1FFF;
      gctx.hy_hash[DICT[i].type_id].map[h] = a;
    }
  }
  stream_init(&gctx.m.seqres_buf); stream_init(&gctx.m.ssbonds); stream_init(&gctx.m.cispeps);
  stream_init(&gctx.m.helices); stream_init(&gctx.m.sheets); stream_init(&gctx.m.sites);
  stream_init(&gctx.m.conects); stream_init(&gctx.m.dbrefs); 
  
  stream_init(&gctx.meta_buffer);
  stream_init(&gctx.cctx.coords);
  cb_init(&gctx.cctx.cb);
  gctx.initialized = 1;
}

static void reset_gctx_for_file() {
  stream_reset(&gctx.m.seqres_buf); stream_reset(&gctx.m.ssbonds); stream_reset(&gctx.m.cispeps);
  stream_reset(&gctx.m.helices); stream_reset(&gctx.m.sheets); stream_reset(&gctx.m.sites);
  stream_reset(&gctx.m.conects); stream_reset(&gctx.m.dbrefs); 
  
  // Keep meta_buffer reset on the outer struct
  stream_reset(&gctx.meta_buffer);
  stream_reset(&gctx.cctx.coords);
  cb_reset(&gctx.cctx.cb);
  
  memset(&gctx.m, 0, sizeof(EncodeMetadata));
  gctx.m.day = 1; gctx.m.month = 1; gctx.m.year = 1973;
  gctx.m.custom_het_idx = 0x81; gctx.m.resolution = -1.0f;
  
  gctx.cctx.current_model_idx = 0;
  gctx.cctx.total_models_minus_1 = 0;
  gctx.cctx.current_c1 = 0; gctx.cctx.current_c2 = 0;
  
  gctx.cctx.st.type = -1;
  memset(gctx.cctx.st.h_present, 0, sizeof(gctx.cctx.st.h_present));
  memset(gctx.cctx.st.hy_present, 0, sizeof(gctx.cctx.st.hy_present));
  for (int k = 0; k < 256; k++) { alt_init(&gctx.cctx.st.h_atoms[k]); alt_init(&gctx.cctx.st.hy_atoms[k]); }
}

// --- Binary Builders ---
static inline uint32_t pack_ec_number(const char* ec_str) {
  int class_val = 0, sub = 0, subsub = 0, serial = 0;
  char temp[32]; strncpy(temp, ec_str, 31); temp[31] = '\0';
  char* t = temp; char* p = strchr(t, '.');
  if (p) { *p = '\0'; class_val = atoi(t); t = p + 1; }
  p = strchr(t, '.'); if (p) { *p = '\0'; sub = atoi(t); t = p + 1; }
  p = strchr(t, '.'); if (p) { *p = '\0'; if (t[0] == '-') subsub = 0; else if (t[0] == 'n') subsub = -atoi(t+1); else subsub = atoi(t); t = p + 1; }
  if (t[0] == '-') serial = 0; else if (t[0] == 'n') serial = -atoi(t+1); else serial = atoi(t);
  uint32_t packed = 0; packed |= ((class_val & 0x07) << 29); packed |= ((sub & 0x1F) << 24); packed |= (((uint32_t)((int8_t)subsub) & 0xFF) << 16); packed |= ((uint32_t)((int16_t)serial) & 0xFFFF);
  return packed;
}

void extract_tlv(char* raw_text, StreamBuffer* out_buf, int is_source, int af_mode) {
  char* tok = strtok(raw_text, ";");
  StreamBuffer mol_buf; stream_init(&mol_buf);
  while (tok) {
    char* colon = strchr(tok, ':');
    if (colon) {
      *colon = '\0'; char* key = tok; char* val = colon + 1;
      trim_string(key); trim_string(val);
      uint8_t token_id = 0; const char** keys = is_source ? SOURCE_KEYS : COMPND_KEYS; int max_keys = is_source ? 36 : 10;
      for (int i = 1; i < max_keys; i++) { if (strcmp(key, keys[i]) == 0) { token_id = i; break; } }
      if (token_id == 1) {
        if (mol_buf.len > 0) { uint16_t mlen = mol_buf.len; stream_write(out_buf, &mlen, 2); stream_write(out_buf, mol_buf.data, mol_buf.len); stream_reset(&mol_buf); }
        uint8_t vlen = strlen(val); stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, val, vlen);
      } else if (token_id > 1) {
        if (!is_source && token_id == COMPND_EC) { uint32_t packed_ec = pack_ec_number(val); uint8_t vlen = 4; stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, &packed_ec, 4); }
        else if (is_source && af_mode && (token_id == SOURCE_ORGANISM_SCIENTIFIC || token_id == SOURCE_ORGANISM_COMMON || token_id == SOURCE_EXPRESSION_SYSTEM || token_id == SOURCE_EXPRESSION_SYSTEM_COMMON)) continue;
        else if (is_source && (token_id == SOURCE_ORGANISM_TAXID || token_id == SOURCE_EXPRESSION_SYSTEM_TAXID)) { uint32_t taxid = atoi(val); uint8_t vlen = 4; stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, &taxid, 4); }
        else { uint8_t vlen = strlen(val); stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, val, vlen); }
      }
    } tok = strtok(NULL, ";");
  }
  if (mol_buf.len > 0) { uint16_t mlen = mol_buf.len; stream_write(out_buf, &mlen, 2); stream_write(out_buf, mol_buf.data, mol_buf.len); }
  stream_free(&mol_buf);
}

static uint64_t build_binary_headers(EncodeMetadata* m, StreamBuffer* meta_buffer, int include_hydrogens, int af_mode, int afdb_mode, ResidueDict* UDICT) {
  StreamBuffer compnd_tlv; stream_init(&compnd_tlv);
  if (strlen(m->raw_compnd) > 0) extract_tlv(m->raw_compnd, &compnd_tlv, 0, af_mode);

  StreamBuffer source_tlv; stream_init(&source_tlv);
  if (strlen(m->raw_source) > 0) extract_tlv(m->raw_source, &source_tlv, 1, af_mode);

  // Restore these required variable declarations:
  uint64_t title_len = strlen(m->title), expdta_len = strlen(m->expdta);
  int has_title = (title_len > 0), has_het = (m->custom_het_idx > 0x81), has_expdta = (expdta_len > 0);
  int has_res = (m->resolution > 0.0f), has_prog = (strlen(m->program) > 0);

  uint64_t meta_mask = 0; uint64_t zero = 0;
  if (af_mode) meta_mask |= (1ULL << META_AF_MODE);
  if (has_res) meta_mask |= (1ULL << META_RESOLUTION);
  if (has_title) meta_mask |= (1ULL << META_TITLE);
  if (has_expdta) meta_mask |= (1ULL << META_EXPDTA);
  if (has_prog) meta_mask |= (1ULL << META_PROGRAM);
  if (has_het) meta_mask |= (1ULL << META_HET_DICT);
  if (m->helices.len > 0) meta_mask |= (1ULL << META_HELIX);
  if (m->sheets.len > 0) meta_mask |= (1ULL << META_SHEET);
  if (m->ssbonds.len > 0) meta_mask |= (1ULL << META_SSBOND);
  if (m->cispeps.len > 0) meta_mask |= (1ULL << META_CISPEP);
  if (m->sites.len > 0) meta_mask |= (1ULL << META_SITE);
  if (include_hydrogens) meta_mask |= (1ULL << META_HYDROGENS);
  if (compnd_tlv.len > 0) meta_mask |= (1ULL << META_COMPND);
  if (source_tlv.len > 0) meta_mask |= (1ULL << META_SOURCE);
  if (m->dbrefs.len > 0) meta_mask |= (1ULL << META_DBREF);
  if (afdb_mode) meta_mask |= (1ULL << META_AFDB);
  if (m->seqres_buf.len > 0 && !afdb_mode) meta_mask |= (1ULL << META_SEQRES);

  stream_write(meta_buffer, &meta_mask, 8);

  // Section 2: Title Section (16-21)
  if (meta_mask & (1ULL << META_TITLE)) { stream_write(meta_buffer, &title_len, 8); stream_write(meta_buffer, m->title, title_len); size_t pad = (8 - (title_len % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_COMPND)) { uint64_t L = compnd_tlv.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, compnd_tlv.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SOURCE)) { uint64_t L = source_tlv.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, source_tlv.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_EXPDTA)) { stream_write(meta_buffer, &expdta_len, 8); stream_write(meta_buffer, m->expdta, expdta_len); size_t pad = (8 - (expdta_len % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_RESOLUTION)) { uint32_t pad=0; stream_write(meta_buffer, &m->resolution, 4); stream_write(meta_buffer, &pad, 4); }
  if (meta_mask & (1ULL << META_PROGRAM)) { uint64_t p_len = 56; stream_write(meta_buffer, &p_len, 8); stream_write(meta_buffer, m->program, 56); }

  // Section 3: Primary Structure Section (32-33)
  if (meta_mask & (1ULL << META_DBREF))  { uint64_t L = m->dbrefs.len;     stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->dbrefs.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SEQRES)) { uint64_t L = m->seqres_buf.len; stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->seqres_buf.data, L); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }

  // Section 4: Heterogen Section (36)
  if (meta_mask & (1ULL << META_HET_DICT)) {
    uint64_t dict_len = 1; for (int i = 0x81; i < m->custom_het_idx; i++) dict_len += 5 + UDICT[i].h_end_norm * 6;
    stream_write(meta_buffer, &dict_len, 8);
    uint8_t custom_count = m->custom_het_idx - 0x81; stream_write(meta_buffer, &custom_count, 1);
    for (int i = 0x81; i < m->custom_het_idx; i++) {
      char r_name[4] = {0}; memcpy(r_name, UDICT[i].resName, 3); r_name[3] = ' '; stream_write(meta_buffer, r_name, 4);
      uint8_t a_count = UDICT[i].h_end_norm; stream_write(meta_buffer, &a_count, 1);
      for (int a = 0; a < a_count; a++) { stream_write(meta_buffer, UDICT[i].h_atoms[a], 4); stream_write(meta_buffer, UDICT[i].h_elements[a], 2); }
    }
    size_t pad = (8 - (dict_len % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad);
  }

  // Section 5: Secondary Structure Section (40-41)
  if (meta_mask & (1ULL << META_HELIX))  { uint64_t L = m->helices.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->helices.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SHEET))  { uint64_t L = m->sheets.len;     stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->sheets.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }

  // Section 6: Connectivity Annotation Section (44-45)
  if (meta_mask & (1ULL << META_SSBOND)) { uint64_t L = m->ssbonds.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->ssbonds.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_CISPEP)) { uint64_t L = m->cispeps.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->cispeps.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }

  // Section 7: Miscellaneous Features Section (48)
  if (meta_mask & (1ULL << META_SITE))   { uint64_t L = m->sites.len;      stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->sites.data, L);      size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }stream_free(&compnd_tlv); stream_free(&source_tlv);
  return meta_mask;
}

// --- Coordinate Handlers ---
static void process_ter(CoordinateContext* ctx, EncodeMetadata* m, int include_hydrogens, int af_mode) {
  ctx->cb.has_ter_raw = 1;
  if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
  if (ctx->cb.num_atoms > 0) {
    unsigned char ch = ctx->cb.c1; 
    ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
    ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
    flush_chain(&ctx->cb, &ctx->coords);
  } else { ctx->cb.has_ter_raw = 0; }
    
  ctx->st.type = -1; 
  memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
  for (int k = 0; k < 256; k++) { alt_init(&ctx->st.h_atoms[k]); alt_init(&ctx->st.hy_atoms[k]); }
}

static void save_model_file(const char* base_out, uint16_t total_m1, uint16_t cur_m1, const char* gid, const char* gclass, int day, int month, int year, StreamBuffer* meta, StreamBuffer* coords, StreamBuffer* conects) {
  char filename[1024];
  if (strcmp(base_out, "-") == 0) strcpy(filename, "-");
  else if (total_m1 > 0) {
    const char* ext = strrchr(base_out, '.');
    if (ext) { int baselen = ext - base_out; snprintf(filename, sizeof(filename), "%.*s_m%d%s", baselen, base_out, cur_m1 + 1, ext); }
    else snprintf(filename, sizeof(filename), "%s_m%d", base_out, cur_m1 + 1);
  } else strcpy(filename, base_out);

  FILE* fout = (strcmp(filename, "-") == 0) ? stdout : fopen(filename, "wb");
  if (!fout) { perror("Failed to open output"); exit(1); }

  uint8_t header[64] = {0};
  memcpy(header, "bPDB", 4); header[4] = 0; header[5] = (uint8_t)day; header[6] = (uint8_t)month;
  header[7] = (year >= 1973 && year <= 2228) ? (uint8_t)(year - 1973) : 0;
  memcpy(header + 8, gid, strlen(gid) > 12 ? 12 : strlen(gid));
  *(uint16_t*)(header + 20) = total_m1; *(uint16_t*)(header + 22) = cur_m1;
  memcpy(header + 24, gclass, strlen(gclass) > 40 ? 40 : strlen(gclass));
  fwrite(header, 1, 64, fout);
  fwrite(meta->data, 1, meta->len, fout);
  fwrite(coords->data, 1, coords->len, fout);
  if (conects->len > 0) fwrite(conects->data, 1, conects->len, fout);

  if (fout != stdout) {
    uint64_t zero = 0; size_t current_size = ftell(fout);
    size_t eof_pad = (8 - (current_size % 8)) % 8;
    if (eof_pad > 0) fwrite(&zero, 1, eof_pad, fout);
    fclose(fout);
  }
}

static void flush_final_to_disk(CoordinateContext* ctx, const char* out_bpdb, EncodeMetadata* m, int include_hydrogens, int af_mode, int afdb_mode) {
  stream_reset(&gctx.meta_buffer);
  build_binary_headers(m, &gctx.meta_buffer, include_hydrogens, af_mode, afdb_mode, gctx.UDICT);
  save_model_file(out_bpdb, ctx->total_models_minus_1, ctx->current_model_idx, 
                  m->pdbID, m->classification, m->day, m->month, m->year, 
                  &gctx.meta_buffer, &ctx->coords, &gctx.m.conects);
}

static void process_endmdl(CoordinateContext* ctx, const char* out_bpdb, EncodeMetadata* m, int include_hydrogens, int af_mode, int afdb_mode) {
  if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
  if (ctx->cb.num_atoms > 0) {
    unsigned char ch = ctx->cb.c1; 
    ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
    ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
    flush_chain(&ctx->cb, &ctx->coords);
  }
  uint64_t end_stream = (0x3FFFULL << 18); stream_write(&ctx->coords, &end_stream, 8);
  flush_final_to_disk(ctx, out_bpdb, m, include_hydrogens, af_mode, afdb_mode);
    
  stream_reset(&ctx->coords);
  ctx->current_model_idx++;
  ctx->st.type = -1; 
  memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
  for (int k = 0; k < 256; k++) { alt_init(&ctx->st.h_atoms[k]); alt_init(&ctx->st.hy_atoms[k]); }
}

static void process_file_boundary(CoordinateContext* ctx, EncodeMetadata* m, int include_hydrogens, int af_mode) {
  if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
  if (ctx->cb.num_atoms > 0) {
    unsigned char ch = ctx->cb.c1; 
    ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
    ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
    flush_chain(&ctx->cb, &ctx->coords);
  }
  uint64_t frag_stream = (0x3FFEULL << 18); stream_write(&ctx->coords, &frag_stream, 8);
  ctx->st.type = -1; 
  memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
  for (int k = 0; k < 256; k++) { alt_init(&ctx->st.h_atoms[k]); alt_init(&ctx->st.hy_atoms[k]); }
}

static void process_eof_flush(CoordinateContext* ctx, const char* out_bpdb, EncodeMetadata* m, int include_hydrogens, int af_mode, int afdb_mode) {
  if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
  if (ctx->cb.num_atoms > 0) {
    unsigned char ch = ctx->cb.c1; 
    ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
    ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
    flush_chain(&ctx->cb, &ctx->coords);
  }
  uint64_t end_stream = (0x3FFFULL << 18); stream_write(&ctx->coords, &end_stream, 8);
  flush_final_to_disk(ctx, out_bpdb, m, include_hydrogens, af_mode, afdb_mode);
}

// --- $O(1)$ Atom Lookups in Main Stream ---
static void process_atom(const char* line, int line_len, CoordinateContext* ctx, EncodeMetadata* m, int include_hydrogens, int af_mode) {
  if (line_len < 54) return;
  char altLoc = line[16], iCode = line[26], c1 = line[20], c2 = line[21];
  int resSeq = ext_int(line, line_len, 23, 4);

  uint32_t res_hash_val = res_hash(line + 17);
  int type = gctx.res_lut[res_hash_val];
    
  // Auto-vivify Custom HETATM
  if (type == 0) {
    if (!include_hydrogens && line[13] == 'H') return;
    if (m->custom_het_idx < 192) {
      type = m->custom_het_idx++;
      memcpy(gctx.UDICT[type].resName, line + 17, 3); gctx.UDICT[type].resName[3] = '\0';
      gctx.res_lut[res_hash_val] = type;
      gctx.UDICT[type].h_end_norm = 0; gctx.UDICT[type].h_end_cterm = 0;
      memset(gctx.h_hash[type].map, 0xFF, 8192);
    } else return;
  }

  int h_idx = -1, hy_idx = -1;
  uint32_t raw_atom_int = *(uint32_t*)(line + 12);
  uint32_t h = (raw_atom_int ^ (raw_atom_int >> 13) ^ (raw_atom_int >> 21)) & 0x1FFF;

  int a_idx = gctx.h_hash[type].map[h];
  if (a_idx != 0xFF && *(uint32_t*)(gctx.UDICT[type].h_atoms[a_idx]) == raw_atom_int) {
    h_idx = a_idx;
  } else if (include_hydrogens) {
    a_idx = gctx.hy_hash[type].map[h];
    if (a_idx != 0xFF && *(uint32_t*)(gctx.UDICT[type].hy_atoms[a_idx]) == raw_atom_int) hy_idx = a_idx;
  }

  if (h_idx == -1 && hy_idx == -1) {
    if (type >= 0x80 && gctx.UDICT[type].h_end_norm < 255) {
      h_idx = gctx.UDICT[type].h_end_norm++; gctx.UDICT[type].h_end_cterm = gctx.UDICT[type].h_end_norm;
      *(uint32_t*)(gctx.UDICT[type].h_atoms[h_idx]) = raw_atom_int; gctx.UDICT[type].h_atoms[h_idx][4] = '\0';
      gctx.UDICT[type].h_elements[h_idx][0] = line_len >= 78 ? line[76] : ' ';
      gctx.UDICT[type].h_elements[h_idx][1] = line_len >= 78 ? line[77] : ' ';
      gctx.UDICT[type].h_elements[h_idx][2] = '\0';
      gctx.h_hash[type].map[h] = h_idx;
    } else return;
  }

  int new_chain = (c1 != ctx->current_c1 || c2 != ctx->current_c2);
  int32_t x_int = fast_parse_coord(line + 30); int32_t y_int = fast_parse_coord(line + 38); int32_t z_int = fast_parse_coord(line + 46);
  int32_t temp_int = fast_parse_f62(line + 60);

  int32_t occ_int = 100; int charge = 0;
  if (!af_mode && line_len >= 60) {
    occ_int = fast_parse_f62(line + 54);
    if (line_len >= 79) {
      char chg = line[78];
      if (chg >= '0' && chg <= '9') { charge = chg - '0'; if (line_len >= 80 && line[79] == '-') charge = -charge; }
      if (charge < -3) charge = -3; else if (charge > 4) charge = 4;
    }
  }

  if (resSeq != ctx->st.resSeq || iCode != ctx->st.iCode || new_chain || ctx->cb.atom_bytes >= 580000 || type != ctx->st.type) {
    int has_oxt = 0;
    if (ctx->st.type > 0 && gctx.UDICT[ctx->st.type].h_end_cterm > gctx.UDICT[ctx->st.type].h_end_norm) {
      if (ctx->st.h_present[gctx.UDICT[ctx->st.type].h_end_cterm - 1]) has_oxt = 1;
    }
    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), has_oxt, include_hydrogens, af_mode);
    if (new_chain || ctx->cb.atom_bytes >= 580000 || (resSeq < ctx->cb.last_resSeq) || (resSeq - ctx->cb.last_resSeq > 100) || has_oxt) {
      if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60) || has_oxt) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = has_oxt;
        flush_chain(&ctx->cb, &ctx->coords);
      }
      if (new_chain) { ctx->current_c1 = c1; ctx->current_c2 = c2; }
      ctx->cb.c1 = c1; ctx->cb.c2 = c2;
    }
    ctx->st.resSeq = resSeq; ctx->st.iCode = iCode; ctx->st.type = type;
    memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
    for (int k = 0; k < 256; k++) { alt_init(&ctx->st.h_atoms[k]); alt_init(&ctx->st.hy_atoms[k]); }
  } else {
    int is_dup = 0;
    if (h_idx != -1 && ctx->st.h_present[h_idx]) { for (int c = 0; c < ctx->st.h_atoms[h_idx].count; c++) { if (ctx->st.h_atoms[h_idx].alt[c] == altLoc) { is_dup = 1; break; } } }
    else if (hy_idx != -1 && ctx->st.hy_present[hy_idx]) { for (int c = 0; c < ctx->st.hy_atoms[hy_idx].count; c++) { if (ctx->st.hy_atoms[hy_idx].alt[c] == altLoc) { is_dup = 1; break; } } }
    if (is_dup) {
      flush_residue(&ctx->cb, &ctx->st, gctx.UDICT, (ctx->cb.num_residues == 0), 0, include_hydrogens, af_mode);
      memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
      for (int k = 0; k < 256; k++) { alt_init(&ctx->st.h_atoms[k]); alt_init(&ctx->st.hy_atoms[k]); }
    }
  }

  if (h_idx != -1) { alt_push(&ctx->st.h_atoms[h_idx], temp_int, x_int, y_int, z_int, occ_int, charge, altLoc); ctx->st.h_present[h_idx] = 1; }
  else if (hy_idx != -1) { alt_push(&ctx->st.hy_atoms[hy_idx], temp_int, x_int, y_int, z_int, occ_int, charge, altLoc); ctx->st.hy_present[hy_idx] = 1; }
}

// --- Metadata Handlers (Bounds Checked) ---
static void parse_header(const char* line, int line_len, EncodeMetadata* m, int af_mode) {
  (void)af_mode;
  ext_str(m->classification, line, line_len, 11, 40);
  char date_str[10]; ext_str(date_str, line, line_len, 51, 9);
  if (strlen(date_str) == 9) {
    char d[3]={0}, mn[4]={0}, y[3]={0}; memcpy(d, date_str, 2); memcpy(mn, date_str+3, 3); memcpy(y, date_str+7, 2);
    m->day = atoi(d); m->month = parse_month(mn); int parsed_y = atoi(y);
    m->year = parsed_y + ((parsed_y < 50) ? 2000 : 1900);
  }
  char temp_id[5]; ext_str(temp_id, line, line_len, 63, 4);
  if (strlen(temp_id) > 0 && strcmp(temp_id, "XXXX") != 0) { strcpy(m->pdbID, temp_id); m->id_from_header = 1; }
}

static void parse_title(const char* line, int line_len, EncodeMetadata* m, int af_mode) {
  char chunk[71]; ext_str_rtrim(chunk, line, line_len, 11, 70); strcat(m->title, chunk);
  if (!m->id_from_header && !af_mode) {
    const char *p1 = strchr(line, '('), *p2 = strchr(p1 ? p1 : line, ')');
    if (p1 && p2 && (p2 > p1 + 4)) { size_t len = p2 - p1 - 1; memcpy(m->pdbID, p1+1, len > 12 ? 12 : len); m->pdbID[len > 12 ? 12 : len] = '\0'; }
  }
}

static void parse_dbref(const char* line, int line_len, EncodeMetadata* m, int af_mode) {
  MetaDBRef db; memset(&db, ' ', sizeof(db));
  if (af_mode) {
    db.chain = ext_char(line, line_len, 13); db.seq_begin = ext_int(line, line_len, 15, 4);
    db.iCode_begin = ext_char(line, line_len, 19); if (db.iCode_begin == '\0') db.iCode_begin = ' ';
    db.seq_end = ext_int(line, line_len, 21, 4); db.iCode_end = ext_char(line, line_len, 25); if (db.iCode_end == '\0') db.iCode_end = ' ';
    char tmp_name[9] = {0}; ext_str(tmp_name, line, line_len, 27, 6); memcpy(db.db_name, tmp_name, strlen(tmp_name));
    char tmp_acc[13] = {0}; ext_str(tmp_acc, line, line_len, 34, 10); memcpy(db.db_accession, tmp_acc, strlen(tmp_acc));
    char tmp_id[13] = {0}; ext_str(tmp_id, line, line_len, 45, 12); memcpy(db.db_id_code, tmp_id, strlen(tmp_id));
    db.db_seq_begin = ext_int(line, line_len, 58, 5); db.db_iCode_begin = ext_char(line, line_len, 63); if (db.db_iCode_begin == '\0') db.db_iCode_begin = ' ';
    db.db_seq_end = ext_int(line, line_len, 65, 5); db.db_iCode_end = ext_char(line, line_len, 70); if (db.db_iCode_end == '\0') db.db_iCode_end = ' ';
  } else {
    db.chain = ext_char(line, line_len, 13); db.seq_begin = ext_int(line, line_len, 15, 4);
    db.iCode_begin = ext_char(line, line_len, 19); if (db.iCode_begin == '\0') db.iCode_begin = ' ';
    db.seq_end = ext_int(line, line_len, 21, 4); db.iCode_end = ext_char(line, line_len, 25); if (db.iCode_end == '\0') db.iCode_end = ' ';
    char tmp_name[9] = {0}; ext_str(tmp_name, line, line_len, 27, 6); memcpy(db.db_name, tmp_name, strlen(tmp_name));
    char tmp_acc[13] = {0}; ext_str(tmp_acc, line, line_len, 34, 8); memcpy(db.db_accession, tmp_acc, strlen(tmp_acc));
    char tmp_id[13] = {0}; ext_str(tmp_id, line, line_len, 43, 12); memcpy(db.db_id_code, tmp_id, strlen(tmp_id));
    db.db_seq_begin = ext_int(line, line_len, 56, 5); db.db_iCode_begin = ext_char(line, line_len, 61); if (db.db_iCode_begin == '\0') db.db_iCode_begin = ' ';
    db.db_seq_end = ext_int(line, line_len, 63, 5); db.db_iCode_end = ext_char(line, line_len, 68); if (db.db_iCode_end == '\0') db.db_iCode_end = ' ';
  }
  stream_write(&m->dbrefs, &db, sizeof(MetaDBRef));
  char tmp_acc_check[13] = {0}; memcpy(tmp_acc_check, db.db_accession, 12); trim_string(tmp_acc_check);
  if (!m->id_from_header && !af_mode && strncmp(db.db_name, "UNP", 3) == 0 && strlen(tmp_acc_check) > 3) { strncpy(m->pdbID, tmp_acc_check, 12); m->pdbID[12] = '\0'; }
}

static void parse_seqres(const char* line, int line_len, EncodeMetadata* m) {
  char chain = ext_char(line, line_len, 12);
  m->has_seqres_chain[(unsigned char)chain] = 1;
  if (chain != m->current_seqres_chain) {
    if (m->current_seqres_count > 0) { 
      uint16_t* count_ptr = (uint16_t*)(m->seqres_buf.data + m->seqres_buf.len - m->current_seqres_count - 2); 
      *count_ptr = m->current_seqres_count; 
    }
    m->current_seqres_chain = chain; m->current_seqres_count = 0;
    MetaSeqresHeader sh = {chain, 0}; 
    stream_write(&m->seqres_buf, &sh, sizeof(MetaSeqresHeader)); // FIX
  }
  for (int i = 0; i < 13; i++) {
    int offset = 19 + i * 4; if (offset + 2 >= line_len) break;
    char res[4] = {0}; memcpy(res, line + offset, 3); char trimmed[4]; strcpy(trimmed, res); trim_string(trimmed);
    if (strlen(trimmed) == 0) continue;
    uint8_t type_id = gctx.res_lut[res_hash(res)];
    stream_write(&m->seqres_buf, &type_id, 1); // FIX
    m->current_seqres_count++;
  }
}

// --- The Refactored Main Encode Pass (Single-Pass MMAP Architecture) ---
void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb, int include_hydrogens, int af_mode, int afdb_mode) {
  init_gctx();

  for (int f = 0; f < num_files; f++) {
    reset_gctx_for_file();

    int fd = open(in_pdbs[f], O_RDONLY);
    if (fd < 0) { perror("Input file error"); continue; }
    
    struct stat st;
    if (fstat(fd, &st) < 0) { close(fd); continue; }
    
    const char* map = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (map == MAP_FAILED) { close(fd); continue; }

    if (af_mode) {
      char temp[128]; strcpy(temp, in_pdbs[f]);
      char* base = strrchr(temp, '/'); base = base ? base + 1 : temp;
      char* af = strstr(base, "AF-"); if (af) memmove(base, af + 3, strlen(af + 3) + 1);
      char* flag = strstr(base, "-F"); if (flag) *flag = '\0';
      char* dot = strrchr(base, '.'); if (dot) *dot = '\0';
      strncpy(gctx.m.pdbID, base, 12); gctx.m.pdbID[12] = '\0';
    }

    const char* ptr = map;
    const char* end = map + st.st_size;

    while (ptr < end) {
      const char* nl = memchr(ptr, '\n', end - ptr);
      int line_len = nl ? (nl - ptr) : (end - ptr);

      if (line_len >= 6) {
        uint64_t tag = (*(uint64_t*)ptr) & 0x0000FFFFFFFFFFFFULL; // LE exactly maps to 6 characters
        switch (tag) {
        case 0x000020204d4f5441ULL: // "ATOM  "
        case 0x00004d5441544548ULL: // "HETATM"
          process_atom(ptr, line_len, &gctx.cctx, &gctx.m, include_hydrogens, af_mode); break;
        case 0x0000202020524554ULL: // "TER   "
          process_ter(&gctx.cctx, &gctx.m, include_hydrogens, af_mode); break;
        case 0x00004c444d444e45ULL: // "ENDMDL"
          process_endmdl(&gctx.cctx, out_bpdb, &gctx.m, include_hydrogens, af_mode, afdb_mode); break;
        case 0x0000524544414548ULL: // "HEADER"
          parse_header(ptr, line_len, &gctx.m, af_mode); break;
        case 0x000020454c544954ULL: // "TITLE "
          parse_title(ptr, line_len, &gctx.m, af_mode); break;
        case 0x0000204645524244ULL: // "DBREF "
          parse_dbref(ptr, line_len, &gctx.m, af_mode); break;
        case 0x0000534552514553ULL: // "SEQRES"
          parse_seqres(ptr, line_len, &gctx.m); break;
        case 0x00004c444d4d554eULL: // "NUMMDL"
          gctx.m.parsed_num_models = ext_int(ptr, line_len, 11, 4); break;
        case 0x0000204c45444f4dULL: // "MODEL "
          gctx.m.counted_models++; break;
        case 0x0000444e504d4f43ULL: // "COMPND"
          { char chunk[71]; ext_str(chunk, ptr, line_len, 11, 70); if (strlen(gctx.m.raw_compnd) > 0) strcat(gctx.m.raw_compnd, " "); strcat(gctx.m.raw_compnd, chunk); } break;
        case 0x0000454352554f53ULL: // "SOURCE"
          { char chunk[70]; ext_str(chunk, ptr, line_len, 11, 69); if (strlen(gctx.m.raw_source) > 0) strcat(gctx.m.raw_source, " "); strcat(gctx.m.raw_source, chunk); } break;
        case 0x0000415444505845ULL: // "EXPDTA"
          { char chunk[70]; ext_str(chunk, ptr, line_len, 11, 69); if (strlen(gctx.m.expdta) > 0) strcat(gctx.m.expdta, " "); strncat(gctx.m.expdta, chunk, 127 - strlen(gctx.m.expdta)); } break;
        case 0x00004b52414d4552ULL: // "REMARK"
          { 
            int rem_num = ext_int(ptr, line_len, 8, 3);
            if (rem_num == 2) { char t[12]; ext_str(t, ptr, line_len, 12, 11); if (strcmp(t, "RESOLUTION.") == 0) { float res = ext_float(ptr, line_len, 24, 7); if (res > 0.0f) gctx.m.resolution = res; } }
            else if (rem_num == 3) { 
              char full[80]; ext_str(full, ptr, line_len, 12, 69); char* prog = strstr(full, "PROGRAM"); char* ver = strstr(full, "VERSION"); char* colon = strchr(full, ':');
              if (prog && colon && strlen(gctx.m.program) == 0) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; strcpy(temp, start); trim_string(temp); strncat(gctx.m.program, temp, 55); }
              else if (ver && colon && strlen(gctx.m.program) > 0 && !strstr(gctx.m.program, "v")) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; strcpy(temp, start); trim_string(temp); strncat(gctx.m.program, " v", 55 - strlen(gctx.m.program)); strncat(gctx.m.program, temp, 55 - strlen(gctx.m.program)); }
            }
          } break;
        case 0x0000444e4f425353ULL: // "SSBOND"
          { MetaSSBond b; memset(&b, 0, sizeof(b)); ext_str(b.res1.resName, ptr, line_len, 12, 3); b.res1.chain = ext_char(ptr, line_len, 16); b.res1.resSeq = ext_int(ptr, line_len, 18, 4); b.res1.iCode = ext_char(ptr, line_len, 22); ext_str(b.res2.resName, ptr, line_len, 26, 3); b.res2.chain = ext_char(ptr, line_len, 30); b.res2.resSeq = ext_int(ptr, line_len, 32, 4); b.res2.iCode = ext_char(ptr, line_len, 36); ext_str(b.sym1, ptr, line_len, 60, 6); ext_str(b.sym2, ptr, line_len, 67, 6); b.length = ext_float(ptr, line_len, 74, 5); stream_write(&gctx.m.ssbonds, &b, sizeof(MetaSSBond)); } break;
        case 0x0000504550534943ULL: // "CISPEP"
          { MetaCispep c; memset(&c, 0, sizeof(c)); ext_str(c.res1.resName, ptr, line_len, 12, 3); c.res1.chain = ext_char(ptr, line_len, 16); c.res1.resSeq = ext_int(ptr, line_len, 18, 4); c.res1.iCode = ext_char(ptr, line_len, 22); ext_str(c.res2.resName, ptr, line_len, 26, 3); c.res2.chain = ext_char(ptr, line_len, 30); c.res2.resSeq = ext_int(ptr, line_len, 32, 4); c.res2.iCode = ext_char(ptr, line_len, 36); c.modNum = ext_int(ptr, line_len, 44, 3); c.angle = ext_float(ptr, line_len, 54, 6); stream_write(&gctx.m.cispeps, &c, sizeof(MetaCispep)); } break;
        case 0x00002058494c4548ULL: // "HELIX "
          { MetaHelix h; memset(&h, 0, sizeof(h)); ext_str(h.helixID, ptr, line_len, 12, 3); ext_str(h.start.resName, ptr, line_len, 16, 3); h.start.chain = ext_char(ptr, line_len, 20); h.start.resSeq = ext_int(ptr, line_len, 22, 4); h.start.iCode = ext_char(ptr, line_len, 26); ext_str(h.end.resName, ptr, line_len, 28, 3); h.end.chain = ext_char(ptr, line_len, 32); h.end.resSeq = ext_int(ptr, line_len, 34, 4); h.end.iCode = ext_char(ptr, line_len, 38); h.helix_class = ext_int(ptr, line_len, 39, 2); ext_str(h.comment, ptr, line_len, 41, 30); h.length = ext_int(ptr, line_len, 72, 5); stream_write(&gctx.m.helices, &h, sizeof(MetaHelix)); } break;
        case 0x0000205445454853ULL: // "SHEET "
          { MetaSheet s; memset(&s, 0, sizeof(s)); s.strand_num = ext_int(ptr, line_len, 8, 3); ext_str(s.sheetID, ptr, line_len, 12, 3); s.num_strands = ext_int(ptr, line_len, 15, 2); ext_str(s.start.resName, ptr, line_len, 18, 3); s.start.chain = ext_char(ptr, line_len, 22); s.start.resSeq = ext_int(ptr, line_len, 23, 4); s.start.iCode = ext_char(ptr, line_len, 27); ext_str(s.end.resName, ptr, line_len, 29, 3); s.end.chain = ext_char(ptr, line_len, 33); s.end.resSeq = ext_int(ptr, line_len, 34, 4); s.end.iCode = ext_char(ptr, line_len, 38); s.sense = ext_int(ptr, line_len, 39, 2); stream_write(&gctx.m.sheets, &s, sizeof(MetaSheet)); } break;
        case 0x0000202045544953ULL: // "SITE  "
          { int seqNum = ext_int(ptr, line_len, 8, 3); if (seqNum == 1) { MetaSiteHeader sh; memset(&sh, 0, sizeof(sh)); ext_str(sh.siteID, ptr, line_len, 12, 3); sh.num_res = ext_int(ptr, line_len, 16, 2); stream_write(&gctx.m.sites, &sh, sizeof(MetaSiteHeader)); } for (int i = 0; i < 4; i++) { int off = 19 + i * 11; char check[4]; ext_str(check, ptr, line_len, off, 3); if (strlen(check) == 0) continue; PdbResId r; memset(&r, 0, sizeof(r)); ext_str(r.resName, ptr, line_len, off, 3); r.chain = ext_char(ptr, line_len, off + 4); r.resSeq = ext_int(ptr, line_len, off + 5, 4); r.iCode = ext_char(ptr, line_len, off + 9); stream_write(&gctx.m.sites, &r, sizeof(PdbResId)); } } break;
        case 0x00005443454e4f43ULL: // "CONECT"
          { int32_t vals[5] = {0}, targets = 0; for (int c = 0; c < 5; c++) { vals[c] = ext_int(ptr, line_len, 7 + c * 5, 5); if (c > 0 && vals[c] > 0) targets++; } if (vals[0] > 0 && targets > 0) { stream_write(&gctx.m.conects, &targets, 4); stream_write(&gctx.m.conects, &vals[0], 4); for (int c = 1; c <= 4; c++) { if (vals[c] > 0) stream_write(&gctx.m.conects, &vals[c], 4); } } } break;
        }
      }
      ptr = nl ? nl + 1 : end;
    }

    if (gctx.m.current_seqres_count > 0) { 
      uint16_t* count_ptr = (uint16_t*)(gctx.m.seqres_buf.data + gctx.m.seqres_buf.len - gctx.m.current_seqres_count - 2); 
      *count_ptr = gctx.m.current_seqres_count; 
    }
    
    uint32_t num_models = (gctx.m.parsed_num_models > 0) ? gctx.m.parsed_num_models : gctx.m.counted_models;
    if (num_models == 0) num_models = 1;
    gctx.cctx.total_models_minus_1 = num_models - 1;

    if (gctx.cctx.coords.len > 0) {
      if (f < num_files - 1) process_file_boundary(&gctx.cctx, &gctx.m, include_hydrogens, af_mode);
      else process_eof_flush(&gctx.cctx, out_bpdb, &gctx.m, include_hydrogens, af_mode, afdb_mode);
    }

    munmap((void*)map, st.st_size);
    close(fd);
  }
}
