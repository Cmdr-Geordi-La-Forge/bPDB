#include "bpdb.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void stream_init(StreamBuffer* buf) {
  buf->cap = 1024 * 1024 * 16; buf->data = (uint8_t*)malloc(buf->cap); buf->len = 0;
}
void stream_write(StreamBuffer* buf, const void* src, size_t size) {
  if (buf->len + size > buf->cap) {
    buf->cap = (buf->cap == 0) ? 1024 : (buf->cap * 2);
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

static inline uint32_t pack_ec_number(const char* ec_str) {
  int class_val = 0, sub = 0, subsub = 0, serial = 0;
  char temp[32]; strncpy(temp, ec_str, 31); temp[31] = '\0';
  char* t = temp; char* p = strchr(t, '.');
  if (p) { *p = '\0'; class_val = atoi(t); t = p + 1; }
  p = strchr(t, '.');
  if (p) { *p = '\0'; sub = atoi(t); t = p + 1; }
  p = strchr(t, '.');
  if (p) {
    *p = '\0';
    if (t[0] == '-') subsub = 0; else if (t[0] == 'n') subsub = -atoi(t+1); else subsub = atoi(t);
    t = p + 1;
  }
  if (t[0] == '-') serial = 0; else if (t[0] == 'n') serial = -atoi(t+1); else serial = atoi(t);

  uint32_t packed = 0;
  packed |= ((class_val & 0x07) << 29);
  packed |= ((sub & 0x1F) << 24);
  packed |= (((uint32_t)((int8_t)subsub) & 0xFF) << 16);
  packed |= ((uint32_t)((int16_t)serial) & 0xFFFF);
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

      uint8_t token_id = 0;
      const char** keys = is_source ? SOURCE_KEYS : COMPND_KEYS;
      int max_keys = is_source ? 36 : 10;
      for (int i = 1; i < max_keys; i++) {
        if (strcmp(key, keys[i]) == 0) { token_id = i; break; }
      }

      if (token_id == 1) {
        if (mol_buf.len > 0) {
          uint16_t mlen = mol_buf.len; stream_write(out_buf, &mlen, 2);
          stream_write(out_buf, mol_buf.data, mol_buf.len); mol_buf.len = 0;
        }
        uint8_t vlen = strlen(val); stream_write(&mol_buf, &token_id, 1);
        stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, val, vlen);
      } else if (token_id > 1) {
        if (!is_source && token_id == COMPND_EC) {
          uint32_t packed_ec = pack_ec_number(val); uint8_t vlen = 4;
          stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1);
          stream_write(&mol_buf, &packed_ec, 4);
        } else if (is_source && af_mode && (token_id == SOURCE_ORGANISM_SCIENTIFIC || token_id == SOURCE_ORGANISM_COMMON || token_id == SOURCE_EXPRESSION_SYSTEM || token_id == SOURCE_EXPRESSION_SYSTEM_COMMON)) {
          continue;
        } else if (is_source && (token_id == SOURCE_ORGANISM_TAXID || token_id == SOURCE_EXPRESSION_SYSTEM_TAXID)) {
          uint32_t taxid = atoi(val); uint8_t vlen = 4;
          stream_write(&mol_buf, &token_id, 1); stream_write(&mol_buf, &vlen, 1);
          stream_write(&mol_buf, &taxid, 4);
        } else {
          uint8_t vlen = strlen(val); stream_write(&mol_buf, &token_id, 1);
          stream_write(&mol_buf, &vlen, 1); stream_write(&mol_buf, val, vlen);
        }
      }
    }
    tok = strtok(NULL, ";");
  }
  if (mol_buf.len > 0) {
    uint16_t mlen = mol_buf.len; stream_write(out_buf, &mlen, 2);
    stream_write(out_buf, mol_buf.data, mol_buf.len);
  }
  stream_free(&mol_buf);
}

typedef struct {
  int32_t *b, *x, *y, *z, *occ;
  int8_t *charge; char *alt;
  int count, cap;
} AltArray;

void alt_init(AltArray* a) {
  a->cap = 2; a->count = 0;
  a->b = (int32_t*)malloc(2 * sizeof(int32_t)); a->x = (int32_t*)malloc(2 * sizeof(int32_t));
  a->y = (int32_t*)malloc(2 * sizeof(int32_t)); a->z = (int32_t*)malloc(2 * sizeof(int32_t));
  a->occ = (int32_t*)malloc(2 * sizeof(int32_t));
  a->charge = (int8_t*)malloc(2 * sizeof(int8_t)); a->alt = (char*)malloc(2 * sizeof(char));
}

void alt_push(AltArray* a, int32_t b, int32_t x, int32_t y, int32_t z, int32_t occ, int8_t charge, char alt) {
  if (a->count == a->cap) {
    a->cap *= 2;
    a->b = (int32_t*)realloc(a->b, a->cap * sizeof(int32_t)); a->x = (int32_t*)realloc(a->x, a->cap * sizeof(int32_t));
    a->y = (int32_t*)realloc(a->y, a->cap * sizeof(int32_t)); a->z = (int32_t*)realloc(a->z, a->cap * sizeof(int32_t));
    a->occ = (int32_t*)realloc(a->occ, a->cap * sizeof(int32_t));
    a->charge = (int8_t*)realloc(a->charge, a->cap * sizeof(int8_t)); a->alt = (char*)realloc(a->alt, a->cap * sizeof(char));
  }
  int i = a->count++;
  a->b[i] = b; a->x[i] = x; a->y[i] = y; a->z[i] = z; a->occ[i] = occ; a->charge[i] = charge; a->alt[i] = alt;
}

void alt_free(AltArray* a) {
  free(a->b); free(a->x); free(a->y); free(a->z); free(a->occ); free(a->charge); free(a->alt);
  a->cap = 0; a->count = 0;
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
  cb->current_b = 0;
  cb->last_type = 0;
}

void cb_free(ChainBuilder* cb) { free(cb->seq_track); free(cb->atom_track); }

void flush_chain(ChainBuilder* cb, StreamBuffer* buf) {
  if (cb->num_atoms == 0) {
    cb->has_ter = 0; cb->has_ter_raw = 0; cb->is_n_term = 0; cb->is_c_term = 0;
    return;
  }
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
  
  cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0;
  cb->has_ter = 0; cb->has_ter_raw = 0; cb->is_n_term = 0; cb->is_c_term = 0; cb->current_b = 0;
  cb->last_type = 0;
}

static inline void split_msb_lsb(int32_t val, int16_t* msb, int16_t* lsb) {
  int32_t base = (val + 4096) >> 13;
  *msb = (int16_t)base;
  *lsb = (int16_t)(val - (base << 13));
}

static inline void encode_alt_array(ChainBuilder* cb, AltArray* arr, int force_jump, int af_mode) {
  for (int a = 0; a < arr->count; a++) {
    int32_t x = arr->x[a], y = arr->y[a], z = arr->z[a], b_int = arr->b[a];
    if (b_int > 4194303) b_int = 4194303;
    if (b_int < 0) b_int = 0;

    if (cb->num_atoms == 0) { cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int; }

    int32_t dx = x - cb->current_x, dy = y - cb->current_y;
    int32_t dz = z - cb->current_z, db = b_int - cb->current_b;

    int alt_flag = (a < arr->count - 1) ? 1 : 0;
    int32_t occ_int = arr->occ[a]; if (occ_int > 100) occ_int = 100; else if (occ_int < 0) occ_int = 0;
    int8_t raw_charge = arr->charge[a];
    uint8_t charge_u = (raw_charge + 3) & 0x07;

    // AFDB Mode utilizes 15-bit deltas (-16384 to 16383). Standard uses 13-bit (-4096 to 4095).
    int limit = af_mode ? 16383 : 4095;
    int needs_jump = (cb->num_atoms == 0) || (force_jump && a == 0) ||
      (dx > limit || dx < -(limit+1) || dy > limit || dy < -(limit+1) ||
       dz > limit || dz < -(limit+1) || db > limit || db < -(limit+1));

    uint8_t* ptr = &cb->atom_track[cb->atom_bytes];
    if (needs_jump) {
      int16_t b_msb, b_lsb, x_msb, x_lsb, y_msb, y_lsb, z_msb, z_lsb;
      // Absolute coordinates are packed in 15-bit halves
      split_msb_lsb(b_int, &b_msb, &b_lsb);
      split_msb_lsb(x, &x_msb, &x_lsb);
      split_msb_lsb(y, &y_msb, &y_lsb);
      split_msb_lsb(z, &z_msb, &z_lsb);

      *(uint16_t*)(ptr)     = (b_msb << 1) | 1;
      *(uint16_t*)(ptr + 2) = (x_msb << 1) | 0;
      *(uint16_t*)(ptr + 4) = (y_msb << 1) | 0;
      *(uint16_t*)(ptr + 6) = (z_msb << 1) | 0;
      
      if (af_mode) {
          *(uint16_t*)(ptr + 8)  = (b_lsb << 1) | 0;
          *(uint16_t*)(ptr + 10) = (x_lsb << 1) | 0;
          *(uint16_t*)(ptr + 12) = (y_lsb << 1) | 0;
          *(uint16_t*)(ptr + 14) = (z_lsb << 1) | 0;
      } else {
          *(uint16_t*)(ptr + 8)  = (b_lsb << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
          *(uint16_t*)(ptr + 10) = (x_lsb << 3) | (((occ_int >> 4) & 3) << 1) | ((charge_u >> 2) & 1);
          *(uint16_t*)(ptr + 12) = (y_lsb << 3) | (((occ_int >> 2) & 3) << 1) | ((charge_u >> 1) & 1);
          *(uint16_t*)(ptr + 14) = (z_lsb << 3) | ((occ_int & 3) << 1) | (charge_u & 1);
      }
      cb->atom_bytes += 16;
    } else {
      if (af_mode) {
          *(uint16_t*)(ptr)     = (db << 1) | 0;
          *(uint16_t*)(ptr + 2) = (dx << 1) | 0;
          *(uint16_t*)(ptr + 4) = (dy << 1) | 0;
          *(uint16_t*)(ptr + 6) = (dz << 1) | 0;
      } else {
          *(uint16_t*)(ptr)     = (db << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
          *(uint16_t*)(ptr + 2) = (dx << 3) | (((occ_int >> 4) & 3) << 1) | ((charge_u >> 2) & 1);
          *(uint16_t*)(ptr + 4) = (dy << 3) | (((occ_int >> 2) & 3) << 1) | ((charge_u >> 1) & 1);
          *(uint16_t*)(ptr + 6) = (dz << 3) | ((occ_int & 3) << 1) | (charge_u & 1);
      }
      cb->atom_bytes += 8;
    }

    cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int;
    cb->num_atoms++;
  }
}

void flush_residue(ChainBuilder* cb, EncoderState* st, ResidueDict* UDICT, int is_first, int is_last, int include_hydrogens, int af_mode) {
  if (st->type <= 0) return;
  if (cb->num_residues == 0) cb->start_resSeq = st->resSeq;

  if (cb->num_residues > 0) {
    int diff = st->resSeq - cb->last_resSeq;
    if (diff == 0) {
      cb->seq_track[cb->seq_bytes++] = 0x7F; cb->num_residues++;
    } else if (diff > 1 && diff <= 100) {
      for (int g = 1; g < diff; g++) { cb->seq_track[cb->seq_bytes++] = 0x00; cb->num_residues++; }
    }
  }

  if (st->iCode != ' ' && st->iCode != 0) {
    cb->seq_track[cb->seq_bytes++] = 0x7C;
    cb->seq_track[cb->seq_bytes++] = st->iCode;
    cb->num_residues++;
  }

  cb->seq_track[cb->seq_bytes++] = st->type;
  cb->last_type = st->type;
  cb->num_residues++;
  int is_hetatm = (st->type >= 0x60);

  int h_start = is_first ? 0 : UDICT[st->type].h_start_norm;
  int h_end   = is_last  ? UDICT[st->type].h_end_cterm : UDICT[st->type].h_end_norm;

  for (int i = h_start; i < h_end; i++) {
    if (!st->h_present[i] || st->h_atoms[i].count == 0) {
      memset(&cb->atom_track[cb->atom_bytes], 0, 8);
      cb->atom_bytes += 8; cb->num_atoms++;
      continue;
    }
    encode_alt_array(cb, &st->h_atoms[i], (is_hetatm && i == h_start), af_mode);
  }

  if (include_hydrogens) {
    int hy_start = is_first ? 0 : UDICT[st->type].hy_start_norm;
    int hy_end   = is_last  ? UDICT[st->type].hy_end_cterm : UDICT[st->type].hy_end_norm;
    for (int i = hy_start; i < hy_end; i++) {
      if (!st->hy_present[i] || st->hy_atoms[i].count == 0) {
        memset(&cb->atom_track[cb->atom_bytes], 0, 8);
        cb->atom_bytes += 8; cb->num_atoms++;
        continue;
      }
      encode_alt_array(cb, &st->hy_atoms[i], 0, af_mode);
    }
  }

  if (!is_hetatm) {
    int anchor_idx = (st->type >= 0x20 && st->type < 0x60) ? 2 : 8;
    if (anchor_idx >= h_start && anchor_idx < h_end && st->h_present[anchor_idx]) {
      cb->current_x = st->h_atoms[anchor_idx].x[0];
      cb->current_y = st->h_atoms[anchor_idx].y[0];
      cb->current_z = st->h_atoms[anchor_idx].z[0];
      cb->current_b = st->h_atoms[anchor_idx].b[0];
    }
  }
  cb->last_resSeq = st->resSeq; cb->last_iCode = st->iCode;
}

static void save_model_file(const char* base_out, uint16_t total_m1, uint16_t cur_m1,
                            const char* gid, const char* gclass, int day, int month, int year,
                            StreamBuffer* meta, StreamBuffer* coords, StreamBuffer* conects) {
  char filename[1024];
  if (strcmp(base_out, "-") == 0) {
    strcpy(filename, "-");
  } else if (total_m1 > 0) {
    const char* ext = strrchr(base_out, '.');
    if (ext) {
      int baselen = ext - base_out;
      snprintf(filename, sizeof(filename), "%.*s_m%d%s", baselen, base_out, cur_m1 + 1, ext);
    } else {
      snprintf(filename, sizeof(filename), "%s_m%d", base_out, cur_m1 + 1);
    }
  } else {
    strcpy(filename, base_out);
  }

  FILE* fout = (strcmp(filename, "-") == 0) ? stdout : fopen(filename, "wb");
  if (!fout) { perror("Failed to open output"); exit(1); }

  uint8_t header[64] = {0};
  memcpy(header, "bPDB", 4); header[4] = 0; header[5] = (uint8_t)day; header[6] = (uint8_t)month;
  header[7] = (year >= 1973 && year <= 2228) ? (uint8_t)(year - 1973) : 0;
  memcpy(header + 8, gid, strlen(gid) > 12 ? 12 : strlen(gid));
  *(uint16_t*)(header + 20) = total_m1;
  *(uint16_t*)(header + 22) = cur_m1;
  memcpy(header + 24, gclass, strlen(gclass) > 40 ? 40 : strlen(gclass));
  fwrite(header, 1, 64, fout);

  fwrite(meta->data, 1, meta->len, fout);
  fwrite(coords->data, 1, coords->len, fout);
  if (conects->len > 0) fwrite(conects->data, 1, conects->len, fout);

  if (fout != stdout) {
    uint64_t zero = 0;
    size_t current_size = ftell(fout);
    size_t eof_pad = (8 - (current_size % 8)) % 8;
    if (eof_pad > 0) fwrite(&zero, 1, eof_pad, fout);
    fclose(fout);
  }
}

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

static void init_encode_meta(EncodeMetadata* m) {
  memset(m, 0, sizeof(EncodeMetadata));
  m->day = 1; m->month = 1; m->year = 1973;
  m->custom_het_idx = 0x81;
  m->resolution = -1.0f;
  stream_init(&m->seqres_buf); stream_init(&m->ssbonds); stream_init(&m->cispeps);
  stream_init(&m->helices); stream_init(&m->sheets); stream_init(&m->sites);
  stream_init(&m->conects); stream_init(&m->dbrefs);
}

static void free_encode_meta(EncodeMetadata* m) {
  stream_free(&m->seqres_buf); stream_free(&m->ssbonds); stream_free(&m->cispeps);
  stream_free(&m->helices); stream_free(&m->sheets); stream_free(&m->sites);
  stream_free(&m->conects); stream_free(&m->dbrefs);
}

// --- Strict Column Parsing Helpers ---

static void ext_str(char* dst, const char* line, int col, int width) {
  int len = strlen(line);
  int start = col - 1;
  if (start >= len) { dst[0] = '\0'; return; }
  int cp = width;
  if (start + cp > len) cp = len - start;
  while(cp > 0 && (line[start+cp-1] == '\n' || line[start+cp-1] == '\r')) cp--;
  if (cp > 0) { memcpy(dst, line+start, cp); dst[cp] = '\0'; } else dst[0] = '\0';
  trim_string(dst);
}

static int ext_int(const char* line, int col, int width) {
  char buf[32]; ext_str(buf, line, col, width);
  return strlen(buf) > 0 ? atoi(buf) : 0;
}

static float ext_float(const char* line, int col, int width) {
  char buf[32]; ext_str(buf, line, col, width);
  return strlen(buf) > 0 ? atof(buf) : 0.0f;
}

static char ext_char(const char* line, int col) {
  int start = col - 1;
  if (start >= (int)strlen(line)) return ' ';
  return line[start] == '\n' || line[start] == '\r' ? ' ' : line[start];
}

// --- PDB Format Dispatch Handlers ---

typedef void (*RecordParser)(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int include_hydrogens);

static void parse_header(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  ext_str(m->classification, line, 11, 40);
  char date_str[10]; ext_str(date_str, line, 51, 9);
  if (strlen(date_str) == 9) {
    char d[3]={0}, mn[4]={0}, y[3]={0};
    memcpy(d, date_str, 2); memcpy(mn, date_str+3, 3); memcpy(y, date_str+7, 2);
    m->day = atoi(d); m->month = parse_month(mn); int parsed_y = atoi(y);
    m->year = parsed_y + ((parsed_y < 50) ? 2000 : 1900);
  }
  char temp_id[5]; ext_str(temp_id, line, 63, 4);
  if (strlen(temp_id) > 0 && strcmp(temp_id, "XXXX") != 0) {
    strcpy(m->pdbID, temp_id); m->id_from_header = 1;
  }
}

static void parse_title(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)inc_h;
  if (!is_first_file) return;
  char chunk[71]; ext_str(chunk, line, 11, 70);
  strcat(m->title, chunk);
  if (!m->id_from_header && !af_mode) {
    // FIX: Add 'const' qualifiers to prevent compiler warnings
    const char *p1 = strchr(line, '('), *p2 = strchr(p1 ? p1 : line, ')');
    if (p1 && p2 && (p2 > p1 + 4)) {
      size_t len = p2 - p1 - 1;
      memcpy(m->pdbID, p1+1, len > 12 ? 12 : len); m->pdbID[len > 12 ? 12 : len] = '\0';
    }
  }
}

static void parse_compnd(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  char chunk[71]; ext_str(chunk, line, 11, 70);
  if (strlen(m->raw_compnd) > 0) strcat(m->raw_compnd, " ");
  strcat(m->raw_compnd, chunk);
}

static void parse_source(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  char chunk[70]; ext_str(chunk, line, 11, 69);
  if (strlen(m->raw_source) > 0) strcat(m->raw_source, " ");
  strcat(m->raw_source, chunk);
}

static void parse_expdta(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  char chunk[70]; ext_str(chunk, line, 11, 69);
  if (strlen(m->expdta) > 0) strcat(m->expdta, " ");
  strncat(m->expdta, chunk, 127 - strlen(m->expdta));
}

static void parse_remark(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  int rem_num = ext_int(line, 8, 3);
  if (rem_num == 2) {
    char t[12]; ext_str(t, line, 12, 11);
    if (strcmp(t, "RESOLUTION.") == 0) {
      float res = ext_float(line, 24, 7);
      if (res > 0.0f) m->resolution = res;
    }
  } else if (rem_num == 3) {
    char full[80]; ext_str(full, line, 12, 69);
    char* prog = strstr(full, "PROGRAM");
    char* ver = strstr(full, "VERSION");
    char* colon = strchr(full, ':');
    if (prog && colon && strlen(m->program) == 0) {
      char temp[64] = {0}; char* start = colon + 1;
      while (*start == ' ' || *start == '\t') start++;
      strcpy(temp, start); trim_string(temp); strncat(m->program, temp, 55);
    } else if (ver && colon && strlen(m->program) > 0 && !strstr(m->program, "v")) {
      char temp[64] = {0}; char* start = colon + 1;
      while (*start == ' ' || *start == '\t') start++;
      strcpy(temp, start); trim_string(temp);
      strncat(m->program, " v", 55 - strlen(m->program));
      strncat(m->program, temp, 55 - strlen(m->program));
    }
  }
}

static void parse_dbref(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)inc_h;
  if (!is_first_file) return;
  MetaDBRef db; memset(&db, ' ', sizeof(db));

  if (af_mode) {
    // AFDB shifted columns for 10-char accessions and larger sequence ranges
    db.chain = ext_char(line, 13);
    db.seq_begin = ext_int(line, 15, 4);
    db.iCode_begin = ' ';
    db.seq_end = ext_int(line, 20, 5);
    db.iCode_end = ' ';
        
    char tmp_name[9] = {0}; ext_str(tmp_name, line, 26, 3); memcpy(db.db_name, tmp_name, strlen(tmp_name));
    char tmp_acc[13] = {0}; ext_str(tmp_acc, line, 33, 10); memcpy(db.db_accession, tmp_acc, strlen(tmp_acc));
    char tmp_id[13] = {0}; ext_str(tmp_id, line, 42, 12); memcpy(db.db_id_code, tmp_id, strlen(tmp_id));
        
    db.db_seq_begin = ext_int(line, 55, 6);
    db.db_iCode_begin = ' ';
    db.db_seq_end = ext_int(line, 61, 6);
    db.db_iCode_end = ' ';
  } else {
    // Standard wwPDB column bounds
    db.chain = ext_char(line, 13);
    db.seq_begin = ext_int(line, 15, 4);
    db.iCode_begin = ext_char(line, 19); if (db.iCode_begin == '\0') db.iCode_begin = ' ';
    db.seq_end = ext_int(line, 21, 4);
    db.iCode_end = ext_char(line, 25); if (db.iCode_end == '\0') db.iCode_end = ' ';
        
    char tmp_name[9] = {0}; ext_str(tmp_name, line, 27, 6); memcpy(db.db_name, tmp_name, strlen(tmp_name));
    char tmp_acc[13] = {0}; ext_str(tmp_acc, line, 34, 8); memcpy(db.db_accession, tmp_acc, strlen(tmp_acc));
    char tmp_id[13] = {0}; ext_str(tmp_id, line, 43, 12); memcpy(db.db_id_code, tmp_id, strlen(tmp_id));
        
    db.db_seq_begin = ext_int(line, 56, 5);
    db.db_iCode_begin = ext_char(line, 61); if (db.db_iCode_begin == '\0') db.db_iCode_begin = ' ';
    db.db_seq_end = ext_int(line, 63, 5);
    db.db_iCode_end = ext_char(line, 68); if (db.db_iCode_end == '\0') db.db_iCode_end = ' ';
  }
    
  stream_write(&m->dbrefs, &db, sizeof(MetaDBRef));
    
  char tmp_acc_check[13] = {0};
  memcpy(tmp_acc_check, db.db_accession, 12); trim_string(tmp_acc_check);
  if (!m->id_from_header && !af_mode && strncmp(db.db_name, "UNP", 3) == 0 && strlen(tmp_acc_check) > 3) {
    strncpy(m->pdbID, tmp_acc_check, 12); m->pdbID[12] = '\0';
  }
}

static void parse_nummdl(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  m->parsed_num_models = ext_int(line, 11, 4);
}

static void parse_model(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)line; (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  m->counted_models++;
}

static void parse_seqres(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  char chain = ext_char(line, 12);
  m->has_seqres_chain[(unsigned char)chain] = 1;
  if (chain != m->current_seqres_chain) {
    if (m->current_seqres_count > 0) {
      uint16_t* count_ptr = (uint16_t*)(m->seqres_buf.data + m->seqres_buf.len - m->current_seqres_count - 2);
      *count_ptr = m->current_seqres_count;
    }
    m->current_seqres_chain = chain; m->current_seqres_count = 0;
    MetaSeqresHeader sh = {chain, 0}; stream_write(&m->seqres_buf, &sh, sizeof(MetaSeqresHeader));
  }
  for (int i = 0; i < 13; i++) {
    int offset = 19 + i * 4;
    if (offset + 2 >= (int)strlen(line)) break;
        
    // Exact copy to preserve formatting for dict lookup
    char res[4] = {0}; memcpy(res, line + offset, 3); 
    char trimmed[4]; strcpy(trimmed, res); trim_string(trimmed);
    if (strlen(trimmed) == 0) continue;

    uint8_t type_id = 0;
    for (int t = 1; t < 256; t++) {
      if ((UDICT[t].h_end_norm > 0 || UDICT[t].h_end_cterm > 0) && strncmp(UDICT[t].resName, res, 3) == 0) {
        type_id = t; break;
      }
    }
    stream_write(&m->seqres_buf, &type_id, 1); m->current_seqres_count++;
  }
}

static void parse_ssbond(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  MetaSSBond b; memset(&b, 0, sizeof(b));
  ext_str(b.res1.resName, line, 12, 3); b.res1.chain = ext_char(line, 16);
  b.res1.resSeq = ext_int(line, 18, 4); b.res1.iCode = ext_char(line, 22);
  ext_str(b.res2.resName, line, 26, 3); b.res2.chain = ext_char(line, 30);
  b.res2.resSeq = ext_int(line, 32, 4); b.res2.iCode = ext_char(line, 36);
  ext_str(b.sym1, line, 60, 6);
  ext_str(b.sym2, line, 67, 6);
  b.length = ext_float(line, 74, 5);
  stream_write(&m->ssbonds, &b, sizeof(MetaSSBond));
}

static void parse_cispep(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  MetaCispep c; memset(&c, 0, sizeof(c));
  ext_str(c.res1.resName, line, 12, 3); c.res1.chain = ext_char(line, 16);
  c.res1.resSeq = ext_int(line, 18, 4); c.res1.iCode = ext_char(line, 22);
  ext_str(c.res2.resName, line, 26, 3); c.res2.chain = ext_char(line, 30);
  c.res2.resSeq = ext_int(line, 32, 4); c.res2.iCode = ext_char(line, 36);
  c.modNum = ext_int(line, 44, 3);
  c.angle = ext_float(line, 54, 6);
  stream_write(&m->cispeps, &c, sizeof(MetaCispep));
}

static void parse_helix(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  MetaHelix h; memset(&h, 0, sizeof(h));
  ext_str(h.helixID, line, 12, 3);
  ext_str(h.start.resName, line, 16, 3); h.start.chain = ext_char(line, 20);
  h.start.resSeq = ext_int(line, 22, 4); h.start.iCode = ext_char(line, 26);
  ext_str(h.end.resName, line, 28, 3); h.end.chain = ext_char(line, 32);
  h.end.resSeq = ext_int(line, 34, 4); h.end.iCode = ext_char(line, 38);
  h.helix_class = ext_int(line, 39, 2);
  ext_str(h.comment, line, 41, 30);
  h.length = ext_int(line, 72, 5);
  stream_write(&m->helices, &h, sizeof(MetaHelix));
}

static void parse_sheet(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  MetaSheet s; memset(&s, 0, sizeof(s));
  s.strand_num = ext_int(line, 8, 3);
  ext_str(s.sheetID, line, 12, 3);
  s.num_strands = ext_int(line, 15, 2);
  ext_str(s.start.resName, line, 18, 3); s.start.chain = ext_char(line, 22);
  s.start.resSeq = ext_int(line, 23, 4); s.start.iCode = ext_char(line, 27);
  ext_str(s.end.resName, line, 29, 3); s.end.chain = ext_char(line, 33);
  s.end.resSeq = ext_int(line, 34, 4); s.end.iCode = ext_char(line, 38);
  s.sense = ext_int(line, 39, 2);
  stream_write(&m->sheets, &s, sizeof(MetaSheet));
}

static void parse_site(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)af_mode; (void)inc_h;
  if (!is_first_file) return;
  int seqNum = ext_int(line, 8, 3);
  if (seqNum == 1) {
    MetaSiteHeader sh; memset(&sh, 0, sizeof(sh));
    ext_str(sh.siteID, line, 12, 3);
    sh.num_res = ext_int(line, 16, 2);
    stream_write(&m->sites, &sh, sizeof(MetaSiteHeader));
  }
  for (int i = 0; i < 4; i++) {
    int off = 19 + i * 11;
    char check[4]; ext_str(check, line, off, 3);
    if (strlen(check) == 0) continue;
    PdbResId r; memset(&r, 0, sizeof(r));
    ext_str(r.resName, line, off, 3); r.chain = ext_char(line, off + 4);
    r.resSeq = ext_int(line, off + 5, 4); r.iCode = ext_char(line, off + 9);
    stream_write(&m->sites, &r, sizeof(PdbResId));
  }
}

static void parse_conect(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)UDICT; (void)is_first_file; (void)af_mode; (void)inc_h;
  int32_t vals[5] = {0}, targets = 0;
  for (int c = 0; c < 5; c++) {
    vals[c] = ext_int(line, 7 + c * 5, 5);
    if (c > 0 && vals[c] > 0) targets++;
  }
  if (vals[0] > 0 && targets > 0) {
    stream_write(&m->conects, &targets, 4);
    stream_write(&m->conects, &vals[0], 4);
    for (int c = 1; c <= 4; c++) {
      if (vals[c] > 0) stream_write(&m->conects, &vals[c], 4);
    }
  }
}

static void parse_atom(const char* line, EncodeMetadata* m, ResidueDict* UDICT, int is_first_file, int af_mode, int inc_h) {
  (void)is_first_file; (void)af_mode;
  if (strlen(line) < 54) return;
    
  // Exact memcpy to preserve spacing
  char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
  char raw_atom[5] = {0}; memcpy(raw_atom, line + 12, 4);
    
  char elem[3] = {0};
  if (strlen(line) >= 78) { elem[0] = line[76]; elem[1] = line[77]; }

  if (!inc_h && (raw_atom[1] == 'H' || raw_atom[0] == 'H' || (elem[0] == ' ' && elem[1] == 'D'))) return;

  int type = -1;
  for (int i = 1; i < 256; i++) {
    if ((UDICT[i].h_end_norm > 0 || UDICT[i].h_end_cterm > 0) && strncmp(UDICT[i].resName, res_name, 3) == 0) {
      type = i; break;
    }
  }
  if (type != -1) {
    int a_idx = -1;
    for (int a = 0; a < UDICT[type].h_end_norm; a++) {
      if (strncmp(UDICT[type].h_atoms[a], raw_atom, 4) == 0) { a_idx = a; break; }
    }
    if (a_idx == -1 && type >= 0x80 && UDICT[type].h_end_norm < 255) {
      a_idx = UDICT[type].h_end_norm++; UDICT[type].h_end_cterm = UDICT[type].h_end_norm;
      memcpy(UDICT[type].h_atoms[a_idx], raw_atom, 4); UDICT[type].h_atoms[a_idx][4] = '\0';
      UDICT[type].h_elements[a_idx][0] = elem[0]; UDICT[type].h_elements[a_idx][1] = elem[1]; UDICT[type].h_elements[a_idx][2] = '\0';
    }
  } else if (m->custom_het_idx < 192) {
    type = m->custom_het_idx++; 
    memcpy(UDICT[type].resName, res_name, 3); UDICT[type].resName[3] = '\0';
    UDICT[type].h_end_norm = 1; UDICT[type].h_end_cterm = 1;
    memcpy(UDICT[type].h_atoms[0], raw_atom, 4); UDICT[type].h_atoms[0][4] = '\0';
    UDICT[type].h_elements[0][0] = elem[0]; UDICT[type].h_elements[0][1] = elem[1]; UDICT[type].h_elements[0][2] = '\0';
  }
}

typedef struct {
  char prefix[7];
  RecordParser parse_func;
} RecordHandler;

static const RecordHandler HANDLERS[] = {
  {"HEADER", parse_header},
  {"TITLE ", parse_title},
  {"COMPND", parse_compnd},
  {"SOURCE", parse_source},
  {"EXPDTA", parse_expdta},
  {"REMARK", parse_remark},
  {"DBREF ", parse_dbref},
  {"NUMMDL", parse_nummdl},
  {"MODEL ", parse_model},
  {"SEQRES", parse_seqres},
  {"SSBOND", parse_ssbond},
  {"CISPEP", parse_cispep},
  {"HELIX ", parse_helix},
  {"SHEET ", parse_sheet},
  {"SITE  ", parse_site},
  {"CONECT", parse_conect},
  {"ATOM  ", parse_atom},
  {"HETATM", parse_atom}
};
static const int NUM_HANDLERS = sizeof(HANDLERS) / sizeof(RecordHandler);

static void scan_pdb_metadata(int num_files, const char** in_pdbs, EncodeMetadata* m, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
  for (int f = 0; f < num_files; f++) {
    FILE* fin = fopen(in_pdbs[f], "r");
    if (!fin) { perror("Input file error"); exit(1); }
    char line[128];
    while (fgets(line, sizeof(line), fin)) {
      for (int i = 0; i < NUM_HANDLERS; i++) {
        if (strncmp(line, HANDLERS[i].prefix, 6) == 0) {
          HANDLERS[i].parse_func(line, m, UDICT, (f == 0), af_mode, include_hydrogens);
          break;
        }
      }
    }
    fclose(fin);
  }
  if (m->current_seqres_count > 0) {
    uint16_t* count_ptr = (uint16_t*)(m->seqres_buf.data + m->seqres_buf.len - m->current_seqres_count - 2);
    *count_ptr = m->current_seqres_count;
  }
}

static uint64_t build_binary_headers(EncodeMetadata* m, StreamBuffer* meta_buffer, int include_hydrogens, int af_mode, int afdb_mode, ResidueDict* UDICT) {
  StreamBuffer compnd_tlv; stream_init(&compnd_tlv);
  if (strlen(m->raw_compnd) > 0) extract_tlv(m->raw_compnd, &compnd_tlv, 0, af_mode);

  StreamBuffer source_tlv; stream_init(&source_tlv);
  if (strlen(m->raw_source) > 0) extract_tlv(m->raw_source, &source_tlv, 1, af_mode);

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

  uint64_t meta_header = (meta_mask << 16);
  stream_write(meta_buffer, &meta_header, 8);

  if (meta_mask & (1ULL << META_RESOLUTION)) { uint32_t pad=0; stream_write(meta_buffer, &m->resolution, 4); stream_write(meta_buffer, &pad, 4); }
  if (meta_mask & (1ULL << META_TITLE)) { stream_write(meta_buffer, &title_len, 8); stream_write(meta_buffer, m->title, title_len); size_t pad = (8 - (title_len % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_EXPDTA)) { stream_write(meta_buffer, &expdta_len, 8); stream_write(meta_buffer, m->expdta, expdta_len); size_t pad = (8 - (expdta_len % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_PROGRAM)) { uint64_t p_len = 56; stream_write(meta_buffer, &p_len, 8); stream_write(meta_buffer, m->program, 56); }
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
  if (meta_mask & (1ULL << META_HELIX))  { uint64_t L = m->helices.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->helices.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SHEET))  { uint64_t L = m->sheets.len;     stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->sheets.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SSBOND)) { uint64_t L = m->ssbonds.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->ssbonds.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_CISPEP)) { uint64_t L = m->cispeps.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->cispeps.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SITE))   { uint64_t L = m->sites.len;      stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->sites.data, L);      size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_COMPND)) { uint64_t L = compnd_tlv.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, compnd_tlv.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SOURCE)) { uint64_t L = source_tlv.len;    stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, source_tlv.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_DBREF))  { uint64_t L = m->dbrefs.len;     stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->dbrefs.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SEQRES)) { uint64_t L = m->seqres_buf.len; stream_write(meta_buffer, &L, 8); stream_write(meta_buffer, m->seqres_buf.data, L); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(meta_buffer, &zero, pad); }

  stream_free(&compnd_tlv); stream_free(&source_tlv);
  return meta_mask;
}

// --- Coordinate Encoding Context & Helpers ---

typedef struct {
  StreamBuffer coords;
  EncoderState st;
  ChainBuilder cb;
  char current_c1;
  char current_c2;
  uint16_t current_model_idx;
  uint16_t total_models_minus_1;
} CoordinateContext;

static void process_ter(CoordinateContext* ctx, EncodeMetadata* m, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
    ctx->cb.has_ter_raw = 1;
    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
    if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
        flush_chain(&ctx->cb, &ctx->coords);
    } else { 
        ctx->cb.has_ter_raw = 0; 
    }
    
    ctx->st.type = -1; 
    memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); 
    memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
    for (int k = 0; k < 256; k++) { ctx->st.h_atoms[k].count = 0; ctx->st.hy_atoms[k].count = 0; }
}

static void process_endmdl(CoordinateContext* ctx, const char* out_bpdb, EncodeMetadata* m, StreamBuffer* meta_buffer, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
    if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
        flush_chain(&ctx->cb, &ctx->coords);
    }
    
    uint64_t end_stream = (0x3FFFULL << 18); 
    stream_write(&ctx->coords, &end_stream, 8);
    save_model_file(out_bpdb, ctx->total_models_minus_1, ctx->current_model_idx, m->pdbID, m->classification, m->day, m->month, m->year, meta_buffer, &ctx->coords, &m->conects);
    
    ctx->coords.len = 0; 
    ctx->current_model_idx++;

    ctx->st.type = -1; 
    memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); 
    memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
    for (int k = 0; k < 256; k++) { ctx->st.h_atoms[k].count = 0; ctx->st.hy_atoms[k].count = 0; }
}

static void process_atom(const char* line, CoordinateContext* ctx, EncodeMetadata* m, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
  if (strlen(line) < 54) return;
    
  char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
    
  char altLoc = ext_char(line, 17), iCode = ext_char(line, 27), c1 = ext_char(line, 21), c2 = ext_char(line, 22);
  int resSeq = ext_int(line, 23, 4);

  char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
  char elem[3] = {0};
  if (strlen(line) >= 78) { elem[0] = line[76]; elem[1] = line[77]; }
    
  if (!include_hydrogens && (clean_atom[1] == 'H' || clean_atom[0] == 'H' || (elem[0] == ' ' && elem[1] == 'D'))) return;

  int type = -1;
  for (int i = 1; i < 256; i++) {
    if ((UDICT[i].h_end_norm > 0 || UDICT[i].h_end_cterm > 0) && strncmp(UDICT[i].resName, res_name, 3) == 0) { 
      type = i; break; 
    }
  }
  if (type == -1) return;

  int new_chain = (c1 != ctx->current_c1 || c2 != ctx->current_c2);

  int32_t x_int = (int32_t)round(ext_float(line, 31, 8) * 1000.0f);
  int32_t y_int = (int32_t)round(ext_float(line, 39, 8) * 1000.0f);
  int32_t z_int = (int32_t)round(ext_float(line, 47, 8) * 1000.0f);
  float temp = ext_float(line, 61, 6);

  float occ = 100.0f;
  int charge = 0;
  if (!af_mode) {
    occ = ext_float(line, 55, 6) * 100.0f;
    char chg = ext_char(line, 79);
    if (chg >= '0' && chg <= '9') { charge = chg - '0'; if (ext_char(line, 80) == '-') charge = -charge; }
    if (charge < -3) charge = -3; else if (charge > 4) charge = 4;
  }

  int h_idx = -1, hy_idx = -1;
  for (int i = 0; i < UDICT[type].h_end_cterm; i++) {
    if (strncmp(clean_atom, UDICT[type].h_atoms[i], 4) == 0) { h_idx = i; break; }
  }
  if (include_hydrogens && h_idx == -1) {
    for (int i = 0; i < UDICT[type].hy_end_cterm; i++) {
      if (strncmp(clean_atom, UDICT[type].hy_atoms[i], 4) == 0) { hy_idx = i; break; }
    }
  }
  if (h_idx == -1 && hy_idx == -1) return;

  if (resSeq != ctx->st.resSeq || iCode != ctx->st.iCode || new_chain || ctx->cb.atom_bytes >= 580000 || type != ctx->st.type) {
    int has_oxt = 0;
    if (ctx->st.type > 0 && UDICT[ctx->st.type].h_end_cterm > UDICT[ctx->st.type].h_end_norm) {
      if (ctx->st.h_present[UDICT[ctx->st.type].h_end_cterm - 1]) has_oxt = 1;
    }

    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), has_oxt, include_hydrogens, af_mode);

    if (new_chain || ctx->cb.atom_bytes >= 580000 || (resSeq < ctx->cb.last_resSeq) || (resSeq - ctx->cb.last_resSeq > 100) || has_oxt) {
      if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        // Only enforce auto-TER if the chain is a standard polymer
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60) || has_oxt) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = has_oxt;
        flush_chain(&ctx->cb, &ctx->coords);
      }
      if (new_chain) { ctx->current_c1 = c1; ctx->current_c2 = c2; }
      ctx->cb.c1 = c1; ctx->cb.c2 = c2;
    }

    ctx->st.resSeq = resSeq; ctx->st.iCode = iCode; ctx->st.type = type;
    memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
    for (int k = 0; k < 256; k++) { ctx->st.h_atoms[k].count = 0; ctx->st.hy_atoms[k].count = 0; }
  } else {
    int is_dup_atom = 0;
    if (h_idx != -1 && ctx->st.h_present[h_idx]) {
      for (int c = 0; c < ctx->st.h_atoms[h_idx].count; c++) { if (ctx->st.h_atoms[h_idx].alt[c] == altLoc) { is_dup_atom = 1; break; } }
    } else if (hy_idx != -1 && ctx->st.hy_present[hy_idx]) {
      for (int c = 0; c < ctx->st.hy_atoms[hy_idx].count; c++) { if (ctx->st.hy_atoms[hy_idx].alt[c] == altLoc) { is_dup_atom = 1; break; } }
    }
    if (is_dup_atom) {
      flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), 0, include_hydrogens, af_mode);
      memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
      for (int k = 0; k < 256; k++) { ctx->st.h_atoms[k].count = 0; ctx->st.hy_atoms[k].count = 0; }
    }
  }

  if (h_idx != -1) {
    alt_push(&ctx->st.h_atoms[h_idx], (int32_t)round(temp * 100.0f), x_int, y_int, z_int, (int32_t)round(occ), charge, altLoc);
    ctx->st.h_present[h_idx] = 1;
  } else if (hy_idx != -1) {
    alt_push(&ctx->st.hy_atoms[hy_idx], (int32_t)round(temp * 100.0f), x_int, y_int, z_int, (int32_t)round(occ), charge, altLoc);
    ctx->st.hy_present[hy_idx] = 1;
  }
}

static void process_file_boundary(CoordinateContext* ctx, EncodeMetadata* m, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
    if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
        flush_chain(&ctx->cb, &ctx->coords);
    }
    uint64_t frag_stream = (0x3FFEULL << 18); stream_write(&ctx->coords, &frag_stream, 8);
    
    ctx->st.type = -1; 
    memset(ctx->st.h_present, 0, sizeof(ctx->st.h_present)); 
    memset(ctx->st.hy_present, 0, sizeof(ctx->st.hy_present));
    for (int k = 0; k < 256; k++) { ctx->st.h_atoms[k].count = 0; ctx->st.hy_atoms[k].count = 0; }
}

static void process_eof_flush(CoordinateContext* ctx, const char* out_bpdb, EncodeMetadata* m, StreamBuffer* meta_buffer, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
    if (ctx->st.type > 0) flush_residue(&ctx->cb, &ctx->st, UDICT, (ctx->cb.num_residues == 0), 1, include_hydrogens, af_mode);
    if (ctx->cb.num_atoms > 0) {
        unsigned char ch = ctx->cb.c1; 
        ctx->cb.has_ter = (ctx->cb.has_ter_raw || (m->has_seqres_chain[ch] && ctx->cb.last_type < 0x60)) ? 1 : 0;
        ctx->cb.is_n_term = 1; ctx->cb.is_c_term = 1;
        flush_chain(&ctx->cb, &ctx->coords);
    }
    uint64_t end_stream = (0x3FFFULL << 18); stream_write(&ctx->coords, &end_stream, 8);
    save_model_file(out_bpdb, ctx->total_models_minus_1, ctx->current_model_idx, m->pdbID, m->classification, m->day, m->month, m->year, meta_buffer, &ctx->coords, &m->conects);
}

// --- The Refactored Main Encoder Stream ---

static void encode_coordinate_stream(int num_files, const char** in_pdbs, const char* out_bpdb, EncodeMetadata* m, StreamBuffer* meta_buffer, ResidueDict* UDICT, int include_hydrogens, int af_mode) {
    CoordinateContext ctx; 
    memset(&ctx, 0, sizeof(CoordinateContext));
    
    stream_init(&ctx.coords);
    ctx.st.type = -1;
    for(int i = 0; i < 256; i++) { 
        alt_init(&ctx.st.h_atoms[i]); 
        alt_init(&ctx.st.hy_atoms[i]); 
    }
    cb_init(&ctx.cb); 

    uint32_t num_models = (m->parsed_num_models > 0) ? m->parsed_num_models : m->counted_models;
    if (num_models == 0) num_models = 1;
    ctx.total_models_minus_1 = num_models - 1;

    for (int f = 0; f < num_files; f++) {
        FILE* fin = fopen(in_pdbs[f], "r");
        char line[128];
        
        while (fgets(line, sizeof(line), fin)) {
            if (strncmp(line, "TER   ", 6) == 0) {
                process_ter(&ctx, m, UDICT, include_hydrogens, af_mode);
            } else if (strncmp(line, "ENDMDL", 6) == 0) {
                process_endmdl(&ctx, out_bpdb, m, meta_buffer, UDICT, include_hydrogens, af_mode);
            } else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
                process_atom(line, &ctx, m, UDICT, include_hydrogens, af_mode);
            }
        }
        fclose(fin);

        if (f < num_files - 1 && ctx.coords.len > 0) {
            process_file_boundary(&ctx, m, UDICT, include_hydrogens, af_mode);
        }
    }

    if (ctx.coords.len > 0) {
        process_eof_flush(&ctx, out_bpdb, m, meta_buffer, UDICT, include_hydrogens, af_mode);
    }

    cb_free(&ctx.cb); 
    stream_free(&ctx.coords);
    for(int i = 0; i < 256; i++) { 
        alt_free(&ctx.st.h_atoms[i]); 
        alt_free(&ctx.st.hy_atoms[i]); 
    }
}

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb, int include_hydrogens, int af_mode, int afdb_mode) {
  EncodeMetadata m; init_encode_meta(&m);

  ResidueDict* UDICT = (ResidueDict*)calloc(256, sizeof(ResidueDict));
  for (size_t i = 0; i < NUM_DICT_RES; i++) UDICT[DICT[i].type_id] = DICT[i];

  scan_pdb_metadata(num_files, in_pdbs, &m, UDICT, include_hydrogens, af_mode);

  if (af_mode && !m.id_from_header) {
    char temp[128]; strcpy(temp, in_pdbs[0]);
    char* base = strrchr(temp, '/'); base = base ? base + 1 : temp;
    char* af = strstr(base, "AF-"); if (af) memmove(base, af + 3, strlen(af + 3) + 1);
    char* f = strstr(base, "-F"); if (f) *f = '\0';
    char* dot = strrchr(base, '.'); if (dot) *dot = '\0';
    strncpy(m.pdbID, base, 12); m.pdbID[12] = '\0';
  }

  StreamBuffer meta_buffer; stream_init(&meta_buffer);
  build_binary_headers(&m, &meta_buffer, include_hydrogens, af_mode, afdb_mode, UDICT);

  encode_coordinate_stream(num_files, in_pdbs, out_bpdb, &m, &meta_buffer, UDICT, include_hydrogens, af_mode);

  stream_free(&meta_buffer);
  free_encode_meta(&m);
  free(UDICT);
}
