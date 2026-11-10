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

      if (token_id == 1) { // MOL_ID boundary
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
} ChainBuilder;

void cb_init(ChainBuilder* cb) {
  cb->seq_track = (uint8_t*)malloc(240000); cb->atom_track = (uint8_t*)malloc(120000 * 8 + 100);
  cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0;
  cb->start_resSeq = 0; cb->has_ter = 0; cb->has_ter_raw = 0;
  cb->is_n_term = 0; cb->is_c_term = 0;
  cb->current_b = 0;
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
}

static inline void split_msb_lsb(int32_t val, int16_t* msb, int16_t* lsb) {
  int32_t base = (val + 4096) >> 13;
  *msb = (int16_t)base;
  *lsb = (int16_t)(val - (base << 13));
}

static inline void encode_alt_array(ChainBuilder* cb, AltArray* arr, int force_jump) {
  for (int a = 0; a < arr->count; a++) {
    int32_t x = arr->x[a], y = arr->y[a], z = arr->z[a], b_int = arr->b[a];
    if (b_int > 4194303) b_int = 4194303;
    if (b_int < 0) b_int = 0;

    if (cb->num_atoms == 0) {
      cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int;
    }

    int32_t dx = x - cb->current_x, dy = y - cb->current_y;
    int32_t dz = z - cb->current_z, db = b_int - cb->current_b;

    int alt_flag = (a < arr->count - 1) ? 1 : 0;
    int32_t occ_int = arr->occ[a]; if (occ_int > 100) occ_int = 100; else if (occ_int < 0) occ_int = 0;
    int8_t raw_charge = arr->charge[a];
    uint8_t charge_u = (raw_charge + 3) & 0x07;

    int needs_jump = (cb->num_atoms == 0) || (force_jump && a == 0) ||
      (dx > 4095 || dx < -4096 || dy > 4095 || dy < -4096 ||
       dz > 4095 || dz < -4096 || db > 4095 || db < -4096);

    uint8_t* ptr = &cb->atom_track[cb->atom_bytes];
    if (needs_jump) {
      int16_t b_msb, b_lsb, x_msb, x_lsb, y_msb, y_lsb, z_msb, z_lsb;
      split_msb_lsb(b_int, &b_msb, &b_lsb);
      split_msb_lsb(x, &x_msb, &x_lsb);
      split_msb_lsb(y, &y_msb, &y_lsb);
      split_msb_lsb(z, &z_msb, &z_lsb);

      *(uint16_t*)(ptr)     = (b_msb << 1) | 1;
      *(uint16_t*)(ptr + 2) = (x_msb << 1) | 0;
      *(uint16_t*)(ptr + 4) = (y_msb << 1) | 0;
      *(uint16_t*)(ptr + 6) = (z_msb << 1) | 0;

      *(uint16_t*)(ptr + 8)  = (b_lsb << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
      *(uint16_t*)(ptr + 10) = (x_lsb << 3) | ((occ_int >> 3) & 7);
      *(uint16_t*)(ptr + 12) = (y_lsb << 3) | (occ_int & 7);
      *(uint16_t*)(ptr + 14) = (z_lsb << 3) | charge_u;

      cb->atom_bytes += 16;
    } else {
      uint16_t b_ctrl = (db << 3) | (((occ_int >> 6) & 1) << 2) | (alt_flag << 1) | 0;
      uint16_t x_ctrl = (dx << 3) | ((occ_int >> 3) & 7);
      uint16_t y_ctrl = (dy << 3) | (occ_int & 7);
      uint16_t z_ctrl = (dz << 3) | charge_u;

      *(uint16_t*)(ptr)     = b_ctrl;
      *(uint16_t*)(ptr + 2) = x_ctrl;
      *(uint16_t*)(ptr + 4) = y_ctrl;
      *(uint16_t*)(ptr + 6) = z_ctrl;
      cb->atom_bytes += 8;
    }

    cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->current_b = b_int;
    cb->num_atoms++;
  }
}

void flush_residue(ChainBuilder* cb, EncoderState* st, ResidueDict* UDICT, int is_first, int is_last, int include_hydrogens) {
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
    encode_alt_array(cb, &st->h_atoms[i], (is_hetatm && i == h_start));
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
      encode_alt_array(cb, &st->hy_atoms[i], 0);
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

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb, int include_hydrogens, int af_mode, int afdb_mode) {
  ResidueDict* UDICT = (ResidueDict*)calloc(256, sizeof(ResidueDict));
  int custom_het_idx = 0x81;
  for (size_t i = 0; i < NUM_DICT_RES; i++) {
    UDICT[DICT[i].type_id] = DICT[i];
  }

  StreamBuffer seqres_buf; stream_init(&seqres_buf);
  char current_seqres_chain = '\0';
  uint16_t current_seqres_count = 0;

  char global_title[4096] = {0}; char global_id[13] = {0}; char global_class[128] = {0};
  char global_expdta[128] = {0}; char global_program[56] = {0};
  char raw_compnd[16384] = {0}; char raw_source[16384] = {0};
  float global_resolution = -1.0f; int day=1, month=1, year=1973, id_from_header = 0;

  StreamBuffer ssbonds; stream_init(&ssbonds); StreamBuffer cispeps; stream_init(&cispeps);
  StreamBuffer helices; stream_init(&helices); StreamBuffer sheets;  stream_init(&sheets);
  StreamBuffer sites;   stream_init(&sites); StreamBuffer conects; stream_init(&conects);
  StreamBuffer dbrefs;  stream_init(&dbrefs);

  int parsed_num_models = 0;
  int counted_models = 0;
  int has_seqres_chain[256] = {0};

  for (int f = 0; f < num_files; f++) {
    FILE* fin = fopen(in_pdbs[f], "r");
    if (!fin) { perror("Input file error"); exit(1); }
    char line[128]; char tbuf[16];
    while (fgets(line, sizeof(line), fin)) {
      if (f == 0 && strncmp(line, "HEADER", 6) == 0) {
        if (strlen(line) >= 50) { memcpy(global_class, line + 10, 40); trim_string(global_class); }
        if (strlen(line) >= 59) {
          char d[3]={0}, m[4]={0}, y[3]={0};
          memcpy(d, line+50, 2); memcpy(m, line+53, 3); memcpy(y, line+57, 2);
          day = atoi(d); month = parse_month(m); int parsed_y = atoi(y);
          year = parsed_y + ((parsed_y < 50) ? 2000 : 1900);
        }
        if (strlen(line) >= 66) { char temp_id[5] = {0}; memcpy(temp_id, line+62, 4); trim_string(temp_id); if (strlen(temp_id) > 0 && strcmp(temp_id, "XXXX") != 0) { strcpy(global_id, temp_id); id_from_header = 1; } }
      } else if (f == 0 && strncmp(line, "TITLE ", 6) == 0) {
        char chunk[71]; memset(chunk, ' ', 70); chunk[70] = '\0';
        int copy_len = strlen(line) - 10;
        while(copy_len > 0 && (line[10+copy_len-1] == '\n' || line[10+copy_len-1] == '\r')) copy_len--;
        if (copy_len > 70) copy_len = 70;
        if (copy_len > 0) memcpy(chunk, line+10, copy_len);
        strcat(global_title, chunk);
        if (!id_from_header && !af_mode) { char *p1 = strchr(line, '('), *p2 = strchr(p1 ? p1 : line, ')'); if (p1 && p2 && (p2 > p1 + 4)) { size_t len = p2 - p1 - 1; memcpy(global_id, p1+1, len > 12 ? 12 : len); global_id[len > 12 ? 12 : len] = '\0'; } }
        if (strstr(global_title, "ALPHAFOLD") || strstr(global_title, "AlphaFold")) af_mode = 1;
      } else if (f == 0 && strncmp(line, "COMPND", 6) == 0) {
        char chunk[72] = {0}; int copy_len = strlen(line) - 10;
        while(copy_len > 0 && (line[10+copy_len-1] == '\n' || line[10+copy_len-1] == '\r')) copy_len--;
        memcpy(chunk, line+10, copy_len > 70 ? 70 : copy_len); trim_string(chunk);
        if (strlen(raw_compnd) > 0) { strcat(raw_compnd, " "); }
        strcat(raw_compnd, chunk);
      } else if (f == 0 && strncmp(line, "SOURCE", 6) == 0) {
        char chunk[72] = {0}; int copy_len = strlen(line) - 10;
        while(copy_len > 0 && (line[10+copy_len-1] == '\n' || line[10+copy_len-1] == '\r')) copy_len--;
        memcpy(chunk, line+10, copy_len > 70 ? 70 : copy_len); trim_string(chunk);
        if (strlen(raw_source) > 0) { strcat(raw_source, " "); }
        strcat(raw_source, chunk);
      } else if (f == 0 && strncmp(line, "EXPDTA", 6) == 0) { char chunk[71] = {0}; size_t c_len = strlen(line) - 10; memcpy(chunk, line+10, c_len > 70 ? 70 : c_len); trim_string(chunk); if (strlen(global_expdta) > 0) strcat(global_expdta, " "); strncat(global_expdta, chunk, 127 - strlen(global_expdta));
      } else if (f == 0 && strncmp(line, "REMARK   2 RESOLUTION.", 22) == 0) { float res = atof(line + 22); if (res > 0.0f) global_resolution = res;
      } else if (f == 0 && strncmp(line, "REMARK   3", 10) == 0) {
        char* prog = strstr(line, "PROGRAM"); char* ver = strstr(line, "VERSION"); char* colon = strchr(line, ':');
        if (prog && colon && strlen(global_program) == 0) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; memcpy(temp, start, strlen(start)); trim_string(temp); strncat(global_program, temp, 55); }
        else if (ver && colon && strlen(global_program) > 0 && !strstr(global_program, "v")) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; memcpy(temp, start, strlen(start)); trim_string(temp); strncat(global_program, " v", 55 - strlen(global_program)); strncat(global_program, temp, 55 - strlen(global_program)); }
      } else if (f == 0 && strncmp(line, "DBREF ", 6) == 0) {
        MetaDBRef db; memset(&db, ' ', sizeof(db)); db.chain = line[12];
        memset(tbuf, 0, 16); memcpy(tbuf, line+14, 4); db.seq_begin = atoi(tbuf); db.iCode_begin = line[18] ? line[18] : ' ';
        memset(tbuf, 0, 16); memcpy(tbuf, line+20, 4); db.seq_end = atoi(tbuf); db.iCode_end = line[24] ? line[24] : ' ';

        char tmp_name[9] = {0}; memcpy(tmp_name, line+26, 5); trim_string(tmp_name); memcpy(db.db_name, tmp_name, 8);
        char tmp_acc[13] = {0}; memcpy(tmp_acc, line+33, 8); trim_string(tmp_acc); memcpy(db.db_accession, tmp_acc, 12);
        char tmp_id[13] = {0}; memcpy(tmp_id, line+42, 12); trim_string(tmp_id); memcpy(db.db_id_code, tmp_id, 12);

        memset(tbuf, 0, 16); memcpy(tbuf, line+55, 5); db.db_seq_begin = atoi(tbuf); db.db_iCode_begin = line[60] ? line[60] : ' ';
        memset(tbuf, 0, 16); memcpy(tbuf, line+62, 5); db.db_seq_end = atoi(tbuf); db.db_iCode_end = line[67] ? line[67] : ' ';
        stream_write(&dbrefs, &db, sizeof(MetaDBRef));
        if (f == 0 && !id_from_header && !af_mode && strncmp(db.db_name, "UNP", 3) == 0 && strlen(tmp_acc) > 3) {
          strncpy(global_id, tmp_acc, 12); global_id[12] = '\0';
        }
      } else if (f == 0 && strncmp(line, "NUMMDL", 6) == 0) {
        parsed_num_models = atoi(line + 10);
      } else if (f == 0 && strncmp(line, "MODEL ", 6) == 0) {
        counted_models++;
      } else if (f == 0 && strncmp(line, "SEQRES", 6) == 0) {
        char chain = line[11];
        has_seqres_chain[(unsigned char)chain] = 1;

        if (chain != current_seqres_chain) {
          if (current_seqres_count > 0) {
            // Update the previous chain's residue count in the buffer
            uint16_t* count_ptr = (uint16_t*)(seqres_buf.data + seqres_buf.len - current_seqres_count - 2);
            *count_ptr = current_seqres_count;
          }
          current_seqres_chain = chain;
          current_seqres_count = 0;
          MetaSeqresHeader sh = {chain, 0};
          stream_write(&seqres_buf, &sh, sizeof(MetaSeqresHeader));
        }

        for (int i = 0; i < 13; i++) {
          int offset = 19 + i * 4;
          if (offset + 2 >= (int)strlen(line)) break;
          char res[4] = {0}; memcpy(res, line + offset, 3); trim_string(res);
          if (strlen(res) == 0) continue;

          uint8_t type_id = 0; // Default to UNK if not found
          for (int t = 1; t < 256; t++) {
            if ((UDICT[t].h_end_norm > 0 || UDICT[t].h_end_cterm > 0) && strncmp(UDICT[t].resName, res, 3) == 0) {
              type_id = t; break;
            }
          }
          stream_write(&seqres_buf, &type_id, 1);
          current_seqres_count++;
        }
      } else if (strncmp(line, "SSBOND", 6) == 0) { MetaSSBond b; memset(&b, 0, sizeof(b)); memcpy(b.res1.resName, line+11, 3); b.res1.chain = line[15]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); b.res1.resSeq = atoi(tbuf); b.res1.iCode = line[21]; memcpy(b.res2.resName, line+25, 3); b.res2.chain = line[29]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); b.res2.resSeq = atoi(tbuf); b.res2.iCode = line[35]; if (strlen(line) >= 65) { char tmp[8] = {0}; memcpy(tmp, line+59, 6); if (sscanf(tmp, "%s", b.sym1) != 1) b.sym1[0] = '\0'; } if (strlen(line) >= 72) { char tmp[8] = {0}; memcpy(tmp, line+66, 6); if (sscanf(tmp, "%s", b.sym2) != 1) b.sym2[0] = '\0'; } if (strlen(line) >= 78) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+73, 5); b.length = atof(tbuf); } stream_write(&ssbonds, &b, sizeof(MetaSSBond));
      } else if (strncmp(line, "CISPEP", 6) == 0) { MetaCispep c; memset(&c, 0, sizeof(c)); memcpy(c.res1.resName, line+11, 3); c.res1.chain = line[15]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); c.res1.resSeq = atoi(tbuf); c.res1.iCode = line[21]; memcpy(c.res2.resName, line+25, 3); c.res2.chain = line[29]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); c.res2.resSeq = atoi(tbuf); c.res2.iCode = line[35]; if (strlen(line) >= 46) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+43, 3); c.modNum = atoi(tbuf); } if (strlen(line) >= 59) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+53, 6); c.angle = atof(tbuf); } stream_write(&cispeps, &c, sizeof(MetaCispep));
      } else if (strncmp(line, "HELIX ", 6) == 0) { MetaHelix h; memset(&h, 0, sizeof(h)); memcpy(h.helixID, line+11, 3); trim_string(h.helixID); memcpy(h.start.resName, line+15, 3); h.start.chain = line[19]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+21, 4); h.start.resSeq = atoi(tbuf); h.start.iCode = line[25]; memcpy(h.end.resName, line+27, 3); h.end.chain = line[31]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); h.end.resSeq = atoi(tbuf); h.end.iCode = line[37]; if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); h.helix_class = atoi(tbuf); } if (strlen(line) >= 70) { memcpy(h.comment, line+40, 30); trim_string(h.comment); } if (strlen(line) >= 76) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+71, 5); h.length = atoi(tbuf); } stream_write(&helices, &h, sizeof(MetaHelix));
      } else if (strncmp(line, "SHEET ", 6) == 0) { MetaSheet s; memset(&s, 0, sizeof(s)); memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+7, 3); s.strand_num = atoi(tbuf); memcpy(s.sheetID, line+11, 3); trim_string(s.sheetID); if (strlen(line) >= 16) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+14, 2); s.num_strands = atoi(tbuf); } memcpy(s.start.resName, line+17, 3); s.start.chain = line[21]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+22, 4); s.start.resSeq = atoi(tbuf); s.start.iCode = line[26]; memcpy(s.end.resName, line+28, 3); s.end.chain = line[32]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); s.end.resSeq = atoi(tbuf); s.end.iCode = line[37]; if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); s.sense = atoi(tbuf); } stream_write(&sheets, &s, sizeof(MetaSheet));
      } else if (strncmp(line, "SITE  ", 6) == 0) {
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+7, 3); int seqNum = atoi(tbuf);
        if (seqNum == 1) {
          MetaSiteHeader sh; memset(&sh, 0, sizeof(sh));
          memcpy(sh.siteID, line+11, 3); trim_string(sh.siteID);
          memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+15, 2); sh.num_res = atoi(tbuf);
          stream_write(&sites, &sh, sizeof(MetaSiteHeader));
        }
        for (uint32_t i=0; i<4; i++) {
          int off = 18 + i*11;
          if (off + 9 < (int)strlen(line)) {
            char check[4] = {0}; memcpy(check, line+off, 3); trim_string(check);
            if (strlen(check) == 0) continue;
            PdbResId r; memset(&r, 0, sizeof(r));
            memcpy(r.resName, line+off, 3); r.chain = line[off+4];
            memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+off+5, 4); r.resSeq = atoi(tbuf);
            r.iCode = line[off+9];
            stream_write(&sites, &r, sizeof(PdbResId));
          }
        }
      } else if (strncmp(line, "CONECT", 6) == 0) {
        int len = strlen(line); int32_t vals[5] = {0}, targets = 0;
        for(int c=0; c<5; c++) { int start = 6 + c*5; if (start < len) { char buf[6] = {0}; int copy_len = len - start; if (copy_len > 5) copy_len = 5; memcpy(buf, line + start, copy_len); vals[c] = atoi(buf); if (c > 0 && vals[c] > 0) targets++; } }
        if (vals[0] > 0 && targets > 0) { stream_write(&conects, &targets, 4); stream_write(&conects, &vals[0], 4); for(int c=1; c<=4; c++) if (vals[c] > 0) stream_write(&conects, &vals[c], 4); }
      } else if (strncmp(line, "ENDMDL", 6) == 0) {
        continue;
      } else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
        char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
        char atom_name[5] = {0}; memcpy(atom_name, line + 12, 4);
        char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
        if (!include_hydrogens && (clean_atom[1] == 'H' || clean_atom[0] == 'H' || (line[76] == ' ' && line[77] == 'D'))) continue;

        int type = -1;
        for (int i = 1; i < 256; i++) {
          if ((UDICT[i].h_end_norm > 0 || UDICT[i].h_end_cterm > 0) && strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; }
        }

        if (type != -1) {
          int a_idx = -1;
          for (int a = 0; a < UDICT[type].h_end_norm; a++) {
            if (strncmp(UDICT[type].h_atoms[a], atom_name, 4) == 0) { a_idx = a; break; }
          }
          if (a_idx == -1 && type >= 0x80 && UDICT[type].h_end_norm < 255) {
            a_idx = UDICT[type].h_end_norm++; UDICT[type].h_end_cterm = UDICT[type].h_end_norm;
            memcpy(UDICT[type].h_atoms[a_idx], atom_name, 4); UDICT[type].h_atoms[a_idx][4] = '\0';
            UDICT[type].h_elements[a_idx][0] = line[76]; UDICT[type].h_elements[a_idx][1] = line[77]; UDICT[type].h_elements[a_idx][2] = '\0';
          }
        } else if (custom_het_idx < 192) {
          type = custom_het_idx++; memcpy(UDICT[type].resName, res_name, 3); UDICT[type].resName[3] = '\0';
          UDICT[type].h_end_norm = 1; UDICT[type].h_end_cterm = 1;
          memcpy(UDICT[type].h_atoms[0], atom_name, 4); UDICT[type].h_atoms[0][4] = '\0';
          UDICT[type].h_elements[0][0] = line[76]; UDICT[type].h_elements[0][1] = line[77]; UDICT[type].h_elements[0][2] = '\0';
        }
      }
    }
    fclose(fin);
  }

  if (current_seqres_count > 0) {
    uint16_t* count_ptr = (uint16_t*)(seqres_buf.data + seqres_buf.len - current_seqres_count - 2);
    *count_ptr = current_seqres_count;
  }

  if (af_mode && !id_from_header) {
    char temp[128]; strcpy(temp, in_pdbs[0]);
    char* base = strrchr(temp, '/'); base = base ? base + 1 : temp;
    char* af = strstr(base, "AF-"); if (af) memmove(base, af + 3, strlen(af + 3) + 1);
    char* f = strstr(base, "-F"); if (f) *f = '\0';
    char* dot = strrchr(base, '.'); if (dot) *dot = '\0';
    strncpy(global_id, base, 12); global_id[12] = '\0';
  }

  StreamBuffer compnd_tlv; stream_init(&compnd_tlv);
  if (strlen(raw_compnd) > 0) extract_tlv(raw_compnd, &compnd_tlv, 0, af_mode);

  StreamBuffer source_tlv; stream_init(&source_tlv);
  if (strlen(raw_source) > 0) extract_tlv(raw_source, &source_tlv, 1, af_mode);

  uint32_t num_models = (parsed_num_models > 0) ? parsed_num_models : counted_models;
  if (num_models == 0) num_models = 1;
  uint16_t total_models_minus_1 = num_models - 1;

  uint64_t title_len = strlen(global_title), expdta_len = strlen(global_expdta);
  int has_title = (title_len > 0), has_het = (custom_het_idx > 0x81), has_expdta = (expdta_len > 0);
  int has_res = (global_resolution > 0.0f), has_prog = (strlen(global_program) > 0);

  uint64_t meta_mask = 0; uint64_t zero = 0;
  if (af_mode) meta_mask |= (1ULL << META_AF_MODE);
  if (has_res) meta_mask |= (1ULL << META_RESOLUTION);
  if (has_title) meta_mask |= (1ULL << META_TITLE);
  if (has_expdta) meta_mask |= (1ULL << META_EXPDTA);
  if (has_prog) meta_mask |= (1ULL << META_PROGRAM);
  if (has_het) meta_mask |= (1ULL << META_HET_DICT);
  if (helices.len > 0) meta_mask |= (1ULL << META_HELIX);
  if (sheets.len > 0) meta_mask |= (1ULL << META_SHEET);
  if (ssbonds.len > 0) meta_mask |= (1ULL << META_SSBOND);
  if (cispeps.len > 0) meta_mask |= (1ULL << META_CISPEP);
  if (sites.len > 0) meta_mask |= (1ULL << META_SITE);
  if (include_hydrogens) meta_mask |= (1ULL << META_HYDROGENS);
  if (compnd_tlv.len > 0) meta_mask |= (1ULL << META_COMPND);
  if (source_tlv.len > 0) meta_mask |= (1ULL << META_SOURCE);
  if (dbrefs.len > 0) meta_mask |= (1ULL << META_DBREF);
  if (afdb_mode) meta_mask |= (1ULL << META_AFDB);
  if (seqres_buf.len > 0 && !afdb_mode) meta_mask |= (1ULL << META_SEQRES);

  StreamBuffer meta_buffer; stream_init(&meta_buffer);
  uint64_t meta_header = 0 | (meta_mask << 16);
  stream_write(&meta_buffer, &meta_header, 8);

  if (meta_mask & (1ULL << META_RESOLUTION)) { uint32_t pad=0; stream_write(&meta_buffer, &global_resolution, 4); stream_write(&meta_buffer, &pad, 4); }
  if (meta_mask & (1ULL << META_TITLE)) { stream_write(&meta_buffer, &title_len, 8); stream_write(&meta_buffer, global_title, title_len); size_t pad = (8 - (title_len % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_EXPDTA)) { stream_write(&meta_buffer, &expdta_len, 8); stream_write(&meta_buffer, global_expdta, expdta_len); size_t pad = (8 - (expdta_len % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_PROGRAM)) { uint64_t p_len = 56; stream_write(&meta_buffer, &p_len, 8); stream_write(&meta_buffer, global_program, 56); }
  if (meta_mask & (1ULL << META_HET_DICT)) {
    uint64_t dict_len = 1; for (int i = 0x81; i < custom_het_idx; i++) dict_len += 5 + UDICT[i].h_end_norm * 6;
    stream_write(&meta_buffer, &dict_len, 8);
    uint8_t custom_count = custom_het_idx - 0x81; stream_write(&meta_buffer, &custom_count, 1);
    for (int i = 0x81; i < custom_het_idx; i++) {
      char r_name[4] = {0}; memcpy(r_name, UDICT[i].resName, 3); r_name[3] = ' '; stream_write(&meta_buffer, r_name, 4);
      uint8_t a_count = UDICT[i].h_end_norm; stream_write(&meta_buffer, &a_count, 1);
      for (int a = 0; a < a_count; a++) { stream_write(&meta_buffer, UDICT[i].h_atoms[a], 4); stream_write(&meta_buffer, UDICT[i].h_elements[a], 2); }
    }
    size_t pad = (8 - (dict_len % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad);
  }
  if (meta_mask & (1ULL << META_HELIX))  { uint64_t L = helices.len;    stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, helices.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SHEET))  { uint64_t L = sheets.len;     stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, sheets.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SSBOND)) { uint64_t L = ssbonds.len;    stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, ssbonds.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_CISPEP)) { uint64_t L = cispeps.len;    stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, cispeps.data, L);    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SITE))   { uint64_t L = sites.len;      stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, sites.data, L);      size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_COMPND)) { uint64_t L = compnd_tlv.len; stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, compnd_tlv.data, L); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SOURCE)) { uint64_t L = source_tlv.len; stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, source_tlv.data, L); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_DBREF))  { uint64_t L = dbrefs.len;     stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, dbrefs.data, L);     size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }
  if (meta_mask & (1ULL << META_SEQRES)) { uint64_t L = seqres_buf.len; stream_write(&meta_buffer, &L, 8); stream_write(&meta_buffer, seqres_buf.data, L); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) stream_write(&meta_buffer, &zero, pad); }

  StreamBuffer coords; stream_init(&coords);
  EncoderState st; memset(&st, 0, sizeof(st)); st.type = -1;
  for(int i=0; i<256; i++) { alt_init(&st.h_atoms[i]); alt_init(&st.hy_atoms[i]); }

  ChainBuilder cb; cb_init(&cb); char current_c1 = '\0', current_c2 = '\0';
  uint16_t current_model_idx = 0;

  for (int f = 0; f < num_files; f++) {
    FILE* fin = fopen(in_pdbs[f], "r");
    char line[128];
    while (fgets(line, sizeof(line), fin)) {
      if (strncmp(line, "TER   ", 6) == 0) {
        cb.has_ter_raw = 1;
        if (st.type > 0) flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), 1, include_hydrogens);
        if (cb.num_atoms > 0) {
          unsigned char ch = cb.c1; cb.has_ter = (cb.has_ter_raw || has_seqres_chain[ch]) ? 1 : 0;
          cb.is_n_term = 1; cb.is_c_term = 1;
          flush_chain(&cb, &coords);
        } else { cb.has_ter_raw = 0; }
        st.type = -1; memset(st.h_present, 0, sizeof(st.h_present)); memset(st.hy_present, 0, sizeof(st.hy_present));
        for (int k = 0; k < 256; k++) { st.h_atoms[k].count = 0; st.hy_atoms[k].count = 0; }
        continue;
      }
      if (strncmp(line, "ENDMDL", 6) == 0) {
        if (st.type > 0) flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), 1, include_hydrogens);
        if (cb.num_atoms > 0) {
          unsigned char ch = cb.c1; cb.has_ter = (cb.has_ter_raw || has_seqres_chain[ch]) ? 1 : 0;
          cb.is_n_term = 1; cb.is_c_term = 1;
          flush_chain(&cb, &coords);
        }
        uint64_t end_stream = (0x3FFFULL << 18); stream_write(&coords, &end_stream, 8);
        save_model_file(out_bpdb, total_models_minus_1, current_model_idx, global_id, global_class, day, month, year, &meta_buffer, &coords, &conects);
        coords.len = 0; current_model_idx++;

        st.type = -1; memset(st.h_present, 0, sizeof(st.h_present)); memset(st.hy_present, 0, sizeof(st.hy_present));
        for (int k = 0; k < 256; k++) { st.h_atoms[k].count = 0; st.hy_atoms[k].count = 0; }
        continue;
      }

      if (strncmp(line, "ATOM  ", 6) != 0 && strncmp(line, "HETATM", 6) != 0) continue;

      char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
      char altLoc = line[16], iCode = line[26], c1 = line[20], c2 = line[21];
      int resSeq = atoi(line + 22);

      char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
      if (!include_hydrogens && (clean_atom[1] == 'H' || clean_atom[0] == 'H' || (line[76] == ' ' && line[77] == 'D'))) continue;

      int type = -1;
      for (int i = 1; i < 256; i++) {
        if ((UDICT[i].h_end_norm > 0 || UDICT[i].h_end_cterm > 0) && strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; }
      }
      if (type == -1) continue;

      int new_chain = (c1 != current_c1 || c2 != current_c2);

      char tbuf2[16];
      memset(tbuf2, 0, sizeof(tbuf2)); memcpy(tbuf2, line + 30, 8); int32_t x_int = (int32_t)round(atof(tbuf2) * 1000.0f);
      memset(tbuf2, 0, sizeof(tbuf2)); memcpy(tbuf2, line + 38, 8); int32_t y_int = (int32_t)round(atof(tbuf2) * 1000.0f);
      memset(tbuf2, 0, sizeof(tbuf2)); memcpy(tbuf2, line + 46, 8); int32_t z_int = (int32_t)round(atof(tbuf2) * 1000.0f);
      memset(tbuf2, 0, sizeof(tbuf2)); memcpy(tbuf2, line + 60, 6); float temp = atof(tbuf2);

      float occ = 100.0f;
      int charge = 0;
      if (!af_mode) {
        memset(tbuf2, 0, sizeof(tbuf2)); memcpy(tbuf2, line + 54, 6); occ = atof(tbuf2) * 100.0f;
        if (line[78] >= '0' && line[78] <= '9') { charge = line[78] - '0'; if (line[79] == '-') charge = -charge; }
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
      if (h_idx == -1 && hy_idx == -1) continue;

      if (resSeq != st.resSeq || iCode != st.iCode || new_chain || cb.atom_bytes >= 580000 || type != st.type) {

        int has_oxt = 0;
        if (st.type > 0 && UDICT[st.type].h_end_cterm > UDICT[st.type].h_end_norm) {
          if (st.h_present[UDICT[st.type].h_end_cterm - 1]) has_oxt = 1;
        }

        if (st.type > 0) flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), has_oxt, include_hydrogens);

        if (new_chain || cb.atom_bytes >= 580000 || (resSeq < cb.last_resSeq) || (resSeq - cb.last_resSeq > 100) || has_oxt) {
          if (cb.num_atoms > 0) {
            unsigned char ch = cb.c1; cb.has_ter = (cb.has_ter_raw || has_seqres_chain[ch] || has_oxt) ? 1 : 0;
            cb.is_n_term = 1; cb.is_c_term = has_oxt;
            flush_chain(&cb, &coords);
          }
          if (new_chain) { current_c1 = c1; current_c2 = c2; }
          cb.c1 = c1; cb.c2 = c2;
        }

        st.resSeq = resSeq; st.iCode = iCode; st.type = type;
        memset(st.h_present, 0, sizeof(st.h_present)); memset(st.hy_present, 0, sizeof(st.hy_present));
        for (int k = 0; k < 256; k++) { st.h_atoms[k].count = 0; st.hy_atoms[k].count = 0; }
      } else {
        int is_dup_atom = 0;
        if (h_idx != -1 && st.h_present[h_idx]) {
          for (int c = 0; c < st.h_atoms[h_idx].count; c++) { if (st.h_atoms[h_idx].alt[c] == altLoc) { is_dup_atom = 1; break; } }
        } else if (hy_idx != -1 && st.hy_present[hy_idx]) {
          for (int c = 0; c < st.hy_atoms[hy_idx].count; c++) { if (st.hy_atoms[hy_idx].alt[c] == altLoc) { is_dup_atom = 1; break; } }
        }
        if (is_dup_atom) {
          flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), 0, include_hydrogens);
          memset(st.h_present, 0, sizeof(st.h_present)); memset(st.hy_present, 0, sizeof(st.hy_present));
          for (int k = 0; k < 256; k++) { st.h_atoms[k].count = 0; st.hy_atoms[k].count = 0; }
        }
      }

      if (h_idx != -1) {
        alt_push(&st.h_atoms[h_idx], (int32_t)round(temp * 100.0f), x_int, y_int, z_int, (int32_t)round(occ), charge, altLoc);
        st.h_present[h_idx] = 1;
      } else if (hy_idx != -1) {
        alt_push(&st.hy_atoms[hy_idx], (int32_t)round(temp * 100.0f), x_int, y_int, z_int, (int32_t)round(occ), charge, altLoc);
        st.hy_present[hy_idx] = 1;
      }
    }
    fclose(fin);

    if (f < num_files - 1 && coords.len > 0) {
      if (st.type > 0) flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), 1, include_hydrogens);
      if (cb.num_atoms > 0) {
        unsigned char ch = cb.c1; cb.has_ter = (cb.has_ter_raw || has_seqres_chain[ch]) ? 1 : 0;
        cb.is_n_term = 1; cb.is_c_term = 1;
        flush_chain(&cb, &coords);
      }
      uint64_t frag_stream = (0x3FFEULL << 18); stream_write(&coords, &frag_stream, 8);
      st.type = -1; memset(st.h_present, 0, sizeof(st.h_present)); memset(st.hy_present, 0, sizeof(st.hy_present));
      for (int k = 0; k < 256; k++) { st.h_atoms[k].count = 0; st.hy_atoms[k].count = 0; }
    }
  }

  if (coords.len > 0) {
    if (st.type > 0) flush_residue(&cb, &st, UDICT, (cb.num_residues == 0), 1, include_hydrogens);
    if (cb.num_atoms > 0) {
      unsigned char ch = cb.c1; cb.has_ter = (cb.has_ter_raw || has_seqres_chain[ch]) ? 1 : 0;
      cb.is_n_term = 1; cb.is_c_term = 1;
      flush_chain(&cb, &coords);
    }
    uint64_t end_stream = (0x3FFFULL << 18); stream_write(&coords, &end_stream, 8);
    save_model_file(out_bpdb, total_models_minus_1, current_model_idx, global_id, global_class, day, month, year, &meta_buffer, &coords, &conects);
  }

  cb_free(&cb);
  stream_free(&coords); stream_free(&conects); stream_free(&ssbonds);
  stream_free(&cispeps); stream_free(&helices); stream_free(&sheets);
  stream_free(&sites); stream_free(&meta_buffer);
  stream_free(&compnd_tlv); stream_free(&source_tlv); stream_free(&dbrefs);
  stream_free(&seqres_buf);

  for(int i=0; i<256; i++) { alt_free(&st.h_atoms[i]); alt_free(&st.hy_atoms[i]); }
  free(UDICT);
}
