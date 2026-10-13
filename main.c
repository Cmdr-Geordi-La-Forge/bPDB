#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "bpdb.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static const char* Z_TO_ELEM[128] = {
  "  ", " H", "HE", "LI", "BE", " B", " C", " N", " O", " F", "NE",
  "NA", "MG", "AL", "SI", " P", " S", "CL", "AR", " K", "CA", "SC",
  "TI", " V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE",
  "AS", "SE", "BR", "KR", "RB", "SR", " Y", "ZR", "NB", "MO", "TC",
  "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE",
  "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
  "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", " W", "RE", "OS",
  "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR",
  "RA", "AC", "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK", "CF",
  "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT",
  "DS", "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG",
  "  ", "  ", "  ", "  ", "  ", "  ", "  ", "  ", "  "
};

static const char* M_STR[] = {"   ", "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                              "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

static const char* MOL_FIELD_NAMES[] = {
  "MOLECULE", "CHAIN", "FRAGMENT", "SYNONYM", "EC", "ENGINEERED", "MUTATION", "OTHER_DETAILS",
  "SYNTHETIC", "ORGANISM_SCIENTIFIC", "ORGANISM_COMMON", "ORGANISM_TAXID", "STRAIN", "VARIANT",
  "CELL_LINE", "ATCC", "ORGAN", "TISSUE", "CELL", "ORGANELLE", "SECRETION", "CELLULAR_LOCATION",
  "PLASMID", "GENE", "EXPRESSION_SYSTEM", "EXPRESSION_SYSTEM_COMMON", "EXPRESSION_SYSTEM_TAXID",
  "EXPRESSION_SYSTEM_STRAIN", "EXPRESSION_SYSTEM_VARIANT", "EXPRESSION_SYSTEM_CELL_LINE",
  "EXPRESSION_SYSTEM_ATCC_NUMBER", "EXPRESSION_SYSTEM_ORGAN", "EXPRESSION_SYSTEM_TISSUE",
  "EXPRESSION_SYSTEM_CELL", "EXPRESSION_SYSTEM_ORGANELLE", "EXPRESSION_SYSTEM_CELLULAR_LOCATION",
  "EXPRESSION_SYSTEM_VECTOR_TYPE", "EXPRESSION_SYSTEM_VECTOR", "EXPRESSION_SYSTEM_PLASMID",
  "EXPRESSION_SYSTEM_GENE", "OTHER_DETAILS"
};

static inline int match_key(const char* key, const char* target) {
  size_t len = strlen(target);
  if (strncmp(key, target, len) != 0) return 0;
  while (key[len] != '\0') {
    if (key[len] != ' ') return 0;
    len++;
  }
  return 1;
}

static const char* JRNL_FIELD_NAMES[] = {
  "AUTH", "TITL", "EDIT", "REF", "PUBL", "REFN", "PMID", "DOI"
};

static inline int get_jrnl_field_id(const char* key) {
  for (int i = 0; i < BPDB_JRNL_FIELD_MAX; i++) {
    if (match_key(key, JRNL_FIELD_NAMES[i])) return i;
  }
  return -1;
}

static void encode_parse_title(const char* line, char* title_buf, int* title_len) {
  size_t tlen = strlen(line);
  if (tlen > 10) {
    size_t copy_len = tlen - 10;
    if (line[10 + copy_len - 1] == '\n') copy_len--;

    int start_idx = 10;
    if (line[8] >= '0' && line[8] <= '9' && line[10] == ' ') {
      start_idx++; copy_len--;
    }

    memcpy(title_buf + *title_len, line + start_idx, copy_len);
    title_buf[*title_len + copy_len] = '\n';
    *title_len += copy_len + 1;
    title_buf[*title_len] = '\0'; // Fix: Force null-termination
  }
}

static void encode_parse_compnd_source(const char* line, char* compnd_buf, int* compnd_len, char* source_buf, int* source_len) {
  int is_source = (line[0] == 'S');
  char* dest = is_source ? source_buf : compnd_buf;
  int* dlen = is_source ? source_len : compnd_len;

  size_t len = strlen(line);
  while(len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) len--;

  if (len > 10) {
    size_t start = 10;
    while(start < len && line[start] == ' ') start++;
    if (start < len) {
      memcpy(dest + *dlen, line + start, len - start);
      *dlen += (len - start);
      dest[*dlen] = ' ';
      *dlen += 1;
      dest[*dlen] = '\0'; // Fix: Force null-termination
    }
  }
}

static void encode_parse_text_pool(const char* line, int offset, char* buf, int* len) {
  size_t tlen = strlen(line);
  while (tlen > 0 && (line[tlen-1] == '\n' || line[tlen-1] == '\r')) tlen--;
  if (tlen > (size_t)offset) {
    size_t start = offset;
    while(start < tlen && line[start] == ' ') start++;
    if (start < tlen) {
      if (*len > 0 && buf[*len-1] != ' ') { buf[*len] = ' '; *len += 1; }
      memcpy(buf + *len, line + start, tlen - start);
      *len += (tlen - start);
    }
  }
  buf[*len] = '\0'; // Fix: Force null-termination
}

static void encode_parse_multiline_block(const char* line, int offset, char* buf, int* len) {
  size_t tlen = strlen(line);
  while (tlen > 0 && (line[tlen-1] == '\n' || line[tlen-1] == '\r')) tlen--;
  if (tlen > (size_t)offset) {
    size_t copy_len = tlen - offset;
    memcpy(buf + *len, line + offset, copy_len);
    *len += copy_len;
  }
  buf[*len] = '\n';
  *len += 1;
  buf[*len] = '\0'; // Fix: Force null-termination
}

static void write_jrnl_multiline(char* out_data, size_t* out_idx, const char* record, const char* sub_record, const char* text, int is_remark) {
  const char* p = text;
  int len = strlen(text);
  int line_num = 1;
  while (len > 0) {
    int chunk = len > 60 ? 60 : len;
    if (chunk == 60) {
      int c = chunk;
      while(c > 0 && p[c-1] != ' ' && p[c-1] != '-') c--;
      if (c > 0) chunk = c;
    }

    char tmp[81];
    memset(tmp, ' ', 80);
    tmp[80] = '\n';

    memcpy(tmp, record, 6);
    if (is_remark) tmp[9] = '1';
    memcpy(tmp + 12, sub_record, 4);

    if (line_num > 1) {
      char num[16];
      snprintf(num, sizeof(num), "%2d", line_num);
      memcpy(tmp + 16, num, 2);
    }

    memcpy(tmp + 19, p, chunk);
    memcpy(out_data + (*out_idx * 81), tmp, 81);

    p += chunk;
    len -= chunk;
    while (len > 0 && *p == ' ') { p++; len--; }
    line_num++;
    (*out_idx)++;
  }
}

uint8_t get_z_num(const char* elem) {
  for (int i = 1; i < 119; i++) {
    if (elem[0] == Z_TO_ELEM[i][0] && elem[1] == Z_TO_ELEM[i][1]) return i;
  }
  return 0;
}

bpdb_date_t parse_pdb_date(const char* str) {
  bpdb_date_t d = {0, 0, 0};
  if (!str || strlen(str) < 9) return d;
  d.day = atoi(str);
  const char* m = str + 3;
  if (strncmp(m,"JAN",3)==0) d.month = 1; else if (strncmp(m,"FEB",3)==0) d.month = 2;
  else if (strncmp(m,"MAR",3)==0) d.month = 3; else if (strncmp(m,"APR",3)==0) d.month = 4;
  else if (strncmp(m,"MAY",3)==0) d.month = 5; else if (strncmp(m,"JUN",3)==0) d.month = 6;
  else if (strncmp(m,"JUL",3)==0) d.month = 7; else if (strncmp(m,"AUG",3)==0) d.month = 8;
  else if (strncmp(m,"SEP",3)==0) d.month = 9; else if (strncmp(m,"OCT",3)==0) d.month = 10;
  else if (strncmp(m,"NOV",3)==0) d.month = 11; else if (strncmp(m,"DEC",3)==0) d.month = 12;
  d.year = atoi(str + 7);
  return d;
}

static inline int32_t parse_pdb_coord_int(const char* str) {
  int32_t val = 0;
  int sign = 1;
  for (int i = 0; i < 8; i++) {
    if (str[i] == '-') sign = -1;
    else if (str[i] >= '0' && str[i] <= '9') {
      val = val * 10 + (str[i] - '0');
    }
  }
  return val * sign;
}

int16_t calc_base(double coord) { return (int16_t)(((int32_t)round(coord * 1000.0)) >> 16); }
uint16_t calc_diff(double coord) { return (uint16_t)(((int32_t)round(coord * 1000.0)) & 0xFFFF); }

double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

size_t write_padded(const void* data, size_t len, FILE* out) {
  fwrite(data, 1, len, out);
  size_t rem = len % 16;
  if (rem != 0) {
    size_t pad = 16 - rem;
    uint8_t zeroes[16] = {0};
    fwrite(zeroes, 1, pad, out);
    return len + pad;
  }
  return len;
}

int get_mol_field_id(const char* key, int is_source) {
  for (int i = 0; i < BPDB_MOL_FIELD_MAX; i++) {
    if (match_key(key, MOL_FIELD_NAMES[i])) {
      if (i == BPDB_MOL_FIELD_COMPND_OTHER || i == BPDB_MOL_FIELD_SOURCE_OTHER) {
        return is_source ? BPDB_MOL_FIELD_SOURCE_OTHER : BPDB_MOL_FIELD_COMPND_OTHER;
      }
      return i;
    }
  }
  return -1;
}

const char* get_mol_field_name(int id) {
  if (id == BPDB_MOL_FIELD_COMPND_OTHER || id == BPDB_MOL_FIELD_SOURCE_OTHER) return "OTHER_DETAILS";
  return MOL_FIELD_NAMES[id];
}

void parse_meta_blob(char* blob, int is_source, char* fields[32][BPDB_MOL_FIELD_MAX]) {
  int current_mol = 0;
  char* token = strtok(blob, ";");

  while(token) {
    while(*token == ' ') token++;
    char* colon = strchr(token, ':');

    if (colon) {
      *colon = '\0';
      char* key = token;
      char* val = colon + 1;
      while(*val == ' ') val++;

      char* end = val + strlen(val) - 1;
      while(end > val && (*end == ' ' || *end == '\n' || *end == '\r')) {
        *end = '\0'; end--;
      }

      if (match_key(key, "MOL_ID")) {
        current_mol = atoi(val) - 1;
      } else if (current_mol >= 0 && current_mol < 32) {
        int field_id = get_mol_field_id(key, is_source);
        if (field_id >= 0) {
          if (fields[current_mol][field_id]) {
            char* old = fields[current_mol][field_id];
            char* merged = malloc(strlen(old) + strlen(val) + 3);
            sprintf(merged, "%s, %s", old, val);
            free(old);
            fields[current_mol][field_id] = merged;
          } else {
            fields[current_mol][field_id] = strdup(val);
          }
        }
      }
    }
    token = strtok(NULL, ";");
  }
}

void write_pdb_multiline(char* out_data, size_t* out_idx, const char* record, const char* text, int* line_num) {
  const char* p = text;
  int len = strlen(text);
  while (len > 0) {
    int is_cont = (*line_num > 1);
    int start_idx = is_cont ? 11 : 10;
    int max_chunk = 80 - start_idx;

    int chunk = len > max_chunk ? max_chunk : len;
    if (chunk == max_chunk) {
      int c = chunk;
      while(c > 0 && p[c-1] != ' ' && p[c-1] != ';') c--;
      if (c > 0) chunk = c;
    }
    char tmp[81];
    memset(tmp, ' ', 80);
    tmp[80] = '\n';

    memcpy(tmp, record, strlen(record));
    if (is_cont) {
      char num[16];
      snprintf(num, sizeof(num), "%2d", *line_num);
      memcpy(tmp + 8, num, 2);
    }
    memcpy(tmp + start_idx, p, chunk);
    memcpy(out_data + (*out_idx * 81), tmp, 81);

    p += chunk;
    len -= chunk;
    while (len > 0 && *p == ' ') { p++; len--; }
    (*line_num)++;
    (*out_idx)++;
  }
}

// --- Encoding Handlers ---
static void encode_parse_header(const char* line, bpdb_header_t* hdr, int* has_hdr) {
  memset(hdr->classification, ' ', 40);
  memset(hdr->idCode, ' ', 12);
  size_t clen = MIN(40, strlen(line) > 10 ? strlen(line+10) : 0);
  if (clen > 0 && line[10+clen-1] == '\n') clen--;
  memcpy(hdr->classification, line + 10, clen);

  if (strlen(line) >= 59) hdr->depDate = parse_pdb_date(line + 50);

  size_t id_len = MIN(12, strlen(line) > 62 ? strlen(line+62) : 0);
  if (id_len > 0 && line[62+id_len-1] == '\n') id_len--;
  memcpy(hdr->idCode, line + 62, id_len);

  *has_hdr = 1;
}

static void encode_flush_metadata(FILE* out, uint32_t* dir_table, uint32_t* current_offset,
                                  bpdb_header_t* hdr, int has_hdr,
                                  char* title_buf, int title_len,
                                  char* compnd_buf, int compnd_len,
                                  char* source_buf, int source_len,
                                  char** text_pools, int* text_lens,
                                  uint32_t nummdl, bpdb_caveat_t* caveat_hdr,
                                  bpdb_id_list_t* obslte_hdr, char obslte_ids[64][12],
                                  bpdb_id_list_t* sprsde_hdr, char sprsde_ids[64][12],
                                  char* jrnl_fields[1][BPDB_JRNL_FIELD_MAX],
                                  char* rem1_fields[64][BPDB_JRNL_FIELD_MAX],
                                  int current_rem1,
                                  char** remarks_pools, int* remarks_lens,
                                  uint32_t resolution_mA, int has_remark2,
                                  bpdb_revdat_rec_t* revdats, char revdat_records[64][128][6], int num_revdats) {
  if (has_hdr) {
    dir_table[DIR_HEADER] = *current_offset;
    *current_offset += write_padded(hdr, sizeof(bpdb_header_t), out);
  }
  if (title_len > 0) {
    dir_table[DIR_TITLE] = *current_offset;
    title_buf[title_len] = '\0';
    *current_offset += write_padded(title_buf, title_len + 1, out);
  }

  char* fields[32][BPDB_MOL_FIELD_MAX] = {0};
  if (compnd_len > 0) parse_meta_blob(compnd_buf, 0, fields);
  if (source_len > 0) parse_meta_blob(source_buf, 1, fields);

  int dirs[] = {DIR_KEYWDS, DIR_EXPDTA, DIR_AUTHOR, DIR_MDLTYP, DIR_CAVEAT};
  for(int i=0; i<5; i++) {
    if (text_lens[i] > 0) {
      dir_table[dirs[i]] = *current_offset;
      if (i == 4) *current_offset += write_padded(caveat_hdr, sizeof(bpdb_caveat_t), out);
      text_pools[i][text_lens[i]] = '\0';
      *current_offset += write_padded(text_pools[i], text_lens[i] + 1, out);
    }
  }

  if (has_remark2) {
      dir_table[DIR_REMARK_2] = *current_offset;
      *current_offset += write_padded(&resolution_mA, sizeof(uint32_t), out);
  }

  int num_remarks = 0;
  for (int i = 0; i < 1000; i++) if (remarks_lens[i] > 0) num_remarks++;
  if (num_remarks > 0) {
      dir_table[DIR_REMARKS] = *current_offset;
      *current_offset += write_padded(&num_remarks, sizeof(uint32_t), out);
      uint32_t* rem_idx = calloc(num_remarks * 2, sizeof(uint32_t));
      uint32_t idx_pos = *current_offset;
      *current_offset += write_padded(rem_idx, num_remarks * 8, out);

      int r = 0;
      for (int i = 0; i < 1000; i++) {
          if (remarks_lens[i] > 0) {
              rem_idx[r*2] = i;
              rem_idx[r*2 + 1] = *current_offset;
              remarks_pools[i][remarks_lens[i]] = '\0';
              *current_offset += write_padded(remarks_pools[i], remarks_lens[i] + 1, out);
              r++;
          }
      }
      size_t save = ftell(out);
      fseek(out, idx_pos, SEEK_SET);
      fwrite(rem_idx, 8, num_remarks, out);
      fseek(out, save, SEEK_SET);
      free(rem_idx);
  }

  if (num_revdats > 0) {
      dir_table[DIR_REVDAT] = *current_offset;
      *current_offset += write_padded(&num_revdats, sizeof(uint32_t), out);
      for (int i = 0; i < num_revdats; i++) {
          *current_offset += write_padded(&revdats[i], sizeof(bpdb_revdat_rec_t), out);
          *current_offset += write_padded(revdat_records[i], revdats[i].num_records * 6, out);
      }
  }

  if (nummdl > 0) {
    dir_table[DIR_NUMMDL] = *current_offset;
    *current_offset += write_padded(&nummdl, sizeof(uint32_t), out);
  }

  if (obslte_hdr->count > 0) {
    dir_table[DIR_OBSLTE] = *current_offset;
    *current_offset += write_padded(obslte_hdr, sizeof(bpdb_id_list_t), out);
    *current_offset += write_padded(obslte_ids, obslte_hdr->count * 12, out);
  }
  if (sprsde_hdr->count > 0) {
    dir_table[DIR_SPRSDE] = *current_offset;
    *current_offset += write_padded(sprsde_hdr, sizeof(bpdb_id_list_t), out);
    *current_offset += write_padded(sprsde_ids, sprsde_hdr->count * 12, out);
  }

  char* dict_maps[2][64][BPDB_JRNL_FIELD_MAX] = {0};
  memcpy(dict_maps[0], jrnl_fields, sizeof(char*) * 1 * BPDB_JRNL_FIELD_MAX);
  memcpy(dict_maps[1], rem1_fields, sizeof(char*) * 64 * BPDB_JRNL_FIELD_MAX);
  int map_counts[2] = {1, current_rem1 + 1};
  int map_dirs[2] = {DIR_JRNL, DIR_REMARK_1};

  for (int d = 0; d < 2; d++) {
    if (map_counts[d] <= 0) continue;
    dir_table[map_dirs[d]] = *current_offset;
    bpdb_mol_meta_header_t jhdr = { map_counts[d], BPDB_JRNL_FIELD_MAX };
    *current_offset += write_padded(&jhdr, sizeof(jhdr), out);

    uint32_t j_offs[64] = {0};
    uint32_t j_offs_pos = *current_offset;
    *current_offset += write_padded(j_offs, sizeof(uint32_t) * map_counts[d], out);

    for(int i = 0; i < map_counts[d]; i++) {
      int has_data = 0;
      for(int j=0; j<BPDB_JRNL_FIELD_MAX; j++) if(dict_maps[d][i][j]) has_data = 1;

      if (has_data) {
        j_offs[i] = *current_offset;
        uint32_t f_offs[BPDB_JRNL_FIELD_MAX] = {0};
        uint32_t f_offs_pos = *current_offset;
        *current_offset += write_padded(f_offs, sizeof(uint32_t) * BPDB_JRNL_FIELD_MAX, out);

        for(int j = 0; j < BPDB_JRNL_FIELD_MAX; j++) {
          if (dict_maps[d][i][j]) {
            f_offs[j] = *current_offset;
            *current_offset += write_padded(dict_maps[d][i][j], strlen(dict_maps[d][i][j]) + 1, out);
            free(dict_maps[d][i][j]); 
          }
        }
        size_t save = ftell(out);
        fseek(out, f_offs_pos, SEEK_SET);
        fwrite(f_offs, sizeof(uint32_t), BPDB_JRNL_FIELD_MAX, out);
        fseek(out, save, SEEK_SET);
      }
    }
    size_t save = ftell(out);
    fseek(out, j_offs_pos, SEEK_SET);
    fwrite(j_offs, sizeof(uint32_t), map_counts[d], out);
    fseek(out, save, SEEK_SET);
  }

  int num_mols = 0;
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < BPDB_MOL_FIELD_MAX; j++) {
      if(fields[i][j]) { num_mols = i + 1; break; }
    }
  }

  if (num_mols > 0) {
    dir_table[DIR_COMPND] = *current_offset;
    bpdb_mol_meta_header_t mhdr = { num_mols, BPDB_MOL_FIELD_MAX };
    *current_offset += write_padded(&mhdr, sizeof(mhdr), out);

    uint32_t mol_offs[32] = {0};
    uint32_t mol_offs_pos = *current_offset;
    *current_offset += write_padded(mol_offs, sizeof(uint32_t) * num_mols, out);

    for(int i = 0; i < num_mols; i++) {
      int has_data = 0;
      for(int j=0; j<BPDB_MOL_FIELD_MAX; j++) if(fields[i][j]) has_data = 1;

      if (has_data) {
        mol_offs[i] = *current_offset;
        uint32_t field_offs[BPDB_MOL_FIELD_MAX] = {0};
        uint32_t field_offs_pos = *current_offset;
        *current_offset += write_padded(field_offs, sizeof(uint32_t) * BPDB_MOL_FIELD_MAX, out);

        for(int j = 0; j < BPDB_MOL_FIELD_MAX; j++) {
          if (fields[i][j]) {
            // Skip text representations of organism and expression names
            if (j == BPDB_MOL_FIELD_ORG_SCIENTIFIC ||
                j == BPDB_MOL_FIELD_ORG_COMMON ||
                j == BPDB_MOL_FIELD_EXP_SYSTEM ||
                j == BPDB_MOL_FIELD_EXP_SYS_COMMON) {
              free(fields[i][j]);
              continue;
            }

            field_offs[j] = *current_offset;
            if (j == BPDB_MOL_FIELD_EC) {
              uint32_t ecs[64];
              ecs[0] = 0;
              char* ec_str_copy = strdup(fields[i][j]);
              char* token = strtok(ec_str_copy, ",");
              while (token && ecs[0] < 63) {
                while (*token == ' ') token++;
                ecs[ecs[0] + 1] = parse_ec_string(token);
                ecs[0]++;
                token = strtok(NULL, ",");
              }
              *current_offset += write_padded(ecs, (ecs[0] + 1) * sizeof(uint32_t), out);
              free(ec_str_copy);
            } else if (j == BPDB_MOL_FIELD_ORG_TAXID || j == BPDB_MOL_FIELD_EXP_SYS_TAXID) {
              uint32_t taxids[64];
              taxids[0] = 0;
              char* tax_str_copy = strdup(fields[i][j]);
              char* token = strtok(tax_str_copy, ",");
              while (token && taxids[0] < 63) {
                while (*token == ' ') token++;
                taxids[taxids[0] + 1] = (uint32_t)atoi(token);
                taxids[0]++;
                token = strtok(NULL, ",");
              }
              *current_offset += write_padded(taxids, (taxids[0] + 1) * sizeof(uint32_t), out);
              free(tax_str_copy);
            } else {
              *current_offset += write_padded(fields[i][j], strlen(fields[i][j]) + 1, out);
            }
            free(fields[i][j]);
          }
        }

        size_t save = ftell(out);
        fseek(out, field_offs_pos, SEEK_SET);
        fwrite(field_offs, sizeof(uint32_t), BPDB_MOL_FIELD_MAX, out);
        fseek(out, save, SEEK_SET);
      }
    }
    size_t save = ftell(out);
    fseek(out, mol_offs_pos, SEEK_SET);
    fwrite(mol_offs, sizeof(uint32_t), num_mols, out);
    fseek(out, save, SEEK_SET);
  }
  dir_table[DIR_COORD] = *current_offset;
}

static void encode_parse_coord(const char* line, bpdb_residue_t* res, bpdb_atom_t* atoms, uint32_t* current_offset, FILE* out) {
  if (strncmp(line, "TER", 3) == 0) {
    if (res->atom_count > 0) {
      *current_offset += write_padded(res, sizeof(bpdb_residue_t), out);
      *current_offset += write_padded(atoms, sizeof(bpdb_atom_t) * res->atom_count, out);
    }
    res->file_flags = 2; res->atom_count = 1;
    bpdb_atom_t ter_atom = {0};
    memcpy(ter_atom.name, "TER ", 4);
    *current_offset += write_padded(res, sizeof(bpdb_residue_t), out);
    *current_offset += write_padded(&ter_atom, sizeof(bpdb_atom_t), out);
    res->atom_count = 0;
  } else {
    size_t len = strlen(line);
    uint16_t r_num = atoi(line + 22);
    char chain = line[21];
    double x = atof(line + 30);
    double y = atof(line + 38);
    double z = atof(line + 46);
    int16_t bx = calc_base(x), by = calc_base(y), bz = calc_base(z);

    if (res->atom_count > 0 && (res->res_num != r_num || res->chain_id != chain ||
                                res->base_x != bx || res->base_y != by || res->base_z != bz)) {
      *current_offset += write_padded(res, sizeof(bpdb_residue_t), out);
      *current_offset += write_padded(atoms, sizeof(bpdb_atom_t) * res->atom_count, out);
      res->atom_count = 0;
    }

    if (res->atom_count == 0) {
      memcpy(res->name, line + 17, 3);
      res->chain_id = chain;
      res->res_num = r_num;
      res->icode = line[26];
      res->file_flags = (line[0] == 'A') ? 0 : 1;
      res->base_x = bx; res->base_y = by; res->base_z = bz;
    }

    bpdb_atom_t* a = &atoms[res->atom_count++];
    memcpy(a->name, line + 12, 4);
    a->dx = calc_diff(x); a->dy = calc_diff(y); a->dz = calc_diff(z);
    a->occupancy = (uint16_t)round(atof(line + 54) * 1000.0);
    a->bfactor = (uint16_t)MIN(65535, round(atof(line + 60) * 100.0));

    char elem[3] = "  ";
    if (len >= 78) { elem[0] = line[76]; elem[1] = line[77]; }
    a->z_num = get_z_num(elem);

    int8_t charge = 0;
    if (len >= 80) {
      if (line[78] >= '0' && line[78] <= '9') {
        charge = line[78] - '0';
        if (line[79] == '-') charge = -charge;
      } else if (line[79] >= '0' && line[79] <= '9') {
        charge = line[79] - '0';
        if (line[78] == '-') charge = -charge;
      }
    }
    uint8_t altloc = (line[16] == ' ') ? 0 : (line[16] % 16);
    a->bitfield = ((charge & 0x0F) << 4) | (altloc & 0x0F);
  }
}

void encode_bpdb(const char* in_path, const char* out_path) {
  FILE* in = fopen(in_path, "r");
  if (!in) { printf("Error: Cannot open %s\n", in_path); return; }
  FILE* out = fopen(out_path, "wb");
  if (!out) { printf("Error: Cannot open %s\n", out_path); fclose(in); return; }

  uint32_t dir_table[64] = {0};
  uint32_t current_offset = 256;
  fseek(out, 256, SEEK_SET);

  char line[128];
  bpdb_header_t hdr = {0};
  int has_hdr = 0;

  char title_buf[4096] = {0};
  int title_len = 0;
  int metadata_flushed = 0;

  char* compnd_buf = calloc(1, 65536);
  char* source_buf = calloc(1, 65536);
  char* jrnl_temp_buf = calloc(1, 65536);
  int compnd_len = 0, source_len = 0;

  char* text_pools[5];
  int text_lens[5] = {0};
  for(int i=0; i<5; i++) text_pools[i] = calloc(1, 65536); 

  char** remarks_pools = calloc(1000, sizeof(char*));
  int remarks_lens[1000] = {0};
  uint32_t resolution_mA = 0;
  int has_remark2 = 0;

  bpdb_revdat_rec_t revdats[64] = {0};
  char revdat_records[64][128][6] = {0};
  int num_revdats = 0;

  char* jrnl_fields[1][BPDB_JRNL_FIELD_MAX] = {0};
  char* rem1_fields[64][BPDB_JRNL_FIELD_MAX] = {0};
  int current_rem1 = -1;

  bpdb_id_list_t obslte_hdr = {0}, sprsde_hdr = {0};
  char obslte_ids[64][12] = {0}, sprsde_ids[64][12] = {0};
  bpdb_caveat_t caveat_hdr = {0};
  uint32_t nummdl = 0;

  bpdb_residue_t res = {0};
  bpdb_atom_t* atoms = malloc(sizeof(bpdb_atom_t) * 65535);
  if (!atoms) { printf("Memory allocation failed\n"); exit(1); }

  while (fgets(line, sizeof(line), in)) {
    if (strncmp(line, "HEADER", 6) == 0) {
      encode_parse_header(line, &hdr, &has_hdr);
    } else if (strncmp(line, "TITLE ", 6) == 0) {
      encode_parse_title(line, title_buf, &title_len);
    } else if (strncmp(line, "COMPND", 6) == 0 || strncmp(line, "SOURCE", 6) == 0) {
      encode_parse_compnd_source(line, compnd_buf, &compnd_len, source_buf, &source_len);
    } else if (strncmp(line, "KEYWDS", 6) == 0) { encode_parse_text_pool(line, 10, text_pools[0], &text_lens[0]);
    } else if (strncmp(line, "EXPDTA", 6) == 0) { encode_parse_text_pool(line, 10, text_pools[1], &text_lens[1]);
    } else if (strncmp(line, "AUTHOR", 6) == 0) { encode_parse_text_pool(line, 10, text_pools[2], &text_lens[2]);
    } else if (strncmp(line, "MDLTYP", 6) == 0) { encode_parse_text_pool(line, 10, text_pools[3], &text_lens[3]);
    } else if (strncmp(line, "NUMMDL", 6) == 0) { nummdl = atoi(line + 10);
    } else if (strncmp(line, "REVDAT", 6) == 0) {
      int modNum = atoi(line + 7);
      int r_idx = -1;
      for (int i=0; i<num_revdats; i++) {
        if (revdats[i].modNum == modNum) { r_idx = i; break; }
      }
      if (r_idx == -1) {
        r_idx = num_revdats++;
        revdats[r_idx].modNum = modNum;
        revdats[r_idx].modDate = parse_pdb_date(line + 13);
        memcpy(revdats[r_idx].modId, line + 23, 4);
        revdats[r_idx].modType = line[31] - '0';
        revdats[r_idx].num_records = 0;
      }
      for (int i = 39; i <= 60; i += 7) {
        if (strlen(line) > (size_t)i+3 && line[i] != ' ') {
          memcpy(revdat_records[r_idx][revdats[r_idx].num_records], line + i, 6);
          revdats[r_idx].num_records++;
        }
      }
    } else if (strncmp(line, "CAVEAT", 6) == 0) {
      if (caveat_hdr.idCode[0] == '\0') memcpy(caveat_hdr.idCode, line + 11, 4);
      encode_parse_text_pool(line, 19, text_pools[4], &text_lens[4]);
    } else if (strncmp(line, "OBSLTE", 6) == 0 || strncmp(line, "SPRSDE", 6) == 0) {
      int is_obs = (line[0] == 'O');
      bpdb_id_list_t* lhdr = is_obs ? &obslte_hdr : &sprsde_hdr;
      char (*ids)[12] = is_obs ? obslte_ids : sprsde_ids;
      if (lhdr->count == 0) {
        lhdr->date = parse_pdb_date(line + 11);
        memcpy(lhdr->idCode, line + 21, 4);
      }
      for (int i = 31; i < 71; i += 5) {
        if (strlen(line) > (size_t)i+3 && line[i] != ' ') {
          memcpy(ids[lhdr->count++], line + i, 4);
        }
      }
    } else if (strncmp(line, "JRNL  ", 6) == 0 || strncmp(line, "REMARK", 6) == 0) {
      int is_rem1 = 0;
      int is_jrnl = (line[0] == 'J');
      int rem_num = -1;

      if (!is_jrnl) {
        rem_num = atoi(line + 7);
        if (rem_num == 1) is_rem1 = 1;
      }

      if (is_jrnl || is_rem1) {
        if (is_rem1 && strncmp(line + 11, "REFERENCE", 9) == 0) {
          current_rem1 = atoi(line + 21) - 1;
          if (current_rem1 >= 64) current_rem1 = 63;
        } else {
          char key[5];
          memcpy(key, line + 12, 4); key[4] = '\0';
          char* k_end = key + 3; while(k_end >= key && *k_end == ' ') { *k_end = '\0'; k_end--; }
          int fid = get_jrnl_field_id(key);
          if (fid >= 0) {
            char** target = is_rem1 ? &rem1_fields[current_rem1][fid] : &jrnl_fields[0][fid];
            int clen = target[0] ? strlen(target[0]) : 0;
            
            jrnl_temp_buf[0] = '\0'; 
            if (target[0]) {
                strcpy(jrnl_temp_buf, target[0]);
                free(target[0]);
            }
            encode_parse_text_pool(line, 19, jrnl_temp_buf, &clen); 
            target[0] = strdup(jrnl_temp_buf);
          }
        }
      } else if (rem_num == 2) {
        if (strstr(line, "RESOLUTION.")) {
          if (strstr(line, "NOT APPLICABLE")) resolution_mA = 0;
          else resolution_mA = (uint32_t)(atof(line + 22) * 1000.0);
          has_remark2 = 1;
        }
      } else if (rem_num >= 0 && rem_num < 1000) {
        if (!remarks_pools[rem_num]) remarks_pools[rem_num] = calloc(1, 65536);
        encode_parse_multiline_block(line, 11, remarks_pools[rem_num], &remarks_lens[rem_num]);
      }
    } else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0 || strncmp(line, "TER", 3) == 0) {
      if (!metadata_flushed) {
        encode_flush_metadata(out, dir_table, &current_offset, &hdr, has_hdr, title_buf, title_len, compnd_buf, compnd_len, source_buf, source_len,
                              text_pools, text_lens, nummdl, &caveat_hdr, &obslte_hdr, obslte_ids, &sprsde_hdr, sprsde_ids, jrnl_fields, rem1_fields, current_rem1,
                              remarks_pools, remarks_lens, resolution_mA, has_remark2, revdats, revdat_records, num_revdats);
        metadata_flushed = 1;
      }
      encode_parse_coord(line, &res, atoms, &current_offset, out);
    }
  }
  if (res.atom_count > 0) {
    write_padded(&res, sizeof(bpdb_residue_t), out);
    write_padded(atoms, sizeof(bpdb_atom_t) * res.atom_count, out);
  }

  fseek(out, 0, SEEK_SET);
  fwrite(dir_table, sizeof(uint32_t), 64, out);

  free(atoms);
  fclose(in); fclose(out);
  free(compnd_buf);
  free(source_buf);
  
  for(int i=0; i<5; i++) free(text_pools[i]); 
  for(int i=0; i<1000; i++) if (remarks_pools[i]) free(remarks_pools[i]); 
  free(remarks_pools);
}

// --- Decoding Handlers ---
static size_t count_decoded_lines(const uint8_t* file_data, const uint32_t* dir_table, size_t file_size) {
  size_t total = 0;
  if (dir_table[DIR_HEADER] != 0) total++;
  if (dir_table[DIR_TITLE] != 0) {
    const char* title = (const char*)(file_data + dir_table[DIR_TITLE]);
    while (*title) { if (*title++ == '\n') total++; }
  }
  if (dir_table[DIR_COMPND] != 0) {
    const bpdb_mol_meta_header_t* mhdr = (const bpdb_mol_meta_header_t*)(file_data + dir_table[DIR_COMPND]);
    size_t padded_hdr_size = (sizeof(bpdb_mol_meta_header_t) + 15) & ~15;
    const uint32_t* mol_offs = (const uint32_t*)((const uint8_t*)mhdr + padded_hdr_size);
    for (int i = 0; i < mhdr->num_molecules; i++) {
      if (mol_offs[i] == 0) continue;
      const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);
      for (int j = 0; j < mhdr->num_fields; j++) {
        if (field_offs[j] != 0) {
          if (j == BPDB_MOL_FIELD_EC) {
            const uint32_t* ec_ptr = (const uint32_t*)(file_data + field_offs[j]);
            total += 2 + ((ec_ptr[0] * 16) / 60);
          } else if (j == BPDB_MOL_FIELD_ORG_TAXID || j == BPDB_MOL_FIELD_EXP_SYS_TAXID) {
            const uint32_t* ptr = (const uint32_t*)(file_data + field_offs[j]);
            total += 8 + ((ptr[0] * 64) / 60);
          } else {
            const char* val = (const char*)(file_data + field_offs[j]);
            total += 2 + (strlen(val) / 60);
          }
        }
      }
    }
  }
  
  int text_dirs[] = {DIR_KEYWDS, DIR_EXPDTA, DIR_AUTHOR, DIR_MDLTYP, DIR_CAVEAT};
  for (int i=0; i<5; i++) {
    if (dir_table[text_dirs[i]] != 0) {
      const char* val = (const char*)(file_data + dir_table[text_dirs[i]] + (i==4 ? sizeof(bpdb_caveat_t) : 0));
      total += 2 + (strlen(val) / 60); 
    }
  }
  if (dir_table[DIR_NUMMDL] != 0) total++;
  if (dir_table[DIR_REMARK_2] != 0) total += 2;
    
  if (dir_table[DIR_OBSLTE] != 0) {
      const bpdb_id_list_t* hdr = (const bpdb_id_list_t*)(file_data + dir_table[DIR_OBSLTE]);
      total += 2 + (hdr->count / 10);
  }
  if (dir_table[DIR_SPRSDE] != 0) {
      const bpdb_id_list_t* hdr = (const bpdb_id_list_t*)(file_data + dir_table[DIR_SPRSDE]);
      total += 2 + (hdr->count / 10);
  }

  if (dir_table[DIR_REMARKS] != 0) {
    uint32_t num_rem = *(const uint32_t*)(file_data + dir_table[DIR_REMARKS]);
    // FIX: Read from +16 to account for the write_padded() 16-byte alignment
    const uint32_t* rem_idx = (const uint32_t*)(file_data + dir_table[DIR_REMARKS] + 16); 
    for (uint32_t i = 0; i < num_rem; i++) {
        const char* val = (const char*)(file_data + rem_idx[i*2 + 1]);
        while (*val) { if (*val++ == '\n') total++; }
    }
  }

  if (dir_table[DIR_REVDAT] != 0) {
    uint32_t num_rev = *(const uint32_t*)(file_data + dir_table[DIR_REVDAT]);
    size_t off = dir_table[DIR_REVDAT] + 16; 
    for (uint32_t i = 0; i < num_rev; i++) {
        const bpdb_revdat_rec_t* rev = (const bpdb_revdat_rec_t*)(file_data + off);
        off += (sizeof(bpdb_revdat_rec_t) + 15) & ~15;
        int lines = (rev->num_records + 3) / 4;
        if (lines == 0) lines = 1;
        total += lines;
        off += (rev->num_records * 6 + 15) & ~15;
    }
  }

  int jrnl_dirs[] = {DIR_JRNL, DIR_REMARK_1};
  for (int d = 0; d < 2; d++) {
    if (dir_table[jrnl_dirs[d]] != 0) {
      const bpdb_mol_meta_header_t* mhdr = (const bpdb_mol_meta_header_t*)(file_data + dir_table[jrnl_dirs[d]]);
      size_t padded_hdr_size = (sizeof(bpdb_mol_meta_header_t) + 15) & ~15;
      const uint32_t* mol_offs = (const uint32_t*)((const uint8_t*)mhdr + padded_hdr_size);
      
      for (int i = 0; i < mhdr->num_molecules; i++) {
        if (mol_offs[i] == 0) continue;
        total += 2; 
        const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);
        for (int j = 0; j < mhdr->num_fields; j++) {
          if (field_offs[j] != 0) {
            const char* val = (const char*)(file_data + field_offs[j]);
            total += 2 + (strlen(val) / 60); 
          }
        }
      }
    }
  }
    
  if (dir_table[DIR_COORD] != 0) {
    size_t offset = dir_table[DIR_COORD];
    while (offset < file_size) {
      const bpdb_residue_t* res = (const bpdb_residue_t*)(file_data + offset);
      offset += sizeof(bpdb_residue_t) + (res->atom_count * sizeof(bpdb_atom_t));
      total += res->atom_count;
    }
  }
  return total;
}

static void decode_header(const uint8_t* file_data, const uint32_t* dir_table, char* out_data, size_t* out_idx) {
  if (dir_table[DIR_HEADER] == 0) return;
  const bpdb_header_t* hdr = (const bpdb_header_t*)(file_data + dir_table[DIR_HEADER]);
  int len = snprintf(out_data + (*out_idx * 81), 81, "HEADER    %-40.40s%02d-%3s-%02d   %-12.12s",
                     hdr->classification, hdr->depDate.day,
                     M_STR[hdr->depDate.month <= 12 ? hdr->depDate.month : 0],
                     hdr->depDate.year, hdr->idCode);
  if (len > 80) len = 80;
  out_data[(*out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
  (*out_idx)++;
}

static void decode_title(const uint8_t* file_data, const uint32_t* dir_table, char* out_data, size_t* out_idx) {
  if (dir_table[DIR_TITLE] == 0) return;
  const char* title = (const char*)(file_data + dir_table[DIR_TITLE]);
  int line_num = 1;
  const char* ptr = title;
  while (*ptr) {
    const char* next_nl = strchr(ptr, '\n');
    if (!next_nl) break;
    int str_len = next_nl - ptr;
    int len = (line_num == 1)
      ? snprintf(out_data + (*out_idx * 81), 81, "TITLE     %.*s", str_len, ptr)
      : snprintf(out_data + (*out_idx * 81), 81, "TITLE   %2d%.*s", line_num, str_len, ptr);

    if (len > 80) len = 80;
    out_data[(*out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
    ptr = next_nl + 1;
    line_num++; (*out_idx)++;
  }
}

static void decode_compnd_source(const uint8_t* file_data, const uint32_t* dir_table, char* out_data, size_t* out_idx, const uint8_t* tax_dict) {
  if (dir_table[DIR_COMPND] == 0) return;
  const bpdb_mol_meta_header_t* mhdr = (const bpdb_mol_meta_header_t*)(file_data + dir_table[DIR_COMPND]);
  size_t padded_hdr_size = (sizeof(bpdb_mol_meta_header_t) + 15) & ~15;
  const uint32_t* mol_offs = (const uint32_t*)((const uint8_t*)mhdr + padded_hdr_size);

  // Find absolute last COMPND field across all molecules
  int last_compnd_mol = -1;
  int last_compnd_field = -1;
  for (int i = 0; i < mhdr->num_molecules; i++) {
    if (mol_offs[i] == 0) continue;
    const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);
    for (int j = 0; j <= BPDB_MOL_FIELD_COMPND_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) {
        last_compnd_mol = i;
        last_compnd_field = j;
      }
    }
  }

  int compnd_line = 1;
  for (int i = 0; i < mhdr->num_molecules; i++) {
    if (mol_offs[i] == 0) continue;
    const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);

    int has_compnd = 0;
    for (int j = 0; j <= BPDB_MOL_FIELD_COMPND_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) has_compnd = 1;
    }
    if (!has_compnd) continue;

    char buf[8192] = {0};
    snprintf(buf, sizeof(buf), "MOL_ID: %d;", i + 1);
    write_pdb_multiline(out_data, out_idx, "COMPND", buf, &compnd_line);

    for (int j = 0; j <= BPDB_MOL_FIELD_COMPND_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) {
        int is_last = (i == last_compnd_mol && j == last_compnd_field);
        const char* fmt = is_last ? "%s: %s" : "%s: %s;";

        if (j == BPDB_MOL_FIELD_EC) {
          const uint32_t* ec_ptr = (const uint32_t*)(file_data + field_offs[j]);
          char ec_str[1024];
          unpack_ec_array(ec_ptr, ec_str);
          snprintf(buf, sizeof(buf), fmt, get_mol_field_name(j), ec_str);
        } else {
          const char* val = (const char*)(file_data + field_offs[j]);
          snprintf(buf, sizeof(buf), fmt, get_mol_field_name(j), val);
        }
        write_pdb_multiline(out_data, out_idx, "COMPND", buf, &compnd_line);
      }
    }
  }

  // Find absolute last SOURCE field across all molecules
  int last_source_mol = -1;
  int last_source_field = -1;
  for (int i = 0; i < mhdr->num_molecules; i++) {
    if (mol_offs[i] == 0) continue;
    const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);
    for (int j = BPDB_MOL_FIELD_SYNTHETIC; j <= BPDB_MOL_FIELD_SOURCE_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) {
        last_source_mol = i;
        last_source_field = j;
      }
    }
  }

  int source_line = 1;
  for (int i = 0; i < mhdr->num_molecules; i++) {
    if (mol_offs[i] == 0) continue;
    const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);

    int has_source = 0;
    for (int j = BPDB_MOL_FIELD_SYNTHETIC; j <= BPDB_MOL_FIELD_SOURCE_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) has_source = 1;
    }
    if (!has_source) continue;

    char buf[8192] = {0};
    snprintf(buf, sizeof(buf), "MOL_ID: %d;", i + 1);
    write_pdb_multiline(out_data, out_idx, "SOURCE", buf, &source_line);

    for (int j = BPDB_MOL_FIELD_SYNTHETIC; j <= BPDB_MOL_FIELD_SOURCE_OTHER; j++) {
      if (j < mhdr->num_fields && field_offs[j] != 0) {
        int is_last = (i == last_source_mol && j == last_source_field);
        const char* fmt = is_last ? "%s: %s" : "%s: %s;";

        if (j == BPDB_MOL_FIELD_ORG_TAXID || j == BPDB_MOL_FIELD_EXP_SYS_TAXID) {
          const uint32_t* tax_ptr = (const uint32_t*)(file_data + field_offs[j]);

          if (tax_dict != NULL) {
            uint32_t max_id = *(const uint32_t*)tax_dict;
            const uint32_t* sci_offs = (const uint32_t*)(tax_dict + 4);
            const uint32_t* com_offs = (const uint32_t*)(tax_dict + 4 + (max_id + 1) * 4);
            const char* str_pool = (const char*)(tax_dict + 4 + (max_id + 1) * 8);

            char sci_buf[2048] = {0};
            char com_buf[2048] = {0};

            for (uint32_t k = 0; k < tax_ptr[0]; k++) {
              uint32_t tid = tax_ptr[k+1];
              if (tid <= max_id) {
                if (sci_offs[tid]) {
                  if (sci_buf[0]) strcat(sci_buf, ", ");
                  strcat(sci_buf, str_pool + sci_offs[tid]);
                }
                if (com_offs[tid]) {
                  if (com_buf[0]) strcat(com_buf, ", ");
                  strcat(com_buf, str_pool + com_offs[tid]);
                }
              }
            }

            // These dynamically injected text strings always get semicolons because TAXID follows them
            if (sci_buf[0]) {
              const char* label = (j == BPDB_MOL_FIELD_ORG_TAXID) ? "ORGANISM_SCIENTIFIC" : "EXPRESSION_SYSTEM";
              snprintf(buf, sizeof(buf), "%s: %s;", label, sci_buf);
              write_pdb_multiline(out_data, out_idx, "SOURCE", buf, &source_line);
            }
            if (com_buf[0]) {
              const char* label = (j == BPDB_MOL_FIELD_ORG_TAXID) ? "ORGANISM_COMMON" : "EXPRESSION_SYSTEM_COMMON";
              snprintf(buf, sizeof(buf), "%s: %s;", label, com_buf);
              write_pdb_multiline(out_data, out_idx, "SOURCE", buf, &source_line);
            }
          }

          char tax_str[1024];
          unpack_taxid_array(tax_ptr, tax_str);
          snprintf(buf, sizeof(buf), fmt, get_mol_field_name(j), tax_str);
          write_pdb_multiline(out_data, out_idx, "SOURCE", buf, &source_line);
        } else {
          const char* val = (const char*)(file_data + field_offs[j]);
          snprintf(buf, sizeof(buf), fmt, get_mol_field_name(j), val);
          write_pdb_multiline(out_data, out_idx, "SOURCE", buf, &source_line);
        }
      }
    }
  }
}

static void decode_coord(const uint8_t* file_data, const uint32_t* dir_table, size_t file_size, char* out_data, size_t* out_idx) {
  if (dir_table[DIR_COORD] == 0) return;
  size_t offset = dir_table[DIR_COORD];
  int serial = 1;
  while (offset < file_size) {
    const bpdb_residue_t* res = (const bpdb_residue_t*)(file_data + offset);
    offset += sizeof(bpdb_residue_t);
    uint8_t type = res->file_flags & 0x03;

    for (int i = 0; i < res->atom_count; i++) {
      const bpdb_atom_t* atom = (const bpdb_atom_t*)(file_data + offset);
      offset += sizeof(bpdb_atom_t);

      if (type == 2) {
        int len = snprintf(out_data + (*out_idx * 81), 81, "TER   %5d      %.3s %c%4d%c",
                           serial++, res->name, res->chain_id, res->res_num, res->icode);
        if (len > 80) len = 80;
        out_data[(*out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
        (*out_idx)++;
        continue;
      }

      int32_t ix, iy, iz;
      decode_coords_int(res, atom, &ix, &iy, &iz);

      char altloc = (atom->bitfield & 0x0F) == 0 ? ' ' : (atom->bitfield & 0x0F) + 0x40;
      char charge_str[3] = "  ";
      int8_t charge = atom->bitfield >> 4;
      if (charge & 8) charge -= 16;

      if (charge != 0) {
        charge_str[0] = '0' + abs(charge);
        charge_str[1] = (charge > 0) ? '+' : '-';
      }

      int len = snprintf(out_data + (*out_idx * 81), 81,
                         "%s%5d %.4s%c%.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                         (type == 0) ? "ATOM  " : "HETATM", serial++, atom->name, altloc,
                         res->name, res->chain_id, res->res_num, res->icode,
                         ix * 0.001, iy * 0.001, iz * 0.001,
                         atom->occupancy / 1000.0f, atom->bfactor / 100.0f,
                         Z_TO_ELEM[atom->z_num], charge_str);

      if (len > 80) len = 80;
      out_data[(*out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
      (*out_idx)++;
    }
  }
}

static void decode_text_pool(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx, const char* record_name) {
    if (offset == 0) return;
    int is_caveat = (strcmp(record_name, "CAVEAT") == 0);
    const char* text = (const char*)(file_data + offset + (is_caveat ? sizeof(bpdb_caveat_t) : 0));
    
    char prefix[8192] = {0};
    if (is_caveat) {
        const bpdb_caveat_t* cvt = (const bpdb_caveat_t*)(file_data + offset);
        // Fix: Prepend space so the ID falls on Column 12 where it belongs
        snprintf(prefix, sizeof(prefix), " %-4.4s    %s", cvt->idCode, text); 
        text = prefix;
    }
    
    int line_num = 1;
    write_pdb_multiline(out_data, out_idx, record_name, text, &line_num);
}

static void decode_multiline_block(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx, const char* record_prefix) {
    if (offset == 0) return;
    const char* text = (const char*)(file_data + offset);
    const char* ptr = text;
    while (*ptr) {
        const char* next_nl = strchr(ptr, '\n');
        if (!next_nl) break;
        int str_len = next_nl - ptr;
        // Fix: Use raw string `%s` instead of padding it `%-11s` which artificially shifted text
        int len = snprintf(out_data + (*out_idx * 81), 81, "%s%.*s", record_prefix, str_len, ptr);
        if (len > 80) len = 80;
        out_data[(*out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
        ptr = next_nl + 1;
        (*out_idx)++;
    }
}

static void decode_id_list(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx, const char* record_name) {
    if (offset == 0) return;
    const bpdb_id_list_t* hdr = (const bpdb_id_list_t*)(file_data + offset);
    const char* ids = (const char*)((const uint8_t*)hdr + sizeof(bpdb_id_list_t));
    
    char buf[8192] = {0};
    char* p = buf;
    // Fix: Prepend space so the Date string falls exactly on Column 12
    p += sprintf(p, " %02d-%3s-%02d %-4.4s", hdr->date.day, M_STR[hdr->date.month], hdr->date.year, hdr->idCode);
    for (int i = 0; i < hdr->count; i++) {
        p += sprintf(p, " %-4.4s", ids + (i * 12));
    }
    int line_num = 1;
    write_pdb_multiline(out_data, out_idx, record_name, buf, &line_num);
}

static void decode_jrnl_dict(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx, int is_remark) {
    if (offset == 0) return;
    const bpdb_mol_meta_header_t* mhdr = (const bpdb_mol_meta_header_t*)(file_data + offset);
    size_t padded_hdr_size = (sizeof(bpdb_mol_meta_header_t) + 15) & ~15;
    const uint32_t* mol_offs = (const uint32_t*)((const uint8_t*)mhdr + padded_hdr_size);
    
    const char* rec_name = is_remark ? "REMARK" : "JRNL  ";
    
    for (int i = 0; i < mhdr->num_molecules; i++) {
        if (mol_offs[i] == 0) continue;
        const uint32_t* field_offs = (const uint32_t*)(file_data + mol_offs[i]);
        
        if (is_remark) {
            int len = snprintf(out_data + (*out_idx * 81), 81, "REMARK   1 REFERENCE %d", i + 1);
            out_data[(*out_idx * 81) + len] = ' ';
            (*out_idx)++;
        }
        
        for (int j = 0; j < mhdr->num_fields; j++) {
            if (field_offs[j] != 0) {
                const char* val = (const char*)(file_data + field_offs[j]);
                write_jrnl_multiline(out_data, out_idx, rec_name, JRNL_FIELD_NAMES[j], val, is_remark);
            }
        }
    }
}

static void decode_remark2(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx) {
    if (offset == 0) return;
    uint32_t res_mA = *(const uint32_t*)(file_data + offset);
    
    char buf[81];
    memset(buf, ' ', 80); buf[80] = '\n';
    memcpy(buf, "REMARK   2", 10);
    memcpy(out_data + (*out_idx * 81), buf, 81);
    (*out_idx)++;
    
    memset(buf, ' ', 80); buf[80] = '\n';
    memcpy(buf, "REMARK   2 RESOLUTION.", 22);
    if (res_mA == 0) {
        memcpy(buf + 23, "NOT APPLICABLE.", 15);
    } else {
        char tmp[32];
        snprintf(tmp, sizeof(tmp), "%7.2f ANGSTROMS.", res_mA / 1000.0f);
        memcpy(buf + 23, tmp, strlen(tmp));
    }
    memcpy(out_data + (*out_idx * 81), buf, 81);
    (*out_idx)++;
}

static void decode_remarks_dict(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx) {
    if (offset == 0) return;
    uint32_t num_rem = *(const uint32_t*)(file_data + offset);
    const uint32_t* rem_idx = (const uint32_t*)(file_data + offset + 16);
    
    for (uint32_t i = 0; i < num_rem; i++) {
        char prefix[16];
        snprintf(prefix, sizeof(prefix), "REMARK %3d ", rem_idx[i*2]);
        decode_multiline_block(file_data, rem_idx[i*2 + 1], out_data, out_idx, prefix);
    }
}

static void decode_revdat(const uint8_t* file_data, uint32_t offset, char* out_data, size_t* out_idx) {
    if (offset == 0) return;
    uint32_t num_rev = *(const uint32_t*)(file_data + offset);
    size_t off = offset + 16;
    for (uint32_t i = 0; i < num_rev; i++) {
        const bpdb_revdat_rec_t* rev = (const bpdb_revdat_rec_t*)(file_data + off);
        off += (sizeof(bpdb_revdat_rec_t) + 15) & ~15;
        const char* recs = (const char*)(file_data + off);
        
        int lines = (rev->num_records + 3) / 4;
        if (lines == 0) lines = 1;
        
        for (int l = 0; l < lines; l++) {
            char buf[81];
            memset(buf, ' ', 80); buf[80] = '\n';
            memcpy(buf, "REVDAT ", 7);
            
            char tmp[16];
            snprintf(tmp, sizeof(tmp), "%3d", rev->modNum);
            memcpy(buf + 7, tmp, 3);
            
            if (l > 0) {
                snprintf(tmp, sizeof(tmp), "%2d", l + 1);
                memcpy(buf + 10, tmp, 2);
            }
            
            if (l == 0) {
                snprintf(tmp, sizeof(tmp), "%02d-%3s-%02d", rev->modDate.day, M_STR[rev->modDate.month], rev->modDate.year);
                memcpy(buf + 13, tmp, 9);
                memcpy(buf + 23, rev->modId, 4);
                buf[31] = '0' + rev->modType;
            }
            
            for (int k = 0; k < 4; k++) {
                int r_idx = l * 4 + k;
                if (r_idx < rev->num_records) {
                    memcpy(buf + 39 + k * 7, recs + r_idx * 6, 6);
                    for (int c = 0; c < 6; c++) {
                        if (buf[39 + k * 7 + c] == '\0') buf[39 + k * 7 + c] = ' ';
                    }
                }
            }
            
            memcpy(out_data + (*out_idx * 81), buf, 81);
            (*out_idx)++;
        }
        off += (rev->num_records * 6 + 15) & ~15;
    }
}

void decode_bpdb(const char* in_path, const char* out_path) {
  int fd_in = open(in_path, O_RDONLY);
  if (fd_in < 0) return;
  struct stat sb;
  fstat(fd_in, &sb);
  if (sb.st_size < 256) { close(fd_in); return; }

  const uint8_t* file_data = (const uint8_t*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd_in, 0);
  if (file_data == MAP_FAILED) { close(fd_in); return; }
  const uint32_t* dir_table = (const uint32_t*)file_data;

  size_t total_lines = count_decoded_lines(file_data, dir_table, sb.st_size);
  if (total_lines == 0) { munmap((void*)file_data, sb.st_size); close(fd_in); return; }

  size_t out_size = total_lines * 81;
  char* out_data = malloc(out_size);
  if (!out_data) { munmap((void*)file_data, sb.st_size); close(fd_in); return; }

  const char pad[] = "                                                                                \n";
  for (size_t i = 0; i < total_lines; i++) {
    memcpy(out_data + (i * 81), pad, 81);
  }

  size_t out_idx = 0;
  decode_header(file_data, dir_table, out_data, &out_idx);
  decode_title(file_data, dir_table, out_data, &out_idx);

  const uint8_t* tax_dict = NULL;
  size_t dict_size = 0;
  int fd_dict = open("taxdict.bin", O_RDONLY);
  if (fd_dict >= 0) {
    struct stat sb_dict;
    fstat(fd_dict, &sb_dict);
    dict_size = sb_dict.st_size;
    tax_dict = (const uint8_t*)mmap(NULL, dict_size, PROT_READ, MAP_PRIVATE, fd_dict, 0);
    if (tax_dict == MAP_FAILED) tax_dict = NULL;
    close(fd_dict);
  }

  decode_compnd_source(file_data, dir_table, out_data, &out_idx, tax_dict);

  if (tax_dict) munmap((void*)tax_dict, dict_size);

  decode_id_list(file_data, dir_table[DIR_OBSLTE], out_data, &out_idx, "OBSLTE");
  decode_id_list(file_data, dir_table[DIR_SPRSDE], out_data, &out_idx, "SPRSDE");
    
  decode_text_pool(file_data, dir_table[DIR_CAVEAT], out_data, &out_idx, "CAVEAT");
  decode_text_pool(file_data, dir_table[DIR_KEYWDS], out_data, &out_idx, "KEYWDS");
  decode_text_pool(file_data, dir_table[DIR_EXPDTA], out_data, &out_idx, "EXPDTA");
  decode_text_pool(file_data, dir_table[DIR_AUTHOR], out_data, &out_idx, "AUTHOR");
  decode_text_pool(file_data, dir_table[DIR_MDLTYP], out_data, &out_idx, "MDLTYP");

  decode_revdat(file_data, dir_table[DIR_REVDAT], out_data, &out_idx);
  decode_remark2(file_data, dir_table[DIR_REMARK_2], out_data, &out_idx);
  decode_remarks_dict(file_data, dir_table[DIR_REMARKS], out_data, &out_idx);

  decode_jrnl_dict(file_data, dir_table[DIR_JRNL], out_data, &out_idx, 0);

  if (dir_table[DIR_NUMMDL] != 0) {
    uint32_t n_mdl = *(const uint32_t*)(file_data + dir_table[DIR_NUMMDL]);
    int len = snprintf(out_data + (out_idx * 81), 81, "NUMMDL    %4u", n_mdl);
    out_data[(out_idx * 81) + len] = ' ';
    out_idx++;
  }

  decode_coord(file_data, dir_table, sb.st_size, out_data, &out_idx);

  FILE* out = fopen(out_path, "wb");
  if (out) {
    fwrite(out_data, 1, out_size, out);
    fclose(out);
  }

  free(out_data);
  munmap((void*)file_data, sb.st_size);
  close(fd_in);
}

// --- Benchmark Runner ---
void run_benchmark(const char* pdb_path, const char* bpdb_path) {
  printf("Running benchmark on %s vs %s...\n", pdb_path, bpdb_path);
  int iterations = 100;

  int64_t sanity_pdb_ma = 0;
  FILE* f_txt = fopen(pdb_path, "r");
  if (!f_txt) return;
  char line[128];
  while (fgets(line, sizeof(line), f_txt)) {
    if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
      sanity_pdb_ma += parse_pdb_coord_int(line + 30);
      sanity_pdb_ma += parse_pdb_coord_int(line + 38);
      sanity_pdb_ma += parse_pdb_coord_int(line + 46);
    }
  }
  fclose(f_txt);

  int64_t sanity_bpdb_ma = 0;
  int fd_bpdb = open(bpdb_path, O_RDONLY);
  if (fd_bpdb < 0) return;
  struct stat sb_bpdb;
  fstat(fd_bpdb, &sb_bpdb);
  if (sb_bpdb.st_size < 256) { close(fd_bpdb); return; }

  uint8_t* bpdb_data = (uint8_t*)mmap(NULL, sb_bpdb.st_size, PROT_READ, MAP_PRIVATE, fd_bpdb, 0);
  if (bpdb_data == MAP_FAILED) { close(fd_bpdb); return; }

  uint32_t* dir_table_bpdb = (uint32_t*)bpdb_data;
  size_t sanity_offset = dir_table_bpdb[DIR_COORD] ? dir_table_bpdb[DIR_COORD] : 256;

  while (sanity_offset < (size_t)sb_bpdb.st_size) {
    const bpdb_residue_t* res = (const bpdb_residue_t*)(bpdb_data + sanity_offset);
    sanity_offset += sizeof(bpdb_residue_t);
    for (int a = 0; a < res->atom_count; a++) {
      const bpdb_atom_t* atom = (const bpdb_atom_t*)(bpdb_data + sanity_offset);
      sanity_offset += sizeof(bpdb_atom_t);
      if ((res->file_flags & 0x03) != 2) {
        int32_t ix, iy, iz;
        decode_coords_int(res, atom, &ix, &iy, &iz);
        sanity_bpdb_ma += (ix + iy + iz);
      }
    }
  }

  printf("--- Sanity Check (mÅ Sum) ---\n");
  printf("Legacy PDB:  %lld mÅ\n", (long long)sanity_pdb_ma);
  printf("bPDB Binary: %lld mÅ\n", (long long)sanity_bpdb_ma);
  printf("Match: %s\n\n", (sanity_pdb_ma == sanity_bpdb_ma) ? "YES" : "NO");

  printf("--- Performance (100 Iterations) ---\n");

  int fd_pdb = open(pdb_path, O_RDONLY);
  if (fd_pdb < 0) return;
  struct stat sb_pdb;
  fstat(fd_pdb, &sb_pdb);
  const char* pdb_data = (const char*)mmap(NULL, sb_pdb.st_size, PROT_READ, MAP_PRIVATE, fd_pdb, 0);
  if (pdb_data == MAP_FAILED) { close(fd_pdb); return; }

  double start_pdb = get_time();
  float sink_pdb = 0;

  for (int i = 0; i < iterations; i++) {
    const char* ptr = pdb_data;
    const char* end = pdb_data + sb_pdb.st_size;

    while (ptr < end) {
      const char* next_nl = memchr(ptr, '\n', end - ptr);
      if (!next_nl) next_nl = end;

      if ((next_nl - ptr) >= 54 && strncmp(ptr, "ATOM  ", 6) == 0) {
        sink_pdb += parse_pdb_coord_int(ptr + 30) * 0.001f;
        sink_pdb += parse_pdb_coord_int(ptr + 38) * 0.001f;
        sink_pdb += parse_pdb_coord_int(ptr + 46) * 0.001f;
      }
      ptr = next_nl + 1;
    }
  }
  double time_pdb = get_time() - start_pdb;
  munmap((void*)pdb_data, sb_pdb.st_size);
  close(fd_pdb);

  double start_bpdb = get_time();
  float sink_bpdb = 0;

  for (int i = 0; i < iterations; i++) {
    size_t offset = dir_table_bpdb[DIR_COORD] ? dir_table_bpdb[DIR_COORD] : 256;
    while (offset < (size_t)sb_bpdb.st_size) {
      const bpdb_residue_t* res = (const bpdb_residue_t*)(bpdb_data + offset);
      offset += sizeof(bpdb_residue_t);
      for (int a = 0; a < res->atom_count; a++) {
        const bpdb_atom_t* atom = (const bpdb_atom_t*)(bpdb_data + offset);
        offset += sizeof(bpdb_atom_t);
        if ((res->file_flags & 0x03) == 0) {
          float fx, fy, fz;
          decode_coords_simd(res, atom, &fx, &fy, &fz);
          sink_bpdb += fx + fy + fz;
        }
      }
    }
  }
  double time_bpdb = get_time() - start_bpdb;
  munmap(bpdb_data, sb_bpdb.st_size);
  close(fd_bpdb);

  printf("Legacy PDB Time: %f seconds (Sink: %.1f)\n", time_pdb, sink_pdb);
  printf("bPDB (mmap) Time: %f seconds (Sink: %.1f)\n", time_bpdb, sink_bpdb);
  printf("Speedup: %.2fx\n", time_pdb / time_bpdb);
}

int main(int argc, char** argv) {
  if (argc > 1 && strcmp(argv[1], "bench") == 0) {
    encode_bpdb("4cha.pdb", "4cha.bpdb");
    decode_bpdb("4cha.bpdb", "4cha_back.pdb");
    run_benchmark("4cha.pdb", "4cha.bpdb");
  } else {
    printf("Usage: %s bench\n", argv[0]);
  }
  return 0;
}
