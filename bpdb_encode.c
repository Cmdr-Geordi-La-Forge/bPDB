#include "bpdb.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

void stream_init(StreamBuffer* buf) {
  buf->cap = 1024 * 1024 * 16;
  buf->data = (uint8_t*)malloc(buf->cap);
  buf->len = 0;
}

void stream_write(StreamBuffer* buf, const void* src, size_t size) {
  if (buf->len + size > buf->cap) {
    buf->cap *= 2;
    buf->data = (uint8_t*)realloc(buf->data, buf->cap);
  }
  memcpy(buf->data + buf->len, src, size);
  buf->len += size;
}

void stream_free(StreamBuffer* buf) { free(buf->data); }

int parse_month(const char* m) {
  for (int i = 1; i <= 12; i++) if (strncmp(m, MONTHS[i], 3) == 0) return i;
  return 1;
}

void trim_string(char* str) {
  int i = strlen(str) - 1;
  while(i >= 0 && isspace(str[i])) { str[i] = '\0'; i--; }
}

int get_res_type(const char* name) {
  for (int i = 0; i < NUM_DICT_RES; i++) {
    if (strncmp(name, DICT[i].resName, 3) == 0) return i;
  }
  return -1;
}

typedef struct {
  int resSeq; char iCode; int type;
  int32_t b_acc[256], x_acc[256], y_acc[256], z_acc[256];
  uint8_t occ_acc[256]; char alt_acc[256]; int present[256];
  int has_oxt; int32_t oxt_b, oxt_x, oxt_y, oxt_z;
  uint8_t oxt_occ; char oxt_alt;
} EncoderState;

int chunk_residues = 0; size_t chunk_header_offset = 0;
int chunk_active = 0; uint8_t chunk_flags = 0;
int32_t chunk_anchor_x = 0; int32_t chunk_anchor_y = 0; int32_t chunk_anchor_z = 0;

void close_chunk(StreamBuffer* buf) {
  if (!chunk_active) return;
  uint16_t count = chunk_residues;
  memcpy(buf->data + chunk_header_offset, &count, 2);
  memcpy(buf->data + chunk_header_offset + 2, &chunk_flags, 1);
  chunk_active = 0; chunk_residues = 0; chunk_flags = 0;
}

void open_chunk(StreamBuffer* buf, uint8_t frag_id, char c1, char c2, const char* file, int32_t ax, int32_t ay, int32_t az) {
  close_chunk(buf);
  chunk_header_offset = buf->len;
  chunk_anchor_x = ax; chunk_anchor_y = ay; chunk_anchor_z = az;

  uint8_t header[32] = {0};
  header[3] = frag_id; header[4] = c1; header[5] = c2;
  size_t flen = strlen(file);
  memcpy(header + 6, file, flen > 14 ? 14 : flen);
  memcpy(header + 20, &ax, 4); memcpy(header + 24, &ay, 4); memcpy(header + 28, &az, 4);

  stream_write(buf, header, 32);
  chunk_active = 1; chunk_residues = 0; chunk_flags = 0;
}

static void write_residue_buffer(StreamBuffer* buf, EncoderState* st, UniDict* UDICT) {
  if (st->type == -1) return;
  chunk_residues++;

  int32_t base_b = 0, base_x = 0, base_y = 0, base_z = 0;
  for (int i = 0; i < UDICT[st->type].num_atoms; i++) {
    if (st->present[i]) { base_b = st->b_acc[i]; base_x = st->x_acc[i]; base_y = st->y_acc[i]; base_z = st->z_acc[i]; break; }
  }

  int32_t c_dx = base_x - chunk_anchor_x, c_dy = base_y - chunk_anchor_y, c_dz = base_z - chunk_anchor_z;

  uint8_t res_block[8] = {
    (uint8_t)st->type, (base_b >> 14) & 0xFF, st->resSeq & 0xFF, (c_dx >> 14) & 0xFF,
    (st->resSeq >> 8) & 0xFF, (c_dy >> 14) & 0xFF, st->iCode, (c_dz >> 14) & 0xFF
  };
  stream_write(buf, res_block, 8);

  int32_t last_b = 0, last_x = 0, last_y = 0, last_z = 0;
  for (int i = 0; i < UDICT[st->type].num_atoms; i++) {
    if (!st->present[i]) {
      uint8_t pad[8] = {0}; stream_write(buf, pad, 8);
      continue;
    }

    uint8_t occ = st->occ_acc[i] & 0x7F;
    uint8_t alt_idx = (st->alt_acc[i] == ' ' || st->alt_acc[i] == 0) ? 0 : 1;
    uint16_t L0, L1, L2, L3;

    if (i == 0) {
      L0 = ((base_b & 0x3FFF) << 2) | (occ & 0x03);
      L1 = ((c_dx & 0x3FFF) << 2) | ((occ >> 2) & 0x03);
      L2 = ((c_dy & 0x3FFF) << 2) | ((occ >> 4) & 0x03);
      L3 = ((c_dz & 0x3FFF) << 2) | ((occ >> 6) & 0x01) | ((alt_idx & 0x01) << 1);
    } else {
      int32_t db = st->b_acc[i] - last_b; int32_t dx = st->x_acc[i] - last_x;
      int32_t dy = st->y_acc[i] - last_y; int32_t dz = st->z_acc[i] - last_z;

      if (db > 8191) db = 8191; else if (db < -8192) db = -8192;
      if (dx > 8191) dx = 8191; else if (dx < -8192) dx = -8192;
      if (dy > 8191) dy = 8191; else if (dy < -8192) dy = -8192;
      if (dz > 8191) dz = 8191; else if (dz < -8192) dz = -8192;

      L0 = ((db & 0x3FFF) << 2) | (occ & 0x03);
      L1 = ((dx & 0x3FFF) << 2) | ((occ >> 2) & 0x03);
      L2 = ((dy & 0x3FFF) << 2) | ((occ >> 4) & 0x03);
      L3 = ((dz & 0x3FFF) << 2) | ((occ >> 6) & 0x01) | ((alt_idx & 0x01) << 1);
    }
    last_b = st->b_acc[i]; last_x = st->x_acc[i]; last_y = st->y_acc[i]; last_z = st->z_acc[i];

    uint8_t atom_block[8] = { L0 & 0xFF, L0 >> 8, L1 & 0xFF, L1 >> 8, L2 & 0xFF, L2 >> 8, L3 & 0xFF, L3 >> 8 };
    stream_write(buf, atom_block, 8);
  }
}

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb) {
  UniDict* UDICT = (UniDict*)calloc(256, sizeof(UniDict));
  int num_udict = NUM_DICT_RES;
  for (int i = 0; i < NUM_DICT_RES; i++) {
    strcpy(UDICT[i].resName, DICT[i].resName); UDICT[i].num_atoms = DICT[i].num_atoms;
    for (int a = 0; a < DICT[i].num_atoms; a++) {
      strcpy(UDICT[i].atom_names[a], DICT[i].atom_names[a]); strcpy(UDICT[i].elements[a], DICT[i].elements[a]);
    }
  }

  char global_title[4096] = {0}; char global_id[128] = {0}; char global_class[128] = {0};
  char global_expdta[128] = {0}; char global_program[56] = {0};
  float global_resolution = -1.0f;
  int day=1, month=1, year=1973, id_from_header = 0;

  StreamBuffer ssbonds; stream_init(&ssbonds);
  StreamBuffer cispeps; stream_init(&cispeps);
  StreamBuffer helices; stream_init(&helices);
  StreamBuffer sheets;  stream_init(&sheets);
  StreamBuffer sites;   stream_init(&sites);

  for (int f = 0; f < num_files; f++) {
    FILE* fin = fopen(in_pdbs[f], "r");
    if (!fin) { perror("Input file error"); exit(1); }
    char line[128];
    char tbuf[16];

    while (fgets(line, sizeof(line), fin)) {
      if (f == 0 && strncmp(line, "HEADER", 6) == 0) {
        if (strlen(line) >= 50) { memcpy(global_class, line + 10, 40); trim_string(global_class); }
        if (strlen(line) >= 59) {
          char d[3]={0}, m[4]={0}, y[3]={0};
          memcpy(d, line+50, 2); memcpy(m, line+53, 3); memcpy(y, line+57, 2);
          day = atoi(d); month = parse_month(m); int parsed_y = atoi(y);
          year = parsed_y + ((parsed_y < 50) ? 2000 : 1900);
        }
        if (strlen(line) >= 66) {
          char temp_id[5] = {0}; memcpy(temp_id, line+62, 4); trim_string(temp_id);
          if (strlen(temp_id) > 0 && strcmp(temp_id, "XXXX") != 0) { strcpy(global_id, temp_id); id_from_header = 1; }
        }
      } else if (f == 0 && strncmp(line, "TITLE ", 6) == 0) {
        if (strlen(line) > 10) {
          char chunk[71] = {0}; size_t c_len = strlen(line) - 10;
          memcpy(chunk, line+10, c_len > 70 ? 70 : c_len); trim_string(chunk);
          if (strlen(global_title) > 0) strcat(global_title, " ");
          strcat(global_title, chunk);
        }
        if (!id_from_header) {
          char *p1 = strchr(line, '('), *p2 = strchr(p1 ? p1 : line, ')');
          if (p1 && p2 && (p2 > p1 + 4)) {
            size_t len = p2 - p1 - 1; memcpy(global_id, p1+1, len > 14 ? 14 : len); global_id[len > 14 ? 14 : len] = '\0';
          }
        }
      } else if (f == 0 && strncmp(line, "EXPDTA", 6) == 0) {
        char chunk[71] = {0}; size_t c_len = strlen(line) - 10;
        memcpy(chunk, line+10, c_len > 70 ? 70 : c_len); trim_string(chunk);
        if (strlen(global_expdta) > 0) strcat(global_expdta, " ");
        strncat(global_expdta, chunk, 127 - strlen(global_expdta));
      } else if (f == 0 && strncmp(line, "REMARK   2 RESOLUTION.", 22) == 0) {
        float res = atof(line + 22);
        if (res > 0.0f) global_resolution = res;
      } else if (f == 0 && strncmp(line, "REMARK   3", 10) == 0) {
        char* prog = strstr(line, "PROGRAM"); char* ver = strstr(line, "VERSION"); char* colon = strchr(line, ':');
        if (prog && colon && strlen(global_program) == 0) {
          char temp[64] = {0}; char* start = colon + 1;
          while (*start == ' ' || *start == '\t') start++;
          memcpy(temp, start, strlen(start)); trim_string(temp); strncat(global_program, temp, 55);
        } else if (ver && colon && strlen(global_program) > 0 && !strstr(global_program, "v")) {
          char temp[64] = {0}; char* start = colon + 1;
          while (*start == ' ' || *start == '\t') start++;
          memcpy(temp, start, strlen(start)); trim_string(temp);
          strncat(global_program, " v", 55 - strlen(global_program)); strncat(global_program, temp, 55 - strlen(global_program));
        }
      } else if (strncmp(line, "SSBOND", 6) == 0) {
        MetaSSBond b; memset(&b, 0, sizeof(b));
        memcpy(b.res1.resName, line+11, 3); b.res1.chain = line[15];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); b.res1.resSeq = atoi(tbuf); b.res1.iCode = line[21];
        memcpy(b.res2.resName, line+25, 3); b.res2.chain = line[29];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); b.res2.resSeq = atoi(tbuf); b.res2.iCode = line[35];
        if (strlen(line) >= 65) {
          char tmp[8] = {0}; memcpy(tmp, line+59, 6);
          if (sscanf(tmp, "%s", b.sym1) != 1) b.sym1[0] = '\0';
        }
        if (strlen(line) >= 72) {
          char tmp[8] = {0}; memcpy(tmp, line+66, 6);
          if (sscanf(tmp, "%s", b.sym2) != 1) b.sym2[0] = '\0';
        }
        if (strlen(line) >= 78) {
          memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+73, 5); b.length = atof(tbuf);
        }
        stream_write(&ssbonds, &b, sizeof(MetaSSBond));
      } else if (strncmp(line, "CISPEP", 6) == 0) {
        MetaCispep c; memset(&c, 0, sizeof(c));
        memcpy(c.res1.resName, line+11, 3); c.res1.chain = line[15];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); c.res1.resSeq = atoi(tbuf); c.res1.iCode = line[21];
        memcpy(c.res2.resName, line+25, 3); c.res2.chain = line[29];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); c.res2.resSeq = atoi(tbuf); c.res2.iCode = line[35];
        if (strlen(line) >= 46) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+43, 3); c.modNum = atoi(tbuf); }
        if (strlen(line) >= 59) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+53, 6); c.angle = atof(tbuf); }
        stream_write(&cispeps, &c, sizeof(MetaCispep));
      } else if (strncmp(line, "HELIX ", 6) == 0) {
        MetaHelix h; memset(&h, 0, sizeof(h));
        memcpy(h.helixID, line+11, 3); trim_string(h.helixID);
        memcpy(h.start.resName, line+15, 3); h.start.chain = line[19];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+21, 4); h.start.resSeq = atoi(tbuf); h.start.iCode = line[25];
        memcpy(h.end.resName, line+27, 3); h.end.chain = line[31];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); h.end.resSeq = atoi(tbuf); h.end.iCode = line[37];
        if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); h.helix_class = atoi(tbuf); }
        if (strlen(line) >= 70) { memcpy(h.comment, line+40, 30); trim_string(h.comment); }
        if (strlen(line) >= 76) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+71, 5); h.length = atoi(tbuf); }
        stream_write(&helices, &h, sizeof(MetaHelix));
      } else if (strncmp(line, "SHEET ", 6) == 0) {
        MetaSheet s; memset(&s, 0, sizeof(s));
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+7, 3); s.strand_num = atoi(tbuf);
        memcpy(s.sheetID, line+11, 3); trim_string(s.sheetID);
        if (strlen(line) >= 16) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+14, 2); s.num_strands = atoi(tbuf); }
        memcpy(s.start.resName, line+17, 3); s.start.chain = line[21];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+22, 4); s.start.resSeq = atoi(tbuf); s.start.iCode = line[26];
        memcpy(s.end.resName, line+28, 3); s.end.chain = line[32];
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); s.end.resSeq = atoi(tbuf); s.end.iCode = line[37];
        if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); s.sense = atoi(tbuf); }
        stream_write(&sheets, &s, sizeof(MetaSheet));
      } else if (strncmp(line, "SITE  ", 6) == 0) {
        MetaSiteHeader sh; memset(&sh, 0, sizeof(sh));
        memcpy(sh.siteID, line+11, 3); trim_string(sh.siteID);
        memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+15, 2); sh.num_res = atoi(tbuf);
        stream_write(&sites, &sh, sizeof(MetaSiteHeader));
        for (uint32_t i=0; i<sh.num_res; i++) {
          PdbResId r; memset(&r, 0, sizeof(r));
          int off = 18 + i*11;
          if (off + 9 < (int)strlen(line)) {
            memcpy(r.resName, line+off, 3); r.chain = line[off+4];
            memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+off+5, 4); r.resSeq = atoi(tbuf);
            r.iCode = line[off+9];
            stream_write(&sites, &r, sizeof(PdbResId));
          }
        }
      } else if (strncmp(line, "DBREF ", 6) == 0) {
        char* unp = strstr(line, "UNP");
        if (unp) {
          char dummy[16], uid[64];
          if (sscanf(unp, "%s %s", dummy, uid) == 2 && strlen(uid) > 3 && f == 0 && !id_from_header) {
            size_t c = strlen(uid); memcpy(global_id, uid, c > 14 ? 14 : c); global_id[c > 14 ? 14 : c] = '\0';
          }
        }
      } else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
        char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
        char atom_name[5] = {0}; memcpy(atom_name, line + 12, 4);
        char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
        if (clean_atom[1] == 'H' || clean_atom[0] == 'H' || strncmp(clean_atom, " OXT", 4) == 0) continue;

        int type = -1;
        for (int i = 0; i < num_udict; i++) { if (strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; } }
        if (type == -1 && num_udict < 255) {
          type = num_udict++;
          memcpy(UDICT[type].resName, res_name, 3); UDICT[type].num_atoms = 0;
        }
        if (type != -1) {
          int a_idx = -1;
          for (int a = 0; a < UDICT[type].num_atoms; a++) { if (strncmp(UDICT[type].atom_names[a], atom_name, 4) == 0) { a_idx = a; break; } }
          if (a_idx == -1 && UDICT[type].num_atoms < 255) {
            a_idx = UDICT[type].num_atoms++;
            memcpy(UDICT[type].atom_names[a_idx], atom_name, 4);
            UDICT[type].elements[a_idx][0] = line[76]; UDICT[type].elements[a_idx][1] = line[77];
          }
        }
      }
    }
    fclose(fin);
  }

  StreamBuffer coords; stream_init(&coords);
  StreamBuffer conects; stream_init(&conects);
  EncoderState st; memset(&st, 0, sizeof(st)); st.type = -1;
  char current_c1 = '\0', current_c2 = '\0';
  int force_new_chunk = 1;

  for (int f = 0; f < num_files; f++) {
    FILE* fin = fopen(in_pdbs[f], "r");
    char line[128], local_id[16] = {0}; uint8_t local_frag = 0;
    const char* f_ptr = strstr(in_pdbs[f], "-F");
    if (f_ptr && isdigit(f_ptr[2])) local_frag = atoi(f_ptr + 2);
    force_new_chunk = 1;

    while (fgets(line, sizeof(line), fin)) {
      if (strncmp(line, "CONECT", 6) == 0) {
        int len = strlen(line); int32_t vals[5] = {0}, targets = 0;
        for(int c=0; c<5; c++) {
          int start = 6 + c*5;
          if (start < len) {
            char buf[6] = {0}; int copy_len = len - start; if (copy_len > 5) copy_len = 5;
            memcpy(buf, line + start, copy_len); vals[c] = atoi(buf);
            if (c > 0 && vals[c] > 0) targets++;
          }
        }
        if (vals[0] > 0 && targets > 0) {
          stream_write(&conects, &targets, 4); stream_write(&conects, &vals[0], 4);
          for(int c=1; c<=4; c++) if (vals[c] > 0) stream_write(&conects, &vals[c], 4);
        }
        continue;
      }

      if (strncmp(line, "TER   ", 6) == 0) {
        write_residue_buffer(&coords, &st, UDICT);
        if (st.has_oxt) {
          chunk_flags |= CTRL_TER_OXT;
          int32_t base_b = 0, base_x = 0, base_y = 0, base_z = 0;
          for (int i = 0; i < UDICT[st.type].num_atoms; i++) {
            if (st.present[i]) { base_b = st.b_acc[i]; base_x = st.x_acc[i]; base_y = st.y_acc[i]; base_z = st.z_acc[i]; break; }
          }
          int32_t db = st.oxt_b - base_b, dx = st.oxt_x - base_x, dy = st.oxt_y - base_y, dz = st.oxt_z - base_z;
          if (db > 8191) db = 8191; else if (db < -8192) db = -8192;
          if (dx > 8191) dx = 8191; else if (dx < -8192) dx = -8192;
          if (dy > 8191) dy = 8191; else if (dy < -8192) dy = -8192;
          if (dz > 8191) dz = 8191; else if (dz < -8192) dz = -8192;
          uint8_t alt = (st.oxt_alt != ' ' && st.oxt_alt != 0) ? 1 : 0;
          uint16_t L0 = ((db & 0x3FFF) << 2) | (st.oxt_occ & 0x03);
          uint16_t L1 = ((dx & 0x3FFF) << 2) | ((st.oxt_occ >> 2) & 0x03);
          uint16_t L2 = ((dy & 0x3FFF) << 2) | ((st.oxt_occ >> 4) & 0x03);
          uint16_t L3 = ((dz & 0x3FFF) << 2) | ((st.oxt_occ >> 6) & 0x01) | ((alt & 0x01) << 1);
          uint8_t oxt_block[8] = { L0 & 0xFF, L0 >> 8, L1 & 0xFF, L1 >> 8, L2 & 0xFF, L2 >> 8, L3 & 0xFF, L3 >> 8 };
          stream_write(&coords, oxt_block, 8);
        } else {
          chunk_flags |= CTRL_TER;
        }
        close_chunk(&coords);
        st.type = -1; st.has_oxt = 0; force_new_chunk = 1;
        continue;
      }

      if (strncmp(line, "ATOM  ", 6) != 0 && strncmp(line, "HETATM", 6) != 0) continue;

      char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
      char atom_name[5] = {0}; memcpy(atom_name, line + 12, 4);
      char altLoc = line[16], iCode = line[26], c1 = line[20], c2 = line[21];
      int resSeq = atoi(line + 22);

      char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
      if (clean_atom[1] == 'H' || clean_atom[0] == 'H') continue;

      int type = -1;
      for (int i = 0; i < num_udict; i++) { if (strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; } }
      if (type == -1) continue;

      if (strlen(local_id) == 0) { size_t g_len = strlen(global_id); memcpy(local_id, global_id, g_len > 14 ? 14 : g_len); }
      int new_chain = (c1 != current_c1 || c2 != current_c2);
      int32_t x_int = (int32_t)round(atof(line + 30) * 1000.0f), y_int = (int32_t)round(atof(line + 38) * 1000.0f), z_int = (int32_t)round(atof(line + 46) * 1000.0f);
      float occ = atof(line + 54), temp = atof(line + 60);

      int bounds_exceeded = 0;
      if (!force_new_chunk && chunk_active) {
        int32_t test_dx = x_int - chunk_anchor_x;
        if (test_dx > 2000000 || test_dx < -2000000) bounds_exceeded = 1;
      }

      if (resSeq != st.resSeq || iCode != st.iCode || new_chain || force_new_chunk || chunk_residues >= 60000 || bounds_exceeded || type != st.type) {
        write_residue_buffer(&coords, &st, UDICT);
        if (new_chain || force_new_chunk || chunk_residues >= 60000 || bounds_exceeded) {
          if (new_chain) { current_c1 = c1; current_c2 = c2; }
          open_chunk(&coords, local_frag, c1, c2, local_id, x_int, y_int, z_int);
          force_new_chunk = 0;
        }
        st.resSeq = resSeq; st.iCode = iCode; st.type = type;
        memset(st.present, 0, sizeof(st.present)); st.has_oxt = 0;
      }

      if (strncmp(clean_atom, " OXT", 4) == 0) {
        st.has_oxt = 1; st.oxt_b = (int32_t)round(temp * 100.0f);
        st.oxt_x = x_int; st.oxt_y = y_int; st.oxt_z = z_int;
        st.oxt_occ = (uint8_t)round(occ * 100.0f); st.oxt_alt = altLoc;
        continue;
      }

      int atom_idx = -1;
      for (int i = 0; i < UDICT[type].num_atoms; i++) { if (strncmp(clean_atom, UDICT[type].atom_names[i], 4) == 0) { atom_idx = i; break; } }
      if (atom_idx == -1 || st.present[atom_idx]) continue;

      st.b_acc[atom_idx] = (int32_t)round(temp * 100.0f);
      st.x_acc[atom_idx] = x_int; st.y_acc[atom_idx] = y_int; st.z_acc[atom_idx] = z_int;
      st.occ_acc[atom_idx] = (uint8_t)round(occ * 100.0f);
      st.alt_acc[atom_idx] = altLoc; st.present[atom_idx] = 1;
    }
    fclose(fin);
  }
  write_residue_buffer(&coords, &st, UDICT);
  close_chunk(&coords);

  FILE* fout = fopen(out_bpdb, "wb");
  uint64_t zero = 0;

  uint8_t header[64] = {0};
  memcpy(header, "bPDB", 4); header[4] = 1; header[5] = (uint8_t)day; header[6] = (uint8_t)month;
  header[7] = (year >= 1973 && year <= 2228) ? (uint8_t)(year - 1973) : 0;
  memcpy(header + 8, global_id, strlen(global_id) > 15 ? 15 : strlen(global_id));
  memcpy(header + 24, global_class, strlen(global_class) > 40 ? 40 : strlen(global_class));
  fwrite(header, 1, 64, fout);

  uint64_t title_len = strlen(global_title), expdta_len = strlen(global_expdta);
  int has_title = (title_len > 0), has_het = (num_udict > NUM_DICT_RES), has_expdta = (expdta_len > 0);
  int has_res = (global_resolution > 0.0f), has_prog = (strlen(global_program) > 0);

  uint16_t num_meta = has_title + has_het + has_expdta + has_res + has_prog;
  if (helices.len > 0) num_meta++;
  if (sheets.len > 0) num_meta++;
  if (ssbonds.len > 0) num_meta++;
  if (cispeps.len > 0) num_meta++;
  if (sites.len > 0) num_meta++;

  fwrite(&num_meta, 1, 2, fout);

  if (has_title) { uint16_t code = META_TITLE; fwrite(&code, 1, 2, fout); }
  if (has_het) { uint16_t code = META_HET_DICT; fwrite(&code, 1, 2, fout); }
  if (has_expdta) { uint16_t code = META_EXPDTA; fwrite(&code, 1, 2, fout); }
  if (has_res) { uint16_t code = META_RESOLUTION; fwrite(&code, 1, 2, fout); }
  if (has_prog) { uint16_t code = META_PROGRAM; fwrite(&code, 1, 2, fout); }
  if (helices.len > 0) { uint16_t code = META_HELIX; fwrite(&code, 1, 2, fout); }
  if (sheets.len > 0) { uint16_t code = META_SHEET; fwrite(&code, 1, 2, fout); }
  if (ssbonds.len > 0) { uint16_t code = META_SSBOND; fwrite(&code, 1, 2, fout); }
  if (cispeps.len > 0) { uint16_t code = META_CISPEP; fwrite(&code, 1, 2, fout); }
  if (sites.len > 0) { uint16_t code = META_SITE; fwrite(&code, 1, 2, fout); }

  size_t cpad = (8 - ((2 + num_meta * 2) % 8)) % 8;
  if (cpad > 0) fwrite(&zero, 1, cpad, fout);

  if (has_title) {
    fwrite(&title_len, 1, 8, fout); fwrite(global_title, 1, title_len, fout);
    size_t tpad = (8 - (title_len % 8)) % 8; if (tpad > 0) fwrite(&zero, 1, tpad, fout);
  }
  if (has_het) {
    uint64_t dict_len = 1;
    for (int i = NUM_DICT_RES; i < num_udict; i++) dict_len += 5 + UDICT[i].num_atoms * 6;
    fwrite(&dict_len, 1, 8, fout);
    uint8_t custom_count = num_udict - NUM_DICT_RES; fwrite(&custom_count, 1, 1, fout);
    for (int i = NUM_DICT_RES; i < num_udict; i++) {
      char r_name[4] = {0}; memcpy(r_name, UDICT[i].resName, 3); r_name[3] = ' '; fwrite(r_name, 1, 4, fout);
      uint8_t a_count = UDICT[i].num_atoms; fwrite(&a_count, 1, 1, fout);
      for (int a = 0; a < a_count; a++) { fwrite(UDICT[i].atom_names[a], 1, 4, fout); fwrite(UDICT[i].elements[a], 1, 2, fout); }
    }
    size_t dpad = (8 - (dict_len % 8)) % 8; if (dpad > 0) fwrite(&zero, 1, dpad, fout);
  }
  if (has_expdta) {
    fwrite(&expdta_len, 1, 8, fout); fwrite(global_expdta, 1, expdta_len, fout);
    size_t pad = (8 - (expdta_len % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout);
  }
  if (has_res) {
    uint64_t r_len = 4; fwrite(&r_len, 1, 8, fout);
    fwrite(&global_resolution, 1, 4, fout);
    uint32_t pad = 0; fwrite(&pad, 1, 4, fout);
  }
  if (has_prog) {
    uint64_t p_len = 56; fwrite(&p_len, 1, 8, fout);
    fwrite(global_program, 1, 56, fout);
  }

  if (helices.len > 0) { uint64_t L = helices.len; fwrite(&L, 1, 8, fout); fwrite(helices.data, 1, L, fout); }
  if (sheets.len > 0)  { uint64_t L = sheets.len; fwrite(&L, 1, 8, fout); fwrite(sheets.data, 1, L, fout); }
  if (ssbonds.len > 0) { uint64_t L = ssbonds.len; fwrite(&L, 1, 8, fout); fwrite(ssbonds.data, 1, L, fout); }
  if (cispeps.len > 0) { uint64_t L = cispeps.len; fwrite(&L, 1, 8, fout); fwrite(cispeps.data, 1, L, fout); }
  if (sites.len > 0) {
    uint64_t L = sites.len; fwrite(&L, 1, 8, fout); fwrite(sites.data, 1, L, fout);
    size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout);
  }

  // 1. Write Coordinates
  fwrite(coords.data, 1, coords.len, fout);

  // 2. Write EOF Marker
  uint8_t eof_marker[32] = {0}; eof_marker[0] = 0xFF; eof_marker[1] = 0xFF;
  fwrite(eof_marker, 1, 32, fout);

  // 3. Write Flat CONECT Stream
  if (conects.len > 0) fwrite(conects.data, 1, conects.len, fout);

  // 4. Pad EOF to 8 bytes
  size_t current_size = ftell(fout);
  size_t eof_pad = (8 - (current_size % 8)) % 8;
  if (eof_pad > 0) fwrite(&zero, 1, eof_pad, fout);

  fclose(fout); stream_free(&coords); stream_free(&conects);
  stream_free(&ssbonds); stream_free(&cispeps); stream_free(&helices); stream_free(&sheets); stream_free(&sites);
  free(UDICT);
}
