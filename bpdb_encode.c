#include "bpdb.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

void stream_init(StreamBuffer* buf) {
    buf->cap = 1024 * 1024 * 16; buf->data = (uint8_t*)malloc(buf->cap); buf->len = 0;
}
void stream_write(StreamBuffer* buf, const void* src, size_t size) {
    if (buf->len + size > buf->cap) { buf->cap *= 2; buf->data = (uint8_t*)realloc(buf->data, buf->cap); }
    memcpy(buf->data + buf->len, src, size); buf->len += size;
}
void stream_free(StreamBuffer* buf) { free(buf->data); }

int parse_month(const char* m) { for (int i = 1; i <= 12; i++) if (strncmp(m, MONTHS[i], 3) == 0) return i; return 1; }
void trim_string(char* str) { int i = strlen(str) - 1; while(i >= 0 && isspace(str[i])) { str[i] = '\0'; i--; } }

typedef struct {
    int resSeq; char iCode; int type; int is_dup;
    int32_t b_acc[256], x_acc[256], y_acc[256], z_acc[256];
    int32_t occ_acc[256]; char alt_acc[256]; int8_t charge_acc[256]; int present[256];
    int has_oxt; int32_t oxt_b, oxt_x, oxt_y, oxt_z; int32_t oxt_occ; char oxt_alt; int8_t oxt_charge;
} EncoderState;

typedef struct {
    uint8_t* seq_track; uint8_t* atom_track;
    int seq_bytes; int atom_bytes; int num_atoms; int num_residues;
    int32_t current_x, current_y, current_z; 
    int32_t current_occ; uint8_t current_charge;
    int start_resSeq; int last_resSeq; char last_iCode;
    char c1, c2; int has_ter;
} ChainBuilder;

void cb_init(ChainBuilder* cb) {
    cb->seq_track = (uint8_t*)malloc(240000); cb->atom_track = (uint8_t*)malloc(120000 * 8 + 100);
    cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0;
    cb->start_resSeq = 0; cb->has_ter = 0;
    cb->current_occ = 100; cb->current_charge = 8; // Offset by +8 
}
void cb_free(ChainBuilder* cb) { free(cb->seq_track); free(cb->atom_track); }

void flush_chain(ChainBuilder* cb, StreamBuffer* buf) {
    if (cb->num_atoms == 0) return;
    
    uint8_t header[8] = {0};
    *(uint16_t*)(header + 0) = cb->num_residues;
    *(uint16_t*)(header + 2) = cb->num_atoms;
    header[4] = cb->c1 | (cb->has_ter ? 0x80 : 0); 
    header[5] = cb->c2;
    *(uint16_t*)(header + 6) = cb->start_resSeq;
    stream_write(buf, header, 8);

    stream_write(buf, cb->seq_track, cb->seq_bytes);
    uint32_t pad = (8 - (cb->seq_bytes & 7)) & 7;
    if (pad > 0) { uint64_t zero = 0; stream_write(buf, &zero, pad); }
    
    stream_write(buf, cb->atom_track, cb->atom_bytes);

    cb->seq_bytes = 0; cb->atom_bytes = 0; cb->num_atoms = 0; cb->num_residues = 0; 
    cb->has_ter = 0; cb->current_occ = 100; cb->current_charge = 8;
}

void flush_residue(ChainBuilder* cb, EncoderState* st, UniDict* UDICT) {
    if (st->type <= 0) return;

    if (cb->num_residues == 0) cb->start_resSeq = st->resSeq;

    if (cb->num_residues > 0 && !st->is_dup) {
        int diff = st->resSeq - cb->last_resSeq;
        if (diff > 1 && diff <= 100) {
            for (int g = 1; g < diff; g++) { cb->seq_track[cb->seq_bytes++] = 0x00; cb->num_residues++; }
        }
    }

    if (st->is_dup) { cb->seq_track[cb->seq_bytes++] = 0x7F; cb->num_residues++; }

    uint8_t meta = st->type & 0x7F;
    if (st->iCode != ' ' && st->iCode != 0) meta |= 0x80;
    cb->seq_track[cb->seq_bytes++] = meta;
    if (meta & 0x80) cb->seq_track[cb->seq_bytes++] = st->iCode;
    cb->num_residues++;

    for (int i = 0; i < UDICT[st->type].num_atoms; i++) {
        if (!st->present[i]) {
            memset(&cb->atom_track[cb->atom_bytes], 0, 8);
            cb->atom_bytes += 8; cb->num_atoms++;
            continue;
        }

        int32_t x = st->x_acc[i], y = st->y_acc[i], z = st->z_acc[i];
        if (cb->num_atoms == 0) { cb->current_x = x; cb->current_y = y; cb->current_z = z; }

        int32_t dx = x - cb->current_x, dy = y - cb->current_y, dz = z - cb->current_z;
        int32_t b_int = st->b_acc[i]; if (b_int > 4194303) b_int = 4194303; if (b_int < 0) b_int = 0;

        char altLoc = st->alt_acc[i];
        int alt_flag = (altLoc != ' ' && altLoc != 0) ? 1 : 0;
        int32_t occ_int = st->occ_acc[i];
        uint8_t charge_u = (st->charge_acc[i] + 8) & 0x0F;

        // Treat HETATMs identical to normal atoms: Only jump if bounds or state change
        int needs_jump = (cb->num_atoms == 0) || 
                         (dx > 32767 || dx < -32768 || dy > 32767 || dy < -32768 || dz > 32767 || dz < -32768) ||
                         (occ_int != cb->current_occ) || (charge_u != cb->current_charge) || 
                         (b_int > 16383);

        uint8_t* ptr = &cb->atom_track[cb->atom_bytes];
        uint16_t b_ctrl = ((b_int & 0x3FFF) << 2) | (alt_flag << 1) | 0; 

        if (needs_jump) {
            ptr[0] = (charge_u << 1) | 1;           // +0 charge + control bit (1)
            ptr[1] = (b_int >> 14) & 0xFF;          // +1 B-anchor MSB
            ptr[2] = (occ_int >> 16) & 0xFF;        // +2 Occ byte 2 (MSB)
            ptr[3] = (x >> 16) & 0xFF;              // +3 X-anchor MSB
            ptr[4] = (occ_int >> 8) & 0xFF;         // +4 Occ byte 1
            ptr[5] = (y >> 16) & 0xFF;              // +5 Y-anchor MSB
            ptr[6] = occ_int & 0xFF;                // +6 Occ byte 0 (LSB)
            ptr[7] = (z >> 16) & 0xFF;              // +7 Z-anchor MSB

            *(uint16_t*)(ptr + 8) = b_ctrl;
            *(uint16_t*)(ptr + 10) = x & 0xFFFF;
            *(uint16_t*)(ptr + 12) = y & 0xFFFF;
            *(uint16_t*)(ptr + 14) = z & 0xFFFF;

            cb->current_occ = occ_int; cb->current_charge = charge_u;
            cb->atom_bytes += 16;
        } else {
            *(uint16_t*)ptr = b_ctrl;
            *(int16_t*)(ptr + 2) = dx;
            *(int16_t*)(ptr + 4) = dy;
            *(int16_t*)(ptr + 6) = dz;
            cb->atom_bytes += 8;
        }

        cb->current_x = x; cb->current_y = y; cb->current_z = z; cb->num_atoms++;
    }
    cb->last_resSeq = st->resSeq; cb->last_iCode = st->iCode;
}

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb) {
    UniDict* UDICT = (UniDict*)calloc(256, sizeof(UniDict));
    int num_udict = 1; 
    for (int i = 0; i < NUM_DICT_RES; i++) {
        strcpy(UDICT[num_udict].resName, DICT[i].resName); 
        UDICT[num_udict].num_atoms = DICT[i].num_atoms;
        for (int a = 0; a < DICT[i].num_atoms; a++) { 
            strcpy(UDICT[num_udict].atom_names[a], DICT[i].atom_names[a]); 
            strcpy(UDICT[num_udict].elements[a], DICT[i].elements[a]); 
        }
        num_udict++;
    }

    char global_title[4096] = {0}; char global_id[128] = {0}; char global_class[128] = {0};
    char global_expdta[128] = {0}; char global_program[56] = {0};
    float global_resolution = -1.0f;
    int day=1, month=1, year=1973, id_from_header = 0;

    int af_mode = 0;
    if (strstr(in_pdbs[0], "AF-") != NULL) af_mode = 1;

    StreamBuffer ssbonds; stream_init(&ssbonds); StreamBuffer cispeps; stream_init(&cispeps);
    StreamBuffer helices; stream_init(&helices); StreamBuffer sheets;  stream_init(&sheets);
    StreamBuffer sites;   stream_init(&sites);

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
                if (strlen(line) > 10) { char chunk[71] = {0}; size_t c_len = strlen(line) - 10; memcpy(chunk, line+10, c_len > 70 ? 70 : c_len); trim_string(chunk); if (strlen(global_title) > 0) strcat(global_title, " "); strcat(global_title, chunk); }
                if (!id_from_header) { char *p1 = strchr(line, '('), *p2 = strchr(p1 ? p1 : line, ')'); if (p1 && p2 && (p2 > p1 + 4)) { size_t len = p2 - p1 - 1; memcpy(global_id, p1+1, len > 14 ? 14 : len); global_id[len > 14 ? 14 : len] = '\0'; } }
                if (strstr(global_title, "ALPHAFOLD") || strstr(global_title, "AlphaFold")) af_mode = 1;
            } else if (f == 0 && strncmp(line, "EXPDTA", 6) == 0) { char chunk[71] = {0}; size_t c_len = strlen(line) - 10; memcpy(chunk, line+10, c_len > 70 ? 70 : c_len); trim_string(chunk); if (strlen(global_expdta) > 0) strcat(global_expdta, " "); strncat(global_expdta, chunk, 127 - strlen(global_expdta));
            } else if (f == 0 && strncmp(line, "REMARK   2 RESOLUTION.", 22) == 0) { float res = atof(line + 22); if (res > 0.0f) global_resolution = res;
            } else if (f == 0 && strncmp(line, "REMARK   3", 10) == 0) {
                char* prog = strstr(line, "PROGRAM"); char* ver = strstr(line, "VERSION"); char* colon = strchr(line, ':');
                if (prog && colon && strlen(global_program) == 0) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; memcpy(temp, start, strlen(start)); trim_string(temp); strncat(global_program, temp, 55); } 
                else if (ver && colon && strlen(global_program) > 0 && !strstr(global_program, "v")) { char temp[64] = {0}; char* start = colon + 1; while (*start == ' ' || *start == '\t') start++; memcpy(temp, start, strlen(start)); trim_string(temp); strncat(global_program, " v", 55 - strlen(global_program)); strncat(global_program, temp, 55 - strlen(global_program)); }
            } else if (strncmp(line, "SSBOND", 6) == 0) { MetaSSBond b; memset(&b, 0, sizeof(b)); memcpy(b.res1.resName, line+11, 3); b.res1.chain = line[15]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); b.res1.resSeq = atoi(tbuf); b.res1.iCode = line[21]; memcpy(b.res2.resName, line+25, 3); b.res2.chain = line[29]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); b.res2.resSeq = atoi(tbuf); b.res2.iCode = line[35]; if (strlen(line) >= 65) { char tmp[8] = {0}; memcpy(tmp, line+59, 6); if (sscanf(tmp, "%s", b.sym1) != 1) b.sym1[0] = '\0'; } if (strlen(line) >= 72) { char tmp[8] = {0}; memcpy(tmp, line+66, 6); if (sscanf(tmp, "%s", b.sym2) != 1) b.sym2[0] = '\0'; } if (strlen(line) >= 78) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+73, 5); b.length = atof(tbuf); } stream_write(&ssbonds, &b, sizeof(MetaSSBond));
            } else if (strncmp(line, "CISPEP", 6) == 0) { MetaCispep c; memset(&c, 0, sizeof(c)); memcpy(c.res1.resName, line+11, 3); c.res1.chain = line[15]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+17, 4); c.res1.resSeq = atoi(tbuf); c.res1.iCode = line[21]; memcpy(c.res2.resName, line+25, 3); c.res2.chain = line[29]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+31, 4); c.res2.resSeq = atoi(tbuf); c.res2.iCode = line[35]; if (strlen(line) >= 46) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+43, 3); c.modNum = atoi(tbuf); } if (strlen(line) >= 59) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+53, 6); c.angle = atof(tbuf); } stream_write(&cispeps, &c, sizeof(MetaCispep));
            } else if (strncmp(line, "HELIX ", 6) == 0) { MetaHelix h; memset(&h, 0, sizeof(h)); memcpy(h.helixID, line+11, 3); trim_string(h.helixID); memcpy(h.start.resName, line+15, 3); h.start.chain = line[19]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+21, 4); h.start.resSeq = atoi(tbuf); h.start.iCode = line[25]; memcpy(h.end.resName, line+27, 3); h.end.chain = line[31]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); h.end.resSeq = atoi(tbuf); h.end.iCode = line[37]; if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); h.helix_class = atoi(tbuf); } if (strlen(line) >= 70) { memcpy(h.comment, line+40, 30); trim_string(h.comment); } if (strlen(line) >= 76) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+71, 5); h.length = atoi(tbuf); } stream_write(&helices, &h, sizeof(MetaHelix));
            } else if (strncmp(line, "SHEET ", 6) == 0) { MetaSheet s; memset(&s, 0, sizeof(s)); memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+7, 3); s.strand_num = atoi(tbuf); memcpy(s.sheetID, line+11, 3); trim_string(s.sheetID); if (strlen(line) >= 16) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+14, 2); s.num_strands = atoi(tbuf); } memcpy(s.start.resName, line+17, 3); s.start.chain = line[21]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+22, 4); s.start.resSeq = atoi(tbuf); s.start.iCode = line[26]; memcpy(s.end.resName, line+28, 3); s.end.chain = line[32]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+33, 4); s.end.resSeq = atoi(tbuf); s.end.iCode = line[37]; if (strlen(line) >= 40) { memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+38, 2); s.sense = atoi(tbuf); } stream_write(&sheets, &s, sizeof(MetaSheet));
            } else if (strncmp(line, "SITE  ", 6) == 0) { MetaSiteHeader sh; memset(&sh, 0, sizeof(sh)); memcpy(sh.siteID, line+11, 3); trim_string(sh.siteID); memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+15, 2); sh.num_res = atoi(tbuf); stream_write(&sites, &sh, sizeof(MetaSiteHeader)); for (uint32_t i=0; i<sh.num_res; i++) { PdbResId r; memset(&r, 0, sizeof(r)); int off = 18 + i*11; if (off + 9 < (int)strlen(line)) { memcpy(r.resName, line+off, 3); r.chain = line[off+4]; memset(tbuf, 0, sizeof(tbuf)); memcpy(tbuf, line+off+5, 4); r.resSeq = atoi(tbuf); r.iCode = line[off+9]; stream_write(&sites, &r, sizeof(PdbResId)); } }
            } else if (strncmp(line, "DBREF ", 6) == 0) { char* unp = strstr(line, "UNP"); if (unp) { char dummy[16], uid[64]; if (sscanf(unp, "%s %s", dummy, uid) == 2 && strlen(uid) > 3 && f == 0 && !id_from_header) { size_t c = strlen(uid); memcpy(global_id, uid, c > 14 ? 14 : c); global_id[c > 14 ? 14 : c] = '\0'; } }
            } else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
                char res_name[4] = {0}; memcpy(res_name, line + 17, 3);
                char atom_name[5] = {0}; memcpy(atom_name, line + 12, 4);
                char clean_atom[5] = {0}; memcpy(clean_atom, line + 12, 4);
                if (clean_atom[1] == 'H' || clean_atom[0] == 'H' || strncmp(clean_atom, " OXT", 4) == 0) continue;
                
                int type = -1;
                for (int i = 1; i < num_udict; i++) { if (strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; } }
                if (type == -1 && num_udict < 127) { type = num_udict++; memcpy(UDICT[type].resName, res_name, 3); UDICT[type].num_atoms = 0; }
                if (type != -1) { 
                    int a_idx = -1; for (int a = 0; a < UDICT[type].num_atoms; a++) { if (strncmp(UDICT[type].atom_names[a], atom_name, 4) == 0) { a_idx = a; break; } }
                    if (a_idx == -1 && UDICT[type].num_atoms < 255) { a_idx = UDICT[type].num_atoms++; memcpy(UDICT[type].atom_names[a_idx], atom_name, 4); UDICT[type].elements[a_idx][0] = line[76]; UDICT[type].elements[a_idx][1] = line[77]; } 
                }
            }
        }
        fclose(fin);
    }

    StreamBuffer coords; stream_init(&coords); StreamBuffer conects; stream_init(&conects);
    EncoderState st; memset(&st, 0, sizeof(st)); st.type = -1;
    ChainBuilder cb; cb_init(&cb); char current_c1 = '\0', current_c2 = '\0';

    for (int f = 0; f < num_files; f++) {
        FILE* fin = fopen(in_pdbs[f], "r");
        char line[128];
        while (fgets(line, sizeof(line), fin)) {
            if (strncmp(line, "CONECT", 6) == 0) {
                int len = strlen(line); int32_t vals[5] = {0}, targets = 0;
                for(int c=0; c<5; c++) { int start = 6 + c*5; if (start < len) { char buf[6] = {0}; int copy_len = len - start; if (copy_len > 5) copy_len = 5; memcpy(buf, line + start, copy_len); vals[c] = atoi(buf); if (c > 0 && vals[c] > 0) targets++; } }
                if (vals[0] > 0 && targets > 0) { stream_write(&conects, &targets, 4); stream_write(&conects, &vals[0], 4); for(int c=1; c<=4; c++) if (vals[c] > 0) stream_write(&conects, &vals[c], 4); }
                continue;
            }

            if (strncmp(line, "TER   ", 6) == 0) {
                if (st.type > 0) flush_residue(&cb, &st, UDICT);
                if (st.has_oxt) {
                    int32_t dx = st.oxt_x - cb.current_x, dy = st.oxt_y - cb.current_y, dz = st.oxt_z - cb.current_z;
                    int32_t b_int = st.oxt_b; if (b_int > 4194303) b_int = 4194303; if (b_int < 0) b_int = 0;
                    
                    char altLoc = st.oxt_alt;
                    int alt_flag = (altLoc != ' ' && altLoc != 0) ? 1 : 0;
                    int32_t occ_int = st.oxt_occ;
                    uint8_t charge_u = (st.oxt_charge + 8) & 0x0F;

                    int needs_jump = (cb.num_atoms == 0) || 
                                     (dx > 32767 || dx < -32768 || dy > 32767 || dy < -32768 || dz > 32767 || dz < -32768) ||
                                     (occ_int != cb.current_occ) || (charge_u != cb.current_charge) || 
                                     (b_int > 16383);

                    uint8_t* ptr = &cb.atom_track[cb.atom_bytes];
                    uint16_t b_ctrl = ((b_int & 0x3FFF) << 2) | (alt_flag << 1) | 0;

                    if (needs_jump) {
                        ptr[0] = (charge_u << 1) | 1;
                        ptr[1] = (b_int >> 14) & 0xFF;
                        ptr[2] = (occ_int >> 16) & 0xFF;
                        ptr[3] = (st.oxt_x >> 16) & 0xFF;
                        ptr[4] = (occ_int >> 8) & 0xFF;
                        ptr[5] = (st.oxt_y >> 16) & 0xFF;
                        ptr[6] = occ_int & 0xFF;
                        ptr[7] = (st.oxt_z >> 16) & 0xFF;
                        
                        *(uint16_t*)(ptr + 8) = b_ctrl;
                        *(uint16_t*)(ptr + 10) = st.oxt_x & 0xFFFF;
                        *(uint16_t*)(ptr + 12) = st.oxt_y & 0xFFFF;
                        *(uint16_t*)(ptr + 14) = st.oxt_z & 0xFFFF;

                        cb.current_occ = occ_int; cb.current_charge = charge_u;
                        cb.atom_bytes += 16;
                    } else {
                        *(uint16_t*)ptr = b_ctrl;
                        *(int16_t*)(ptr + 2) = dx; *(int16_t*)(ptr + 4) = dy; *(int16_t*)(ptr + 6) = dz;
                        cb.atom_bytes += 8;
                    }
                    cb.current_x = st.oxt_x; cb.current_y = st.oxt_y; cb.current_z = st.oxt_z; cb.num_atoms++;
                }
                cb.has_ter = 1; 
                flush_chain(&cb, &coords);
                st.type = -1; st.has_oxt = 0; st.is_dup = 0;
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
            for (int i = 1; i < num_udict; i++) { if (strncmp(UDICT[i].resName, res_name, 3) == 0) { type = i; break; } }
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
                if (charge < -8) charge = -8; else if (charge > 7) charge = 7;
            }

            int atom_idx = -1;
            for (int i = 0; i < UDICT[type].num_atoms; i++) { if (strncmp(clean_atom, UDICT[type].atom_names[i], 4) == 0) { atom_idx = i; break; } }
            if (atom_idx == -1 && strncmp(clean_atom, " OXT", 4) != 0) continue;

            if (resSeq != st.resSeq || iCode != st.iCode || new_chain || cb.atom_bytes >= 580000 || type != st.type) {
                if (st.type > 0) flush_residue(&cb, &st, UDICT);
                if (new_chain || cb.atom_bytes >= 580000 || (resSeq < cb.last_resSeq) || (resSeq - cb.last_resSeq > 100)) {
                    if (cb.num_atoms > 0) flush_chain(&cb, &coords);
                    if (new_chain) { current_c1 = c1; current_c2 = c2; }
                    cb.c1 = c1; cb.c2 = c2; 
                }
                st.resSeq = resSeq; st.iCode = iCode; st.type = type;
                memset(st.present, 0, sizeof(st.present)); st.has_oxt = 0; st.is_dup = 0;
            } else if (atom_idx != -1 && st.present[atom_idx]) {
                flush_residue(&cb, &st, UDICT);
                st.is_dup = 1; memset(st.present, 0, sizeof(st.present)); st.has_oxt = 0;
            }

            if (strncmp(clean_atom, " OXT", 4) == 0) {
                st.has_oxt = 1; st.oxt_b = (int32_t)round(temp * 100.0f); st.oxt_x = x_int; st.oxt_y = y_int; st.oxt_z = z_int;
                st.oxt_occ = (int32_t)round(occ); st.oxt_alt = altLoc; st.oxt_charge = charge;
                continue;
            }

            st.b_acc[atom_idx] = (int32_t)round(temp * 100.0f); st.x_acc[atom_idx] = x_int; st.y_acc[atom_idx] = y_int; st.z_acc[atom_idx] = z_int;
            st.occ_acc[atom_idx] = (int32_t)round(occ); st.alt_acc[atom_idx] = altLoc; st.charge_acc[atom_idx] = charge; st.present[atom_idx] = 1;
        }
        fclose(fin);
    }
    
    if (st.type > 0) flush_residue(&cb, &st, UDICT);
    if (st.has_oxt) {
        int32_t dx = st.oxt_x - cb.current_x, dy = st.oxt_y - cb.current_y, dz = st.oxt_z - cb.current_z;
        int32_t b_int = st.oxt_b; if (b_int > 4194303) b_int = 4194303; if (b_int < 0) b_int = 0;
        
        char altLoc = st.oxt_alt;
        int alt_flag = (altLoc != ' ' && altLoc != 0) ? 1 : 0;
        int32_t occ_int = st.oxt_occ;
        uint8_t charge_u = (st.oxt_charge + 8) & 0x0F;

        int needs_jump = (cb.num_atoms == 0) || 
                         (dx > 32767 || dx < -32768 || dy > 32767 || dy < -32768 || dz > 32767 || dz < -32768) ||
                         (occ_int != cb.current_occ) || (charge_u != cb.current_charge) || 
                         (b_int > 16383);

        uint8_t* ptr = &cb.atom_track[cb.atom_bytes];
        uint16_t b_ctrl = ((b_int & 0x3FFF) << 2) | (alt_flag << 1) | 0;

        if (needs_jump) {
            ptr[0] = (charge_u << 1) | 1;
            ptr[1] = (b_int >> 14) & 0xFF;
            ptr[2] = (occ_int >> 16) & 0xFF;
            ptr[3] = (st.oxt_x >> 16) & 0xFF;
            ptr[4] = (occ_int >> 8) & 0xFF;
            ptr[5] = (st.oxt_y >> 16) & 0xFF;
            ptr[6] = occ_int & 0xFF;
            ptr[7] = (st.oxt_z >> 16) & 0xFF;
            
            *(uint16_t*)(ptr + 8) = b_ctrl;
            *(uint16_t*)(ptr + 10) = st.oxt_x & 0xFFFF;
            *(uint16_t*)(ptr + 12) = st.oxt_y & 0xFFFF;
            *(uint16_t*)(ptr + 14) = st.oxt_z & 0xFFFF;

            cb.current_occ = occ_int; cb.current_charge = charge_u;
            cb.atom_bytes += 16;
        } else {
            *(uint16_t*)ptr = b_ctrl;
            *(int16_t*)(ptr + 2) = dx; *(int16_t*)(ptr + 4) = dy; *(int16_t*)(ptr + 6) = dz;
            cb.atom_bytes += 8;
        }
        cb.num_atoms++;
    }
    if (st.type > 0 && st.type < 38) cb.has_ter = 1; 
    flush_chain(&cb, &coords);
    cb_free(&cb);

    FILE* fout = fopen(out_bpdb, "wb");
    uint64_t zero = 0;
    
    uint8_t header[64] = {0};
    memcpy(header, "bPDB", 4); header[4] = 5; header[5] = (uint8_t)day; header[6] = (uint8_t)month;
    header[7] = (year >= 1973 && year <= 2228) ? (uint8_t)(year - 1973) : 0;
    memcpy(header + 8, global_id, strlen(global_id) > 15 ? 15 : strlen(global_id));
    memcpy(header + 24, global_class, strlen(global_class) > 40 ? 40 : strlen(global_class));
    fwrite(header, 1, 64, fout);

    uint64_t title_len = strlen(global_title), expdta_len = strlen(global_expdta);
    int has_title = (title_len > 0), has_het = (num_udict > 1 + NUM_DICT_RES), has_expdta = (expdta_len > 0);
    int has_res = (global_resolution > 0.0f), has_prog = (strlen(global_program) > 0);
    
    uint64_t meta_mask = 0;
    if (af_mode) meta_mask |= (1ULL << 0);
    if (has_res) meta_mask |= (1ULL << 1);
    if (has_title) meta_mask |= (1ULL << 2);
    if (has_expdta) meta_mask |= (1ULL << 3);
    if (has_prog) meta_mask |= (1ULL << 4);
    if (has_het) meta_mask |= (1ULL << 5);
    if (helices.len > 0) meta_mask |= (1ULL << 6);
    if (sheets.len > 0) meta_mask |= (1ULL << 7);
    if (ssbonds.len > 0) meta_mask |= (1ULL << 8);
    if (cispeps.len > 0) meta_mask |= (1ULL << 9);
    if (sites.len > 0) meta_mask |= (1ULL << 10);

    uint64_t meta_header = 5 | (meta_mask << 16);
    fwrite(&meta_header, 1, 8, fout);

    if (meta_mask & (1ULL << 1)) { uint32_t pad=0; fwrite(&global_resolution, 1, 4, fout); fwrite(&pad, 1, 4, fout); }
    if (meta_mask & (1ULL << 2)) { fwrite(&title_len, 1, 8, fout); fwrite(global_title, 1, title_len, fout); size_t pad = (8 - (title_len % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 3)) { fwrite(&expdta_len, 1, 8, fout); fwrite(global_expdta, 1, expdta_len, fout); size_t pad = (8 - (expdta_len % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 4)) { uint64_t p_len = 56; fwrite(&p_len, 1, 8, fout); fwrite(global_program, 1, 56, fout); }
    if (meta_mask & (1ULL << 5)) {
        uint64_t dict_len = 1; 
        for (int i = NUM_DICT_RES + 1; i < num_udict; i++) dict_len += 5 + UDICT[i].num_atoms * 6;
        fwrite(&dict_len, 1, 8, fout); 
        uint8_t custom_count = num_udict - NUM_DICT_RES - 1; 
        fwrite(&custom_count, 1, 1, fout);
        for (int i = NUM_DICT_RES + 1; i < num_udict; i++) {
            char r_name[4] = {0}; memcpy(r_name, UDICT[i].resName, 3); r_name[3] = ' '; fwrite(r_name, 1, 4, fout);
            uint8_t a_count = UDICT[i].num_atoms; fwrite(&a_count, 1, 1, fout);
            for (int a = 0; a < a_count; a++) { fwrite(UDICT[i].atom_names[a], 1, 4, fout); fwrite(UDICT[i].elements[a], 1, 2, fout); }
        }
        size_t pad = (8 - (dict_len % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout);
    }
    
    if (meta_mask & (1ULL << 6)) { uint64_t L = helices.len; fwrite(&L, 1, 8, fout); fwrite(helices.data, 1, L, fout); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 7)) { uint64_t L = sheets.len;  fwrite(&L, 1, 8, fout); fwrite(sheets.data, 1, L, fout); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 8)) { uint64_t L = ssbonds.len; fwrite(&L, 1, 8, fout); fwrite(ssbonds.data, 1, L, fout); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 9)) { uint64_t L = cispeps.len; fwrite(&L, 1, 8, fout); fwrite(cispeps.data, 1, L, fout); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }
    if (meta_mask & (1ULL << 10)) { uint64_t L = sites.len; fwrite(&L, 1, 8, fout); fwrite(sites.data, 1, L, fout); size_t pad = (8 - (L % 8)) % 8; if (pad > 0) fwrite(&zero, 1, pad, fout); }

    fwrite(coords.data, 1, coords.len, fout);
    
    uint8_t eof_marker[8] = {0}; 
    eof_marker[0] = 0xFF; eof_marker[1] = 0xFF; 
    fwrite(eof_marker, 1, 8, fout);
    
    if (conects.len > 0) fwrite(conects.data, 1, conects.len, fout);
    
    size_t current_size = ftell(fout);
    size_t eof_pad = (8 - (current_size % 8)) % 8; 
    if (eof_pad > 0) fwrite(&zero, 1, eof_pad, fout);

    fclose(fout); 
    stream_free(&coords); stream_free(&conects); stream_free(&ssbonds); 
    stream_free(&cispeps); stream_free(&helices); stream_free(&sheets); stream_free(&sites);
    free(UDICT);
}
