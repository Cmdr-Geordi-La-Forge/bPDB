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

static inline uint32_t pack_ec(int c, int sc, int ssc, int ser) {
    return ((c & 0x07) << 29) | ((sc & 0x1F) << 24) | 
           (((uint32_t)(ssc & 0xFF)) << 16) | ((uint32_t)(ser & 0xFFFF));
}

uint32_t parse_ec_string(const char* str) {
    int p[4] = {0, 0, 0, 0}, idx = 0;
    while (*str && idx < 4) {
        if (*str == '-') { p[idx++] = 0; str++; } 
        else if (*str == 'n') { p[idx++] = -atoi(str + 1); while (*str && *str != '.') str++; } 
        else if (*str >= '0' && *str <= '9') { p[idx++] = atoi(str); while (*str >= '0' && *str <= '9') str++; } 
        else str++;
    }
    return pack_ec(p[0], p[1], p[2], p[3]);
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

// --- PDB to bPDB ---
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
    
    bpdb_residue_t res = {0};
    bpdb_atom_t* atoms = malloc(sizeof(bpdb_atom_t) * 65535); // Heap allocation prevents stack overflow
    if (!atoms) { printf("Memory allocation failed\n"); exit(1); }

    while (fgets(line, sizeof(line), in)) {
        if (strncmp(line, "HEADER", 6) == 0) {
            memset(hdr.classification, ' ', 40);
            memset(hdr.idCode, ' ', 12);
            size_t clen = MIN(40, strlen(line) > 10 ? strlen(line+10) : 0);
            if (clen > 0 && line[10+clen-1] == '\n') clen--;
            memcpy(hdr.classification, line + 10, clen);
            
            if (strlen(line) >= 59) hdr.depDate = parse_pdb_date(line + 50);
            
            size_t id_len = MIN(12, strlen(line) > 62 ? strlen(line+62) : 0);
            if (id_len > 0 && line[62+id_len-1] == '\n') id_len--;
            memcpy(hdr.idCode, line + 62, id_len);
            
            has_hdr = 1;
        } 
        else if (strncmp(line, "TITLE ", 6) == 0) {
            size_t tlen = strlen(line);
            if (tlen > 10) {
                size_t copy_len = tlen - 10;
                if (line[10 + copy_len - 1] == '\n') copy_len--;
                
                // Trims the extra continuation space to prevent double-spacing
                int start_idx = 10;
                if (line[8] >= '0' && line[8] <= '9' && line[10] == ' ') {
                    start_idx++; copy_len--;
                }
                
                memcpy(title_buf + title_len, line + start_idx, copy_len);
                title_buf[title_len + copy_len] = '\n'; 
                title_len += copy_len + 1;
            }
        }
        else if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0 || strncmp(line, "TER", 3) == 0) {
            
            if (!metadata_flushed) {
                if (has_hdr) {
                    dir_table[DIR_HEADER] = current_offset;
                    current_offset += write_padded(&hdr, sizeof(bpdb_header_t), out);
                }
                if (title_len > 0) {
                    dir_table[DIR_TITLE] = current_offset;
                    title_buf[title_len] = '\0';
                    current_offset += write_padded(title_buf, title_len + 1, out);
                }
                dir_table[DIR_COORD] = current_offset;
                metadata_flushed = 1;
            }

            if (strncmp(line, "TER", 3) == 0) {
                if (res.atom_count > 0) {
                    current_offset += write_padded(&res, sizeof(bpdb_residue_t), out);
                    current_offset += write_padded(atoms, sizeof(bpdb_atom_t) * res.atom_count, out);
                }
                res.file_flags = 2; res.atom_count = 1;
                bpdb_atom_t ter_atom = {0}; 
                memcpy(ter_atom.name, "TER ", 4);
                current_offset += write_padded(&res, sizeof(bpdb_residue_t), out);
                current_offset += write_padded(&ter_atom, sizeof(bpdb_atom_t), out);
                res.atom_count = 0;
            } else {
                size_t len = strlen(line);
                uint16_t r_num = atoi(line + 22);
                char chain = line[21];
                double x = atof(line + 30);
                double y = atof(line + 38);
                double z = atof(line + 46);
                int16_t bx = calc_base(x), by = calc_base(y), bz = calc_base(z);

                if (res.atom_count > 0 && (res.res_num != r_num || res.chain_id != chain ||
                    res.base_x != bx || res.base_y != by || res.base_z != bz)) {
                    current_offset += write_padded(&res, sizeof(bpdb_residue_t), out);
                    current_offset += write_padded(atoms, sizeof(bpdb_atom_t) * res.atom_count, out);
                    res.atom_count = 0;
                }

                if (res.atom_count == 0) {
                    memcpy(res.name, line + 17, 3);
                    res.chain_id = chain;
                    res.res_num = r_num;
                    res.icode = line[26];
                    res.file_flags = (line[0] == 'A') ? 0 : 1;
                    res.base_x = bx; res.base_y = by; res.base_z = bz;
                }

                bpdb_atom_t* a = &atoms[res.atom_count++];
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
    }
    if (res.atom_count > 0) {
        write_padded(&res, sizeof(bpdb_residue_t), out);
        write_padded(atoms, sizeof(bpdb_atom_t) * res.atom_count, out);
    }

    fseek(out, 0, SEEK_SET);
    fwrite(dir_table, sizeof(uint32_t), 64, out);

    free(atoms);
    fclose(in); fclose(out);
}

// --- bPDB to PDB (Back-Converter) ---
void decode_bpdb(const char* in_path, const char* out_path) {
    int fd_in = open(in_path, O_RDONLY);
    if (fd_in < 0) return;
    struct stat sb;
    fstat(fd_in, &sb);
    if (sb.st_size < 256) { close(fd_in); return; }

    const uint8_t* file_data = (const uint8_t*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd_in, 0);
    if (file_data == MAP_FAILED) { close(fd_in); return; }
    const uint32_t* dir_table = (const uint32_t*)file_data;

    // 1. Fast pre-pass to calculate exact line count
    size_t total_lines = 0;
    if (dir_table[DIR_HEADER] != 0) total_lines++;
    if (dir_table[DIR_TITLE] != 0) {
        const char* title = (const char*)(file_data + dir_table[DIR_TITLE]);
        while (*title) { if (*title++ == '\n') total_lines++; }
    }
    if (dir_table[DIR_COORD] != 0) {
        size_t offset = dir_table[DIR_COORD];
        while (offset < (size_t)sb.st_size) {
            const bpdb_residue_t* res = (const bpdb_residue_t*)(file_data + offset);
            offset += sizeof(bpdb_residue_t) + (res->atom_count * sizeof(bpdb_atom_t));
            total_lines += res->atom_count;
        }
    }

    if (total_lines == 0) { munmap((void*)file_data, sb.st_size); close(fd_in); return; }

    // 2. RAM Buffer Allocation (Immune to mapping errors)
    size_t out_size = total_lines * 81;
    char* out_data = malloc(out_size);
    if (!out_data) { munmap((void*)file_data, sb.st_size); close(fd_in); return; }

    // 3. Pre-fill entirely with PDB 81-byte padding pattern
    const char pad[81] = "                                                                                \n";
    for (size_t i = 0; i < total_lines; i++) {
        memcpy(out_data + (i * 81), pad, 81);
    }

    // 4. Format directly to mapped output
    size_t out_idx = 0;
    
    if (dir_table[DIR_HEADER] != 0) {
        const bpdb_header_t* hdr = (const bpdb_header_t*)(file_data + dir_table[DIR_HEADER]);
        int len = snprintf(out_data + (out_idx * 81), 81, "HEADER    %-40.40s%02d-%3s-%02d   %-12.12s", 
                hdr->classification, hdr->depDate.day, 
                M_STR[hdr->depDate.month <= 12 ? hdr->depDate.month : 0], 
                hdr->depDate.year, hdr->idCode);
        if (len > 80) len = 80;
        out_data[(out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
        out_idx++;
    }

    if (dir_table[DIR_TITLE] != 0) {
        const char* title = (const char*)(file_data + dir_table[DIR_TITLE]);
        int line_num = 1;
        const char* ptr = title;
        while (*ptr) {
            const char* next_nl = strchr(ptr, '\n');
            if (!next_nl) break;
            int str_len = next_nl - ptr;
            int len = (line_num == 1) 
                ? snprintf(out_data + (out_idx * 81), 81, "TITLE     %.*s", str_len, ptr)
                : snprintf(out_data + (out_idx * 81), 81, "TITLE   %2d%.*s", line_num, str_len, ptr);
            
            if (len > 80) len = 80;
            out_data[(out_idx * 81) + len] = (len == 80) ? '\n' : ' '; 
            ptr = next_nl + 1;
            line_num++; out_idx++;
        }
    }

    if (dir_table[DIR_COORD] != 0) {
        size_t offset = dir_table[DIR_COORD];
        int serial = 1;
        while (offset < (size_t)sb.st_size) {
            const bpdb_residue_t* res = (const bpdb_residue_t*)(file_data + offset);
            offset += sizeof(bpdb_residue_t);
            uint8_t type = res->file_flags & 0x03;

            for (int i = 0; i < res->atom_count; i++) {
                const bpdb_atom_t* atom = (const bpdb_atom_t*)(file_data + offset);
                offset += sizeof(bpdb_atom_t);

                if (type == 2) {
                    int len = snprintf(out_data + (out_idx * 81), 81, "TER   %5d      %.3s %c%4d%c", 
                            serial++, res->name, res->chain_id, res->res_num, res->icode);
                    if (len > 80) len = 80;
                    out_data[(out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
                    out_idx++;
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

                int len = snprintf(out_data + (out_idx * 81), 81, 
                        "%s%5d %.4s%c%.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                        (type == 0) ? "ATOM  " : "HETATM", serial++, atom->name, altloc, 
                        res->name, res->chain_id, res->res_num, res->icode,
                        ix * 0.001, iy * 0.001, iz * 0.001,
                        atom->occupancy / 1000.0f, atom->bfactor / 100.0f,
                        Z_TO_ELEM[atom->z_num], charge_str);
                
                if (len > 80) len = 80;
                out_data[(out_idx * 81) + len] = (len == 80) ? '\n' : ' ';
                out_idx++;
            }
        }
    }

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
