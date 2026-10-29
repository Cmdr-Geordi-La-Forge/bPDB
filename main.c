#include "bpdb.h"
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

static inline int req_bits(int16_t v) {
    for (int b = 1; b <= 16; b++) {
        int max_val = (1 << (b - 1)) - 1;
        int min_val = -(1 << (b - 1));
        if (v >= min_val && v <= max_val) return b;
    }
    return 16;
}

static inline int req_bits_u(uint16_t v) {
    for (int b = 1; b <= 16; b++) {
        if (v < (1U << b)) return b;
    }
    return 16;
}

void stat_bpdb(const char* filepath) {
    FILE* fb = fopen(filepath, "rb");
    if (!fb) { perror("Failed to open bpdb"); return; }
    fseek(fb, 0, SEEK_END); size_t size = ftell(fb); fseek(fb, 0, SEEK_SET);
    uint8_t* buffer = (uint8_t*)malloc(size); fread(buffer, 1, size, fb); fclose(fb);

    printf("\n=================================================\n");
    printf(" bPDB Format Verbatim Statistics\n");
    printf("=================================================\n");
    printf("Total File Size      : %zu bytes\n", size);
    if (size < 64 + 8) { free(buffer); return; }
    printf("Fixed Global Header  : 64 bytes\n");

    size_t offset = 64;
    uint64_t meta_header = *(uint64_t*)(buffer + offset); offset += 8;
    uint64_t mask = meta_header >> 16;
    
    int num_meta = 0;
    for (int i = 0; i < 48; i++) {
        if (mask & (1ULL << i)) num_meta++;
    }

    static uint8_t base_natoms[256] = {0};
    static int dict_init = 0;
    if (!dict_init) {
        for (int i = 0; i < NUM_DICT_RES; i++) base_natoms[DICT[i].type_id] = DICT[i].num_atoms;
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
            int u = 0x81 + c; ptr += 4;
            r_natoms[u] = p_het_dict[ptr++]; ptr += r_natoms[u] * 6;
        }
        dict_scan_off += dict_len + ((8 - (dict_len % 8)) % 8);
    }
    if (mask & (1ULL << 6)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 7)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 8)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 9)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 10)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
    if (mask & (1ULL << 11)) { dict_scan_off += 8; }

    printf("Metadata Block (TLV) : %zu bytes (Sections: %d)\n", dict_scan_off - offset + 8, num_meta);
    offset = dict_scan_off;

    int chains = 0; size_t seq_bytes_total = 0; size_t atom_bytes_total = 0;
    int empty_atoms = 0; int present_atoms = 0; int jumps = 0; int standard_blocks = 0;
    int present_ATOM = 0; int present_HETATM = 0;
    int models = 1; int total_residues = 0;
    uint64_t dx_bits = 0, dy_bits = 0, dz_bits = 0, b_bits = 0;
    int max_dx = 0, max_dy = 0, max_dz = 0, max_b = 0;

    while (offset < size) {
        uint16_t num_residues = *(uint16_t*)(buffer + offset);
        if (num_residues == 0xFFFE) {
            models++;
            offset += 2;
            continue;
        }
        if (num_residues == 0xFFFF) {
            offset += 8; 
            break;
        }
        chains++;
        uint16_t num_atoms = *(uint16_t*)(buffer + offset + 2);
        offset += 8;

        int seq_bytes = 0; uint32_t expected_atoms = 0;

        const uint8_t* seq_ptr = buffer + offset;
        for (uint32_t r = 0; r < num_residues; r++) {
            uint8_t type = seq_ptr[seq_bytes++];
            if (type == 0x7C) {
                seq_bytes++;
            } else if (type == 0x7E || type == 0x7D) {
                expected_atoms++;
            } else if (type != 0x00 && type != 0x7F) {
                expected_atoms += r_natoms[type & 0xFF];
                total_residues++;
            }
        }
        size_t padded_seq = seq_bytes + ((8 - (seq_bytes & 7)) & 7);
        seq_bytes_total += padded_seq;
        offset += padded_seq;

        const uint8_t* atom_ptr = buffer + offset;
        
        int atoms_left = 0; uint32_t parsed_base_atoms = 0;
        int s_idx = 0; uint8_t last_valid_type = 0;
        
        for (uint32_t a = 0; a < num_atoms; a++) {
            if (atoms_left == 0) {
                if (parsed_base_atoms < expected_atoms) {
                    uint8_t meta; int current_has_5p = 0, current_has_3p = 0;
                    do {
                        if (s_idx >= seq_bytes) break;
                        meta = seq_ptr[s_idx++];
                        if (meta == 0x7E) current_has_5p = 1;
                        else if (meta == 0x7D) current_has_3p = 1;
                        else if (meta == 0x7C) s_idx++;
                    } while (meta == 0x00 || meta == 0x7F || meta == 0x7E || meta == 0x7D || meta == 0x7C);
                    
                    if (s_idx <= seq_bytes) {
                        last_valid_type = meta & 0xFF;
                        atoms_left = r_natoms[last_valid_type] + current_has_5p + current_has_3p;
                    } else {
                        atoms_left = 1; 
                    }
                } else {
                    atoms_left = 1;
                }
            }

            int is_hetatm = (last_valid_type >= 0x80);
            int alt_flag = 0;

            if (*(uint64_t*)(atom_ptr) != 0) {
                present_atoms++;
                if (is_hetatm) present_HETATM++;
                else present_ATOM++;

                if ((atom_ptr)[0] & 0x01) { 
                    jumps++;
                    uint16_t b_ctrl = *(uint16_t*)(atom_ptr + 8);
                    alt_flag = (b_ctrl >> 1) & 0x01;
                    atom_ptr += 16;
                    atom_bytes_total += 16;
                } else { 
                    standard_blocks++;
                    uint16_t b_ctrl = *(uint16_t*)atom_ptr;
                    alt_flag = (b_ctrl >> 1) & 0x01;
                    
                    int bb = req_bits_u(b_ctrl >> 2);
                    int bx = req_bits(*(int16_t*)(atom_ptr + 2));
                    int by = req_bits(*(int16_t*)(atom_ptr + 4));
                    int bz = req_bits(*(int16_t*)(atom_ptr + 6));
                    
                    b_bits += bb; dx_bits += bx; dy_bits += by; dz_bits += bz;
                    if (bb > max_b) max_b = bb;
                    if (bx > max_dx) max_dx = bx;
                    if (by > max_dy) max_dy = by;
                    if (bz > max_dz) max_dz = bz;
                    
                    atom_ptr += 8;
                    atom_bytes_total += 8;
                }
            } else {
                empty_atoms++;
                atom_ptr += 8;
                atom_bytes_total += 8;
            }

            if (alt_flag == 0) {
                atoms_left--;
                parsed_base_atoms++;
            }
        }
        offset = (size_t)(atom_ptr - buffer);
    }
    
    size_t conect_start = offset;
    while (offset + 4 <= size) {
        int32_t N = *(int32_t*)(buffer + offset);
        if (N <= 0 || N > 4 || offset + 8 + (N << 2) > size) break;
        offset += 8 + (N << 2);
    }
    
    printf("-------------------------------------------------\n");
    printf("Coordinate Stream    : %d Models, %d Chains, %d Residues\n", models, chains, total_residues);
    printf(" -> Chain Headers    : %d bytes\n", chains * 8);
    printf(" -> Sequence Tracks  : %zu bytes\n", seq_bytes_total);
    printf(" -> Atom Tracks      : %zu bytes\n", atom_bytes_total);
    printf("    * Present Atoms  : %d\n", present_atoms);
    printf("      - ATOM         : %d\n", present_ATOM);
    printf("      - HETATM       : %d\n", present_HETATM);
    printf("    * Block Types\n");
    printf("      - Standard     : %d (8B)\n", standard_blocks);
    if (standard_blocks > 0) {
        printf("        > Avg Bits   : B: %.1f | X: %.1f | Y: %.1f | Z: %.1f\n", 
               (double)b_bits / standard_blocks,
               (double)dx_bits / standard_blocks, 
               (double)dy_bits / standard_blocks, 
               (double)dz_bits / standard_blocks);
        printf("        > Max Bits   : B: %d   | X: %d   | Y: %d   | Z: %d\n", 
               max_b, max_dx, max_dy, max_dz);
    }
    printf("      - Jump Anchors : %d (16B)\n", jumps);
    printf("    * Empty Slots    : %d (%d bytes)\n", empty_atoms, empty_atoms * 8);
    printf("CONECT Stream        : %zu bytes\n", offset - conect_start);
    printf("=================================================\n\n");
    
    free(buffer);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printf("Usage Encode: %s encode <in1.pdb> [in2.pdb...] <out.bpdb>\n", argv[0]);
        printf("Usage Test  : %s test <in.bpdb> <out.pdb>\n", argv[0]);
        printf("Usage Stat  : %s stat <in.bpdb>\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "encode") == 0) {
        const char* out_bpdb = argv[argc - 1];
        int num_files = argc - 3;
        const char** in_pdbs = (const char**)&argv[2];
        encode_bpdb(num_files, in_pdbs, out_bpdb);
        if (strcmp(out_bpdb, "-") != 0) {
            printf("[+] Encoded %d files to %s\n", num_files, out_bpdb);
        }
    } else if (strcmp(argv[1], "test") == 0) {
        FILE* fb = NULL;
        size_t size = 0;
        uint8_t* buffer = NULL;

        if (strcmp(argv[2], "-") == 0) {
            fb = stdin;
            size_t cap = 1024 * 1024;
            buffer = (uint8_t*)malloc(cap);
            size_t bytes;
            while ((bytes = fread(buffer + size, 1, cap - size, fb)) > 0) {
                size += bytes;
                if (size == cap) { cap *= 2; buffer = (uint8_t*)realloc(buffer, cap); }
            }
        } else {
            fb = fopen(argv[2], "rb");
            if (!fb) { perror("Failed to open bpdb"); return 1; }
            fseek(fb, 0, SEEK_END); size = ftell(fb); fseek(fb, 0, SEEK_SET);
            buffer = (uint8_t*)malloc(size); fread(buffer, 1, size, fb); fclose(fb);
        }

        FILE* fout = (strcmp(argv[3], "-") == 0) ? stdout : fopen(argv[3], "w");
        if (!fout) { perror("Failed to open output"); free(buffer); return 1; }

        decode_bpdb(buffer, size, fout); 

        if (fout != stdout) {
            fclose(fout);
            printf("[+] Decoded to %s for diffing.\n", argv[3]);
            struct timespec start, end;
            clock_gettime(CLOCK_MONOTONIC, &start);
            int n = 100;
            for (int i = 0; i < n*1000; i++) decode_bpdb(buffer, size, NULL); 
            clock_gettime(CLOCK_MONOTONIC, &end);
            
            double t_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
            printf("[+] Parsing took %.2f s total (%i,000 files, %.3f us/run)\n", t_taken/1000, n, t_taken/n);
        }
        free(buffer);
    } else if (strcmp(argv[1], "stat") == 0) {
        stat_bpdb(argv[2]);
    }
    return 0;
}
