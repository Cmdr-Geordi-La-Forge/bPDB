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

  uint16_t total_models_minus_1 = *(uint16_t*)(buffer + 20);
  uint16_t current_model_minus_1 = *(uint16_t*)(buffer + 22);

  size_t offset = 64;
  uint64_t meta_header = *(uint64_t*)(buffer + offset); offset += 8;
  uint64_t mask = meta_header >> 16;

  int num_meta = 0;
  for (int i = 0; i < 48; i++) {
    if (mask & (1ULL << i)) num_meta++;
  }

  static ResidueDict UDICT[256];
  static int dict_init = 0;
  if (!dict_init) {
    memset(UDICT, 0, sizeof(UDICT));
    for (size_t i = 0; i < NUM_DICT_RES; i++) {
      UDICT[DICT[i].type_id] = DICT[i];
    }
    dict_init = 1;
  }

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
      UDICT[u].h_end_norm = p_het_dict[ptr++];
      UDICT[u].h_end_cterm = UDICT[u].h_end_norm;
      ptr += UDICT[u].h_end_norm * 6;
    }
    dict_scan_off += dict_len + ((8 - (dict_len % 8)) % 8);
  }
  if (mask & (1ULL << 6)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << 7)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << 8)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << 9)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }
  if (mask & (1ULL << 10)) { uint64_t L = *(uint64_t*)(buffer + dict_scan_off); dict_scan_off += 8 + L + ((8 - (L % 8)) % 8); }

  printf("Metadata Block (TLV) : %zu bytes (Sections: %d)\n", dict_scan_off - offset + 8, num_meta);
  offset = dict_scan_off;

  int chains = 0; size_t seq_bytes_total = 0; size_t atom_bytes_total = 0;
  int empty_atoms = 0; int present_atoms = 0; int jumps = 0; int standard_blocks = 0;
  int present_ATOM = 0; int present_HETATM = 0;
  int models = total_models_minus_1 + 1; int total_residues = 0;
  uint64_t dx_bits = 0, dy_bits = 0, dz_bits = 0, b_bits = 0;
  int max_dx = 0, max_dy = 0, max_dz = 0, max_b = 0;

  size_t conect_start = size;

  while (offset < size) {
    uint16_t num_residues = *(uint16_t*)(buffer + offset);
    if (num_residues == 0xFFFE) {
      offset += 8;
      continue;
    }
    if (num_residues == 0xFFFF) {
      offset += 8;
      conect_start = offset;
      break;
    }
    chains++;
    
    uint8_t c1_raw = buffer[offset + 4];
    uint8_t c2_raw = buffer[offset + 5];
    int is_n_term = (c1_raw & 0x80) != 0;
    int is_c_term = (c2_raw & 0x80) != 0;
    int has_hydrogens = (mask & (1ULL << 11)) ? 1 : 0;
    
    uint16_t num_atoms = *(uint16_t*)(buffer + offset + 2);
    offset += 8;

    int seq_bytes = 0; uint32_t expected_atoms = 0;
    const uint8_t* seq_ptr = buffer + offset;
    int first_seen = 0;
    
    for (uint32_t r = 0; r < num_residues; r++) {
      uint8_t type = seq_ptr[seq_bytes++];
      if (type == 0x7C) {
        seq_bytes++;
      } else if (type != 0x00 && type != 0x7F) {
        int is_first = (!first_seen) && is_n_term;
        int is_last  = (r == (uint32_t)(num_residues - 1)) && is_c_term;
        first_seen = 1;
        
        expected_atoms += (is_last ? UDICT[type].h_end_cterm : UDICT[type].h_end_norm) - (is_first ? 0 : UDICT[type].h_start_norm);
        if (has_hydrogens) {
          expected_atoms += (is_last ? UDICT[type].hy_end_cterm : UDICT[type].hy_end_norm) - (is_first ? 0 : UDICT[type].hy_start_norm);
        }
        total_residues++;
      }
    }
    size_t padded_seq = seq_bytes + ((8 - (seq_bytes & 7)) & 7);
    seq_bytes_total += padded_seq;
    offset += padded_seq;

    const uint8_t* atom_ptr = buffer + offset;
    int atoms_left = 0; uint32_t parsed_base_atoms = 0;
    int s_idx = 0; uint8_t last_valid_type = 0;
    int first_decoded = 0;

    for (uint32_t a = 0; a < num_atoms; a++) {
      if (atoms_left == 0) {
        if (parsed_base_atoms < expected_atoms) {
          uint8_t meta;
          do {
            if (s_idx >= seq_bytes) break;
            meta = seq_ptr[s_idx++];
            if (meta == 0x7C) { s_idx++; }
          } while (meta == 0x00 || meta == 0x7F || meta == 0x7C);

          if (s_idx <= seq_bytes) {
            last_valid_type = meta & 0xFF;
            
            int is_first = (!first_decoded) && is_n_term;
            int is_last  = (s_idx == seq_bytes) && is_c_term;
            first_decoded = 1;
            
            int h_start = is_first ? 0 : UDICT[last_valid_type].h_start_norm;
            int h_end   = is_last  ? UDICT[last_valid_type].h_end_cterm : UDICT[last_valid_type].h_end_norm;
            
            int hy_start = 0, hy_end = 0;
            if (has_hydrogens) {
                hy_start = is_first ? 0 : UDICT[last_valid_type].hy_start_norm;
                hy_end   = is_last  ? UDICT[last_valid_type].hy_end_cterm : UDICT[last_valid_type].hy_end_norm;
            }
            
            atoms_left = (h_end - h_start) + (hy_end - hy_start);
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

  printf("-------------------------------------------------\n");
  printf("Model Tracker        : Model %d of %d\n", current_model_minus_1 + 1, models);
  printf("Coordinate Stream    : %d Chains, %d Residues\n", chains, total_residues);
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

  if (conect_start < size) {
    printf("CONECT Stream        : %zu bytes\n", size - conect_start);
  }
  printf("=================================================\n\n");

  free(buffer);
}

int main(int argc, char** argv) {
  if (argc < 3) {
    printf("Usage Encode: %s encode [-h|--hydrogens] <in1.pdb> [in2.pdb...] <out.bpdb>\n", argv[0]);
    printf("Usage Test  : %s test <in1.bpdb> [in2.bpdb...] <out.pdb>\n", argv[0]);
    printf("Usage Stat  : %s stat <in.bpdb>\n", argv[0]);
    return 1;
  }

  if (strcmp(argv[1], "encode") == 0) {
    int include_hydrogens = 0;
    int file_start = 2;
    if (strcmp(argv[2], "-h") == 0 || strcmp(argv[2], "--hydrogens") == 0) {
      include_hydrogens = 1;
      file_start = 3;
    }
    const char* out_bpdb = argv[argc - 1];
    int num_files = argc - file_start - 1;
    const char** in_pdbs = (const char**)&argv[file_start];
    encode_bpdb(num_files, in_pdbs, out_bpdb, include_hydrogens);
    if (strcmp(out_bpdb, "-") != 0) {
      printf("[+] Encoded %d files to %s\n", num_files, out_bpdb);
    }
  } else if (strcmp(argv[1], "test") == 0) {
    const char* out_pdb = argv[argc - 1];
    int num_files = argc - 3;
    const char** in_bpdbs = (const char**)&argv[2];

    uint8_t** buffers = malloc(num_files * sizeof(uint8_t*));
    size_t* sizes = malloc(num_files * sizeof(size_t));

    for (int i = 0; i < num_files; i++) {
      FILE* fb = NULL;
      sizes[i] = 0;
      buffers[i] = NULL;

      if (strcmp(in_bpdbs[i], "-") == 0) {
        fb = stdin;
        size_t cap = 1024 * 1024;
        buffers[i] = (uint8_t*)malloc(cap);
        size_t bytes;
        while ((bytes = fread(buffers[i] + sizes[i], 1, cap - sizes[i], fb)) > 0) {
          sizes[i] += bytes;
          if (sizes[i] == cap) { cap *= 2; buffers[i] = (uint8_t*)realloc(buffers[i], cap); }
        }
      } else {
        fb = fopen(in_bpdbs[i], "rb");
        if (!fb) { perror("Failed to open bpdb"); continue; }
        fseek(fb, 0, SEEK_END); sizes[i] = ftell(fb); fseek(fb, 0, SEEK_SET);
        buffers[i] = (uint8_t*)malloc(sizes[i]); fread(buffers[i], 1, sizes[i], fb); fclose(fb);
      }
    }

    FILE* fout = (strcmp(out_pdb, "-") == 0) ? stdout : fopen(out_pdb, "w");
    if (!fout) { perror("Failed to open output"); return 1; }

    for (int i = 0; i < num_files; i++) {
      if (sizes[i] > 0) {
        decode_bpdb(buffers[i], sizes[i], fout, (i == 0), (i == num_files - 1));
      }
    }

    if (fout != stdout) {
      fclose(fout);
      printf("[+] Decoded %d files to %s for diffing.\n", num_files, out_pdb);
      struct timespec start, end;
      clock_gettime(CLOCK_MONOTONIC, &start);
      int n = 100;
      for (int j = 0; j < n * 1000; j++) {
        for (int i = 0; i < num_files; i++) {
          if (sizes[i] > 0) decode_bpdb(buffers[i], sizes[i], NULL, (i == 0), (i == num_files - 1));
        }
      }
      clock_gettime(CLOCK_MONOTONIC, &end);
      double t_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
      printf("[+] Parsing took %.2f s total (%i,000 files, %.3f us/run)\n", t_taken/1000, n, t_taken/n);
    }

    for (int i = 0; i < num_files; i++) free(buffers[i]);
    free(buffers); free(sizes);
  } else if (strcmp(argv[1], "stat") == 0) {
    stat_bpdb(argv[2]);
  }
  return 0;
}
