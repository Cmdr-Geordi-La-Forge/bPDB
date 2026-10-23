#include "bpdb.h"
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

void stat_bpdb(const char* filepath) {
    FILE* fb = fopen(filepath, "rb");
    if (!fb) { perror("Failed to open bpdb"); return; }
    fseek(fb, 0, SEEK_END); size_t size = ftell(fb); fseek(fb, 0, SEEK_SET);
    uint8_t* buffer = (uint8_t*)malloc(size); fread(buffer, 1, size, fb); fclose(fb);

    printf("\n=================================================\n");
    printf(" bPDB Format Verbatim Statistics\n");
    printf("=================================================\n");
    printf("Total File Size      : %zu bytes\n", size);
    if (size < 64) { free(buffer); return; }
    printf("Fixed Global Header  : 64 bytes\n");

    size_t offset = 64;
    uint16_t num_meta = *(uint16_t*)(buffer + offset); offset += 2;
    offset += num_meta * 2; offset += (8 - (offset % 8)) % 8;
    
    for (int i = 0; i < num_meta; i++) {
        uint64_t sec_len = *(uint64_t*)(buffer + offset);
        offset += 8 + sec_len + ((8 - (sec_len % 8)) % 8);
    }
    printf("Metadata Block (TLV) : %zu bytes (Sections: %d)\n", offset - 64, num_meta);

    int chunks = 0; size_t seq_bytes_total = 0; size_t atom_bytes_total = 0;
    int empty_atoms = 0; int present_atoms = 0;

    while (offset < size) {
        uint32_t num_atoms = *(uint32_t*)(buffer + offset);
        uint16_t num_residues = *(uint16_t*)(buffer + offset + 4);
        if (num_residues == 0xFFFF) {
            printf("-------------------------------------------------\n");
            printf("Coordinate Stream    : %d Chunks created\n", chunks);
            printf(" -> Chunk Headers    : %d bytes (Overhead)\n", chunks * 32);
            printf(" -> Sequence Tracks  : %zu bytes\n", seq_bytes_total);
            printf(" -> Atom Tracks      : %zu bytes\n", atom_bytes_total);
            printf("    * Present Atoms  : %d (%d bytes)\n", present_atoms, present_atoms * 8);
            printf("    * Empty Slots    : %d (%d bytes)\n", empty_atoms, empty_atoms * 8);
            
            offset += 32; size_t conect_start = offset;
            while (offset + 4 <= size) {
                int32_t N = *(int32_t*)(buffer + offset);
                if (N <= 0 || N > 4 || offset + 8 + (N << 2) > size) break;
                offset += 8 + (N << 2);
            }
            printf("CONECT Stream        : %zu bytes\n", offset - conect_start);
            break;
        }
        chunks++;
        uint32_t seq_bytes = *(uint32_t*)(buffer + offset + 28);
        size_t padded_seq = seq_bytes + ((8 - (seq_bytes % 8)) % 8);
        seq_bytes_total += padded_seq;
        offset += 32 + padded_seq;

        const uint8_t* atom_ptr = buffer + offset;
        for (uint32_t a = 0; a < num_atoms; a++) {
            if (*(uint64_t*)(atom_ptr + (a << 3)) == 0) empty_atoms++;
            else present_atoms++;
        }
        atom_bytes_total += num_atoms * 8;
        offset += num_atoms * 8;
    }
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
        printf("[+] Encoded %d files to %s\n", num_files, out_bpdb);
    } else if (strcmp(argv[1], "test") == 0) {
        FILE* fb = fopen(argv[2], "rb");
        if (!fb) { perror("Failed to open bpdb"); return 1; }
        
        fseek(fb, 0, SEEK_END); size_t size = ftell(fb); fseek(fb, 0, SEEK_SET);
        uint8_t* buffer = (uint8_t*)malloc(size); fread(buffer, 1, size, fb); fclose(fb);

        FILE* fout = fopen(argv[3], "w");
        decode_bpdb(buffer, size, fout); fclose(fout);
        printf("[+] Decoded to %s for diffing.\n", argv[3]);

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        int n = 1000;
        for (int i = 0; i < n*1000; i++) decode_bpdb(buffer, size, NULL); 
        clock_gettime(CLOCK_MONOTONIC, &end);
        
        double t_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
        printf("[+] Parsing took %.2f s total (%.3f us/run)\n", t_taken/n, t_taken/n);
        free(buffer);
    } else if (strcmp(argv[1], "stat") == 0) {
        stat_bpdb(argv[2]);
    }
    return 0;
}
