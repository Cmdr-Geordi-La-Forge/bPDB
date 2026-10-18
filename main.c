#include "bpdb.h"
#include <string.h>
#include <time.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    if (argc < 4) {
        printf("Usage Encode: %s encode <in1.pdb> [in2.pdb...] <out.bpdb>\n", argv[0]);
        printf("Usage Test  : %s test <in.bpdb> <out.pdb>\n", argv[0]);
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
        
        fseek(fb, 0, SEEK_END);
        size_t size = ftell(fb);
        fseek(fb, 0, SEEK_SET);
        uint8_t* buffer = (uint8_t*)malloc(size);
        fread(buffer, 1, size, fb);
        fclose(fb);

        FILE* fout = fopen(argv[3], "w");
        decode_bpdb(buffer, size, fout);
        fclose(fout);
        printf("[+] Decoded to %s for diffing.\n", argv[3]);

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i = 0; i < 1000; i++) decode_bpdb(buffer, size, NULL); 
        clock_gettime(CLOCK_MONOTONIC, &end);
        
        double t_taken = (end.tv_sec - start.tv_sec) * 1e3 + (end.tv_nsec - start.tv_nsec) / 1e6;
        printf("[+] Parsing took %.2f ms total (%.3f us/run)\n", t_taken, t_taken * 1000.0 / 1000);
        
        free(buffer);
    }
    return 0;
}
