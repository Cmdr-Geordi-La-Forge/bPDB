CC = gcc
CFLAGS = -g -O3 -Wall -Wextra -march=native -flto -falign-loops=32

.PHONY: all clean test_rcsb test_af bench

all: 4cha_back.pdb AF-A0A011QK89-F1_back.pdb

taxdict_bin.o: taxdict.bin
	ld -r -b binary $< -o $@

bpdb: main.o bpdb_dict.o bpdb_encode.o bpdb_decode.o taxdict_bin.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c bpdb.h
	$(CC) $(CFLAGS) -c $< -o $@

# --- RCSB Structural Validation ---
4cha.pdb:
	@wget -q -O 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb" || curl -s -o 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb"

4cha.bpdb: bpdb 4cha.pdb
	@./bpdb encode 4cha.pdb 4cha.bpdb

4cha_back.pdb: bpdb 4cha.bpdb
	@./bpdb decode 4cha.bpdb 4cha_back.pdb
	@./bpdb bench 4cha.bpdb

# --- AlphaFold DB Validation ---
AF-A0A011QK89-F1.pdb:
	@wget -q -O AF-A0A011QK89-F1.pdb "https://alphafold.ebi.ac.uk/files/AF-A0A011QK89-F1-model_v6.pdb" || curl -s -o AF-A0A011QK89-F1.pdb "https://alphafold.ebi.ac.uk/files/AF-A0A011QK89-F1-model_v6.pdb"

AF-A0A011QK89-F1.bpdb: bpdb AF-A0A011QK89-F1.pdb
	@./bpdb encode --afdb AF-A0A011QK89-F1.pdb AF-A0A011QK89-F1.bpdb

AF-A0A011QK89-F1_back.pdb: bpdb AF-A0A011QK89-F1.bpdb
	@./bpdb decode AF-A0A011QK89-F1.bpdb AF-A0A011QK89-F1_back.pdb
	@./bpdb bench AF-A0A011QK89-F1.bpdb

# --- Standalone Benchmark Target ---
bench: 4cha.bpdb AF-A0A011QK89-F1.bpdb
	@echo "--- Benchmarking RCSB ---"
	@./bpdb bench 4cha.bpdb
	@echo "--- Benchmarking AFDB ---"
	@./bpdb bench AF-A0A011QK89-F1.bpdb

clean:
	rm -f bpdb *.o *.bpdb *_back.pdb 4cha.pdb AF-A0A011QK89-F1.pdb
