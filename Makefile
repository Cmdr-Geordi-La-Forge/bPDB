CC = gcc
CFLAGS = -g -O3 -Wall -Wextra -march=native -flto -falign-loops=32

.PHONY: all clean test_rcsb test_af

all: 4cha_back.pdb AF-A0A011QK89-F1_back.pdb

bpdb: main.o bpdb_dict.o bpdb_encode.o bpdb_decode.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c bpdb.h
	$(CC) $(CFLAGS) -c $< -o $@

# --- RCSB Structural Validation ---
4cha.pdb:
	@wget -q -O 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb" || curl -s -o 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb"

4cha.bpdb: bpdb 4cha.pdb
	@./bpdb encode 4cha.pdb 4cha.bpdb

4cha_back.pdb: bpdb 4cha.bpdb
	@./bpdb test 4cha.bpdb 4cha_back.pdb

# --- AlphaFold DB Validation ---
AF-A0A011QK89-F1.pdb:
	@wget -q -O AF-A0A011QK89-F1.pdb "https://alphafold.ebi.ac.uk/files/AF-A0A011QK89-F1-model_v6.pdb" || curl -s -o AF-A0A011QK89-F1.pdb "https://alphafold.ebi.ac.uk/files/AF-A0A011QK89-F1-model_v6.pdb"

AF-A0A011QK89-F1.bpdb: bpdb AF-A0A011QK89-F1.pdb
	@./bpdb encode AF-A0A011QK89-F1.pdb AF-A0A011QK89-F1.bpdb

AF-A0A011QK89-F1_back.pdb: bpdb AF-A0A011QK89-F1.bpdb
	@./bpdb test AF-A0A011QK89-F1.bpdb AF-A0A011QK89-F1_back.pdb

clean:
	rm -f bpdb *.o *.bpdb *_back.pdb 4cha.pdb AF-A0A011QK89-F1.pdb
