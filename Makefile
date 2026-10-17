CC = gcc
CFLAGS = -O3 -Wall -Wextra -march=native

.PHONY: all clean

all: 4cha_back.pdb

bpdb: bpdb.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

4cha.pdb:
	@echo "Fetching 4CHA.pdb (if not present)..."
	@wget -q -O 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb" || curl -s -o 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb"

4cha.bpdb: bpdb 4cha.pdb
	@echo "Encoding 4cha.pdb to 4cha.bpdb..."
	@./bpdb encode 4cha.pdb 4cha.bpdb

4cha_back.pdb: bpdb 4cha.bpdb
	@echo "Decoding 4cha.bpdb to 4cha_back.pdb and running benchmark..."
	@./bpdb test 4cha.bpdb 4cha_back.pdb
	@echo ""
	@echo "✅ Done! 4cha.bpdb and 4cha_back.pdb are now in the folder."
	@echo "Note: A direct 'diff 4cha.pdb 4cha_back.pdb' will show differences because hydrogens and HETATMs were intentionally dropped."

clean:
	rm -f bpdb 4cha.bpdb 4cha_back.pdb 4cha_canonical.bpdb 4cha_canonical.pdb
