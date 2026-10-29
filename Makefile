CC = gcc
CFLAGS = -O3 -Wall -Wextra -march=native -flto -falign-loops=32
# CFLAGS = -g -O0 -Wall -Wextra

.PHONY: all clean

all: 4cha_back.pdb

bpdb: main.o bpdb_dict.o bpdb_encode.o bpdb_decode.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c bpdb.h
	$(CC) $(CFLAGS) -c $< -o $@

4cha.pdb:
	@wget -q -O 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb" || curl -s -o 4cha.pdb "https://files.rcsb.org/download/4CHA.pdb"

4cha.bpdb: bpdb 4cha.pdb
	@./bpdb encode 4cha.pdb 4cha.bpdb

4cha_back.pdb: bpdb 4cha.bpdb
	@./bpdb test 4cha.bpdb 4cha_back.pdb

clean:
	rm -f bpdb *.o *.bpdb 4cha_back.pdb
