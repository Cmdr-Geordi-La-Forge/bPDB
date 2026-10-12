CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra
LDFLAGS = -lm

TARGET = bpdb_bench
SRC = main.c

all: $(TARGET)

$(TARGET): $(SRC) bpdb.h
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

test: $(TARGET) 4cha.pdb
	./$(TARGET) bench

4cha.pdb:
	curl -s https://files.rcsb.org/download/4CHA.pdb -o 4cha.pdb

clean:
	rm -f $(TARGET) 4cha.pdb 4cha.bpdb 4cha_back.pdb

.PHONY: all test clean
