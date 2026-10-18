#ifndef BPDB_H
#define BPDB_H

#include <stdint.h>
#include <stdio.h>

#define CTRL_TER        0x04
#define CTRL_TER_OXT    0x08

#define META_TITLE      0x0001
#define META_HET_DICT   0x0002
#define META_EXPDTA     0x0003
#define META_RESOLUTION 0x0004
#define META_PROGRAM    0x0005
#define META_SSBOND     0x0006
#define META_CISPEP     0x0007
#define META_HELIX      0x0008
#define META_SHEET      0x0009
#define META_SITE       0x000A

#define NUM_DICT_RES    38

typedef struct {
  char resName[4];
  int num_atoms;
  char atom_names[24][5];
  char elements[24][3];
} ResidueDict;

typedef struct {
  char resName[4];
  int num_atoms;
  char atom_names[256][5];
  char elements[256][3];
} UniDict;

// 8-byte structural anchor
typedef struct {
  char resName[4];
  char chain;
  char iCode;
  uint16_t resSeq;
} PdbResId;

typedef struct {
  PdbResId res1;
  PdbResId res2;
  char sym1[8];
  char sym2[8];
  float length;
} MetaSSBond;

typedef struct {
  PdbResId res1;
  PdbResId res2;
  int32_t modNum;
  float angle;
} MetaCispep;

typedef struct {
  PdbResId start;
  PdbResId end;
  char helixID[4];
  uint8_t helix_class;
  uint8_t pad;
  uint16_t length;
  char comment[32];
} MetaHelix;

typedef struct {
  PdbResId start;
  PdbResId end;
  char sheetID[4];
  int8_t num_strands;
  int8_t strand_num;
  int8_t sense;
  uint8_t pad;
} MetaSheet;

typedef struct {
  char siteID[4];
  uint32_t num_res;
} MetaSiteHeader;

extern const ResidueDict DICT[NUM_DICT_RES];
extern const char* MONTHS[];

typedef struct {
  uint8_t* data;
  size_t len;
  size_t cap;
} StreamBuffer;

void stream_init(StreamBuffer* buf);
void stream_write(StreamBuffer* buf, const void* src, size_t size);
void stream_free(StreamBuffer* buf);

int parse_month(const char* m);
void trim_string(char* str);

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb);
void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out);

#endif // BPDB_H
