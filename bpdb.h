#ifndef BPDB_H
#define BPDB_H

#include <stdint.h>
#include <stdio.h>

// --- Metadata Bitmask Indices ---
// Engine & Format Flags (Bits 0-15 reserved for bPDB core)
#define META_AF_MODE    0
#define META_AFDB       1
#define META_HYDROGENS  2

// Section 2: Title Section (Bits 16-31)
#define META_TITLE      16
#define META_COMPND     17
#define META_SOURCE     18
#define META_EXPDTA     19
#define META_RESOLUTION 20
#define META_PROGRAM    21

// Section 3: Primary Structure Section (Bits 32-35)
#define META_DBREF      32
#define META_SEQRES     33

// Section 4: Heterogen Section (Bits 36-39)
#define META_HET_DICT   36

// Section 5: Secondary Structure Section (Bits 40-43)
#define META_HELIX      40
#define META_SHEET      41

// Section 6: Connectivity Annotation Section (Bits 44-47)
#define META_SSBOND     44
#define META_CISPEP     45

// Section 7: Miscellaneous Features Section (Bits 48-51)
#define META_SITE       48

// --- Metadata Enums ---
enum CompndToken { 
  COMPND_MOL_ID=1, COMPND_MOLECULE, COMPND_CHAIN, COMPND_FRAGMENT, 
  COMPND_SYNONYM, COMPND_EC, COMPND_ENGINEERED, COMPND_MUTATION, COMPND_OTHER_DETAILS 
};

enum SourceToken { 
  SOURCE_MOL_ID=1, SOURCE_SYNTHETIC, SOURCE_FRAGMENT, SOURCE_ORGANISM_SCIENTIFIC, 
  SOURCE_ORGANISM_COMMON, SOURCE_ORGANISM_TAXID, SOURCE_STRAIN, SOURCE_VARIANT, 
  SOURCE_CELL_LINE, SOURCE_ATCC, SOURCE_ORGAN, SOURCE_TISSUE, SOURCE_CELL, 
  SOURCE_ORGANELLE, SOURCE_SECRETION, SOURCE_CELLULAR_LOCATION, SOURCE_PLASMID, 
  SOURCE_GENE, SOURCE_EXPRESSION_SYSTEM, SOURCE_EXPRESSION_SYSTEM_COMMON, 
  SOURCE_EXPRESSION_SYSTEM_TAXID, SOURCE_EXPRESSION_SYSTEM_STRAIN, 
  SOURCE_EXPRESSION_SYSTEM_VARIANT, SOURCE_EXPRESSION_SYSTEM_CELL_LINE, 
  SOURCE_EXPRESSION_SYSTEM_ATCC_NUMBER, SOURCE_EXPRESSION_SYSTEM_ORGAN, 
  SOURCE_EXPRESSION_SYSTEM_TISSUE, SOURCE_EXPRESSION_SYSTEM_CELL, 
  SOURCE_EXPRESSION_SYSTEM_ORGANELLE, SOURCE_EXPRESSION_SYSTEM_CELLULAR_LOCATION, 
  SOURCE_EXPRESSION_SYSTEM_VECTOR_TYPE, SOURCE_EXPRESSION_SYSTEM_VECTOR, 
  SOURCE_EXPRESSION_SYSTEM_PLASMID, SOURCE_EXPRESSION_SYSTEM_GENE, SOURCE_OTHER_DETAILS}
  ;

// --- SEQRES Packed Header ---
typedef struct {
  char chain;
  uint16_t num_res;
} __attribute__((packed)) MetaSeqresHeader;

typedef struct {
  uint8_t type_id;
  char resName[4];
  
  // Heavy Atom Topology
  uint8_t h_start_norm;
  uint8_t h_end_norm;
  uint8_t h_end_cterm;
  char h_atoms[256][5];
  char h_elements[256][3];

  // Hydrogen Topology
  uint8_t hy_start_norm;
  uint8_t hy_end_norm;
  uint8_t hy_end_cterm;
  char hy_atoms[256][5];
  char hy_elements[256][3];
} ResidueDict;

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
  uint32_t pad;
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

extern const uint8_t _binary_taxdict_bin_start[];

// --- DBREF Packed Struct ---
typedef struct {
  char chain;             // 1 byte
  char db_name[8];        // 8 bytes (e.g., "UNP")
  char db_accession[12];  // 12 bytes
  char db_id_code[12];    // 12 bytes
  int32_t seq_begin;      // 4 bytes
  int32_t seq_end;        // 4 bytes
  int32_t db_seq_begin;   // 4 bytes
  int32_t db_seq_end;     // 4 bytes
  char iCode_begin;       // 1 byte
  char iCode_end;         // 1 byte
  char db_iCode_begin;    // 1 byte
  char db_iCode_end;      // 1 byte
} __attribute__((packed)) MetaDBRef; // 52 bytes

// --- NCBI Taxonomy Dictionary (Zero-Copy) ---
typedef struct {
    uint32_t max_id;
    const uint32_t* sci_offs;
    const uint32_t* com_offs;
    const char* string_pool;
} TaxDict;

static inline TaxDict load_taxdict(void) {
    TaxDict dict;
    const uint8_t* ptr = _binary_taxdict_bin_start;
    dict.max_id = *(const uint32_t*)ptr; ptr += 4;
    dict.sci_offs = (const uint32_t*)ptr; ptr += (dict.max_id + 1) * 4;
    dict.com_offs = (const uint32_t*)ptr; ptr += (dict.max_id + 1) * 4;
    dict.string_pool = (const char*)ptr;
    return dict;
}

static inline const char* get_tax_name(const TaxDict* dict, uint32_t taxid, int is_common) {
    if (taxid > dict->max_id) return NULL;
    uint32_t offset = is_common ? dict->com_offs[taxid] : dict->sci_offs[taxid];
    if (offset == 0 && taxid != 0) return NULL; 
    return dict->string_pool + offset;
}

extern const ResidueDict DICT[];
extern const size_t NUM_DICT_RES;
extern const char* MONTHS[];

typedef struct {
  uint8_t* data;
  size_t len;
  size_t cap;
} StreamBuffer;

// --- Internal Contexts for Refactoring ---

typedef struct {
    uint64_t mask;
    int af_mode;
    int has_hydrogens;
    float resolution;

    const char* title;      uint64_t title_len;
    const char* expdta;     uint64_t expdta_len;
    const char* program;
    
    const uint8_t* compnd;  uint64_t compnd_len;
    const uint8_t* source;  uint64_t source_len;
    const uint8_t* seqres;  uint64_t seqres_len;
    const uint8_t* dbrefs;  uint64_t num_dbrefs;
    const uint8_t* ssbonds; uint64_t num_ssbonds;
    const uint8_t* cispeps; uint64_t num_cispeps;
    const uint8_t* helices; uint64_t num_helices;
    const uint8_t* sheets;  uint64_t num_sheets;
    const uint8_t* sites;   uint64_t len_sites;

    size_t payload_start_offset;
    size_t conect_offset;
} DecodeMetadata;

typedef struct {
    const uint8_t* buffer;
    size_t size;
    FILE* out;
    int is_first_model;
    int is_last_model;

    int day, month, year;
    char pdbID[13];
    char classification[41];
    uint16_t total_models_minus_1;
    uint16_t current_model_minus_1;
    
    size_t total_out_lines;
} DecodeContext;

void stream_init(StreamBuffer* buf);
void stream_write(StreamBuffer* buf, const void* src, size_t size);
void stream_free(StreamBuffer* buf);

int parse_month(const char* m);
void trim_string(char* str);

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb, int include_hydrogens, int af_mode, int afdb_mode);
void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out, int is_first_model, int is_last_model);

#endif // BPDB_H
