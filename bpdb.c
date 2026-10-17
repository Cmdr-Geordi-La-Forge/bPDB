#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#if defined(__AVX2__) || defined(__SSE4_1__)
#include <immintrin.h>
#elif defined(__aarch64__) || defined(__ARM_NEON)
#include <arm_neon.h>
#endif

// -----------------------------------------------------------------------------
// 1. Dictionary and Data Structures
// -----------------------------------------------------------------------------
typedef struct {
    char resName[4];
    int num_atoms;
    char atom_names[24][5];
    char elements[24][3];
} ResidueDict;

const ResidueDict DICT[37] = {
    {"ALA", 5,  {" N  "," CA "," C  "," O  "," CB "}, {" N"," C"," C"," O"," C"}},
    {"ARG", 11, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2"}, {" N"," C"," C"," O"," C"," C"," C"," N"," C"," N"," N"}},
    {"ASN", 8,  {" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2"}, {" N"," C"," C"," O"," C"," C"," O"," N"}},
    {"ASP", 8,  {" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2"}, {" N"," C"," C"," O"," C"," C"," O"," O"}},
    {"CYS", 6,  {" N  "," CA "," C  "," O  "," CB "," SG "}, {" N"," C"," C"," O"," C"," S"}},
    {"GLN", 9,  {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2"}, {" N"," C"," C"," O"," C"," C"," C"," O"," N"}},
    {"GLU", 9,  {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2"}, {" N"," C"," C"," O"," C"," C"," C"," O"," O"}},
    {"GLY", 4,  {" N  "," CA "," C  "," O  "}, {" N"," C"," C"," O"}},
    {"HIS", 10, {" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2"}, {" N"," C"," C"," O"," C"," C"," N"," C"," C"," N"}},
    {"ILE", 8,  {" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1"}, {" N"," C"," C"," O"," C"," C"," C"," C"}},
    {"LEU", 8,  {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"}, {" N"," C"," C"," O"," C"," C"," C"," C"}},
    {"LYS", 9,  {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ "}, {" N"," C"," C"," O"," C"," C"," C"," C"," N"}},
    {"MET", 8,  {" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE "}, {" N"," C"," C"," O"," C"," C"," S"," C"}},
    {"PHE", 11, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "}, {" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," C"}},
    {"PRO", 7,  {" N  "," CA "," C  "," O  "," CB "," CG "," CD "}, {" N"," C"," C"," O"," C"," C"," C"}},
    {"SER", 6,  {" N  "," CA "," C  "," O  "," CB "," OG "}, {" N"," C"," C"," O"," C"," O"}},
    {"THR", 7,  {" N  "," CA "," C  "," O  "," CB "," OG1"," CG2"}, {" N"," C"," C"," O"," C"," O"," C"}},
    {"TRP", 14, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"}, {" N"," C"," C"," O"," C"," C"," C"," C"," N"," C"," C"," C"," C"," C"}},
    {"TYR", 12, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH "}, {" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," C"," O"}},
    {"VAL", 7,  {" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"}, {" N"," C"," C"," O"," C"," C"," C"}},
    {"SEC", 6,  {" N  "," CA "," C  "," O  "," CB "," SE "}, {" N"," C"," C"," O"," C","SE"}},
    {"PYL", 22, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ "," C2 "," C3 "," C4 "," C5 "," C6 "," O2 "," CA2"," CB2"," CG2"," CD2"," CE2"," CE3"," O22"}, {" N"," C"," C"," O"," C"," C"," C"," C"," N"," C"," C"," C"," C"," C"," O"," C"," C"," C"," C"," C"," C"," O"}},
    {"ASX", 9,  {" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2"," ND2"}, {" N"," C"," C"," O"," C"," C"," O"," O"," N"}},
    {"GLX", 10, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2"," NE2"}, {" N"," C"," C"," O"," C"," C"," C"," O"," O"," N"}},
    {"UNK", 14, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," CZ "," CH "," X1 "," X2 "," X3 "," X4 "}, {" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," X"," X"," X"," X"}},
    {"  A", 22, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N6 "," N1 "," C2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," N"," N"," C"," N"," C"}},
    {"  C", 20, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," O"," N"," C"," N"," C"," C"}},
    {"  G", 23, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," N"," C"}},
    {"  U", 20, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," O"," N"," C"," O"," C"," C"}},
    {"  I", 22, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," C"}},
    {" DA", 21, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," N6 "," N1 "," C2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," N"," N"," C"," N"," C"}},
    {" DC", 19, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," N"," C"," C"}},
    {" DG", 22, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," N"," C"}},
    {" DT", 20, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C7 "," C6 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," O"," C"," C"," C"}},
    {" DU", 19, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," O"," C"," C"}},
    {" DI", 21, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," C"}},
    {"  N", 24, {" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N6 "," N1 "," C2 "," N2 "," N3 "," C4 "}, {" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," N"," C"," N"," N"," C"}}
};

const char* MONTHS[] = {"", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

int parse_month(const char* m) {
    for (int i = 1; i <= 12; i++) if (strncmp(m, MONTHS[i], 3) == 0) return i;
    return 1;
}

int get_res_type(const char* name) {
    for (int i = 0; i < 37; i++) {
        if (strncmp(name, DICT[i].resName, 3) == 0) return i;
    }
    return -1;
}

void trim_string(char* str) {
    int i = strlen(str) - 1;
    while(i >= 0 && isspace(str[i])) { str[i] = '\0'; i--; }
}

// -----------------------------------------------------------------------------
// 2. Linear Stream Memory Buffer
// -----------------------------------------------------------------------------
typedef struct {
    uint8_t* data;
    size_t len;
    size_t cap;
} StreamBuffer;

void stream_init(StreamBuffer* buf) {
    buf->cap = 1024 * 1024 * 16;
    buf->data = (uint8_t*)malloc(buf->cap);
    buf->len = 0;
}

void stream_write(StreamBuffer* buf, const void* src, size_t size) {
    if (buf->len + size > buf->cap) {
        buf->cap *= 2;
        buf->data = (uint8_t*)realloc(buf->data, buf->cap);
    }
    memcpy(buf->data + buf->len, src, size);
    buf->len += size;
}

void stream_free(StreamBuffer* buf) {
    free(buf->data);
}

void write_padded_string(FILE* f, const char* str) {
    uint64_t len = strlen(str);
    fwrite(&len, 1, 8, f);
    if (len > 0) fwrite(str, 1, len, f);
    size_t pad = (8 - (len % 8)) % 8;
    if (pad > 0) {
        uint64_t zero = 0;
        fwrite(&zero, 1, pad, f);
    }
}

// -----------------------------------------------------------------------------
// 3. Encoder
// -----------------------------------------------------------------------------
typedef struct {
    int resSeq;
    char iCode;
    int type;
    int32_t b_acc[24], x_acc[24], y_acc[24], z_acc[24];
    uint8_t occ_acc[24];
    char alt_acc[24];
    int8_t charge_acc[24];
    int present[24];
    
    int has_oxt;
    int32_t oxt_b, oxt_x, oxt_y, oxt_z;
    uint8_t oxt_occ;
    char oxt_alt;
    int8_t oxt_charge;
} EncoderState;

int chunk_residues = 0;
size_t chunk_header_offset = 0;
int chunk_active = 0;
uint8_t chunk_flags = 0;

#define CTRL_TER        0x04
#define CTRL_TER_OXT    0x08

void close_chunk(StreamBuffer* buf) {
    if (!chunk_active) return;
    uint16_t count = chunk_residues;
    memcpy(buf->data + chunk_header_offset, &count, 2);
    memcpy(buf->data + chunk_header_offset + 2, &chunk_flags, 1);
    chunk_active = 0;
    chunk_residues = 0;
    chunk_flags = 0;
}

void open_chunk(StreamBuffer* buf, uint8_t frag_id, char c1, char c2, const char* file) {
    close_chunk(buf);
    chunk_header_offset = buf->len;
    
    uint8_t header[16] = {0};
    header[3] = frag_id;
    header[4] = c1;
    header[5] = c2;
    size_t flen = strlen(file);
    memcpy(header + 6, file, flen > 10 ? 10 : flen);
    
    stream_write(buf, header, 16);
    chunk_active = 1;
    chunk_residues = 0;
    chunk_flags = 0;
}

static void write_residue_buffer(StreamBuffer* buf, EncoderState* st) {
    if (st->type == -1) return;
    chunk_residues++;
    
    int32_t base_b = 0, base_x = 0, base_y = 0, base_z = 0;
    for (int i = 0; i < DICT[st->type].num_atoms; i++) {
        if (st->present[i]) {
            base_b = st->b_acc[i]; base_x = st->x_acc[i]; base_y = st->y_acc[i]; base_z = st->z_acc[i];
            break;
        }
    }

    uint8_t res_block[8] = {0};
    res_block[0] = (uint8_t)st->type;
    res_block[1] = (base_b >> 13) & 0xFF;
    res_block[2] = st->resSeq & 0xFF;
    res_block[3] = (base_x >> 13) & 0xFF;
    res_block[4] = (st->resSeq >> 8) & 0xFF;
    res_block[5] = (base_y >> 13) & 0xFF; 
    res_block[6] = st->iCode;
    res_block[7] = (base_z >> 13) & 0xFF;
    stream_write(buf, res_block, 8);

    int32_t last_b = 0, last_x = 0, last_y = 0, last_z = 0;
    for (int i = 0; i < DICT[st->type].num_atoms; i++) {
        uint8_t atom_block[8] = {0};
        uint8_t occ = st->present[i] ? (st->occ_acc[i] & 0x7F) : 0;
        char alt = st->alt_acc[i];
        int8_t charge = st->charge_acc[i];
        
        uint8_t alt_idx = (alt == ' ' || alt == 0) ? 0 : (alt - 'A' + 1);
        if (alt_idx > 7) alt_idx = 7;
        uint8_t alt_charge = (alt_idx << 3) | (charge & 0x07);
        
        uint16_t L0, L1, L2, L3;
        if (i == 0) { 
            L0 = ((base_b & 0x1FFF) << 3) | (occ & 0x07);
            L1 = ((base_x & 0x1FFF) << 3) | ((occ >> 3) & 0x07);
            L2 = ((base_y & 0x1FFF) << 3) | (((occ >> 6) & 0x01) | ((alt_charge & 0x01) << 1) | ((alt_charge & 0x02) << 1));
            L3 = ((base_z & 0x1FFF) << 3) | ((alt_charge >> 2) & 0x07);
            last_b = base_b; last_x = base_x; last_y = base_y; last_z = base_z;
        } else {
            if (!st->present[i]) {
                stream_write(buf, atom_block, 8); 
                continue;
            }
            
            int32_t db = st->b_acc[i] - last_b; int32_t dx = st->x_acc[i] - last_x;
            int32_t dy = st->y_acc[i] - last_y; int32_t dz = st->z_acc[i] - last_z;
            
            if (db > 4095) db = 4095; else if (db < -4096) db = -4096;
            if (dx > 4095) dx = 4095; else if (dx < -4096) dx = -4096;
            if (dy > 4095) dy = 4095; else if (dy < -4096) dy = -4096;
            if (dz > 4095) dz = 4095; else if (dz < -4096) dz = -4096;
            
            L0 = ((db & 0x1FFF) << 3) | (occ & 0x07);
            L1 = ((dx & 0x1FFF) << 3) | ((occ >> 3) & 0x07);
            L2 = ((dy & 0x1FFF) << 3) | (((occ >> 6) & 0x01) | ((alt_charge & 0x01) << 1) | ((alt_charge & 0x02) << 1));
            L3 = ((dz & 0x1FFF) << 3) | ((alt_charge >> 2) & 0x07);
            
            last_b = st->b_acc[i]; last_x = st->x_acc[i]; last_y = st->y_acc[i]; last_z = st->z_acc[i];
        }

        atom_block[0] = L0 & 0xFF; atom_block[1] = L0 >> 8;
        atom_block[2] = L1 & 0xFF; atom_block[3] = L1 >> 8;
        atom_block[4] = L2 & 0xFF; atom_block[5] = L2 >> 8;
        atom_block[6] = L3 & 0xFF; atom_block[7] = L3 >> 8;
        stream_write(buf, atom_block, 8);
    }
}

void encode_bpdb(int num_files, const char** in_pdbs, const char* out_bpdb) {
    StreamBuffer coords;
    stream_init(&coords);

    char global_title[4096] = {0};
    char global_id[128] = {0};
    char global_class[128] = {0};
    int day=1, month=1, year=1973;
    int id_from_header = 0;

    EncoderState st;
    memset(&st, 0, sizeof(st));
    st.type = -1;

    char current_c1 = '\0', current_c2 = '\0';
    int force_new_chunk = 1;

    for (int f = 0; f < num_files; f++) {
        FILE* fin = fopen(in_pdbs[f], "r");
        if (!fin) { perror("Input file error"); exit(1); }

        char line[128];
        char local_id[16] = {0};
        uint8_t local_frag = 0;

        const char* f_ptr = strstr(in_pdbs[f], "-F");
        if (f_ptr && isdigit(f_ptr[2])) local_frag = atoi(f_ptr + 2);

        force_new_chunk = 1;

        while (fgets(line, sizeof(line), fin)) {
            if (f == 0 && strncmp(line, "HEADER", 6) == 0) {
                size_t len = strlen(line);
                if (len >= 50) {
                    memcpy(global_class, line + 10, 40);
                    trim_string(global_class);
                }
                if (len >= 59) {
                    char d[3]={0}, m[4]={0}, y[3]={0};
                    memcpy(d, line+50, 2); memcpy(m, line+53, 3); memcpy(y, line+57, 2);
                    day = atoi(d); month = parse_month(m); int parsed_y = atoi(y);
                    year = parsed_y + ((parsed_y < 50) ? 2000 : 1900);
                }
                if (len >= 66) {
                    char temp_id[5] = {0};
                    memcpy(temp_id, line+62, 4);
                    trim_string(temp_id);
                    if (strlen(temp_id) > 0 && strcmp(temp_id, "XXXX") != 0) {
                        strcpy(global_id, temp_id);
                        id_from_header = 1;
                    }
                }
            }
            if (f == 0 && strncmp(line, "TITLE ", 6) == 0) {
                char chunk[71] = {0};
                size_t line_len = strlen(line);
                if (line_len > 10) {
                    size_t c_len = line_len - 10;
                    if (c_len > 70) c_len = 70;
                    memcpy(chunk, line+10, c_len);
                    trim_string(chunk);
                    if (strlen(global_title) > 0) strcat(global_title, " ");
                    strcat(global_title, chunk);
                }
                
                if (!id_from_header) {
                    char* p1 = strchr(line, '(');
                    char* p2 = strchr(p1 ? p1 : line, ')');
                    if (p1 && p2 && (p2 > p1 + 4)) {
                        size_t len = p2 - p1 - 1;
                        if (len > 15) len = 15;
                        memcpy(global_id, p1+1, len); global_id[len] = '\0';
                    }
                }
            }
            if (strncmp(line, "DBREF ", 6) == 0) {
                char* unp = strstr(line, "UNP");
                if (unp) {
                    char dummy_str[16], uid[64];
                    if (sscanf(unp, "%s %s", dummy_str, uid) == 2 && strlen(uid) > 3) {
                        if (f == 0 && !id_from_header) {
                            size_t c = strlen(uid);
                            memcpy(global_id, uid, c > 15 ? 15 : c);
                            global_id[c > 15 ? 15 : c] = '\0';
                        }
                        size_t l = strlen(uid);
                        memcpy(local_id, uid, l > 10 ? 10 : l);
                    }
                }
            }

            if (strncmp(line, "TER   ", 6) == 0) {
                write_residue_buffer(&coords, &st);
                
                if (st.has_oxt) {
                    chunk_flags |= CTRL_TER_OXT;
                    int32_t base_b = 0, base_x = 0, base_y = 0, base_z = 0;
                    for (int i = 0; i < DICT[st.type].num_atoms; i++) {
                        if (st.present[i]) {
                            base_b = st.b_acc[i]; base_x = st.x_acc[i]; base_y = st.y_acc[i]; base_z = st.z_acc[i];
                            break;
                        }
                    }
                    int32_t db = st.oxt_b - base_b; int32_t dx = st.oxt_x - base_x;
                    int32_t dy = st.oxt_y - base_y; int32_t dz = st.oxt_z - base_z;
                    
                    if (db > 4095) db = 4095; else if (db < -4096) db = -4096;
                    if (dx > 4095) dx = 4095; else if (dx < -4096) dx = -4096;
                    if (dy > 4095) dy = 4095; else if (dy < -4096) dy = -4096;
                    if (dz > 4095) dz = 4095; else if (dz < -4096) dz = -4096;

                    uint8_t occ = st.oxt_occ & 0x7F;
                    uint8_t alt = (st.oxt_alt != ' ' && st.oxt_alt != 0) ? 1 : 0;
                    uint8_t charge = st.oxt_charge & 0x0F;
                    
                    uint8_t alt_idx = (alt == ' ' || alt == 0) ? 0 : (alt - 'A' + 1);
                    if (alt_idx > 7) alt_idx = 7;
                    uint8_t alt_charge = (alt_idx << 3) | (charge & 0x07);
                    
                    uint16_t L0 = ((db & 0x1FFF) << 3) | (occ & 0x07);
                    uint16_t L1 = ((dx & 0x1FFF) << 3) | ((occ >> 3) & 0x07);
                    uint16_t L2 = ((dy & 0x1FFF) << 3) | (((occ >> 6) & 0x01) | ((alt_charge & 0x01) << 1) | ((alt_charge & 0x02) << 1));
                    uint16_t L3 = ((dz & 0x1FFF) << 3) | ((alt_charge >> 2) & 0x07);
                    
                    uint8_t oxt_block[8];
                    oxt_block[0] = L0 & 0xFF; oxt_block[1] = L0 >> 8;
                    oxt_block[2] = L1 & 0xFF; oxt_block[3] = L1 >> 8;
                    oxt_block[4] = L2 & 0xFF; oxt_block[5] = L2 >> 8;
                    oxt_block[6] = L3 & 0xFF; oxt_block[7] = L3 >> 8;
                    stream_write(&coords, oxt_block, 8);
                } else {
                    chunk_flags |= CTRL_TER;
                }
                
                close_chunk(&coords);
                st.type = -1;
                st.has_oxt = 0;
                force_new_chunk = 1;
                continue;
            }
            
            if (strncmp(line, "ATOM  ", 6) != 0) continue;

            char atom_name[5] = {0}, res_name[4] = {0};
            memcpy(atom_name, line + 12, 4);
            memcpy(res_name, line + 17, 3);
            char altLoc = line[16];
            char iCode = line[26];
            int resSeq = atoi(line + 22);
            
            char c1 = line[20];
            char c2 = line[21];

            char clean_atom[5] = {0};
            memcpy(clean_atom, line + 12, 4);
            if (clean_atom[1] == 'H' || clean_atom[0] == 'H') continue;

            int type = get_res_type(res_name);
            if (type == -1) continue;

            if (strlen(local_id) == 0) {
                size_t g_len = strlen(global_id);
                memcpy(local_id, global_id, g_len > 10 ? 10 : g_len);
            }

            int new_chain = (c1 != current_c1 || c2 != current_c2);

            if (resSeq != st.resSeq || iCode != st.iCode || new_chain || force_new_chunk || chunk_residues >= 60000 || type != st.type) {
                write_residue_buffer(&coords, &st);
                
                if (new_chain || force_new_chunk || chunk_residues >= 60000) {
                    if (new_chain) { current_c1 = c1; current_c2 = c2; }
                    open_chunk(&coords, local_frag, c1, c2, local_id);
                    force_new_chunk = 0;
                }
                
                st.resSeq = resSeq;
                st.iCode = iCode;
                st.type = type;
                memset(st.present, 0, sizeof(st.present));
                st.has_oxt = 0;
            }

            float x = atof(line + 30), y = atof(line + 38), z = atof(line + 46);
            float occ = atof(line + 54), temp = atof(line + 60);

            int charge = 0;
            if (line[78] >= '0' && line[78] <= '9') {
                charge = line[78] - '0';
                if (line[79] == '-') charge = -charge;
            }
            if (charge < -8) charge = -8; else if (charge > 7) charge = 7;

            if (strncmp(clean_atom, " OXT", 4) == 0) {
                st.has_oxt = 1;
                st.oxt_b = (int32_t)round(temp * 100.0f);
                st.oxt_x = (int32_t)round(x * 1000.0f);
                st.oxt_y = (int32_t)round(y * 1000.0f);
                st.oxt_z = (int32_t)round(z * 1000.0f);
                st.oxt_occ = (uint8_t)round(occ * 100.0f);
                st.oxt_alt = altLoc;
                st.oxt_charge = charge;
                continue;
            }

            int atom_idx = -1;
            for (int i = 0; i < DICT[type].num_atoms; i++) {
                if (strncmp(clean_atom, DICT[type].atom_names[i], 4) == 0) { atom_idx = i; break; }
            }
            if (atom_idx == -1 || st.present[atom_idx]) continue;

            st.b_acc[atom_idx] = (int32_t)round(temp * 100.0f);
            st.x_acc[atom_idx] = (int32_t)round(x * 1000.0f);
            st.y_acc[atom_idx] = (int32_t)round(y * 1000.0f);
            st.z_acc[atom_idx] = (int32_t)round(z * 1000.0f);
            st.occ_acc[atom_idx] = (uint8_t)round(occ * 100.0f);
            st.alt_acc[atom_idx] = altLoc;
            st.charge_acc[atom_idx] = charge;
            st.present[atom_idx] = 1;
        }
        fclose(fin);
    }
    write_residue_buffer(&coords, &st);
    close_chunk(&coords);

    FILE* fout = fopen(out_bpdb, "wb");
    uint64_t zero = 0;
    
    // 1. Write Fixed Header (8 Bytes)
    uint8_t header[8] = {'b','P','D','B', 1, (uint8_t)day, (uint8_t)month, 
                         (year >= 1973 && year <= 2228) ? (uint8_t)(year - 1973) : 0};
    fwrite(header, 1, 8, fout);

    // 2. Mandatory Blocks (ID, Class, Title)
    write_padded_string(fout, global_id);
    write_padded_string(fout, global_class);
    write_padded_string(fout, global_title);

    // 3. Extensible Metadata Section (TLV) - Empty for now
    uint16_t num_meta = 0;
    fwrite(&num_meta, 1, 2, fout);
    uint16_t m_pad = (8 - (2 % 8)) % 8;
    if (m_pad > 0) fwrite(&zero, 1, m_pad, fout);

    // 4. Coordinate Stream
    fwrite(coords.data, 1, coords.len, fout);

    fclose(fout);
    stream_free(&coords);
}

// -----------------------------------------------------------------------------
// 4. High-Performance Decoder
// -----------------------------------------------------------------------------

#if defined(__AVX2__) || defined(__SSE4_1__)
    typedef __m128i accum_t;
    #define GET_SIMD_B(acc) _mm_extract_epi32((acc), 0)
    #define GET_SIMD_X(acc) _mm_extract_epi32((acc), 1)
    #define GET_SIMD_Y(acc) _mm_extract_epi32((acc), 2)
    #define GET_SIMD_Z(acc) _mm_extract_epi32((acc), 3)
#elif defined(__aarch64__) || defined(__ARM_NEON)
    typedef int32x4_t accum_t;
    #define GET_SIMD_B(acc) vgetq_lane_s32((acc), 0)
    #define GET_SIMD_X(acc) vgetq_lane_s32((acc), 1)
    #define GET_SIMD_Y(acc) vgetq_lane_s32((acc), 2)
    #define GET_SIMD_Z(acc) vgetq_lane_s32((acc), 3)
#else
    typedef struct { int32_t b, x, y, z; } accum_t;
    #define GET_SIMD_B(acc) ((acc).b)
    #define GET_SIMD_X(acc) ((acc).x)
    #define GET_SIMD_Y(acc) ((acc).y)
    #define GET_SIMD_Z(acc) ((acc).z)
#endif

// Global sink to prevent dead-code elimination during benchmark timing
volatile int32_t benchmark_sink = 0;

void decode_bpdb(const uint8_t* buffer, size_t size, FILE* out) {
    if (size < 8) return;
    
    int day = buffer[5], month = buffer[6];
    int year = 1973 + buffer[7];

    size_t offset = 8;
    
    // Read Mandatory Strings
    uint64_t id_len = *(uint64_t*)(buffer + offset); offset += 8;
    const char* p_idCode = (const char*)(buffer + offset);
    offset += id_len; offset += (8 - (id_len % 8)) % 8;

    uint64_t class_len = *(uint64_t*)(buffer + offset); offset += 8;
    const char* p_class = (const char*)(buffer + offset);
    offset += class_len; offset += (8 - (class_len % 8)) % 8;

    uint64_t title_len = *(uint64_t*)(buffer + offset); offset += 8;
    const char* p_title = (const char*)(buffer + offset);
    offset += title_len; offset += (8 - (title_len % 8)) % 8;

    // Read TLV Meta
    uint16_t num_meta = *(uint16_t*)(buffer + offset); offset += 2;
    offset += (8 - (2 % 8)) % 8;

    for (int i = 0; i < num_meta; i++) {
        offset += 8; // skip type code
        uint64_t section_len = *(uint64_t*)(buffer + offset); offset += 8;
        offset += section_len;
        offset += (8 - (section_len % 8)) % 8;
    }

    if (!out) {
        // Benchmark mode: just run through decoding without disk output
        while (offset < size) {
            uint16_t num_residues = *(uint16_t*)(buffer + offset);
            uint8_t flags         = buffer[offset + 2];
            offset += 16;
            
            int32_t last_a0_b = 0, last_a0_x = 0, last_a0_y = 0, last_a0_z = 0;
            
            for (uint32_t r = 0; r < num_residues; r++) {
                const uint8_t* res_ptr = buffer + offset;
                int res_type = res_ptr[0];
                int num_atoms = DICT[res_type].num_atoms;
                offset += 8;

                const uint8_t* a1_ptr = buffer + offset;
                accum_t accum;
                
                int32_t b_m = (int8_t)res_ptr[1]; int32_t x_m = (int8_t)res_ptr[3];
                int32_t y_m = (int8_t)res_ptr[5]; int32_t z_m = (int8_t)res_ptr[7];
                int32_t L0 = a1_ptr[0] | (a1_ptr[1] << 8); int32_t L1 = a1_ptr[2] | (a1_ptr[3] << 8);
                int32_t L2 = a1_ptr[4] | (a1_ptr[5] << 8); int32_t L3 = a1_ptr[6] | (a1_ptr[7] << 8);
                int32_t b_l = (L0 >> 3) & 0x1FFF; int32_t x_l = (L1 >> 3) & 0x1FFF;
                int32_t y_l = (L2 >> 3) & 0x1FFF; int32_t z_l = (L3 >> 3) & 0x1FFF;
                
#if defined(__AVX2__) || defined(__SSE4_1__)
                accum = _mm_set_epi32((z_m << 13) | z_l, (y_m << 13) | y_l, (x_m << 13) | x_l, (b_m << 13) | b_l);
#elif defined(__aarch64__) || defined(__ARM_NEON)
                int32x4_t setup = {(b_m << 13) | b_l, (x_m << 13) | x_l, (y_m << 13) | y_l, (z_m << 13) | z_l}; accum = setup;
#else
                accum.b = (b_m << 13) | b_l; accum.x = (x_m << 13) | x_l; accum.y = (y_m << 13) | y_l; accum.z = (z_m << 13) | z_l;
#endif

                last_a0_b = GET_SIMD_B(accum); last_a0_x = GET_SIMD_X(accum);
                last_a0_y = GET_SIMD_Y(accum); last_a0_z = GET_SIMD_Z(accum);

                offset += 8;
                for (int i = 1; i < num_atoms; i++) {
                    const uint8_t* an_ptr = buffer + offset;
#if defined(__AVX2__) || defined(__SSE4_1__)
                    __m128i an_vec = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr));
                    __m128i delta = _mm_srai_epi32(an_vec, 3);
                    accum = _mm_add_epi32(accum, delta); 
#elif defined(__aarch64__) || defined(__ARM_NEON)
                    int16x4_t an_vec_16 = vreinterpret_s16_u8(vld1_u8(an_ptr));
                    int32x4_t an_vec_32 = vmovl_s16(an_vec_16);
                    int32x4_t delta = vshrq_n_s32(an_vec_32, 3);
                    accum = vaddq_s32(accum, delta);
#else
                    int32_t aL0 = an_ptr[0] | (an_ptr[1] << 8); int32_t aL1 = an_ptr[2] | (an_ptr[3] << 8);
                    int32_t aL2 = an_ptr[4] | (an_ptr[5] << 8); int32_t aL3 = an_ptr[6] | (an_ptr[7] << 8);
                    accum.b += ((int16_t)aL0) >> 3; accum.x += ((int16_t)aL1) >> 3;
                    accum.y += ((int16_t)aL2) >> 3; accum.z += ((int16_t)aL3) >> 3;
#endif
                    offset += 8;
                }
                benchmark_sink += GET_SIMD_X(accum);
            }
            if (flags & CTRL_TER_OXT) {
                // Ensure the variables are "used" to prevent warnings
                benchmark_sink += last_a0_b + last_a0_x + last_a0_y + last_a0_z;
                offset += 8;
            }
        }
        return;
    }

    // Phase 1: Pre-calculate exact number of lines needed for the output
    size_t out_lines = 0;
    if (id_len > 0 || class_len > 0 || title_len > 0) {
        out_lines++; // HEADER
        if (title_len > 0) out_lines += (title_len + 69) / 70;
    }

    size_t scan_off = offset;
    while (scan_off < size) {
        uint16_t num_residues = *(uint16_t*)(buffer + scan_off);
        uint8_t flags = buffer[scan_off + 2];
        scan_off += 16;
        for (uint32_t r = 0; r < num_residues; r++) {
            int res_type = buffer[scan_off];
            int num_atoms = DICT[res_type].num_atoms;
            scan_off += 8;
            for (int i = 0; i < num_atoms; i++) {
                uint8_t occ_byte = buffer[scan_off + 6]; 
                if ((occ_byte & 0x01) || buffer[scan_off] != 0) out_lines++;
                scan_off += 8;
            }
        }
        if (flags & CTRL_TER_OXT) { out_lines++; scan_off += 8; }
        if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) out_lines++;
    }

    // Phase 2: Memory-mapped format
    size_t out_size = out_lines * 81;
    char* out_data = (char*)malloc(out_size);
    memset(out_data, ' ', out_size);
    for (size_t i = 0; i < out_lines; i++) out_data[i * 81 + 80] = '\n';

    char* w_ptr = out_data;
    char tbuf[128];

    if (id_len > 0 || class_len > 0 || title_len > 0) {
        memcpy(w_ptr, "HEADER    ", 10);
        if (class_len > 0) memcpy(w_ptr + 10, p_class, class_len > 40 ? 40 : class_len);
        
        int dlen = snprintf(tbuf, 128, "%02d-%3s-%02d", day, MONTHS[month], year % 100);
        memcpy(w_ptr + 50, tbuf, dlen);
        
        if (id_len > 0) memcpy(w_ptr + 62, p_idCode, id_len > 18 ? 18 : id_len);
        w_ptr += 81;
        
        if (title_len > 0) {
            for (uint64_t i = 0; i < title_len; i += 70) {
                memcpy(w_ptr, "TITLE   ", 8);
                if (i > 0) {
                    int ln = snprintf(tbuf, 128, "%2d", (int)(i/70)+1);
                    memcpy(w_ptr + 8, tbuf, ln);
                }
                size_t copy_len = title_len - i;
                if (copy_len > 70) copy_len = 70;
                memcpy(w_ptr + 10, p_title + i, copy_len);
                w_ptr += 81;
            }
        }
    }

    int atom_serial = 1;

    while (offset < size) {
        uint16_t num_residues = *(uint16_t*)(buffer + offset);
        uint8_t flags         = buffer[offset + 2];
        char c1               = buffer[offset + 4];
        char c2               = buffer[offset + 5];
        offset += 16;
        
        if (c1 == 0) c1 = ' ';
        if (c2 == 0) c2 = ' ';

        int32_t last_a0_b = 0, last_a0_x = 0, last_a0_y = 0, last_a0_z = 0;

        for (uint32_t r = 0; r < num_residues; r++) {
            const uint8_t* res_ptr = buffer + offset;
            int res_type = res_ptr[0];
            int resSeq = res_ptr[2] | (res_ptr[4] << 8);
            char iCode = res_ptr[6];
            int num_atoms = DICT[res_type].num_atoms;
            offset += 8;

            const uint8_t* a1_ptr = buffer + offset;
            accum_t accum;
            
            int32_t b_m = (int8_t)res_ptr[1]; int32_t x_m = (int8_t)res_ptr[3];
            int32_t y_m = (int8_t)res_ptr[5]; int32_t z_m = (int8_t)res_ptr[7];
            int32_t L0 = a1_ptr[0] | (a1_ptr[1] << 8); int32_t L1 = a1_ptr[2] | (a1_ptr[3] << 8);
            int32_t L2 = a1_ptr[4] | (a1_ptr[5] << 8); int32_t L3 = a1_ptr[6] | (a1_ptr[7] << 8);
            int32_t b_l = (L0 >> 3) & 0x1FFF; int32_t x_l = (L1 >> 3) & 0x1FFF;
            int32_t y_l = (L2 >> 3) & 0x1FFF; int32_t z_l = (L3 >> 3) & 0x1FFF;
            
#if defined(__AVX2__) || defined(__SSE4_1__)
            accum = _mm_set_epi32((z_m << 13) | z_l, (y_m << 13) | y_l, (x_m << 13) | x_l, (b_m << 13) | b_l);
#elif defined(__aarch64__) || defined(__ARM_NEON)
            int32x4_t setup = {(b_m << 13) | b_l, (x_m << 13) | x_l, (y_m << 13) | y_l, (z_m << 13) | z_l}; accum = setup;
#else
            accum.b = (b_m << 13) | b_l; accum.x = (x_m << 13) | x_l; accum.y = (y_m << 13) | y_l; accum.z = (z_m << 13) | z_l;
#endif
            last_a0_b = GET_SIMD_B(accum); last_a0_x = GET_SIMD_X(accum);
            last_a0_y = GET_SIMD_Y(accum); last_a0_z = GET_SIMD_Z(accum);

            uint8_t occ = (L0 & 0x07) | ((L1 & 0x07) << 3) | ((L2 & 0x01) << 6);
            if (occ > 0) {
                uint8_t alt_charge = ((L2 >> 1) & 0x03) | ((L3 & 0x07) << 2);
                int8_t charge = alt_charge & 0x0F;
                if (charge & 0x08) charge |= 0xF0;
                uint8_t alt_idx = (alt_charge >> 4) & 0x07;
                char altLoc = alt_idx == 0 ? ' ' : 'A' + alt_idx - 1;
                char chg_str[3] = "  ";
                if (charge > 0) { chg_str[0] = '0' + charge; chg_str[1] = '+'; } 
                else if (charge < 0) { chg_str[0] = '0' - charge; chg_str[1] = '-'; }
                
                int len = snprintf(tbuf, 128, "ATOM  %5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                        atom_serial++, DICT[res_type].atom_names[0], altLoc, DICT[res_type].resName,
                        c1, c2, resSeq, iCode == 0 ? ' ' : iCode,
                        last_a0_x / 1000.0f, last_a0_y / 1000.0f, last_a0_z / 1000.0f, occ / 100.0f, last_a0_b / 100.0f, 
                        DICT[res_type].elements[0], chg_str);
                memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; 
                w_ptr += 81;
            }
            offset += 8;

            for (int i = 1; i < num_atoms; i++) {
                const uint8_t* an_ptr = buffer + offset;
#if defined(__AVX2__) || defined(__SSE4_1__)
                __m128i an_vec = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)an_ptr));
                __m128i delta = _mm_srai_epi32(an_vec, 3);
                accum = _mm_add_epi32(accum, delta); 
#elif defined(__aarch64__) || defined(__ARM_NEON)
                int16x4_t an_vec_16 = vreinterpret_s16_u8(vld1_u8(an_ptr));
                int32x4_t an_vec_32 = vmovl_s16(an_vec_16);
                int32x4_t delta = vshrq_n_s32(an_vec_32, 3);
                accum = vaddq_s32(accum, delta);
#else
                int32_t aL0 = an_ptr[0] | (an_ptr[1] << 8); int32_t aL1 = an_ptr[2] | (an_ptr[3] << 8);
                int32_t aL2 = an_ptr[4] | (an_ptr[5] << 8); int32_t aL3 = an_ptr[6] | (an_ptr[7] << 8);
                accum.b += ((int16_t)aL0) >> 3; accum.x += ((int16_t)aL1) >> 3;
                accum.y += ((int16_t)aL2) >> 3; accum.z += ((int16_t)aL3) >> 3;
#endif
                
                uint16_t aL0 = an_ptr[0] | (an_ptr[1] << 8); uint16_t aL1 = an_ptr[2] | (an_ptr[3] << 8);
                uint16_t aL2 = an_ptr[4] | (an_ptr[5] << 8); uint16_t aL3 = an_ptr[6] | (an_ptr[7] << 8);
                occ = (aL0 & 0x07) | ((aL1 & 0x07) << 3) | ((aL2 & 0x01) << 6);
                
                if (occ > 0) {
                    uint8_t alt_charge = ((aL2 >> 1) & 0x03) | ((aL3 & 0x07) << 2);
                    int8_t charge = alt_charge & 0x0F;
                    if (charge & 0x08) charge |= 0xF0;
                    uint8_t alt_idx = (alt_charge >> 4) & 0x07;
                    char altLoc = alt_idx == 0 ? ' ' : 'A' + alt_idx - 1;
                    char chg_str[3] = "  ";
                    if (charge > 0) { chg_str[0] = '0' + charge; chg_str[1] = '+'; } 
                    else if (charge < 0) { chg_str[0] = '0' - charge; chg_str[1] = '-'; }
                    
                    int32_t v_b = GET_SIMD_B(accum); int32_t v_x = GET_SIMD_X(accum);
                    int32_t v_y = GET_SIMD_Y(accum); int32_t v_z = GET_SIMD_Z(accum);
                    
                    int len = snprintf(tbuf, 128, "ATOM  %5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                            atom_serial++, DICT[res_type].atom_names[i], altLoc, DICT[res_type].resName,
                            c1, c2, resSeq, iCode == 0 ? ' ' : iCode,
                            v_x / 1000.0f, v_y / 1000.0f, v_z / 1000.0f, occ / 100.0f, v_b / 100.0f, 
                            DICT[res_type].elements[i], chg_str);
                    memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; 
                    w_ptr += 81;
                }
                offset += 8;
            }
        }

        if (flags & CTRL_TER_OXT) {
            const uint8_t* oxt_ptr = buffer + offset;
            uint16_t mL0 = oxt_ptr[0] | (oxt_ptr[1] << 8);
            uint16_t mL1 = oxt_ptr[2] | (oxt_ptr[3] << 8);
            uint16_t mL2 = oxt_ptr[4] | (oxt_ptr[5] << 8);
            uint16_t mL3 = oxt_ptr[6] | (oxt_ptr[7] << 8);
            
            uint8_t occ = (mL0 & 0x07) | ((mL1 & 0x07) << 3) | ((mL2 & 0x01) << 6);
            if (occ > 0) {
                int32_t db = ((int16_t)mL0) >> 3; int32_t dx = ((int16_t)mL1) >> 3;
                int32_t dy = ((int16_t)mL2) >> 3; int32_t dz = ((int16_t)mL3) >> 3;
                
                int32_t oxt_b = last_a0_b + db; int32_t oxt_x = last_a0_x + dx;
                int32_t oxt_y = last_a0_y + dy; int32_t oxt_z = last_a0_z + dz;
                
                uint8_t alt_charge = ((mL2 >> 1) & 0x03) | ((mL3 & 0x07) << 2);
                int8_t charge = alt_charge & 0x0F;
                if (charge & 0x08) charge |= 0xF0;
                uint8_t alt_idx = (alt_charge >> 4) & 0x07;
                char altLoc = alt_idx == 0 ? ' ' : 'A' + alt_idx - 1;
                char chg_str[3] = "  ";
                if (charge > 0) { chg_str[0] = '0' + charge; chg_str[1] = '+'; } 
                else if (charge < 0) { chg_str[0] = '0' - charge; chg_str[1] = '-'; }
                
                int len = snprintf(tbuf, 128, "ATOM  %5d %-4s%c%3s%c%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
                        atom_serial++, " OXT", altLoc, DICT[0].resName, c1, c2, 0, ' ',
                        oxt_x / 1000.0f, oxt_y / 1000.0f, oxt_z / 1000.0f, occ / 100.0f, oxt_b / 100.0f, " O", chg_str);
                memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; 
                w_ptr += 81;
            }
            offset += 8;
        }

        if ((flags & CTRL_TER) || (flags & CTRL_TER_OXT)) {
            int len = snprintf(tbuf, 128, "TER   %5d      %3s%c%c%4d%c", 
                                atom_serial++, DICT[0].resName, c1, c2, 0, ' ');
            memcpy(w_ptr, tbuf, len); w_ptr[len] = ' '; 
            w_ptr += 81;
        }
    }
    
    fwrite(out_data, 1, out_size, out);
    free(out_data);
}

// -----------------------------------------------------------------------------
// 5. Test Runner & I/O
// -----------------------------------------------------------------------------
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
