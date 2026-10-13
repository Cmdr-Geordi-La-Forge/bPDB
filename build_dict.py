import struct

sci_dict = {}
gb_common = {}  
reg_common = {} 
max_id = 0

def clean_common_name(name):
    """Strip adjectives frequently omitted by the PDB."""
    if name.startswith("DOMESTIC "):
        return name[9:]
    if name.startswith("WILD "):
        return name[5:]
    return name

with open('names.dmp', 'r') as f:
    for line in f:
        parts = [p.strip() for p in line.split('|')]
        if len(parts) < 4: continue
        
        tid = int(parts[0])
        name = parts[1].upper()
        name_class = parts[3]
        
        if tid > max_id: 
            max_id = tid
            
        if name_class == 'scientific name':
            sci_dict[tid] = name.encode('ascii', 'ignore') + b'\0'
        elif name_class == 'genbank common name':
            if tid not in gb_common: gb_common[tid] = []
            gb_common[tid].append(name)
        elif 'common name' in name_class:
            if tid not in reg_common: reg_common[tid] = []
            reg_common[tid].append(name)

com_dict = {}
for tid in range(max_id + 1):
    if tid in gb_common:
        best_name = min(gb_common[tid], key=len)
        best_name = clean_common_name(best_name)
        com_dict[tid] = best_name.encode('ascii', 'ignore') + b'\0'
        
    elif tid in reg_common:
        best_name = min(reg_common[tid], key=len)
        best_name = clean_common_name(best_name)
        com_dict[tid] = best_name.encode('ascii', 'ignore') + b'\0'

sci_offs, com_offs = [0] * (max_id + 1), [0] * (max_id + 1)
string_pool = bytearray(b'\0') 

for tid in range(max_id + 1):
    if tid in sci_dict:
        sci_offs[tid] = len(string_pool)
        string_pool += sci_dict[tid]
    if tid in com_dict:
        com_offs[tid] = len(string_pool)
        string_pool += com_dict[tid]

with open('taxdict.bin', 'wb') as f:
    f.write(struct.pack('<I', max_id))
    f.write(struct.pack(f'<{max_id+1}I', *sci_offs))
    f.write(struct.pack(f'<{max_id+1}I', *com_offs))
    f.write(string_pool)
    print(f"Built taxdict.bin (Max ID: {max_id}, Size: {f.tell() / 1024 / 1024:.2f} MB)")
