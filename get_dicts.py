import gzip
import os

FILE_PATH = "components.cif.gz"

AA_MAP = {
    "ALA": 1, "ARG": 2, "ASN": 3, "ASP": 4, "CYS": 5, "GLN": 6, "GLU": 7, "GLY": 8, "HIS": 9, "ILE": 10,
    "LEU": 11, "LYS": 12, "MET": 13, "PHE": 14, "PRO": 15, "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20,
    "SEC": 21, "PYL": 22, "ASX": 23, "GLX": 24, "UNK": 25
}
NA_MAP = {
    "A": 1, "C": 2, "G": 3, "U": 4, "I": 5,
    "DA": 6, "DC": 7, "DG": 8, "DT": 9, "DU": 10, "DI": 11, "N": 12
}
SOLV_MAP = {"HOH": 1}

TARGETS = set(list(AA_MAP.keys()) + list(NA_MAP.keys()) + list(SOLV_MAP.keys()))

NOMENCLATURE_MAP = {"H5'1": "H5'", "H5'2": "H5''", "H2'1": "H2'", "H2'2": "H2''"}

MANUAL_EXTRAS = {
    "PYL": ([" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ "," C2 "," C3 "," C4 "," C5 "," C6 "," O2 "," CA2"," CB2"," CG2"," CD2"," CE2"," CE3"," O22"," H  "," HA "],
            [" N"," C"," C"," O"," C"," C"," C"," C"," N"," C"," C"," C"," C"," C"," O"," C"," C"," C"," C"," C"," C"," O"," H"," H"]),
    "ASX": ([" N  "," CA "," C  "," O  "," CB "," CG "," XD1"," XD2"," H  "," HA "," HB2"," HB3"],
            [" N"," C"," C"," O"," C"," C"," X"," X"," H"," H"," H"," H"]),
    "GLX": ([" N  "," CA "," C  "," O  "," CB "," CG "," CD "," XE1"," XE2"," H  "," HA "," HB2"," HB3"," HG2"," HG3"],
            [" N"," C"," C"," O"," C"," C"," C"," X"," X"," H"," H"," H"," H"," H"," H"]),
    "UNK": ([" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," CZ "," CH "," X1 "," X2 "," X3 "," X4 "," H  "," HA "],
            [" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," X"," X"," X"," X"," H"," H"]),
    "N":   ([" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N6 "," N1 "," C2 "," N2 "," N3 "," C4 "," H5'","H5''"," H4'"," H3'"," H2'","HO2'"],
            [" P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," N"," C"," N"," N"," C"," H"," H"," H"," H"," H"," H"]),
    "HOH": ([" O  ", " H1 ", " H2 "], [" O", " H", " H"])
}

def format_atom_name(name, element):
    name = name.strip('"')
    name = NOMENCLATURE_MAP.get(name, name)
    element = element.strip('"')

    if len(element) == 2:
        # Two-letter elements (Se, Fe, Ca) start in column 13 (no leading space)
        return f"{name:<4}"
    else:
        # One-letter elements start in column 14, unless 4-chars long (like HXT)
        return f" {name:<3}" if len(name) < 4 else f"{name:<4}"

def partition_atoms(atoms, elements, res_type):
    heavy = []
    hydro = []

    if res_type == "AA":
        n_term_heavy = set()
        c_term_heavy = {"OXT"}
        n_term_hydro = {"H1", "H2", "H3"}
        c_term_hydro = {"HXT"}
    elif res_type == "NA":
        n_term_heavy = {"OP3"}
        c_term_heavy = {"O3T"}
        n_term_hydro = {"HOP3", "HOP2", "HO5'"}
        c_term_hydro = {"HO3'"}
    else:
        n_term_heavy, c_term_heavy, n_term_hydro, c_term_hydro = set(), set(), set(), set()

    for a, e in zip(atoms, elements):
        astrip = a.strip()
        if "H" in e or astrip in n_term_hydro or astrip in c_term_hydro:
            hydro.append((a, e))
        else:
            heavy.append((a, e))

    def sort_group(grp, n_set, c_set):
        n_part = [x for x in grp if x[0].strip() in n_set]
        c_part = [x for x in grp if x[0].strip() in c_set]
        core_part = [x for x in grp if x[0].strip() not in n_set and x[0].strip() not in c_set]

        if "HO5'" in n_set:
            n_part = [x for x in n_part if x[0].strip() != "HO5'"] + [x for x in n_part if x[0].strip() == "HO5'"]

        return n_part, core_part, c_part

    h_n, h_core, h_c = sort_group(heavy, n_term_heavy, c_term_heavy)
    hy_n, hy_core, hy_c = sort_group(hydro, n_term_hydro, c_term_hydro)

    heavy_sorted = h_n + h_core + h_c
    hydro_sorted = hy_n + hy_core + hy_c

    return (heavy_sorted, len(h_n), len(h_n)+len(h_core), len(heavy_sorted),
            hydro_sorted, len(hy_n), len(hy_n)+len(hy_core), len(hydro_sorted))

def generate_dict():
    residues = {t: {"atoms": [], "elements": []} for t in TARGETS}

    if os.path.exists(FILE_PATH):
        with gzip.open(FILE_PATH, 'rt') as f:
            in_atom_loop = False; headers = []
            for line in f:
                line = line.strip()
                if not line: continue
                if line == "loop_":
                    in_atom_loop = False; headers = []; continue
                if line.startswith("_chem_comp_atom."):
                    in_atom_loop = True; headers.append(line); continue

                if in_atom_loop and not line.startswith("_"):
                    tokens = line.split()
                    if len(tokens) >= len(headers):
                        try:
                            c_id = tokens[headers.index("_chem_comp_atom.comp_id")].strip('"')
                            if c_id in TARGETS:
                                raw_atom = tokens[headers.index("_chem_comp_atom.atom_id")]
                                raw_elem = tokens[headers.index("_chem_comp_atom.type_symbol")]

                                a_id = format_atom_name(raw_atom, raw_elem)

                                e_sym = raw_elem.strip('"').upper()
                                e_sym = f" {e_sym}" if len(e_sym) == 1 else e_sym[:2]

                                residues[c_id]["atoms"].append(a_id)
                                residues[c_id]["elements"].append(e_sym)
                        except ValueError:
                            in_atom_loop = False
                elif in_atom_loop and line.startswith("_"):
                    in_atom_loop = False

    for k, v in MANUAL_EXTRAS.items():
        if k not in residues or not residues[k]["atoms"]:
            residues[k] = {"atoms": v[0], "elements": v[1]}

    # --- INJECT POLYMER TERMINAL ATOMS ---
    for aa in AA_MAP:
        if aa in residues:
            atoms = residues[aa]["atoms"]
            elements = residues[aa]["elements"]

            if " H2 " in atoms:
                h2_idx = atoms.index(" H2 ")
                if " H1 " not in atoms:
                    atoms.insert(h2_idx, " H1 ")
                    elements.insert(h2_idx, " H")

                h2_idx = atoms.index(" H2 ")
                if " H3 " not in atoms:
                    atoms.insert(h2_idx + 1, " H3 ")
                    elements.insert(h2_idx + 1, " H")
            else:
                for atom, elem in zip([" H1 ", " H2 ", " H3 "], [" H", " H", " H"]):
                    if atom not in atoms:
                        atoms.append(atom)
                        elements.append(elem)

            for atom, elem in zip([" HXT", " OXT"], [" H", " O"]):
                if atom not in atoms:
                    atoms.append(atom)
                    elements.append(elem)

    for na in NA_MAP:
        if na in residues:
            atoms = residues[na]["atoms"]
            elements = residues[na]["elements"]

            for atom, elem in zip([" OP3", " O3T"], [" O", " O"]):
                if atom not in atoms:
                    atoms.append(atom)
                    elements.append(elem)

            for atom, elem in zip(["HOP3", "HOP2", "HO5'", "HO3'"], [" H", " H", " H", " H"]):
                if atom not in atoms:
                    atoms.append(atom)
                    elements.append(elem)

    def add_residue(base_id, name, atoms, elements, res_type):
        h_sorted, h_sn, h_en, h_ec, hy_sorted, hy_sn, hy_en, hy_ec = partition_atoms(atoms, elements, res_type)

        fmt_name = f'"{name:>3}"' if name in NA_MAP else f'"{name:<3}"'

        h_a_str = '","'.join([x[0] for x in h_sorted])
        h_e_str = '","'.join([x[1] for x in h_sorted])
        hy_a_str = '","'.join([x[0] for x in hy_sorted])
        hy_e_str = '","'.join([x[1] for x in hy_sorted])

        h_a_fmt = f'{{"{h_a_str}"}}' if h_a_str else "{}"
        h_e_fmt = f'{{"{h_e_str}"}}' if h_e_str else "{}"
        hy_a_fmt = f'{{"{hy_a_str}"}}' if hy_a_str else "{}"
        hy_e_fmt = f'{{"{hy_e_str}"}}' if hy_e_str else "{}"

        return f'    {{{base_id:#04x}, {fmt_name}, {h_sn}, {h_en}, {h_ec}, {h_a_fmt}, {h_e_fmt}, {hy_sn}, {hy_en}, {hy_ec}, {hy_a_fmt}, {hy_e_fmt}}},'

    lines = ["    // 0x01 - 0x1F: Nucleic Acids"]
    for res, idx in sorted(NA_MAP.items(), key=lambda x: x[1]):
        lines.append(add_residue(0x00 + idx, res, residues[res]["atoms"], residues[res]["elements"], "NA"))

    lines.append("\n    // 0x20 - 0x3F: Amino Acids")
    for res, idx in sorted(AA_MAP.items(), key=lambda x: x[1]):
        lines.append(add_residue(0x20 + idx, res, residues[res]["atoms"], residues[res]["elements"], "AA"))

    lines.append("\n    // 0x60 - 0x7F: Solvents")
    for res, idx in sorted(SOLV_MAP.items(), key=lambda x: x[1]):
        lines.append(add_residue(0x60 + idx, res, residues[res]["atoms"], residues[res]["elements"], "SOLV"))

    out = '#include "bpdb.h"\n\n'
    out += 'const ResidueDict DICT[] = {\n'
    out += "\n".join(filter(None, lines))
    out += '\n};\n\n'
    out += 'const size_t NUM_DICT_RES = sizeof(DICT) / sizeof(DICT[0]);\n'
    out += 'const char* MONTHS[] = {"", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};\n'

    with open("bpdb_dict.c", "w") as f:
        f.write(out)

if __name__ == "__main__":
    generate_dict()
