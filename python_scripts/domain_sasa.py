"""
Calculate SASA and BSA using Bio.PDB ShrakeRupley.
Extracts chains from the complex to compute monomer SASA.
CAMSAP = chain A
Tubulin alpha = chains B, C
Tubulin beta  = chains D, E
BSA = (SASA_CAMSAP_isolated + SASA_TUBULIN_isolated) - SASA_COMPLEX
"""

import os # Operating system
import sys
import glob # Finds all CIF files in directory
import csv # Writes/reads csv files
from Bio.PDB import MMCIFParser, is_aa # CIF file parser
from Bio.PDB.SASA import ShrakeRupley # Calculates SASA
from collections import defaultdict
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from copy import deepcopy


# CONFIGURATION

CAMSAP_CHAIN = 'A'
TUBULIN_ALPHA = ['B', 'C']
TUBULIN_BETA  = ['D', 'E']
ALL_TUBULIN = TUBULIN_ALPHA + TUBULIN_BETA

# File name suffix

CIF_SUFFIX = '_tubulin_complex_model'

# DIRECTORIES

INPUT_DIR = "cifs"
OUTPUT_CSV = "sasa_bsa_results.csv"
DOMAIN_OUTPUT_CSV = "domain_sasa_summary.csv"
MAPPING_FILE = "domain_ranges.csv"

# BOUNDARIES
# Using UniProt annotations

boundaries = {
    "CAMSAP1" : {
    "n_terminus" :	(1,	215), # Start amino acid to start of CH domain
    "ch" : (216, 331), # CH domain
    "ch_cc1_linker" : (332, 872), # Linker region between CH and CC1
    "cc1" : (873, 909), # CC1 domain
    "cc1_cc2_linker" :	(910, 1015), # Linker region between CC1 and CC2
    "cc2" : (1016, 1048), # CC2 domain
    "cc2_cc3_linker" :	(1049, 1290), # Linker region between CC2 and CC3
    "cc3" :	(1291, 1343), # CC3 domain
    "cc" : (873, 1343), # Full coiled-coil region
    "cc3_ckk_linker" :	(1344, 1462), # Linker region between CC3 and CKK domain
    "ckk" :	(1463, 1597), # CKK domain
    "c_terminus" :	(1598, 1602), # End of CKK domain to end amino acid
    },
    "CAMSAP2" : {
    "n_terminus" :	(1,	221),
    "ch" : (222, 335),
    "ch_cc1_linker" : (336, 755),
    "cc1" : (756, 793),
    "cc1_cc2_linker" :	(794, 886),
    "cc2" : (887, 926),
    "cc2_cc3_linker" :	(927, 1165),
    "cc3" :	(1166, 1238),
    "cc" : (756, 1238),
    "cc3_ckk_linker" :	(1239, 1348),
    "ckk" :	(1349, 1483),
    "c_terminus" :	(1484, 1489),
    },
    "CAMSAP3" : {
    "n_terminus" :	(1,	202),
    "ch" : (203, 312),
    "ch_cc1_linker" : (313, 593),
    "cc1" : (594, 628),
    "cc1_cc2_linker" :	(629, 695),
    "cc2" : (696, 729),
    "cc2_cc3_linker" :	(730, 895),
    "cc3" :	(896, 936),
    "cc" : (594, 936),
    "cc3_ckk_linker" :	(937, 1108),
    "ckk" :	(1109, 1243),
    "c_terminus" :	(1244, 1249),
    }
}

# HELPER FUNCTIONS

def get_family(filename):
    """
    Get paralog name from filename
    """
    f = os.path.basename(filename).lower()
    if 'camsap1' in f:
        return 'CAMSAP1'
    elif 'camsap2' in f:
        return 'CAMSAP2'
    if 'camsap3' in f:
        return 'CAMSAP3'
    return None

def load_domain_mapping(mapping_file):
    """
    Load domain_ranges.csv into a dictionary of pattern -> list of (domain, start, end).
    """
    mapping = {}
    if not os.path.exists(mapping_file):
        print(f"WARNING: {mapping_file} not found. Domain SASA will be skipped.")
        return mapping

    with open(mapping_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            pattern = row[0].strip().lower()
            domain = row[1].strip()
            try:
                start = int(row[2])
                end = int(row[3])
            except ValueError:
                continue
            if pattern not in mapping:
                mapping[pattern] = []
            mapping[pattern].append((domain, start, end))

    print(f"Loaded {len(mapping)} construct mappings from {mapping_file}")
    return mapping

def get_domain_ranges(filename, mapping):
    """
    Get domain ranges from a cif file by an exact match against the construct pattern,
    with '.cif' and 'tubulin_complex_model' removed.
    """
    base = os.path.basename(filename)
    if base.lower().endswith('.cif'):
        base = base[:-4]
    base_lower = base.lower()

    pattern = base.lower()
    if pattern.endswith(CIF_SUFFIX):
        pattern = pattern[:-len(CIF_SUFFIX)]

    if pattern in mapping:
        return mapping[pattern]
    if base_lower in mapping:
        return mapping[base_lower]

    print(f"WARNING! No exact domain mapping found for pattern '{pattern}'"
          f"(from {base}). Domain SASA skipped for this file")

    return None

def compute_sasa_for_structure(structure):
    """
    Runs ShrakeRupley on the whole structure, and returns the total SASA summed
    across all residues.
    """
    sr = ShrakeRupley()
    sr.compute(structure, level='R')
    total = 0.0
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue) and hasattr(residue, 'sasa'):
                    total += residue.sasa
    return total

def get_per_residue_sasa(structure, chain_id):
    """
    Reads back per-residue SASA for one chain in a structure that was
    computed in the previous function.
    """
    sasa_per_res = {}
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if is_aa(residue) and hasattr(residue, 'sasa'):
                        sasa_per_res[residue.id[1]] = residue.sasa
    return sasa_per_res

def extract_chains_to_structure(structure, chain_ids):
    """
    Extract specified chains into a new Structure object.
    """
    # Create new structure

    new_struct = Structure.Structure('extracted')
    new_model = Model.Model(0)
    new_struct.add(new_model)

    for model in structure:
        for chain in model:
            if chain.id in chain_ids:
                new_chain = Chain.Chain(chain.id)
                for residue in chain:
                    if is_aa(residue):
                        new_chain.add(residue.copy())
                new_model.add(new_chain)

    return new_struct

def calculate_sasa_bsa(cif_file):
    """
    This function returns sasa_complex, sasa_camsap, sasa_tubulin, bsa,
    structure, struct_camsap.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('complex', cif_file)

    # SASA for complex (all chains)

    sasa_complex = compute_sasa_for_structure(structure)

    # SASA for isolated CAMSAP chain (chain A only)

    struct_camsap = extract_chains_to_structure(structure, [CAMSAP_CHAIN])
    sasa_camsap = compute_sasa_for_structure(struct_camsap)

    # SASA for isolated tubulin chains (chains B, C, D, E)

    struct_tubulin = extract_chains_to_structure(structure, ALL_TUBULIN)
    sasa_tubulin = compute_sasa_for_structure(struct_tubulin)

    # BSA

    bsa = (sasa_camsap + sasa_tubulin) - sasa_complex

    return sasa_complex, sasa_camsap, sasa_tubulin, bsa, structure, struct_camsap

def compute_domain_sasa_from_mapping(cif_file, structure, struct_camsap, mapping):
    """
    Compute per domain sasa and bsa using domain_ranges.csv. Works for any
    construct - full length or fragmented - due to domain_ranges.csv giving
    renumbered domain boundaries based on AlphaFold3 renumbering of residues

    'structure' refers to the full complex as a whole, SASA already computed
    in the context of being in a complex with tubulin

    'struct_camsap' refers to the isolated camsap chain A alone

    SASA is not recalculated in this function, it only reads back values
    computed in calculate_sasa_bsa
    """
    # Get domain ranges for this CIF

    ranges = get_domain_ranges(cif_file, mapping)
    if not ranges:
        return [], None

    # Build domain map: domain -> (start, end)

    domain_map = {domain: (start, end) for domain, start, end in ranges}

    # Get per residue SASA for CAMSAP in full complex
    # and isolated
    
    sasa_per_res_complex = get_per_residue_sasa(structure, CAMSAP_CHAIN)
    sasa_per_res_isolated = get_per_residue_sasa(struct_camsap, CAMSAP_CHAIN)

    # Determine family for output

    family = get_family(cif_file)
    if family is None:
        family = "Unknown"

    # Group SASA by domain for both states

    domain_sasa_complex = defaultdict(lambda: {'sum': 0.0, 'count': 0})
    domain_sasa_isolated = defaultdict(lambda: {'sum': 0.0, 'count': 0})

    # Helper to assign residue to domain

    def assign_domain(resnum):
        for domain, (start, end) in domain_map.items():
            if start <= resnum <= end:
                return domain
        return None

    for resnum, sasa in sasa_per_res_complex.items():
        domain = assign_domain(resnum)
        if domain is not None:
            domain_sasa_complex[domain]['sum'] += sasa
            domain_sasa_complex[domain]['count'] += 1

    for resnum, sasa in sasa_per_res_isolated.items():
        domain = assign_domain(resnum)
        if domain is not None:
            domain_sasa_isolated[domain]['sum'] += sasa
            domain_sasa_isolated[domain]['count'] += 1

    # Combine all domains

    all_domains = set(domain_sasa_complex.keys()) | set(domain_sasa_isolated.keys())

    domain_entries = []
    for domain in all_domains:
        sasa_comp = domain_sasa_complex.get(domain, {'sum': 0.0})['sum']
        sasa_isol = domain_sasa_isolated.get(domain, {'sum': 0.0})['sum']
        bsa_domain = sasa_isol - sasa_comp

        domain_entries.append({
            'domain': domain,
            'sasa_sum': sasa_comp,
            'bsa_domain': bsa_domain
        })

    return domain_entries, family

# MAIN

def main():
    if not os.path.isdir(INPUT_DIR):
        print(f"ERROR: Input folder '{INPUT_DIR}' not found.")
        sys.exit(1)

    cif_files = glob.glob(os.path.join(INPUT_DIR, "*.cif"))
    if not cif_files:
        print(f"No CIF files found in '{INPUT_DIR}'")
        return

    print(f"Found {len(cif_files)} CIF files.\n")

    mapping = load_domain_mapping(MAPPING_FILE)

    all_results = []
    domain_results = []

    for cif_path in cif_files:
        base = os.path.splitext(os.path.basename(cif_path))[0]
        print(f"Processing: {base}")

        try:
            sasa_complex, sasa_camsap, sasa_tubulin, bsa, structure, struct_camsap = calculate_sasa_bsa(cif_path)
            domain_entries, family = compute_domain_sasa_from_mapping(cif_path, structure, struct_camsap, mapping)
        except Exception as e:
            print(f"ERROR: {e}")
            continue

        print(f"  SASA Complex: {sasa_complex:.2f}")
        print(f"  SASA CAMSAP (isolated): {sasa_camsap:.2f}")
        print(f"  SASA Tubulin (isolated): {sasa_tubulin:.2f}")
        print(f"  BSA:          {bsa:.2f}")

        all_results.append({
            'file': base,
            'sasa_complex': sasa_complex,
            'sasa_camsap': sasa_camsap,
            'sasa_tubulin': sasa_tubulin,
            'bsa': bsa,
        })

        for entry in domain_entries:
            domain_results.append({
                'file':base,
                'family':family,
                'domain':entry['domain'],
                'sasa_domain':entry['sasa_sum'],
                'bsa_domain':entry['bsa_domain']
                })

    if all_results:
        with open(OUTPUT_CSV, 'w', newline='') as f:
            fieldnames = ['file', 'sasa_complex', 'sasa_camsap', 'sasa_tubulin', 'bsa']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)

        print(f"Whole protein summary saved to {OUTPUT_CSV}")

        if domain_results:
            with open(DOMAIN_OUTPUT_CSV, 'w', newline='') as f:
                fieldnames = ['file', 'family', 'domain', 'sasa_domain', 'bsa_domain']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(domain_results)

            print(f"FINISHED! Summary file saved to {DOMAIN_OUTPUT_CSV}")
        else:
            print("No domain data produced")

if __name__ == "__main__":
    main()