# GENERATE 48 ALPHAFOLD JSON FILES FROM CONSTRUCTS

# This script reads construct sequences from a FASTA file
# and creates 48 individual files for AlphaFold3
# Each JSON file will contain ->
# One CAMSAP construct (A)
# 2 alpha chains (B,C)
# 2 beta chains (D, E)
# Ligands (GTP, GDP, Mg2+)

# IMPORT LIBRARIES

from Bio import SeqIO # For reading FASTA files
import json # For writing JSON files
import os # For creating dictionaries
import re # Clean tubulin sequences

# TUBULIN SEQUENCES

# These are the human tubulin sequences used in all models (retrieved from UniProt)
# >sp|P68363|TBA1B_HUMAN Tubulin alpha-1B
# >sp|P07437|TBB5_HUMAN Tubulin beta chain
# Hardcoded here so chains are consistent for all JSON files

ALPHA_TUBULIN = re.sub(r"\s+", "",
"""MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGK
HVPRAVFVDLEPTVIDEVRTGTYRQLFHPEQLITGKEDAANNYARGHYTIGKEIIDLVLD
RIRKLADQCTGLQGFLVFHSFGGGTGSGFTSLLMERLSVDYGKKSKLEFSIYPAPQVSTA
VVEPYNSILTTHTTLEHSDCAFMVDNEAIYDICRRNLDIERPTYTNLNRLISQIVSSITA
SLRFDGALNVDLTEFQTNLVPYPRIHFPLATYAPVISAEKAYHEQLSVAEITNACFEPAN
QMVKCDPRHGKYMACCLLYRGDVVPKDVNAAIATIKTKRSIQFVDWCPTGFKVGINYQPP
TVVPGGDLAKVQRAVCMLSNTTAIAEAWARLDHKFDLMYAKRAFVHWYVGEGMEEGEFSE
AREDMAALEKDYEEVGVDSVEGEGEEEGEEY
""")
BETA_TUBULIN = re.sub(r"\s+", "",
"""MREIVHIQAGQCGNQIGAKFWEVISDEHGIDPTGTYHGDSDLQLDRISVYYNEATGGKYV
PRAILVDLEPGTMDSVRSGPFGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELVDSVLDVV
RKEAESCDCLQGFQLTHSLGGGTGSGMGTLLISKIREEYPDRIMNTFSVVPSPKVSDTVV
EPYNATLSVHQLVENTDETYCIDNEALYDICFRTLKLTTPTYGDLNHLVSATMSGVTTCL
RFPGQLNADLRKLAVNMVPFPRLHFFMPGFAPLTSRGSQQYRALTVPELTQQVFDAKNMM
AACDPRHGRYLTVAAVFRGRMSMKEVDEQMLNVQNKNSSYFVEWIPNNVKTAVCDIPPRG
LKMAVTFIGNSTAIQELFKRISEQFTAMFRRKAFLHWYTGEGMDEMEFTEAESNMNDLVS
EYQQYQDATAEEEEDFGEEAEEEA
""")

# LIGAND SETTINGS

# INCLUDE_LIGANDS can be set to False for ligands to be omitted from model

INCLUDE_LIGANDS = True

# Chain F and G -> GTP (2 copies)
# Chain H and I -> GDP (2 copies)
# Chain J, K, L, M -> Mg2+ (4 copies)

LIGANDS = [
    {"ligand":{"id":["F","G"], "ccdCodes": ["GTP"]}},
    {"ligand":{"id":["H","I"], "ccdCodes": ["GDP"]}},
    {"ligand":{"id":["J","K","L","M"], "ccdCodes": ["MG"]}},
]

# LOAD CONSTRUCTS FROM FASTA FILE

def load_constructs(fasta_file):
    """
    Load construct sequences from a fasta file.

    ARGS:
    fasta_file (str) -> path to a fasta file (e.g., "all_constructs.fasta").

    Returns -> dict -> keys are construct IDs and values are the sequences as strings.
    
    """

    constructs = {} # Store in dictionary

    # SeqIO.parse() reads a FASTA file and returns

    for record in SeqIO.parse(fasta_file, "fasta"):

        # record.id is the FASTA header (e.g., "CAMSAP1")
        # str(record.seq) converts the sequence to a string

        constructs[record.id] = str(record.seq)

    return constructs

def generate_json(protein, construct_type, sequence, output_dir="json_files"):
    """
    Generate an AlphaFold input file - JSON file - for a construct.

    ARGS:
    protein -> str -> protein name (e.g, "CAMSAP2")
    construct_type -> str -> construct type (e.g., "ch", "full", "delCH")
    sequence -> str -> the CAMSAP sequence for this construct (retrieved from UniProt)
    output_dir -> str -> directory to save JSON files ("json_files")

    RETURNS:
    A path to the created JSON file.
    """

    # Create the output directory

    os.makedirs(output_dir, exist_ok=True)

    # Create unique name for the constructs JSON file
    # Example -> ("CAMSAP1_full_tubulin_complex")

    name = f"{protein}_{construct_type}_tubulin_complex"

    # Build the sequences array
    # CAMSAP construct first
    # Followed by alpha and beta tubulin
    # And then the ligands (if using)

    sequences = [
        {
            "protein": {
                "id": "A",
                "sequence": sequence 
            }
        },
        {
            "protein": {
                "id": ["B","C"],
                "sequence": ALPHA_TUBULIN
            }
        },
        {
            "protein": {
                "id": ["D","E"],
                "sequence": BETA_TUBULIN
            }
        }
    ]

    if INCLUDE_LIGANDS:
        sequences.extend(LIGANDS)

    json_data = {
        "name" : name,
        "sequences" : sequences,
        "modelSeeds" : [1],
        "dialect" : "alphafold3",
        "version" : 1
    }

    # Writes the JSON to a file

    filename = f"{output_dir}/{protein}_{construct_type}.json"
    with open (filename, "w") as f:
        json.dump(json_data, f, indent=2) # indent 2 makes the file human readable

    return filename

def main():
    """
    Loads constructs from FASTA file and generates a JSON file for each one
    which is needed for AlphaFold3 submission.

    -> Reads all constructs from FASTA file
    -> For all constructs, splits the ID to get protein name and construct type
    -> Generates a JSON file for each construct

    """

    # Load constructs from FASTA file

    constructs = load_constructs("cc_constructs.fasta")

    # Loop through each construct
    # .items() returns key -> construct_id ("CAMSAP1_ch")
    # and value -> construct sequence

    for construct_id, sequence in constructs.items():
        parts = construct_id.split("_") # .split() splits the string at an underscore so ID is split into protein and construct type
        protein = parts[0] # Protein is first element
        construct_type = "_".join(parts[1:]) # Takes everything after first element (construct type)
        generate_json(protein, construct_type, sequence) # Creates JSON file

    print(f"Finshed! Generated {len(constructs)} JSON files!")

if __name__ == "__main__":
    main()