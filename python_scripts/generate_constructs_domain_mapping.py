# This script generates all CAMSAP construct FASTA file and a matching domain_ranges.csv
# that maps each constructs residue position (starting from 1) to domain names
# This will be used when calculating contacts/sasa per domain in a full-length construct

# IMPORT PYTHON LIBRARIES

from Bio import SeqIO # Reading/writing FASTA files
from Bio.SeqRecord import SeqRecord # Creates sequence records
from Bio.Seq import Seq # Handles biological sequences
import re # Clean sequence
import csv # Write a csv file

# DEFINE PROTEIN SEQUENCES

# This dictionary stores full-length amino acid protein sequences for each CAMSAP.
# The keys are protein names (e.g., "CAMSAP3") and the values are actual protein sequence (strings of amino acids).

sequences = {
    "CAMSAP1" : re.sub(r"\s+", "",
"""MVDASGRAAAEGWRKMEAPPDGAADLVPLDRYDAARAKIAANLQWICAKAYGRDNIPEDL
RDPFYVDQYEQEHIKPPVIKLLLSSELYCRVCSLILKGDQVAALQGHQSVIQALSRKGIY
VMESDDTPVTESDLSRAPIKMSAHMAMVDALMMAYTVEMISIEKVVASVKRFSTFSASKE
LPYDLEDAMVFWINKVNLKMREITEKEVKLKQQLLESPAHQKVRYRREHLSARQSPYFPL
LEDLMRDGSDGAALLAVIHYYCPEQMKLDDICLKEVTSMADSLYNIRLLREFSNEYLNKC
FYLTLEDMLYAPLVLKPNVMVFIAELFWWFENVKPDFVQPRDVQELKDAKTVLHQKSSRP
PVPISNATKRSFLGSPAAGTLAELQPPVQLPAEGCHRHYLHPEEPEYLGKGTAAFSPSHP
LLPLRQKQQKSIQGEDIPDQRHRSNSLTRVDGQPRGAAIAWPEKKTRPASQPTPFALHHA
ASCEVDPSSGDSISLARSISKDSLASNIVNLTPQNQPHPTATKSHGKSLLSNVSIEDEEE
ELVAIVRADVVPQQADPEFPRASPRALGLTANARSPQGQLDTSESKPDSFFLEPLMPAVL
KPAKEKQVITKEDERGEGRPRSIVSRRPSEGPQPLVRRKMTGSRDLNRTFTPIPCSEFPM
GIDPTETGPLSVETAGEVCGGPLALGGFDPFPQGPSTDGFFLHVGRADEDTEGRLYVSCS
KSPNSHDSEPWTLLRQDSDSDVVDIEEAEHDFMGEAHPVVFSRYIGEEESAKLQEDMKVK
EHEDKDDASGRSSPCLSTASQMSSVSMASGSVKMTSFAERKLQRLNSCETKSSTSSSQKT
TPDASESCPAPLTTWRQKREQSPSQHGKDPASLLASELVQLHMQLEEKRRAIEAQKKKME
ALSARQRLKLGKAAFLHVVKKGKAEAAPPLRPEHFAKEYSQHNGEDCGDAVSKTEDFLVK
EEQREELLHEPQDVDKESLAFAQQHKAKDPVALHELERNKVISAALLEDTVGEVVDVNEC
DLSIEKLNETISTLQQAILKISQQQEQLLMKSPTVPVPGSKNNSQDHKVKAPVHFVEPLS
PTGVAGHRKAPRLGQGRNSRSGRPAELKVPKDRPQGSSRSKTPTPSVETLPHLRPFPASS
HPRTPTDPGLDSALEPSGDPHGKCLFDSYRLHDESNQRTLTLSSSKDANILSEQMSLKEV
LDASVKEVGSSSSDVSGKESVPVEEPLRSRASLIEVDLSDLKAPDEDGELVSLDGSADLV
SEGDQKPGVGFFFKDEQKAEDELAKKRAAFLLKQQRKAEEARVRKQQLEAEVELKRDEAR
RKAEEDRVRKEEEKARRELIKQEYLRRKQQQILEEQGLGKPKSKPKKPRPKSVHREESCS
DSGTKCSSTPDNLSRTQSGSSLSLASAATTEPESVHSGGTPSQRVESMEALPILSRNPSR
STDRDWETASAASSLASVAEYTGPKLFKEPSSKSNKPIIHNAISHCCLAGKVNEPHKNSI
LEELEKCDANHYIILFRDAGCQFRALYCYYPDTEEIYKLTGTGPKNITKKMIDKLYKYSS
DRKQFNLIPAKTMSVSVDALTIHNHLWQPKRPAVPKKAQTRK"""),
    "CAMSAP2" : re.sub(r"\s+", "",
"""MGDAADPREMRKTFIVPAIKPFDHYDFSRAKIACNLAWLVAKAFGTENVPEELQEPFYTD
QYDQEHIKPPVVNLLLSAELYCRAGSLILKSDAAKPLLGHDAVIQALAQKGLYVTDQEKL
VTERDLHKKPIQMSAHLAMIDTLMMAYTVEMVSIEKVIACAQQYSAFFQATDLPYDIEDA
VMYWINKVNEHLKDIMEQEQKLKEHHTVEAPGGQKSPSKWFWKLVPARYRKEQTLLKQLP
CIPLVENLLKDGTDGCALAALIHFYCPDVVRLEDICLKETMSLADSLYNLQLIQEFCQEY
LNQCCHFTLEDMLYAASSIKSNYLVFMAELFWWFEVVKPSFVQPRVVRPQGAEPVKDMPS
IPVLNAAKRNVLDSSSDFPSSGEGATFTQSHHHLPSRYSRPQAHSSASGGIRRSSSMSYV
DGFIGTWPKEKRSSVHGVSFDISFDKEDSVQRSTPNRGITRSISNEGLTLNNSHVSKHIR
KNLSFKPINGEEEAESIEEELNIDSHSDLKSCVPLNTNELNSNENIHYKLPNGALQNRIL
LDEFGNQIETPSIEEALQIIHDTEKSPHTPQPDQIANGFFLHSQEMSILNSNIKLNQSSP
DNVTDTKGALSPITDNTEVDTGIHVPSEDIPETMDEDSSLRDYTVSLDSDMDDASKFLQD
YDIRTGNTREALSPCPSTVSTKSQPGSSASSSSGVKMTSFAEQKFRKLNHTDGKSSGSSS
QKTTPEGSELNIPHVVAWAQIPEETGLPQGRDTTQLLASEMVHLRMKLEEKRRAIEAQKK
KMEAAFTKQRQKMGRTAFLTVVKKKGDGISPLREEAAGAEDEKVYTDRAKEKESQKTDGQ
RSKSLADIKESMENPQAKWLKSPTTPIDPEKQWNLASPSEETLNEGEILEYTKSIEKLNS
SLHFLQQEMQRLSLQQEMLMQMREQQSWVISPPQPSPQKQIRDFKPSKQAGLSSAIAPFS
SDSPRPTHPSPQSSNRKSASFSVKSQRTPRPNELKITPLNRTLTPPRSVDSLPRLRRFSP
SQVPIQTRSFVCFGDDGEPQLKESKPKEEVKKEELESKGTLEQRGHNPEEKEIKPFESTV
SEVLSLPVTETVCLTPNEDQLNQPTEPPPKPVFPPTAPKNVNLIEVSLSDLKPPEKADVP
VEKYDGESDKEQFDDDQKVCCGFFFKDDQKAENDMAMKRAALLEKRLRREKETQLRKQQL
EAEMEHKKEETRRKTEEERQKKEDERARREFIRQEYMRRKQLKLMEDMDTVIKPRPQVVK
QKKQRPKSIHRDHIESPKTPIKGPPVSSLSLASLNTGDNESVHSGKRTPRSESVEGFLSP
SRCGSRNGEKDWENASTTSSVASGTEYTGPKLYKEPSAKSNKHIIQNALAHCCLAGKVNE
GQKKKILEEMEKSDANNFLILFRDSGCQFRSLYTYCPETEEINKLTGIGPKSITKKMIEG
LYKYNSDRKQFSHIPAKTLSASVDAITIHSHLWQTKRPVTPKKLLPTKA"""),
    "CAMSAP3" : re.sub(r"\s+", "",
"""MVEAAPPGPGPLRRTFLVPEIKSLDQYDFSRAKAAASLAWVLRAAFGGAEHVPPELWEPF
YTDQYAQEHVKPPVTRLLLSAELYCRAWRQALPQLETPPNPSALLALLARRGTVPALPER
PVREADLRHQPILMGAHLAVIDALMAAFAFEWTKTLPGPLALTSLEHKLLFWVDTTVRRL
QEKTEQEAAQRASPAAPADGAAPAQPSIRYRKDRVVARRAPCFPTVTSLQDLASGAALAA
TIHCYCPQLLRLEEVCLKDPMSVADSLYNLQLVQDFCASRLPRGCPLSLEDLLYVPPPLK
VNLVVMLAELFMCFEVLKPDFVQVKDLPDGHAASPRGTEASPPQNNSGSSSPVFTFRHPL
LSSGGPQSPLRGSTGSLKSSPSMSHMEALGKAWNRQLSRPLSQAVSFSTPFGLDSDVDVV
MGDPVLLRSVSSDSLGPPRPAPARTPTQPPPEPGDLPTIEEALQIIHSAEPRLLPDGAAD
GSFYLHSPEGPSKPSLASPYLPEGTSKPLSDRPTKAPVYMPHPETPSKPSPCLVGEASKP
PAPSEGSPKAVASSPAATNSEVKMTSFAERKKQLVKAEAEAGAGSPTSTPAPPEALSSEM
SELSARLEEKRRAIEAQKRRIEAIFAKHRQRLGKSAFLQVQPREASGEAEAEAEEADSGP
VPGGERPAGEGQGEPTSRPKAVTFSPDLGPVPHEGLGEYNRAVSKLSAALSSLQRDMQRL
TDQQQRLLAPPEAPGSAPPPAAWVIPGPTTGPKAASPSPARRVPATRRSPGPGPSQSPRS
PKHTRPAELRLAPLTRVLTPPHDVDSLPHLRKFSPSQVPVQTRSSILLAEETPPEEPAAR
PGLIEIPLGSLADPAAEDEGDGSPAGAEDSLEEEASSEGEPRVGLGFFYKDEDKPEDEMA
QKRASLLERQQRRAEEARRRKQWQEVEKEQRREEAARLAQEEAPGPAPLVSAVPMATPAP
AARAPAEEEVGPRKGDFTRQEYERRAQLKLMDDLDKVLRPRAAGSGGPGRGGRRATRPRS
GCCDDSALARSPARGLLGSRLSKIYSQSTLSLSTVANEAHNNLGVKRPTSRAPSPSGLMS
PSRLPGSRERDWENGSNASSPASVPEYTGPRLYKEPSAKSNKFIIHNALSHCCLAGKVNE
PQKNRILEEIEKSKANHFLILFRDSSCQFRALYTLSGETEELSRLAGYGPRTVTPAMVEG
IYKYNSDRKRFTQIPAKTMSMSVDAFTIQGHLWQGKKPTTPKKGGGTPK""")
}

LINKER = 'GGGGSGGGGS'

# DEFINE DOMAIN BOUNDRIES

# This dictionary stores the start and end boundries of regions for each CAMSAP.
# The numbers are 1-based. Each key is a protein name (e.g., "CAMSAP2"), and the value
# is another dictionary. The inner dictionary maps region names (e.g., "CKK") to a tuple
# of start, end positions.

boundries = {
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

def slice_region(seq, start, end):
    """
    Extracts a region from a protein using 1-based boundries.

    PARAMETERS:

    seq(str) -> full protein sequence
    start(int) -> start residue
    end(int) -> end residue

    RETURNS:

    str -> extracted region sequence
    """

    # Covert 1-based (UniProt) to 0-based indexing (python)

    return seq[start-1:end]

def get_length(b, domain):
    """
    Returns length of a domain (end - start + 1)
    """
    if domain not in b:
        return 0
    s, e = b[domain]
    return e - s + 1

def get_domain_order(protein, ctype):
    """
    Returns a list of (domain_name, length) for the given construct type.
    The order matches the order in generate_constructs.
    """
    b = boundries[protein]  # your "boundries" dictionary

    # Full-length
    if ctype == "full":
        return [
            ("n_terminus", get_length(b, "n_terminus")),
            ("ch", get_length(b, "ch")),
            ("ch_cc1_linker", get_length(b, "ch_cc1_linker")),
            ("cc1", get_length(b, "cc1")),
            ("cc1_cc2_linker", get_length(b, "cc1_cc2_linker")),
            ("cc2", get_length(b, "cc2")),
            ("cc2_cc3_linker", get_length(b, "cc2_cc3_linker")),
            ("cc3", get_length(b, "cc3")),
            ("cc3_ckk_linker", get_length(b, "cc3_ckk_linker")),
            ("ckk", get_length(b, "ckk")),
            ("c_terminus", get_length(b, "c_terminus")),
        ]

    # Individual domains and linkers

    if ctype in ["ch", "cc", "ckk", "cc1", "cc2", "cc3", "n_terminus",
                 "ch_cc1_linker", "cc1_cc2_linker", "cc2_cc3_linker",
                 "cc3_ckk_linker", "c_terminus"]:
        if ctype in b:
            return [(ctype, get_length(b, ctype))]

    # Fused constructs
    if ctype == "ch_ckk":
        ch_len = get_length(b, "ch")
        ckk_len = get_length(b, "ckk")
        linker_len = len(LINKER)  # "GGGGSGGGGS"
        return [("ch", ch_len), ("linker", linker_len), ("ckk", ckk_len)]

    if ctype == "cc_ckk":
        cc_len = get_length(b, "cc")
        linker_len = get_length(b, "cc3_ckk_linker")
        ckk_len = get_length(b, "ckk")
        return [("cc", cc_len), ("cc3_ckk_linker", linker_len), ("ckk", ckk_len)]

    if ctype == "ch_cc":
        ch_len = get_length(b, "ch")
        linker_len = get_length(b, "ch_cc1_linker")
        cc_len = get_length(b, "cc")
        return [("ch", ch_len), ("ch_cc1_linker", linker_len), ("cc", cc_len)]

    # Note to add: these following constructs will only be used in the final analysis if
    # the server decides to work and I can actually model them

    if ctype == "cc3_d2_ckk":
        cc3_len = get_length(b, "cc3")
        linker_len = get_length(b, "cc3_ckk_linker")
        ckk_len = get_length(b, "ckk")
        return [("cc3", cc3_len), ("cc3_ckk_linker", linker_len), ("ckk", ckk_len)]

    if ctype == "cc2_mbd_cc3_d2_ckk":
        cc2_len = get_length(b, "cc2")
        mbd_len = get_length(b, "cc2_cc3_linker")
        cc3_len = get_length(b, "cc3")
        d2_len = get_length(b, "cc3_ckk_linker")
        ckk_len = get_length(b, "ckk")
        return [("cc2", cc2_len), ("cc2_cc3_linker", mbd_len), ("cc3", cc3_len), ("cc3_ckk_linker", d2_len), ("ckk", ckk_len)]

    if ctype == "cc3_flex_ckk":
        cc3_len = get_length(b, "cc3")
        linker_len = len(LINKER)
        ckk_len = get_length(b, "ckk")
        return [("cc3", cc3_len), ("linker", linker_len), ("ckk", ckk_len)]

    if ctype == "cc2_flex_cc3_flex_ckk":
        cc2_len = get_length(b, "cc2")
        linker_len = len(LINKER)
        cc3_len = get_length(b, "cc3")
        ckk_len = get_length(b, "ckk")
        return [("cc2", cc2_len), ("linker", linker_len), ("cc3", cc3_len), ("linker", linker_len), ("ckk", ckk_len)]

    if ctype == "cc2_flex_cc3_d2_ckk":
        cc2_len = get_length(b, "cc2")
        linker_len = len(LINKER)
        cc3_len = get_length(b, "cc3")
        d2_len = get_length(b, "cc3_ckk_linker")
        ckk_len = get_length(b, "ckk")
        return [("cc2", cc2_len), ("linker", linker_len), ("cc3", cc3_len), ("cc3_ckk_linker", d2_len), ("ckk", ckk_len)]

    if ctype == "cc2_mbd_cc3_flex_ckk":
        cc2_len = get_length(b, "cc2")
        mbd_len = get_length(b, "cc2_cc3_linker")
        cc3_len = get_length(b, "cc3")
        linker_len = len(LINKER)
        ckk_len = get_length(b, "ckk")
        return [("cc2", cc2_len), ("cc2_cc3_linker", mbd_len), ("cc3", cc3_len), ("linker", linker_len), ("ckk", ckk_len)]

    if ctype == "cc2_mbd_cc3_d2":
        cc2_len = get_length(b, "cc2")
        mbd_len = get_length(b, "cc2_cc3_linker")
        cc3_len = get_length(b, "cc3")
        d2_len = get_length(b, "cc3_ckk_linker")
        return [("cc2", cc2_len), ("cc2_cc3_linker", mbd_len), ("cc3", cc3_len), ("cc3_ckk_linker", d2_len)]

    # Deletions

    if ctype == "delCH":
        return [
            ("n_terminus", get_length(b, "n_terminus")),
            ("ch_cc1_linker", get_length(b, "ch_cc1_linker")),
            ("cc1", get_length(b, "cc1")),
            ("cc1_cc2_linker", get_length(b, "cc1_cc2_linker")),
            ("cc2", get_length(b, "cc2")),
            ("cc2_cc3_linker", get_length(b, "cc2_cc3_linker")),
            ("cc3", get_length(b, "cc3")),
            ("cc3_ckk_linker", get_length(b, "cc3_ckk_linker")),
            ("ckk", get_length(b, "ckk")),
            ("c_terminus", get_length(b, "c_terminus")),
        ]

    if ctype == "delCC":
        return [
            ("n_terminus", get_length(b, "n_terminus")),
            ("ch", get_length(b, "ch")),
            ("ch_cc1_linker", get_length(b, "ch_cc1_linker")),
            ("cc3_ckk_linker", get_length(b, "cc3_ckk_linker")),
            ("ckk", get_length(b, "ckk")),
            ("c_terminus", get_length(b, "c_terminus")),
        ]

    if ctype == "delCKK":
        return [
            ("n_terminus", get_length(b, "n_terminus")),
            ("ch", get_length(b, "ch")),
            ("ch_cc1_linker", get_length(b, "ch_cc1_linker")),
            ("cc1", get_length(b, "cc1")),
            ("cc1_cc2_linker", get_length(b, "cc1_cc2_linker")),
            ("cc2", get_length(b, "cc2")),
            ("cc2_cc3_linker", get_length(b, "cc2_cc3_linker")),
            ("cc3", get_length(b, "cc3")),
            ("cc3_ckk_linker", get_length(b, "cc3_ckk_linker")),
            ("c_terminus", get_length(b, "c_terminus")),
        ]

    # If we get here, construct type unknown

    print(f"WARNING: Unknown construct type '{ctype}' for {protein}")
    return []

def generate_domain_ranges(constructs, output_file="domain_ranges.csv"):
    """
    Generate domain_ranges.csv directly from the constructs list.
    Each construct's domains are renumbered from 1 to match how 
    AlphaFold3 will number that chain in the output CIF
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['# Domain ranges from generate_constructs.py'])
        writer.writerow(['# Format: pattern,domain,start,end'])

        for rec in constructs:

            # Get the construct name with full description

            name = rec.id
            seq = str(rec.seq)

            # The construct name tells us which protein and which type
            # e.g., "CAMSAP1_full" -> protein = CAMSAP1, ctype = full

            parts = name.split('_')
            protein = parts[0]
            ctype = '_'.join(parts[1:])

            # Get the domain order for this construct type

            domain_lengths = get_domain_order(protein, ctype)
            if not domain_lengths:
                continue

            # Assign ranges (renumbered from 1)

            cum = 1
            for domain, length in domain_lengths:
                if length == 0:
                    continue
                end = cum + length - 1

                # Use lowercase pattern to match CIF filenames

                pattern = name.lower()
                writer.writerow([pattern, domain, cum, end])
                cum += length

    print(f"Domain ranges saved to {output_file}")

def generate_constructs(sequences, boundries, linker = "GGGGSGGGGS"):
    """
    Generate all AlphaFold constructs for each CAMSAP

    This functions creates:
    1. Full-length proteins
    2. Individual domains (CH, CC, CKK, CC1, CC2, CC3)
    3. Linker regions (e.g., CH-CC1 linker)
    4. Fused constructs with linker (CH+CKK with linker, CC+CKK, CH+CC)
    5. Deletion constructs (-CH, -CC, CKK)

    PARAMETERS:

    sequences(dict) -> protein name & full sequence
    boundries(dict) -> protein name & region (start, end)
    linker(str) -> flexible linker region for CH+CKK construct (GGGGSGGGGS)

    RETURNS:

    list -> BioPython SeqRecord objects containing sequence & ID & description
    """

    # List for all SeqRecord objects

    constructs = []

    # Loop through each protein

    for protein, seq in sequences.items():

        # boundries for a specific protein

        b = boundries[protein]

        # FULL-LENGTH CONSTRUCTS (3)

        constructs.append(
            SeqRecord(
                Seq(seq), # Coverts str to BioPython Seq object
                id=f"{protein}_full", # Creates unique identifier (e.g., CAMSAP1_full)
                description=f"{protein} full-length" # Creates description
            )
        )

        # INDIVIDUAL DOMAINS (18)
        # Extract individual domains for each CAMSAP

        for region in ["ch", "cc", "ckk", "cc1", "cc2", "cc3"]:
            if region in b: # Checks if this region exists
                start, end = b[region]
                region_seq = slice_region(seq, start, end)
                constructs.append(
                    SeqRecord(
                        Seq(region_seq),
                        id=f"{protein}_{region}",
                        description=f"{protein} {region}"
                    )
                )

        # LINKER REGIONS (18)
        # Extract all linker regions individually

        for region in ["n_terminus", "ch_cc1_linker", "cc1_cc2_linker",
                       "cc2_cc3_linker", "cc3_ckk_linker", "c_terminus"]:
            if region in b: # Checks if this region exists
                start, end = b[region]
                region_seq = slice_region(seq, start, end)
                constructs.append(
                    SeqRecord(
                        Seq(region_seq),
                        id=f"{protein}_{region}",
                        description=f"{protein} {region}"
                    )
                )

       # CH + CKK WITH FLEXIBLE LINKER (GGGGSGGGGS) (3)
       # Fuses CH and CKK with a flexible linker

        if "ch" in b and "ckk" in b:
            ch_seq = slice_region(seq, b["ch"][0], b["ch"][1])
            ckk_seq = slice_region(seq,b["ckk"][0], b["ckk"][1])
            fused = ch_seq + linker + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_ch_ckk",
                    description=f"{protein} ch + linker + ckk"
                )
            )

        # COILED COIL REGION WITH CKK DOMAIN AND NATIVE LINKER (3)
        # Fuses coiled coil domain with ckk domain using native cc3_ckk linker

        if "cc" in b and "cc3_ckk_linker" in b and "ckk" in b:
            cc_seq = slice_region(seq, b["cc"][0], b["cc"][1])
            cc3_ckk_linker_seq = slice_region(seq, b["cc3_ckk_linker"][0], b["cc3_ckk_linker"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc_seq + cc3_ckk_linker_seq + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc_ckk",
                    description=f"{protein} cc + native linker + ckk"
                )
            )

        # CH DOMAIN WITH COILED COIL REGION AND NATIVE LINKER (3)
        # Fused ch domain with coiled coil region using native ch_cc1 linker

        if "ch" in b and "ch_cc1_linker" in b and "cc" in b:
            ch_seq = slice_region(seq, b["ch"][0], b["ch"][1])
            ch_cc1_linker_seq = slice_region(seq, b["ch_cc1_linker"][0], b["ch_cc1_linker"][1])
            cc_seq = slice_region(seq, b["cc"][0], b["cc"][1])
            fused = ch_seq + ch_cc1_linker_seq + cc_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_ch_cc",
                    description=f"{protein} ch + native linker + cc"
                )
            )

        # DELETION CONSTRUCTS (9)
        # Removes specific domains from full-length protein sequence

        # Removes the ch domain

        if "ch" in b:
            ch_start, ch_end = b["ch"] # Residues before and after ch domain
            del_ch = seq[:ch_start-1] + seq[ch_end:]
            constructs.append(
                SeqRecord(
                    Seq(del_ch),
                    id=f"{protein}_delCH",
                    description=f"{protein} with ch deleted"
                )
            )

        # Removes the coiled coil region

        if "cc" in b:
            cc_start, cc_end = b["cc"]
            del_cc = seq[:cc_start-1] + seq[cc_end:]
            constructs.append(
                SeqRecord(
                    Seq(del_cc),
                    id=f"{protein}_delCC",
                    description=f"{protein} without cc domain"
                )
            )

        # Removes the ckk domain

        if "ckk" in b:
            ckk_start, ckk_end = b["ckk"]
            del_ckk = seq[:ckk_start-1] + seq[ckk_end:]
            constructs.append(
                SeqRecord(
                    Seq(del_ckk),
                    id=f"{protein}_delCKK",
                    description=f"{protein} without ckk domain"
                )
            )

        # MBD + D2 constructs (server permitted)
        # C - Terminal constructs

        # cc3 + D2 + CKK

        if "cc3" in b and "cc3_ckk_linker" in b and "ckk" in b:
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            d2_seq = slice_region(seq, b["cc3_ckk_linker"][0], b["cc3_ckk_linker"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc3_seq + d2_seq + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc3_d2_ckk",
                    description=f"{protein} cc3 + d2 + ckk"
                )
            )

        # cc2 + mbd + cc3 + d2 + ckk

        if "cc2" in b and "cc2_cc3_linker" in b and "cc3" in b and "cc3_ckk_linker" in b and "ckk" in b:
            cc2_seq = slice_region(seq, b["cc2"][0], b["cc2"][1])
            mbd_seq = slice_region(seq, b["cc2_cc3_linker"][0], b["cc2_cc3_linker"][1])
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            d2_seq = slice_region(seq, b["cc3_ckk_linker"][0], b["cc3_ckk_linker"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc2_seq + mbd_seq + cc3_seq + d2_seq + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc2_mbd_cc3_d2_ckk",
                    description=f"{protein} cc2 + mbd + cc3 + d2 + ckk"
                )
            )

        # cc3 + flexible linker (GGGGSGGGGS) + ckk

        if "cc3" in b and "ckk" in b:
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc3_seq + linker + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc3_flex_ckk",
                    description=f"{protein} cc3 + flex + ckk"
                )
            )

        # cc2 + flex linker + cc3 + flex linker + ckk

        if "cc2" in b and "cc3" in b and "ckk" in b:
            cc2_seq = slice_region(seq, b["cc2"][0], b["cc2"][1])
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc2_seq + linker + cc3_seq + linker + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc2_flex_cc3_flex_ckk",
                    description=f"{protein} cc2 + flex + cc3 + flex + ckk"
                )
            )

        # cc2 + flex linker + cc3 + d2 + ckk

        if "cc2" in b and "cc3" in b and "cc3_ckk_linker" in b and "ckk" in b:
            cc2_seq = slice_region(seq, b["cc2"][0], b["cc2"][1])
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            d2_seq = slice_region(seq, b["cc3_ckk_linker"][0], b["cc3_ckk_linker"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc2_seq + linker + cc3_seq + d2_seq + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc2_flex_cc3_d2_ckk",
                    description=f"{protein} cc2 + flex + cc3 + d2 + ckk"
                )
            )

        # cc2 + mbd + cc3 + flex linker + ckk

        if "cc2" in b and "cc2_cc3_linker" in b and "cc3" in b and "ckk" in b:
            cc2_seq = slice_region(seq, b["cc2"][0], b["cc2"][1])
            mbd_seq = slice_region(seq, b["cc2_cc3_linker"][0], b["cc2_cc3_linker"][1])
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            ckk_seq = slice_region(seq, b["ckk"][0], b["ckk"][1])
            fused = cc2_seq + mbd_seq + cc3_seq + linker + ckk_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc2_mbd_cc3_flex_ckk",
                    description=f"{protein} cc2 + mbd + cc3 + flex + ckk"
                )
            )

         # cc2 + mbd + cc3 + d2 (without ckk)

        if "cc2" in b and "cc2_cc3_linker" in b and "cc3" in b and "cc3_ckk_linker" in b:
            cc2_seq = slice_region(seq, b["cc2"][0], b["cc2"][1])
            mbd_seq = slice_region(seq, b["cc2_cc3_linker"][0], b["cc2_cc3_linker"][1])
            cc3_seq = slice_region(seq, b["cc3"][0], b["cc3"][1])
            d2_seq = slice_region(seq, b["cc3_ckk_linker"][0], b["cc3_ckk_linker"][1])
            fused = cc2_seq + mbd_seq + cc3_seq + d2_seq
            constructs.append(
                SeqRecord(
                    Seq(fused),
                    id=f"{protein}_cc2_mbd_cc3_d2",
                    description=f"{protein} cc2 + mbd + cc3 + d2"
                )
            )
    return constructs

def save_to_fasta(constructs, output_file = "all_constructs.fasta"):
    """
    Write all constructs to a FASTA file

    FASTA format -> 
    >ID Description
    Sequence line 1
    Sequence line 2...

    PARAMETERS:

    constructs(list) -> list of SeqRecord objects
    output_file (str) -> output filename ("all_constructs.fasta")

    """

    # Use SeqIO.write() to write all constructs to a specific file titled "all_constructs.fasta"
    # "fasta" specifies output format

    SeqIO.write(constructs, output_file, "fasta")
    print(f"Saved {len(constructs)} constructs to {output_file}")

def save_constructs_list(constructs, output_file = "constructs_list.csv"):
    """
    Save a summary of all constructs to an output file

    CSV columns -> ID, Description, Length

    PARAMETERS:

    Constructs (list) -> list of SeqRecord objects
    output_file (str) -> output filename ("constructs_list.csv")

    """
    
    # Open the file in write mode with newline

    with open(output_file, 'w', newline= '') as f:
        writer = csv.writer(f) # Create a csv writer
        writer.writerow(["ID","Description","Length"]) # Write header row
        for rec in constructs:
            writer.writerow([rec.id, rec.description, len(rec.seq)]) # Write one row per construct
    
    print(f"Saved construct list to {output_file}")

def main():

    constructs = generate_constructs(sequences, boundries, linker= "GGGGSGGGGS")
    save_to_fasta(constructs, "all_constructs.fasta")
    save_constructs_list(constructs, "constructs_list.csv")

    generate_domain_ranges(constructs, "domain_ranges.csv")
    print(f"Total: {len(constructs)} constructs!")

if __name__ == "__main__":
    main()