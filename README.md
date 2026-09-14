CAMSAP-tubulin interactions analysis

Contains Python and R scripts and data for an MSc thesis
investigating how non-CKK regions contribute to CAMSAP-tubulin
interactions

This repository contains custom scripts and data used to analyse 
AlphaFold3 structural predictions of CAMSAP1, CAMSAP2, and CAMSAP3
in a complex with two alpha/beta-tubulin heterodimers. The analysis
examines contacts, buried surface area (BSA), SASA, pLDDT confidence
scoring and sequence conservation across constructs.

CONTENTS

Python scripts:

generate_json.py -> Generates AlphaFold3 JSON input files
domain_sasa.py -> Calculates BSA and SASA for whole complex & per domain
calculate_contacts.py -> Calculates contacts & pLDDT score for whole constuct, per domain and per-residue
generate_constructs_domain_mapping.py -> Generates constructs given domain boundaries & maps domain boundaries to each construct
per_residue_contact_plot.py -> Calculates percent identity per position from a MAFFT alignment & generates conservation plot

R scripts:

R scripts were used purely for data visualisation and figure creation

Data:

per_residue_summary -> Folder containing per-residue contacts and pLDDT scores for all contsructs
contacts_summary -> CSV file containing contacts summary for all constructs
domain_sasa_summary -> CSV file containing SASA/BSA summary per domain in all constructs
per_domain_summary -> CSV file containing contacts summary per domain in all constructs
sasa_bsa_results -> CSV file containing SASA/BSA for full construct(s)

REQUIRES:
Python 3.x
Biopython
R with ggplot2

This script needs to be run in the same directory in which the CIF files are contained.

NOTES

CAMSAP is modelled as chain A
Alpha tubulin is modelled as chain B, C
Beta tubulin is modelled as chain D, E
SASA and BSA are calculated used Shrake-Rupley as implemented by BioPython

Author: Parice Lott
