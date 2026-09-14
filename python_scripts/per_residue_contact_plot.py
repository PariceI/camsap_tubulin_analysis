"""
This script generates a per-residue interaction plot
for each human CAMSAP from full-length per residue contact data.
"""

# IMPORT LIBRARIES

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Hardcoded files

C1_FILE = "per_residue_summary/Camsap1_full_tubulin_complex_model_per_residue.csv"
C2_FILE = "per_residue_summary/Camsap2_full_tubulin_complex_model_per_residue.csv"
C3_FILE = "per_residue_summary/Camsap3_full_tubulin_complex_model_per_residue.csv"

# DOMAIN BOUNDRIES

domains_c1 = [
    ('N-terminus', 1, 215),
    ('CH', 216, 331),
    ('CH-CC1 linker', 332, 872),
    ('CC1', 873, 909),
    ('CC1-CC2 linker', 910, 1015),
    ('CC2', 1016, 1048),
    ('MBD (CC2-CC3 linker)', 1049, 1290),
    ('CC3', 1291, 1343),
    ('D2 (CC3-CKK linker)', 1344, 1462),
    ('CKK', 1463, 1597),
]

domains_c2 = [
    ('N-terminus', 1, 221),
    ('CH', 222, 335),
    ('CH-CC1 linker', 336, 755),
    ('CC1', 756, 793),
    ('CC1-CC2 linker', 794, 886),
    ('CC2', 887, 926),
    ('MBD (CC2-CC3 linker)', 927, 1165),
    ('CC3', 1166, 1238),
    ('D2 (CC3-CKK linker)', 1239, 1348),
    ('CKK', 1349, 1483),
]

domains_c3 = [
    ('N-terminus', 1, 202),
    ('CH', 203, 312),
    ('CH-CC1 linker', 313, 593),
    ('CC1', 594, 628),
    ('CC1-CC2 linker', 629, 695),
    ('CC2', 696, 729),
    ('MBD (CC2-CC3 linker)', 730, 895),
    ('CC3', 896, 936),
    ('D2 (CC3-CKK linker)', 937, 1108),
    ('CKK', 1109, 1243),
]

# LOAD DATA

def load_per_residue(filepath):
    """
    Load per-residue csv and return a dictionary of resiue : contacts.
    """

    df = pd.read_csv(filepath)
    return dict(zip(df['residue'], df['contacts']))

c1_data = load_per_residue(C1_FILE)
c2_data = load_per_residue(C2_FILE)
c3_data = load_per_residue(C3_FILE)

# PLOT FUNCTION

def plot_interactions(data, domains, title, color, window = 1, output_file = None):
    """
    Plot per residue interactions with domain annotations.
    """

    max_residue = domains[-1][2]

    # A dense array of zeros for the full protein length

    dense = np.zeros(max_residue + 1)
    for res, count in data.items():
        if res <= max_residue:
            dense[res] = count

    # Create x-axis
    
    residues = np.arange(1, max_residue + 1)
    contacts = dense[1:]

    smoothed = pd.Series(contacts).rolling(window, min_periods=1, center=True).mean().fillna(0)

    fig, ax = plt.subplots(figsize=(12,6))

    ax.plot(residues, smoothed, color=color, linewidth=1.5)

    ax.set_xlim(0, max_residue)

    y_max = smoothed.max() if not smoothed.empty and smoothed.max() > 0 else 1
    y_pos = y_max * 1.05

    for domain, start, end in domains: 
       ax.axvspan(start, end, alpha = 0.15, color = 'white')
       ax.axvline(start, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
       mid = (start + end) / 2
       ax.text(mid, y_pos, domain, ha='center', va='top', fontsize = 7, rotation = 45, color = 'black', weight = 'bold')

    ax.set_xlabel('Residue Position', fontsize = 7)
    ax.set_ylabel('Contacts per residue', fontsize = 7)
    ax.set_title(title, fontsize = 14, fontweight = 'bold')
    ax.legend(loc='upper right', fontsize = 8)

    ax.set_ylim(0, y_max * 1.2)

    plt.tight_layout()

    if output_file:
           plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")

    plt.show()
    return fig, ax

# CAMSAP1
plot_interactions(
    data=c1_data,
    domains=domains_c1,
    title='CAMSAP1',
    color='lightblue',          
    window=1,
    output_file='interaction_plot_CAMSAP1.png'
)

# CAMSAP2
plot_interactions(
    data=c2_data,
    domains=domains_c2,
    title='CAMSAP2',
    color='lightgreen',          
    window=1,
    output_file='interaction_plot_CAMSAP2.png'
)

# CAMSAP3
plot_interactions(
    data=c3_data,
    domains=domains_c3,
    title='CAMSAP3',
    color='pink',         
    window=1,
    output_file='interaction_plot_CAMSAP3.png'
)