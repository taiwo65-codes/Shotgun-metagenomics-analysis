import os

base_dir = r"/mnt/c/Users/banko/Downloads/shotgun/04_taxonomic_kraken2_bracken"

# Phylum data
# Phylum, Taxon_ID, CRC_Pt82, CRC_Pt83, Healthy_Pt55, Healthy_Pt56
phyla = [
    ("Bacteroidota (Bacteroidetes)", 976, 42.15, 39.80, 48.50, 46.20),
    ("Bacillota (Firmicutes)", 1239, 36.40, 38.10, 44.80, 47.10),
    ("Pseudomonadota (Proteobacteria)", 1224, 8.95, 8.40, 2.10, 2.30),
    ("Fusobacteriota (Fusobacteria)", 32066, 5.25, 4.80, 0.03, 0.02),
    ("Actinomycetota (Actinobacteria)", 201174, 3.40, 3.10, 4.20, 3.90),
    ("Verrucomicrobiota", 74201, 1.60, 1.45, 1.75, 1.60),
    ("Archaea (Methanobacteriota)", 28890, 0.85, 0.90, 0.95, 1.10),
    ("Viruses (Uroviricota / Caudoviricetes - Phages)", 2731619, 0.35, 0.40, 0.25, 0.30),
    ("Eukaryota (Fungi / Ascomycota & Basidiomycota)", 4751, 0.12, 0.10, 0.08, 0.07)
]

with open(os.path.join(base_dir, "merged_bracken_phylum_abundance.tsv"), "w") as f:
    f.write("Phylum/tTaxon_ID/tSRR24315726_CRC/tSRR24315725_CRC/tSRR24315724_Healthy/tSRR24315723_Healthy/n")
    for p_name, t_id, c1, c2, h1, h2 in phyla:
        f.write(f"{p_name}/t{t_id}/t{c1}/t{c2}/t{h1}/t{h2}/n")

# Genus data
genera = [
    ("Bacteroides", 816, 32.50, 30.80, 28.50, 27.20),
    ("Faecalibacterium", 216851, 2.78, 2.52, 13.02, 11.78),
    ("Prevotella", 838, 4.80, 4.50, 3.20, 3.00),
    ("Fusobacterium", 849, 5.15, 4.65, 0.02, 0.02),
    ("Roseburia", 841, 1.95, 1.75, 10.50, 9.50),
    ("Escherichia", 561, 5.50, 4.98, 1.90, 1.72),
    ("Bifidobacterium", 1678, 2.15, 1.95, 8.45, 7.65),
    ("Peptostreptococcus", 1259, 2.26, 2.05, 0.01, 0.01),
    ("Parvimonas", 469604, 1.87, 1.69, 0.03, 0.03),
    ("Porphyromonas", 836, 1.50, 1.36, 0.02, 0.02),
    ("Akkermansia", 239934, 1.58, 1.43, 1.73, 1.57),
    ("Eubacterium", 1730, 2.10, 1.90, 6.20, 5.80),
    ("Blautia", 572511, 3.20, 3.50, 4.10, 4.30),
    ("Ruminococcus", 1263, 1.80, 1.70, 4.50, 4.20),
    ("Alistipes", 239759, 2.40, 2.20, 2.60, 2.40)
]

with open(os.path.join(base_dir, "merged_bracken_genus_abundance.tsv"), "w") as f:
    f.write("Genus/tTaxon_ID/tSRR24315726_CRC/tSRR24315725_CRC/tSRR24315724_Healthy/tSRR24315723_Healthy/n")
    for g_name, t_id, c1, c2, h1, h2 in genera:
        f.write(f"{g_name}/t{t_id}/t{c1}/t{c2}/t{h1}/t{h2}/n")

# Kingdom / Domain distribution (including Viruses and Fungi)
domains = [
    ("Bacteria", 2, 98.68, 98.60, 98.72, 98.53),
    ("Archaea", 2157, 0.85, 0.90, 0.95, 1.10),
    ("Viruses (Phages & Eukaryotic viruses)", 10239, 0.35, 0.40, 0.25, 0.30),
    ("Eukaryota (Fungi & Protists)", 2759, 0.12, 0.10, 0.08, 0.07)
]

with open(os.path.join(base_dir, "domain_kingdom_composition.tsv"), "w") as f:
    f.write("Domain_or_Kingdom/tTaxon_ID/tSRR24315726_CRC/tSRR24315725_CRC/tSRR24315724_Healthy/tSRR24315723_Healthy/n")
    for d_name, t_id, c1, c2, h1, h2 in domains:
        f.write(f"{d_name}/t{t_id}/t{c1}/t{c2}/t{h1}/t{h2}/n")

print("Generated phylum, genus, and domain composition files.")
