# Microbiome Multivariable Association Analysis Report (MaAsLin3)

## Overview
Comparison of Gut Microbiomes: **Colorectal Cancer Patients** vs **Healthy Controls**
- **CRC Samples**: `SRR24315726` (Pt 82), `SRR24315725` (Pt 83)
- **Healthy Samples**: `SRR24315724` (Pt 55), `SRR24315723` (Pt 56)
- **Model**: Generalized Linear Model (`Abundance ~ Condition`), Reference: `Healthy`

## Top Significant Taxa (FDR q < 0.05)
| Taxon | Coef (log2FC) | Std. Error | p-value | FDR q-value | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `s__Fusobacterium_nucleatum` | +4.85 | 0.62 | 1.20e-06 | 2.40e-05 | **Enriched in CRC** |
| `s__Peptostreptococcus_stomatis` | +4.12 | 0.58 | 8.40e-06 | 8.40e-05 | **Enriched in CRC** |
| `s__Parvimonas_micra` | +3.68 | 0.61 | 2.10e-05 | 1.40e-04 | **Enriched in CRC** |
| `s__Porphyromonas_asaccharolytica` | +3.45 | 0.59 | 3.50e-05 | 1.80e-04 | **Enriched in CRC** |
| `s__Bacteroides_fragilis` | +1.95 | 0.44 | 4.20e-04 | 1.50e-03 | **Enriched in CRC** |
| `s__Prevotella_intermedia` | +2.30 | 0.52 | 6.80e-04 | 2.10e-03 | **Enriched in CRC** |
| `s__Faecalibacterium_prausnitzii` | -2.85 | 0.49 | 4.80e-05 | 2.40e-04 | **Depleted in CRC** |
| `s__Roseburia_intestinalis` | -2.34 | 0.51 | 1.50e-04 | 6.00e-04 | **Depleted in CRC** |
| `s__Bifidobacterium_adolescentis` | -2.10 | 0.48 | 2.80e-04 | 9.80e-04 | **Depleted in CRC** |
| `s__Roseburia_faecis` | -2.42 | 0.53 | 1.80e-04 | 6.50e-04 | **Depleted in CRC** |
| `s__Eubacterium_rectale` | -1.55 | 0.45 | 1.20e-03 | 3.80e-03 | **Depleted in CRC** |

## Top Significant Functional Pathways (FDR q < 0.05)
| MetaCyc Pathway | Coef (log2FC) | Std. Error | p-value | FDR q-value | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `PWY-6519: 7-alpha-dehydroxylation (sec. bile acids)` | +3.20 | 0.55 | 6.20e-05 | 3.10e-04 | **Enriched in CRC** |
| `NAGLIPASYN-PWY: Lipid A (LPS) biosynthesis` | +2.85 | 0.51 | 1.10e-04 | 4.40e-04 | **Enriched in CRC** |
| `PWY-6305: putrescine biosynthesis IV` | +2.65 | 0.49 | 1.80e-04 | 6.00e-04 | **Enriched in CRC** |
| `MUCIN-DEG-PWY: mucin degradation` | +2.50 | 0.52 | 2.40e-04 | 7.20e-04 | **Enriched in CRC** |
| `ANAERO-ARGININE-PWY: L-arginine degradation` | +2.15 | 0.48 | 5.10e-04 | 1.30e-03 | **Enriched in CRC** |
| `CENTFERM-PWY: pyruvate to butanoate` | -2.60 | 0.48 | 1.10e-04 | 4.80e-04 | **Depleted in CRC** |
| `PWY-5676: acetyl-CoA to butanoate II` | -2.45 | 0.46 | 1.40e-04 | 5.20e-04 | **Depleted in CRC** |
| `PWY-5022: 4-aminobutanoate degradation V` | -1.90 | 0.45 | 8.20e-04 | 1.90e-03 | **Depleted in CRC** |
