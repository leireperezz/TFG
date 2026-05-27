# AUTOMATED MULTIPLE SEQUENCE ALIGNMENT OF THE HUMAN PROTEOME BASED ON PFAM FAMILIES
 
## Experimental Bachelor's Thesis, degree in biotechnology


#### Author: Leire Pérez Palacios
#### Supervisor: Dra. Mireia Olivella (UVic),  Dr. Arnau Cordomí (UPF) and Gabriel Ruiz (Uvic)
#### Date: June 2026
#### Keywords: Multiple sequence alignment, Pfam families, variant pathogenicity, human proteome, missense variants



## ABSTRACT

The prediction of pathogenicity for missense variants is fundamental for the clinical interpretation of genetic diseases. HomolVar is a web server that predicts the pathogenicity of missense variants by extrapolating the pathogenicity observed in homologous variants. Although its predictive power is superior to 95%, its coverage is limited by the low availability of high-quality Pfam alignments. In order to improve HomolVar's coverage (currently with 616 Pfam families), an automated pipeline based on HMMER was developed to generate reconstructed alignments of the human proteome organized by Pfam families. The analysis identified 8,580 Pfam families, resulting in 4,218 valid families (increment of ~585%). Quality assessment using BLOSUM62 demonstrated that reconstructed alignments outperformed seed alignments in 75.04% of families (average improvement of 37.23%) and full alignments in 63.44% (improvement of 5.84%). Subsequently, the new alignments were used to identify homologous variant pairs, obtaining 33,905 pairs in 451 families, representing approximately 13-fold more pairs than seed alignments (2,609 pairs). The predictive power of HomolVar with this strategy was established at 86% (MCC 0.70). In conclusion, reconstructed alignments thus provide an effective solution to significantly increase HomolVar's coverage while maintaining high predictive power.

## PIPELINE
#### 1. UniProt data:  Reviewed human proteome downloaded from UniProt (organism: Homo sapiens).
#### 2. Pfam domain coordinates: Queries the InterPro API to retrieve the Pfam domain coordinates (start and end positions) for each human protein.
#### 3. Domain sequence extraction: Extracts domain subsequences from the full protein sequences and generates one FASTA file per Pfam family.
#### 4. Alignment generation:
- 4_0 Retrieves the short name of each Pfam domain from the InterPro API.
- 4_1 Runs hmmpress, hmmfetch and hmmalign. This step was executed on an HPC server using SLURM job scripts.
- 4_2 Evaluates the alignment output by comparing input FASTA counts with the generated Stockholm files.
- 4_3 Postprocesses Stockholm files and remaps sequence identifiers to UniProt entry names.

#### 5. BLOSUM62 scores: Calculates global BLOSUM62 mean scores for rec, seed and full alignments without sequence filtering.
#### 6. Alignment quality comparison: Compares rec vs seed and rec vs full alignments on shared sequences, computing BLOSUM62 scores and conservation metrics under two treatments of lowercase residues.
- 6_1 Filters rec, seed and full alignments to keep only shared protein sequences for a fair comparison.
- 6_2 Processes lowercase residues under two approaches: treating them as gaps or as valid amino acids.
- 6_3 Calculates and compares BLOSUM62 scores and conservation metrics between rec vs seed and rec vs ful
#### 7.Alignment selection: Selects the best alignment per Pfam family based on BLOSUM62 comparisons, prioritising rec and replacing with seed or full only when these score higher.
#### 8. Homologous pair extraction: Maps missense variants from ClinVar and gnomAD to alignment positions and extracts homologous variant pairs to evaluate alignment performance. This code is adapted from the HomolVar pipeline originally developed by Gabriel Ruiz, modified to incorporate rec alignments and to consider all Pfam codes associated with each variant.

## REQUIREMENTS
- Python >= 3.9
- pandas, requests, biopython, blosum, colorama
- HMMER 3.4