# maf

This contains scripts and data for analyzing mutation annotation format (MAF) files, particularly in the context of The Cancer Genome Atlas (TCGA) data.

Key Components:

download_TCGA_clinical_biospecimen_data.R: Script for downloading clinical and biospecimen data from TCGA.

TCGA_clinical_analysis.Rmd: R Markdown document for performing clinical analyses on TCGA data.

g4_enrichment.R: Script for analyzing G-quadruplex (G4) enrichment in genomic data.

get_seq_from_maf.R: Script to extract sequences from MAF files.

intersect_maf_bed.R: Script to find intersections between MAF files and BED files with G4 region, useful for genomic region analyses.

g4_maf.Rdata: R data file, containing processed MAF data related to G-quadruplex analyses.

These tools facilitate the processing and analysis of MAF files, enabling users to conduct clinical analyses, sequence extractions, and genomic region intersections, particularly in cancer genomics research
