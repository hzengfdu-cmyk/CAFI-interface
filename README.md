# CAFI Interface Analysis Code for Bladder Cancer

This repository contains the R and Python scripts used for the main computational analyses of the study. The scripts are organised by analysis module and document the methods and key parameters used in the manuscript.

CAFI refers to the spatial interface between cancer-associated fibroblast-associated stroma and invadopodiahi epithelial cells.

## Analysis scripts

| File | Analysis module |
|---|---|
| `01_scRNA_QC_BBKNN.py` | Single-cell quality control, Scrublet, BBKNN integration and Leiden clustering |
| `02_malignant_identification.py` | Malignant epithelial-cell identification with inferCNVpy and scMalignantFinder |
| `03_invadopodia_scoring.R` | Invadopodia scoring, cell states and sample groups |
| `04_single_cell_DEG_GSEA.R` | Single-cell differential expression and GSEA |
| `05_SuperCell_pySCENIC.R` | SuperCell construction and pySCENIC parameters |
| `06_RCTD.R` | Spatial-transcriptomics deconvolution with RCTD |
| `07_CAFI_interface.R` | CAFI-interface identification and proximal/distal regions |
| `08_BANKSY.R` | BANKSY tissue domains and the tumour invasive front |
| `09_spatial_pseudobulk_TIDE.R` | Spatial pseudobulk analysis, CAFI signatures and TIDE inputs |
| `10_CellChat.R` | CellChat communication analysis and probability-weighted Ro/e |
| `11_LR_four_layer_scoring.R` | Four-layer ligand–receptor evidence scoring |
| `12_bulk_integration_scoring.R` | Bulk-cohort integration, ComBat and signature scoring |
| `13_MFP_TIDE.py` | MFP microenvironment subtyping and bulk TIDE analysis |
| `14_survival_Cox.R` | Kaplan–Meier and Cox survival analyses |
| `15_genomic_analysis.R` | Mutation and copy-number analyses |
| `16_mIF_segmentation.py` | Multiplex immunofluorescence cell segmentation |
| `17_mIF_phenotyping_interface.R` | mIF phenotyping, CAFI interfaces and density analyses |
| `18_mIF_distance_analysis.R` | mIF distance binning and correlation analyses |

## Software

The analyses used R 4.4.1 and Python 3.11.6.

- R: Seurat 5.1.0, spacexr 2.2.1, BANKSY 1.0.0, DESeq2 1.44.0, muscat 1.18.0, CellChat 1.6.1, SuperCell 1.0.1 and IOBR 0.99.0.
- Python: Scanpy 1.10.1, AnnData 0.10.8, BBKNN 1.6.0, infercnvpy 0.6.0, pySCENIC 0.12.1, Cellpose 4.2.1.1 and TIDEpy 1.3.9.

## Data and usage

The scripts use `data/` for input files and `results/` for output files. Expression matrices, clinical data, images and intermediate objects are not included in this repository. Dataset and gene-set sources are described in Supplementary Tables 2, 3 and 6.

The scripts document the analysis workflow and key parameters. Local input objects, metadata fields, image channels and software environments may require adaptation before analysis.
