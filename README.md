# Ligand-receptor-LR-interaction-analysis
This repository contains code and workflows for performing ligand–receptor (LR) interaction analysis on Xenium spatial transcriptomics data using stLearn

Ligand–Receptor Interaction Analysis for Xenium Data Using stLearn
The pipeline integrates spatial coordinates with gene expression profiles to infer cell–cell communication within the tissue microenvironment, enabling the identification of spatially informed signaling interactions between neighboring cells.

# Key Features:

Spatially Informed LR Inference: Identifies ligand–receptor interactions by incorporating spatial proximity and gene expression from Xenium single-cell data.

Xenium-Specific Processing: Optimized for high-resolution Xenium spatial transcriptomics output, including cell segmentation and spatial coordinates.

Cell-Type–Aware Communication: Infers signaling interactions between annotated cell types or spatial domains.

Spatial Neighborhood Analysis: Detects local signaling niches and microenvironment-driven interactions.

Visualization: Spatial interaction maps, LR network plots, heatmaps, and dot plots to explore signaling patterns across tissue sections.

Integration with stLearn Framework: Leverages stLearn’s spatial graph learning and neighborhood modeling capabilities.

# Requirements:

Python (for stLearn)

stLearn

Scanpy / AnnData

Xenium spatial transcriptomics output

<img width="4800" height="4800" alt="813_scc_CCI_LR_chord_plot_all" src="https://github.com/user-attachments/assets/0c01b91c-4c78-4be9-9a39-cc72a14dcd1f" />


<img width="9600" height="3000" alt="813_scc_CCI_LR" src="https://github.com/user-attachments/assets/db8988f6-3c8f-45b9-a650-4edc38b12c03" />


<img width="15000" height="4800" alt="813_scc_Detailed_Cell_Annotations_clustering_gridding" src="https://github.com/user-attachments/assets/cf08748b-fc54-48b7-977d-75b09712cf5f" />

<img width="15000" height="4800" alt="813_scc_Detailed_Cell_Annotations_gridding" src="https://github.com/user-attachments/assets/bff1d22e-e30d-4c2c-840e-d612f99f48c6" />


