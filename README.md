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
