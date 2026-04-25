## Phylogenetic Analysis Pipeline

An automated bioinformatics workflow for fetching sequences from NCBI, performing Multiple Sequence Alignment (MSA), concatenating markers, and inferring trees using both Maximum Likelihood (RAxML-NG) and Bayesian Inference (MrBayes).

## Workflow Overview

Fetch: Downloads sequences from NCBI Entrez via accession numbers.

Clean: Standardizes FASTA headers and maps accessions to specific strain identifiers.

Align: Aligns sequences using MAFFT.

Concatenate: (Optional) Merges multiple genetic markers into a single supermatrix.

Model Selection: Determines the best-fit substitution model using ModelTest-NG.

Inference:

  ML Tree: Runs RAxML-NG with automatic bootstrap convergence.

  Bayesian Tree: Runs MrBayes with optimized convergence criteria.

Visualization: Generates publication-quality Phylograms and Cladograms in PDF and SVG formats.

## Environment & Requirements

This pipeline is optimized and tested for the following software and library versions:

System Dependencies
Tool            Version     Description
Python          3.14.3      Core programming language
Grep            2.32.5      Text-processing utility
MAFFT           v7.526      Multiple Sequence Alignment
Modeltest-ng    v0.1.7      Evolutionary model selection
RAxML-NG        v1.2.0      Maximum Likelihood tree inference
MrBayes         3.2.7       Bayesian Inference of Phylogeny

Python Libraries
Library         Version     Use Case
biopython       1.87        Sequence manipulation and I/O
toytree         3.0.11      Tree manipulation and rooting
toyplot         2.0.0       Canvas rendering and PDF/SVG export
requests        Latest      API calls to NCBI Entrez


## Usage

The pipeline is executed through the central run.sh script.

1. Single Marker Analysis (e.g., ITS)

To run a specific marker, pass the path to the text file containing your accession numbers:
Bash

    bash scripts/run.sh ITS

2. Multi-marker Concatenated Analysis

To run a concatenated analysis of all markers present in the importante/ folder:
Bash

    bash scripts/run.sh concatenated

## Project Structure

scripts/run.sh: The main entry point for all analyses.

scripts/pipeline.py: The engine coordinating Step 5 through Step 7.

scripts/fetch_sequences.py: Handles NCBI Entrez downloading.

scripts/simple_fasta.py: Cleans headers and renames sequences based on a strain dictionary.

scripts/Seqconcact.py: Concatenates individual gene alignments.

scripts/desenhar_arvore.py: Tree visualization for single markers.

scripts/desenhar_arvore_concatenated.py: Tree visualization for multi-marker matrices.

## Outputs

The pipeline creates a specific directory for each run (e.g., arvore_ITS/ or arvore_concatenated/) containing:

.fasta: Raw and aligned sequences.

raxml.raxml.support: ML tree with bootstrap values.

  *.mrbayes.nex.con.tre: Bayesian consensus tree.

  .pdf / .svg: High-resolution tree visualizations (Phylogram and Cladogram).
