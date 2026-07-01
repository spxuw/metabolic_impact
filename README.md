# Microbiome Perturbation Modeling of Plasma Metabolic Impact

This repository provides code for a perturbation-based framework to model how changes in gut microbial community structure may propagate to host plasma metabolite profiles. The workflow uses species-level microbiome profiles from curatedMetagenomicData to train a community reconstruction model, performs in silico single- and pairwise-species perturbations, predicts the resulting post-perturbation microbiome composition, and projects these predicted compositions onto a microbiome–metabolite association matrix to estimate metabolic impact.

## Overview

The goal of this project is to move beyond differential abundance analysis and quantify how depletion of specific microbial species, or combinations of species, may alter the predicted plasma metabolome. The framework is designed as a data-driven approximation of microbiome–metabolite coupling rather than a fully mechanistic metabolic reconstruction model.

The analysis consists of four main steps:

1. **Train a community reconstruction model**  
   A continuous neural ODE (cNODE)-based model is trained to predict species composition from binary species assemblage patterns using species-level abundance profiles from curatedMetagenomicData.

2. **Perform in silico perturbations**  
   Individual species, or pairs of species, are removed from the observed community assemblage. The trained model is then used to predict the restructured post-perturbation microbial community.

3. **Project microbiome composition to plasma metabolites**  
   Predicted microbiome compositions are multiplied by a reference microbiome–metabolite association matrix to infer plasma metabolite profiles before and after perturbation.

4. **Quantify metabolic impact**  
   The metabolic impact of a perturbation is quantified as the displacement between the original and perturbed predicted metabolite profiles. Larger values indicate stronger predicted perturbation of the plasma metabolome.

## Conceptual Workflow

```text
Species abundance data
        |
        v
Binary species assemblage representation
        |
        v
Train cNODE model to reconstruct species composition
        |
        v
In silico perturbation
(single-species or pairwise-species removal)
        |
        v
Predicted post-perturbation microbiome composition
        |
        v
Projection using microbiome–metabolite association matrix
        |
        v
Predicted plasma metabolite shift
        |
        v
Metabolic impact score and downstream analyses
```

## Input Data

### Microbiome data

The model is trained using species-level gut microbiome profiles derived from curatedMetagenomicData. These profiles are transformed into binary assemblage representations indicating species presence or absence, while retaining quantitative abundance information for model training and evaluation.

### Microbiome–metabolite association matrix

The microbiome–metabolite matrix represents data-driven associations between microbial taxa and plasma metabolites. This matrix is used to project microbiome composition into predicted metabolite space. Because this matrix reflects statistical coupling rather than direct biochemical causality, results should be interpreted as predicted metabolite-level consequences of microbiome perturbation, not as fully mechanistic flux estimates.

## Main Outputs

### Predicted post-perturbation microbiome composition

For each perturbation event, the trained cNODE model generates a predicted microbial composition after removal of a target species or species pair.

### Metabolic impact score

For each perturbation event, the predicted metabolite profile after perturbation is compared with the original predicted metabolite profile. The resulting distance is used as a global metabolic impact score.

### Genus/species-level impact profiles

Metabolic impact scores can be summarized across samples to identify taxa with disproportionate predicted effects on the plasma metabolome.

### Metabolite vulnerability

Metabolites can be ranked by their average sensitivity to microbial perturbations, identifying plasma metabolites that are especially vulnerable to changes in gut microbial community structure.

### Pairwise perturbation effects

Pairwise species-removal simulations allow comparison between observed double-removal impact and the additive expectation from corresponding single-species removals. This analysis is used to evaluate nonlinear and partially redundant perturbation effects.

## Example Analyses

The framework supports several downstream analyses:

- ranking high-impact microbial taxa;
- identifying metabolites vulnerable to microbial perturbation;
- comparing metabolic impact between healthy and disease-associated microbiomes;
- evaluating interaction-aware versus interaction-free perturbation models;
- testing robustness across microbiome–metabolite association thresholds;
- decomposing global metabolic impact into metabolite-class or metabolite-subclass contributions;
- validating predicted genus–metabolite relationships in independent dietary intervention datasets.

## Interpretation

This framework is intended to provide a scalable systems-level approximation of how microbiome disruption may translate into plasma metabolic consequences. It does not infer direct biochemical causality for individual genus–metabolite links. Instead, it uses population-scale microbiome–metabolite coupling to estimate how perturbations to microbial community structure may propagate to host metabolic profiles.

This distinction is important: the model is most useful for prioritizing taxa, metabolites, and perturbation patterns for downstream validation, rather than serving as a replacement for genome-scale metabolic models, isotope tracing, or experimental microbial knockouts.

## Repository Structure

A suggested repository structure is:

```text
.
├── README.md
├── data/
│   ├── microbiome_profiles/
│   ├── assemblage_matrices/
│   └── microbiome_metabolite_matrix/
├── scripts/
│   ├── 01_prepare_data.py
│   ├── 02_train_cnode.py
│   ├── 03_single_species_perturbation.py
│   ├── 04_pairwise_species_perturbation.py
│   ├── 05_project_to_metabolites.py
│   └── 06_downstream_analysis.R
├── results/
│   ├── predicted_compositions/
│   ├── perturbation_scores/
│   └── figures/
└── environment.yml
```

## Minimal Usage

### 1. Prepare microbiome assemblage data

```bash
python scripts/01_prepare_data.py \
  --input data/microbiome_profiles/species_abundance.tsv \
  --output data/assemblage_matrices/
```

### 2. Train the cNODE model

```bash
python scripts/02_train_cnode.py \
  --assemblage data/assemblage_matrices/binary_species_matrix.tsv \
  --abundance data/microbiome_profiles/species_abundance.tsv \
  --output results/model/
```

### 3. Run single-species perturbations

```bash
python scripts/03_single_species_perturbation.py \
  --model results/model/ \
  --assemblage data/assemblage_matrices/binary_species_matrix.tsv \
  --output results/predicted_compositions/single_species/
```

### 4. Run pairwise-species perturbations

```bash
python scripts/04_pairwise_species_perturbation.py \
  --model results/model/ \
  --assemblage data/assemblage_matrices/binary_species_matrix.tsv \
  --species-list data/selected_species.txt \
  --output results/predicted_compositions/pairwise_species/
```

### 5. Project predicted microbiome compositions to metabolite space

```bash
python scripts/05_project_to_metabolites.py \
  --composition results/predicted_compositions/single_species/ \
  --association-matrix data/microbiome_metabolite_matrix/assoc_mat.tsv \
  --output results/perturbation_scores/
```

## Notes on Robustness

Recommended sensitivity analyses include:

- comparing Euclidean and cosine distance-based metabolic impact scores;
- repeating analyses across multiple thresholds of the microbiome–metabolite association matrix;
- comparing interaction-aware perturbations with interaction-free and randomized baselines;
- evaluating single-species and pairwise-species perturbation consistency;
- validating predicted genus–metabolite relationships using independent intervention datasets when available.

## Citation

If you use this code or framework, please cite the associated manuscript:

```text
[Manuscript citation to be added]
```

## License

Please specify the license for this repository, for example MIT, Apache-2.0, or GPL-3.0.
