# Microbiome Perturbation Modeling of Plasma Metabolic Impact

This repository provides code for a perturbation-based framework to model how changes in gut microbial community structure may propagate to host plasma metabolite profiles. The workflow uses genus-level microbiome profiles from curatedMetagenomicData to train a community reconstruction model, performs in silico single- and pairwise-species perturbations, predicts the resulting post-perturbation microbiome composition, and projects these predicted compositions onto a microbiome–metabolite association matrix to estimate metabolic impact.

## Overview

The goal of this project is to move beyond differential abundance analysis and quantify how depletion of specific microbial species, or combinations of species, may alter the predicted plasma metabolome. The framework is designed as a data-driven approximation of microbiome–metabolite coupling rather than a fully mechanistic metabolic reconstruction model.

The analysis consists of four main steps:

1. **Train a community reconstruction model**  
   A continuous neural ODE (cNODE)-based model is trained to predict species composition from binary species assemblage patterns using genus-level abundance profiles from curatedMetagenomicData.

2. **Perform in silico perturbations**  
   Individual genera, or pairs of genus, are removed from the observed community assemblage. The trained model is then used to predict the restructured post-perturbation microbial community.

3. **Project microbiome composition to plasma metabolites**  
   Predicted microbiome compositions are multiplied by a reference microbiome–metabolite association matrix to infer plasma metabolite profiles before and after perturbation.

4. **Quantify metabolic impact**  
   The metabolic impact of a perturbation is quantified as the displacement between the original and perturbed predicted metabolite profiles. Larger values indicate stronger predicted perturbation of the plasma metabolome.

## Conceptual Workflow

```text
Species abundance data
        |
Binary species assemblage representation
        |
Train cNODE model to reconstruct species composition
        |
In silico perturbation
(single-genus or pairwise-genura removal)
        |
Predicted post-perturbation microbiome composition
        |
Projection using microbiome–metabolite association matrix
        |
Predicted plasma metabolite shift
        |
Metabolic impact score and downstream analyses
```

## Input Data

### Microbiome data

The model is trained using genus-level gut microbiome profiles derived from curatedMetagenomicData (version 3.16.1). These profiles are transformed into binary assemblage representations indicating genus presence or absence, while retaining quantitative abundance information for model training and evaluation.

To predict the post-perturbation microbiome composition, we removed each genus in each sample (e.g., by setting z_i = 0, if we need to remove a present genus i in a sample, and z is binary genus assemblage).

### Microbiome–metabolite association matrix

The microbiome–metabolite matrix represents data-driven associations between microbial taxa and plasma metabolites (source paper: An online atlas of human plasma metabolite signatures of gut microbiome composition). This matrix is used to project microbiome composition into predicted metabolite space. Because this matrix reflects statistical coupling rather than direct biochemical causality, results should be interpreted as predicted metabolite-level consequences of microbiome perturbation, not as fully mechanistic flux estimates.

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



## Minimal Usage

```bash
python python Keystone_cNODE_metabolite.py --dataset "KeyPlasma"
The input data in folder KeyPlasma includes: Z_train.csv, P_train.csv, Z_trainr.csv and Z_testr.csv. This will generate two key outputs: qtst1.csv (predicted post-perturbation microbime comppsition of healthy samples) and qtst2.csv (predicted post-perturbation microbime comppsition of diseased samples).
```
