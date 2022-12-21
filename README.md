# Measuring Inter-Brain Synchrony: Methods and Pitfalls

[![DOI](https://zenodo.org/badge/580886256.svg)](https://zenodo.org/badge/latestdoi/580886256)

Marten de Vries

> Collecting EEG data for two participants simultaneously during a task (i.e.,
hyperscanning) allows us to study their social interaction. Of particular
interest is their inter-brain synchrony (IBS), i.e. how functionally similar
their neural oscillations are. Using a full IBS analysis of a tacit coordination
experiment and simulations, we study the effect of different methodological
choices in an IBS measurement pipeline. Three ways to quantify IBS are studied:
the phase-locking value (PLV), the circular correlation (CCorr) and the
imaginary part of coherency (ImagCoh). Each measures functional similarity in
its own way and has its advantages and disadvantages. We find the CCorr measure
to be less stable than the others, but still recommend its use along with the
PLV measure because of its natural interpretation and good performance on
simulated data. We present a robust version of the circular correlation measure
and make recommendations for how to best perform an IBS analysis.

This repository contains the code of my final project for the Computational
Cognitive Science study at the University of Groningen.

## What's in here?

- ``code/prototype.m``; The main MATLAB file. Implements the time shuffle and
  dyad shuffle permutation tests.
- ``code/calculateMeasures.m``: Contains MATLAB implementations of the phase
  locking value, circular correlation, imaginary part of coherency and robust
  circular correlation (Algorithm 2) measures.
- ``stats/prediction.Rmd``; Python code that trains and evaluates classifiers.
  R code that plots the prediction results.
- ``stats/new_stats.Rmd``; all other R code. Includes statistics, visualization,
  and simulation code. Some highlights that might be useful for other projects
  include:
  - line 75, a function to plot two phase components against each other, showing
    the circularity of phase angles using repetition. Used to generate e.g.
    Figure 6.
  - line 725, a function that performs a permutation test. Nearby is also code
    to perform the Benjami-Hochberg FDR procedure.
  - line 837, an R implementation of the circular correlation coefficient.
  - line 849, an R implementation of the phase locking value.
  - line 853, an R implementation of the imaginary part of coherency.
  - line 916, an implementation of Algorithm 1, which generates example signals
    that minimize some criterium. Can be used to e.g. find example signals for
    specific measure values. See the simulation chapter for more details.
- ``stats/results/``; mostly graphs used in the thesis or final presentation.
- ``thesis/``; the LaTeX files of the thesis.
- ``presentation/``; the slides of the final presentation.
- ``slides/``; slides with progress reports and questions that came up
  throughout this project.
- ``extra/``; side projects.


## What's not in here?

All data necessary to re-run my code. Some intermediate data files used to
generate plots have been included in ``stats/results/``, but most original files
are too big, potentially privacy-sensitive and not mine to share.

Also parts of the MATLAB code base are missing. E.g. code to pre-process EEG
data and load the data for the inter-brain synchrony analysis. The reason for
this is that that code is not mine to share, as it was not written by me.
