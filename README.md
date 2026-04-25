# SDSS2507 Alpha/Gamma Orientation Relationship EBSD Analysis

This repository contains MATLAB/MTEX scripts, raw EBSD `.ang` files, and generated figure outputs for analyzing alpha/gamma orientation relationships in LPBF-fabricated super duplex stainless steel 2507.

The workflow evaluates crystallographic orientation relationships between ferrite/alpha and austenite/gamma phases and generates clean EBSD maps, deviation histograms, and publication-ready outputs.

## Repository structure

```text
SDSS2507-alphaGamma-OR-EBSD/
├── SDSS2507_paper_alphaGamma_OR_withDeviationHists_cleanMaps.m
├── README.md
├── data/
│   └── raw/
│       ├── 01_AS.ang
│       ├── 02_SR400.ang
│       ├── 03_SR450.ang
│       ├── 04_SR500.ang
│       ├── 05_SR550.ang
│       └── 06_SA1100.ang
└── figures/
    └── exported figure files
