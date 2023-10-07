# Trimer-based TCR embedding and Repertoire Functional Units

This repository provides the codebase and associated files to implement the RFU method. Functions are designed to run in R console. 

## Input files
RFU method is compatible with the input files for [GIANA](https://github.com/s175573/GIANA). Users should follow the instructions of this repository to generate input files directly from TCR-seq samples produced from Adaptive Biotechnology. 

## Usage
In the R console, first load the associated data files:

`load('trimerMDSfit_small.Rdata')`

This command creates a variable named `fit` in the environment. 

`load('km5000noMax.Rdata')`

This command creates a variable named  `km5000noMax`.

