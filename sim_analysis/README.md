# Overview

This is a project that was used to process the MC simulation for $K^*(892)$ in Au+Au at $\sqrt{s_{NN}} = 200$ GeV. The project consist of macros located in src/ that can be compiled by running 

```sh
make all
```

The output of compilation are executables .exe that can be run

Parameters for different macros are located in structure Par that can be found in a header file with the same name as macros file located in lib/ directory. Some macros require output files from other macros located in analysis directory.

# Requirements

- [ROOT6](https://root.cern/) or newer compiled GNU GCC++17 or newer
- [ProgressBar](https://github.com/Sergeyir/ProgressBar)
- MC TTrees need to be deposited in ../data/Run7AuAu

MC TTrees can be found in simulation output directory specified in AN in analysis organisation section

# Macros
- HeatMapper.cpp - processes MC TTrees into heatmaps
- EmbAnalyze.cpp - processes MC TTrees into histograms from which embedding can be calculated
- SingleAnalyze.cpp - processes MC TTrees into histograms from which $\epsilon_{reg}$ can be calculated
- WidthAnalyze.cpp - processses MC TTrees into histograms from which width of GDF can be calculated
- PairAnalyze.cpp - processes MC TTrees into histograms from which reconstruction efficiency of $K^*(892)$ can be calculated
