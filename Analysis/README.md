# Overview

This is a project that was used to process the real data events and simulation for $K^*(892)$ in Au+Au at $\sqrt{s_{NN}} = 200$ GeV. The project consist of macros located in src/ that can be run with [ROOT CLING](https://root.cern/cling/). Parameters for different macros are located in structure Par that can be found either in the macros file or header file with the same name as macros file located in lib/ directory. Some macros require output files from other macros and some macros need to be run with different parameters to output all results.

# Requirements

- [ROOT6](https://root.cern/) or newer compiled GNU GCC++17 or newer
- Data files from taxi 19062 need to be placed in ../data/Run7AuAu directory
- Data files from taxi 19079 need to be placed in ../data/taxi/Run7AuAu/19079 directory

Taxi directories with data files are presented in AN in analysis organisation section

# Macros
- ADCEff.cc - adc cut efficiency caclulation
- InteractCanv.cc - GUI heatmaps figucial cutter
- MS.cc - drawing charged hadrons $m^2$ distribution and extracting yields
- GoodRunsSelector.cc - chi2/ndf and average to fit run by run DC runs selector
- FitSpecctra.cpp - charged hadrons spectra approximation for MC $p_T$ reweight
- DM.cc - drawing heatmaps and acceptance uncertainty calculation
- RegEff.cc - calculation of $\epsilon_{reg}$ for EMCal
- M2Eff.cc - calculation of $\epsilon_{m^2}$ and $\epsilon_{id}$
- Width.cc - calculation of $\sigma_{GDF}$
- PariEff.cc - calculation of reconstruction efficiency of $K^*(892)$
- InvM.cc - calculation of invariant mass distribution, calculation of invariant $p_T$ spectra, calculation of systematical uncertainties
- GuiInvM.cc - $K^*(892)$ approximation tweaking GUI tool for bad approximations
- DrawSpectra.cc - drawing invariant $p_T$ spectra
- DrawRAB.cc - drawing RAB
- DrawRCP.cc - drawing RCP
