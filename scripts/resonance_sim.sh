#!/bin/bash

for file in $@; do
   ./bin/AnalyzeSimResonance $file
   ./bin/AnalyzeSimResonance $file -s ptscale 0.995
   ./bin/AnalyzeSimResonance $file -s ptscale 1.005
   ./bin/AnalyzeSimResonance $file -s acceptance -1
   ./bin/AnalyzeSimResonance $file -s acceptance 1
   ./bin/AnalyzeSimResonance $file -s cuts -1
   ./bin/AnalyzeSimResonance $file -s cuts 1
   ./bin/EstimateRecEffOfResonance $file
done
