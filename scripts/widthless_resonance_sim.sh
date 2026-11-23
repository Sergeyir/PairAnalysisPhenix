#!/bin/bash

for file in $@; do
   ./bin/AnalyzeSimWidthlessResonance $file
   ./bin/EstimateGaussianBroadening $file
done
