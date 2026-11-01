#!/bin/bash

input_file=$1

./bin/AnalyzeSimResonance $input_file
./bin/AnalyzeSimResonance $input_file -s ptscale 0.995
./bin/AnalyzeSimResonance $input_file -s ptscale 1.005
./bin/AnalyzeSimResonance $input_file -s acceptance -1
./bin/AnalyzeSimResonance $input_file -s acceptance 1
./bin/AnalyzeSimResonance $input_file -s cuts -1
./bin/AnalyzeSimResonance $input_file -s cuts 1
