#!/bin/bash

input_file=$1

./bin/AnalyzeSimResonance $input_file
./bin/AnalyzeSimResonance $input_file ptscale 0.995
./bin/AnalyzeSimResonance $input_file ptscale 1.005
./bin/AnalyzeSimResonance $input_file acceptance -1
./bin/AnalyzeSimResonance $input_file acceptance 1
./bin/AnalyzeSimResonance $input_file cuts -1
./bin/AnalyzeSimResonance $input_file cuts 1
