#!/bin/bash

input_file='ccsm3_sir_de_input.dat'
output_file='ccsm3_sir_de_output.dat'
output_dir='1D_DE_test'

experiments="fig10bare\
 fig10bare1sd\
 fig10bare2sd\
 fig10bareminus1sd\
 fig10bareminus2sd\
 fig10pond\
 fig10pond1sd\
 fig10pondminus1sd\
 fig11bare15mminus1sd\
 fig11bare15mminus2sd\
 fig11bare15mplus1sd\
 fig11bare15mplus2sd\
 fig12mp1sd\
 fig12mpminus1sd\
 fig7bareice\
 fig7ponded\
 fig8bare15m\
 fig8bare1m\
 fig9pond"

for experiment in $experiments
do
    exp_input_file="1D_DE/ccsm3_sir_de_input_$experiment.dat"
    exp_output_file="$output_dir/ccsm3_sir_de_output_$experiment.dat"
    echo "Running $exp_input_file -> $exp_output_file"
    cp $exp_input_file $input_file
    ./a.out
    cp $output_file $exp_output_file
done

ls -ltr $output_dir
