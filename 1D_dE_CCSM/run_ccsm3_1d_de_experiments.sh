!#/bin/bash

input_file='ccsm3_sir_de_input.dat'
output_dir='1D_DE_test'

cp 1D_DE/ccsm3_sir_de_input_fig10bare.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10bare.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10bare1sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10bare1sd.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10bare2sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10bare2sd.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10bareminus1sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10bareminus1sd.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10bareminus2sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10bareminus2sd.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10pond.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10pond.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10pond1sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10pond1sd.dat"

cp 1D_DE/ccsm3_sir_de_input_fig10pondminus1sd.dat $input_file
./a.out
cp ccsm3_sir_de_output.dat "$output_dir/ccsm3_sir_de_output_fig10pondminus1sd.dat"
