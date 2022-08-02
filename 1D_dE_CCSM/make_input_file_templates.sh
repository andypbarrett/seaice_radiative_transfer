!#/bin/bash

template_file='ccsm3_sir_de_input.dat.sav'

for fp in `ls 1D_DE/ccsm3_sir_de_output_*.dat`
do
    dirpath=`dirname $fp`
    stem=`basename $fp .dat`
    tag=${stem##*_}
    input_file="$dirpath/ccsm3_sir_de_input_$tag.dat"
    echo "Copying $fp to $input_file"
    cp $template_file $input_file
done
