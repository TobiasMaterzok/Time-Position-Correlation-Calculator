#!/bin/bash

# interpolate the input files to the same update frequency
function prepare_file(){
    awk '{if((substr($1,0,1) != "#") && (substr($1,0,1) != "@")) print $0;}' "$1" > tmp_b.xvg
    awk -v "tj=$2" '{if($1%tj==0) print $0}' tmp_b.xvg > "$1"
    rm tmp_b.xvg
}

UPDATE_FREQUENCY_POSITION_PROPERTY=20

# prepare the input files
prepare_file run_pullx.xvg "$UPDATE_FREQUENCY_POSITION_PROPERTY"
prepare_file energy.xvg "$UPDATE_FREQUENCY_POSITION_PROPERTY"

# calculate the correlation
python time-pos-corr-calc.py -expl 1 -ncol_d 1 -ncol_p 1 -bw 0.33 -f energy.xvg -x run_pullx.xvg -o energy_position.xvg
