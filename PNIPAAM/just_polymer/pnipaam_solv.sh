#!/bin/bash

scripts_dir=/home/tpascal/scripts/

#make the original box a bit smaller
${scripts_dir}/updateBGFBox.pl -b pnipaam.bgf -c "15 15 120 90 90 90" -s test.bgf

#create an array of 12 P(Nipaam) chains
${scripts_dir}/replicate.pl -b test.bgf -d "4 3 1" -s 12_pnipaam_dry.bgf

#add solvent with lower density to avoid clashes in initial structure
${scripts_dir}/addSolvent.pl -i 12_pnipaam_dry.bgf -f "DREIDING spcew" -n "density:0.9" -w spc -s 12_pnipaam_solv.bgf