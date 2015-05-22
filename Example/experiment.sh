#!/bin/bash

# This is a small script to generate a few avalanches and plot their
# distributions. Refer to Readme.md for more details


##
# The following lines will generate $navalanches avalanches with
# several avg. branching ratios (refer to paper).
navalanches=20000

# Subcritical case:
./size_av.sh -f experiment_s08.dat -s 0.80 -n $navalanches

# Critical case
./size_av.sh -f experiment_s10.dat -s 1.00 -n $navalanches


##
# Splitting avalanche sizes and times

cp $PGH_gen_programs/SplitAvalanches.py ./
python SplitAvalanches.py                      ## Change here to pypy if you want better performance!
rm SplitAvalanches.py


##
# Next we construct the distributions of both avalanches 

# Compiling the current version of histo_create
gfortran $PGH_gen_programs/histo_create.f90 -o histo_create

file='Size_Exp_SizeAvalanches_ErdosRenyi_s0.80..K10..N100000.data'
./histo_create linear $file 2000 pdf_experiment_s08.dat
./histo_create log-cumulative $file 150 cdf_experiment_s08.dat

file='Size_Exp_SizeAvalanches_ErdosRenyi_s1.00..K10..N100000.data'
./histo_create linear $file 5000 pdf_experiment_s10.dat
./histo_create log-cumulative $file 150 cdf_experiment_s10.dat


##
# Plotting 
gnuplot plot.gp


# Removing histo_create from here
rm -f histo_create
