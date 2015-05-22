#!/bin/bash

# This is a small script to generate a few avalanches and plot their
# distributions. Refer to Readme.md for more details


##
# The following lines will generate $navalanches avalanches with
# several avg. branching ratios (refer to paper).
navalanches=100

# Subcritical case:
./size_av.sh -f experiment_s08.dat -s 0.80 -n $navalanches

# Critical case
./size_av.sh -f experiment_s10.dat -s 1.00 -n $navalanches



##
# Next we construct the distributions of both avalanches 

# Compiling the current version of histo_create
gfortran $PGH_gen_programs/histo_create.f90 -o histo_create

./histo_create linear experiment_s08.dat 1000 pdf_experiment_s08.dat
./histo_create log-cumulative experiment_s08.dat 150 cdf_experiment_s08.dat

./histo_create linear experiment_s10.dat 1000 pdf_experiment_s10.dat
./histo_create log-cumulative experiment_s10.dat 150 cdf_experiment_s10.dat


##
# Plotting 



# Removing histo_create from here
rm -f histo_create
