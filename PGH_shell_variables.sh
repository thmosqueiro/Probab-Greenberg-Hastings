#
# This shell script defines a set of variables to make this
# PGH-project more customizable.
# 
# Run this script before to define these variables.
# 

# Getting current location
unset PGH_root
PGH_root=$(pwd)
export PGH_root

# Defining paths
unset PGH_core_files
export PGH_core_files=$PGH_root'/core'

unset PGH_gen_programs
export PGH_gen_programs=$PGH_root'/gen_programs'
