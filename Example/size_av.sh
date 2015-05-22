#!/bin/sh

rm -f size_av *~

nneurons=100000
ntrials=300
k=10
Tsimul=10000
sigma="1.00d0"
optim="-O3"
fbox=""
file="saving.DAT"
log_file="size_avalanches.log"


dataprog=$( date +"%m-%d-%y"_"%H:%M" )
# echo $dataprog

set -- `getopt N:n:T:f:k:s:o:x $*`
[ $# -lt 1 ] && exit 1
while [ $# -gt 0 ]
do
    case "$1" in
	-N)     nneurons=$2; shift;;
	-n)	ntrials=$2; shift;;
	-T)	Tsimul=$2; shift;;
	-f)	file=$2; shift;;
	-k)     k=$2; shift;;
	-s)	sigma=$2; shift;;
	-o)     optim="-O${2}"; echo "using ${optim} optimization."; shift;;
	-x)	fbox="-fbounds-check"; echo "using bounds check" ; shift;;
	*)	break;;
    esac
    shift
done

#echo $sigma

gfortran $fbox -c $PGH_core_files/*.f90
gfortran -c sa.f90
gfortran $optim -o size_av *.o
chmod +x size_av

rm -f *.o *.~

time ./size_av $RANDOM $nneurons $ntrials $sigma $k $Tsimul $file $dataprog $log_file $pl

mv fort.1 Exp_SizeAvalanches_ErdosRenyi_s$sigma..K$k..N$nneurons..$dataprog.log

mv $file Exp_SizeAvalanches_ErdosRenyi_s$sigma..K$k..N$nneurons..$dataprog.data

rm -f *.o *.~ size_av
