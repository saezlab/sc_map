#!/bin/sh
# Resources, e.g. a total time of 15 hours...
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
# Ensures that the Linux environment for the job is the same as the one we're working in:
#PBS -V
# stderr redirection
#PBS -e ./results/panc8/test.err
# stdout redirection
#PBS -o ./results/panc8/test.log

source ~/.bashrc
conda activate scmap

# Change directory to the project
if [ ! -z ${PBS_O_WORKDIR+x} ];then
  cd ${PBS_O_WORKDIR};
fi


Rscript ./src/sc_mapping.R --S1 ./data/panc8/fluidigmc1.rds --S2 ./data/panc8/smartseq2.rds --METHOD $METHOD --METRIC $METRIC --N $NGENES --OUTDIR ./results/panc8

