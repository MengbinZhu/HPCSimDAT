#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N fj_data
#$ -pe impi 1
#$ -V
#$ -q q2.q
#source /disk1/soft/intel/impi/3.2.0.011/bin64/mpivars.sh
#mpirun -r ssh -np 1 /disk1/home/ljp/fengjie/attr/lorenz/localr/proba/a1.out
#mpirun -r ssh -np 1 /disk1/home/ljp/fengjie/QG_EPS/a.out
#mpirun -r ssh -np 1 /disk1/home/ljp/fengjie/assess/NLLE/a1.out
  ifort -o fengjie.out qgcont.f90
  nohup ./fengjie.out >fj_data.out &
  echo "qsub ok"
#mpirun -r ssh -np 1 /disk1/home/ljp/fengjie/assess/analysis/fengjie.out
