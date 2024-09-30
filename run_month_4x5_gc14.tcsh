#!/bin/tcsh
#PBS -P q90
#PBS -q normal
#PBS -N month_test
#PBS -l walltime=12:00:00
#PBS -l wd
#PBS -l storage=gdata/m19+scratch/m19
#PBS -l ncpus=8
#PBS -l mem=40000MB
#PBS -o ./out.run
#PBS -e ./out.run

limit stacksize unlimited
limit coredumpsize unlimited
limit memorylocked unlimited
setenv OMP_NUM_THREADS $PBS_NCPUS
setenv OMP_STACKSIZE 2240000

  #-------------------------------------------------------------
  # Run the new file
  #-------------------------------------------------------------
  ./gcclassic > log.geos

  #-------------------------------------------------------------
  # Check the exit status
  #-------------------------------------------------------------
  set errstat=$?
  if ($errstat != 0) then
      # A brief nap so PBS kills us in normal termination
      # Prefer to be killed by PBS if PBS detected some resource
      # excess
      sleep 5
      echo "Job returned an error status $errstat - stopping job sequence."
      echo "There is now an error file at:"
      echo `date`
      exit $errstat
  endif


echo "Finished job "

