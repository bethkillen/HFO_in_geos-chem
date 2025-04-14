#!/bin/tcsh
#PBS -P q90
#PBS -q normal
#PBS -N multi-month
#PBS -l walltime=13:00:00
#PBS -l wd
#PBS -l storage=gdata/m19
#PBS -l ncpus=16
#PBS -l mem=25000MB
#PBS -o ./out.run
#PBS -e ./err.run

limit stacksize unlimited
limit coredumpsize unlimited
limit memorylocked unlimited
setenv OMP_NUM_THREADS $PBS_NCPUS
setenv OMP_PROC_BIND TRUE
setenv OMP_STACKSIZE 2240000

# =============================================================================
#  Daisy-chaining PBS tcsh script (adapted from an NCI script):
#  Resubmits this script after running to run the each month
#
#  USAGE:
#  qsub -v NJOB=$NJOB,NMONTHS=$NMONTHS,STARTYEAR=$STARTYEAR,STARTMONTH=$STARTMONTH run_multimonth.tcsh
#
#  WHERE:
#  $NJOB       is the current job number (defaults to 1 if not included)
#  $NMONTHS    is total number of months you plan to run
#  $STARTYEAR  is the first year you plan to run
#  $STARTMONTH is the first month you plan to run
#  
#  EXAMPLE:
#  qsub -v NJOB=1,NMONTHS=3,STARTYEAR=2015,STARTMONTH=4 run_multimonth.tcsh
#
#  This example would run 3 months sequentially, starting from April 2015.
#
#  REQUIRES:
#  changedate.pl (perl script)
#  input.geos-template (copy of input.geos with all settings as you plan to run them)
#
# =============================================================================

#-------------------------------------------------------------
# Setup variables
#-------------------------------------------------------------
set RUNSCRIPT=run_multimonth.tcsh
set INPUTFILE=geoschem_config.yml-template

# If NJOB is undefined, this must be the first one
if (! $?NJOB) then
   set NJOB=1
endif

#-------------------------------------------------------------
# Update input.geos file with correct dates
#-------------------------------------------------------------
./changedate.pl $INPUTFILE $STARTYEAR $STARTMONTH

#-------------------------------------------------------------
# Set up the log file and run the job
#-------------------------------------------------------------
set logmonth = `printf "%02d" $STARTMONTH`
echo "Running GEOS-Chem for $STARTYEAR $logmonth"
./gcclassic > log.geos_$STARTYEAR$logmonth

#-------------------------------------------------------------
# Check the exit status
#-------------------------------------------------------------
set errstat=$?
if ($errstat != 0) then
    # A brief nap so PBS kills us in normal termination
    # Prefer to be killed by PBS if PBS detected some resource
    # excess
    sleep 5
    echo "Job number $NJOB returned an error status $errstat - stopping job sequence."
    echo "There is now an error file at:"
    echo `date`
    exit $errstat
endif

# Output info to log - may not work properly on Gadi?
#/opt/pbs/default/bin/pbs_rusage $PBS_JOBID | grep 'Usage\|JobId\|Service' | tr -d '\n' | sed 's/$/\n/' >> ./Out_All_Runs.txt

#-------------------------------------------------------------
# Are we in an incomplete job sequence - more jobs to run ?
#-------------------------------------------------------------
if ( $NJOB < $NMONTHS ) then

    #-------------------------------------------------------------
    # Increment counter, years and months and submit the next job:
    #-------------------------------------------------------------
    @ njob=$NJOB
    @ njob++
    set NJOB=$njob
    echo "Submitting job number $NJOB in sequence of $NMONTHS jobs"

    if ( $STARTMONTH < 12 ) then
      @ smonth=$STARTMONTH
      @ smonth++
      set STARTMONTH=$smonth
    else
      @ syear=$STARTYEAR
      @ syear++
      set STARTYEAR=$syear
      set STARTMONTH=01
    endif

    echo "  for new year month: $STARTYEAR $STARTMONTH"

    qsub -v NMONTHS=$NMONTHS,NJOB=$NJOB,STARTYEAR=$STARTYEAR,STARTMONTH=$STARTMONTH $RUNSCRIPT

#-------------------------------------------------------------
# No more jobs
#-------------------------------------------------------------
else
    echo "Finished last job in sequence of $NMONTHS jobs"
endif
