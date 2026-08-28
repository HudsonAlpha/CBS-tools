#!/bin/bash
#SBATCH --job-name="discovery"			# name of this job
#SBATCH -p normal						# name of the partition (queue)
#SBATCH -N 1 							# number of nodes
#SBATCH --ntasks=1						# only really important for arrays
#SBATCH --cpus-per-task=64				# number of cores/tasks in this job
##SBATCH -t 						    # time allocated for this job days-hours:mins:sec
#SBATCH --mem=300G
##SBATCH --array=0-n					# job array, n is number of jobs
#SBATCH --output=CBS_discovery_%J.out	# standard output, exact path best
#SBATCH --error=CBS_discovery_%J.err	# optional error file
######################################################################################################################################
### Environment, modules, path variables
######################################################################################################################################

# path to meryl
export PATH=$PATH:/YOUR_PATH/meryl-1.4.1/bin
# path to CBS tools
export PATH=$PATH:/YOUR_PATH/CBS_tools/git_src/

######################################################################################################################################
### File and parameter inputs
######################################################################################################################################
#put slash at the end if this is a directory
FASTA_DIR=

# list prefixes in the arrays below - sometimes files aren't named with M/F identifiers for grepping, so this is best
 
MALE_PREFIXES=()
FEMALE_PREFIXES=()

# full path to file name or integers separated by commas
COMBOS=''

OUTPUT=''

######################################################################################################################################
### Start job
######################################################################################################################################

MALES=("${MALE_PREFIXES[@]/#/${FASTA_DIR}}")
FEMALES=("${FEMALE_PREFIXES[@]/#/${FASTA_DIR}}")
MEM_GB=$(echo $SLURM_MEM_PER_NODE | awk '{print $1/1024}')

date "+%c jobID: $SLURM_JOBID   mem_detected: $MEM_GB"

#default run
CBS_Discovery.sh -k 21 -m ${MALES[@]} -f ${FEMALES[@]} -t $SLURM_CPUS_PER_TASK -g $MEM_GB -o ${OUTPUT}

# high diversity check
#CBS_Discovery.sh -k 21 -m ${MALES[@]} -f ${FEMALES[@]} -t $SLURM_CPUS_PER_TASK -g $MEM_GB -o ${OUTPUT} -c ${COMBOS} 


date "+%c --- Finished!"

# END 
