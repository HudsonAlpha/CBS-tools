#!/bin/bash
#SBATCH --job-name="phase"			  # name of this job
#SBATCH -p normal					        # name of the partition (queue)
#SBATCH -N 1 						          # number of nodes
#SBATCH --ntasks=1					      # only really important for arrays
#SBATCH --cpus-per-task=4			    # number of cores/tasks in this job
##SBATCH --array=1					      # job array, n is number of jobs and number of genome references for mapping if choosing to array
#SBATCH --output=CBS_phase_%J.out	# standard output, exact path best, use %J_%a if using arrays
#SBATCH --error=CBS_phase_%J.err	# optional error file, exact path best, use %J_%a if using arrays
######################################################################################################################################
### Environment, modules, path variables
######################################################################################################################################

# path to CBS tools
export PATH=$PATH:YOUR_PATH/CBS_tools/git_src/

# load conda environment
conda activate CBS-Phase

######################################################################################################################################
### File and parameter inputs
######################################################################################################################################

# Reference genome assemblies for mapping
GENOMES=()

# Path to CBS-Discovery results - make sure you have write permissions here!
DISCOVERY_OUT=

# Kmer size
KMER=

# Sex-specific k-mer option
#SEX_KMERS=''

# Output dir
#OUTPUT=

######################################################################################################################################
### Start job
######################################################################################################################################

# If choosing to array instead of listing them all in the same job, use this:
#GENOME=${GENOMES[${SLURM_ARRAY_TASK_ID}-1]}
#CBS_Phase.sh -d ${DISCOVERY_OUT} -g ${GENOME} -t ${SLURM_CPUS_PER_TASK} -k ${KMER}

# using defaults with plotting
CBS_Phase.sh -d ${DISCOVERY_OUT} -g ${GENOMES[@]} -t ${SLURM_CPUS_PER_TASK} -k ${KMER} -H

date "+%c --- Finished!"

# END
