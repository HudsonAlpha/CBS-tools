#!/bin/bash
#SBATCH --job-name="phaser"			# name of this job
#SBATCH -p normal					# name of the partition (queue)
#SBATCH -N 1 						# number of nodes
#SBATCH --ntasks=1					# only really important for arrays
#SBATCH --cpus-per-task=24			# number of cores/tasks in this job
##SBATCH -t 						# time allocated for this job days-hours:mins:sec
#SBATCH --mem=50G
##SBATCH --array=1					# job array, n is number of jobs and number of genome references for mapping if choosing to array
#SBATCH --output=CBS_phaser_%J.out	# standard output, exact path best, use %J_%a if using arrays
#SBATCH --error=CBS_phaser_%J.err	# optional error file, exact path best, use %J_%a if using arrays
######################################################################################################################################
### Environment, modules, path variables
######################################################################################################################################

# path to CBS tools
export PATH=$PATH:YOUR_PATH/CBS_tools/git_src/

# load the R environment
source /cluster/home/lwhitt/.bashrc
conda activate CBS_Phaser

# modules or conda environment initialization here:
#module load cluster/bwa/0.7.17
#module load samtools/1.19.2
#module load cluster/bedtools/2.28.0 

######################################################################################################################################
### File and parameter inputs
######################################################################################################################################

# Reference genome assemblies for mapping
GENOMES=()

# Path to previous step of the pipeline: CBS_meryl_pipeline - make sure you have write permissions here!
DISCOVERY_OUT=

# Sex-specific k-mer option
SEX_KMERS=''

# Kmer size
#KMER=

# Output dir
#OUTPUT=

######################################################################################################################################
### Start job
######################################################################################################################################

# If choosing to array instead of listing them all in the same job, use this:
#GENOME=${GENOMES[${SLURM_ARRAY_TASK_ID}-1]}
#CBS_Phaser.sh -d ${DISCOVERY_OUT} -g ${GENOME} -t ${SLURM_CPUS_PER_TASK} -s ${SEX_KMERS} -o ${OUTPUT}

# using default kmer size and output
CBS_Phaser_plottingopts_included.sh -d ${DISCOVERY_OUT} -g ${GENOMES[@]} -t ${SLURM_CPUS_PER_TASK} -s ${SEX_KMERS} -o ${OUTPUT} -H

date "+%c --- Finished!"

# END
