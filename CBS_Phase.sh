#!/bin/bash
#
# requires samtools and bwa-mem, R plotting is optional
################################################################################
#                           CBS tools: CBS Phase.                             #
#         Sex chromosome phasing between haplotypes of genome assembly         #
################################################################################
#                       0. Reading command-line input                          #
################################################################################

#Optional
output='./CBS_Phase_out'
t_num=1
ref_genomes=()
redo=0
sex_kmers=''
histogram_plots=0
bin_size=1000000
min_size=2000000
plot_prefix="kmer_hist"

print_usage() {
	echo "Usage: $0 [Options...] -d {CBS DISCOVERY DIR} -g {REF GENOME FASTAS}"
	echo
	echo "To display this screen use -h for help"
	echo
	echo "Required:"
	echo "	-d  Directory name of CBS Discovery results, provide the full path"
	echo
	echo "  -g  Reference genome assembly fasta(s) to align k-mers against"
	echo "			Multiple haplotypes or genomes can be listed separated by spaces."
	echo
	echo "	-k	K-mer size"

	echo "Mapping Options:"
	echo "	-s  [Male|Female] Sex-specific k-mers to map (case sensitive)"
	echo "			Ideally set to the heterogametic sex determined with Discovery."
	echo "			If results were inconclusive, Phase can map both sets of sex-specific k-mers by omitting this option (DEFAULT)."
	echo
	echo "	-t	Number of threads to use, DEFAULT: 1"
	echo 
	echo
	echo "Plotting Options:"
	echo "	-H	Quick histogram plots to check k-mer peaks for phasing sex-linked contigs between haplotypes"
	echo
	echo "	-b	Bin-size used in ggcoverage histogram plots, DEFAULT: 1000000"
	echo
	echo "	-m	Minimum scaffold or contig size to map to in reference assembly, DEFAULT: 2000000"
	echo
	echo "	-p	Prefix to add to plot file names, DEFAULT: kmer_hist"
	echo
	echo "General Options:"
	echo "	-o  Output directory path and name, DEFAULT: ./CBS_Phase_out"
	echo
	echo "	-r	Redo/overwrite results in output directory"
	echo "			CBS Phase's default behavior won't overwrite files of the same name in output directory."
	echo "			If additional k-mer lists are added to an existing Phase run, it will skip existing outputs and only run on new inputs." 
	echo "			If an error occured during the job, it would be best to use -r or delete the suspect files manually and rerun without -r."
}

date "+%c Starting CBS Phase"

while getopts ':d:g:s:t:k:b:m:p:o:rHh' flag
do
    case "${flag}" in
		h) print_usage && exit 1 ;;
		d) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -d requires an argument" && print_usage && exit 1; fi
		   cbs_meryl_out="${OPTARG}" ;;
		g) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -g requires an argument" && print_usage && exit 1; fi
			ref_genomes+=("${OPTARG}")
		   while [[ ${!OPTIND} && ${!OPTIND} != -* ]]; do ref_genomes+=("${!OPTIND}") && ((OPTIND++)); done ;;
		s) if [ "$OPTARG" != "Male" ] && [ "$OPTARG" != "Female" ]; then 
				echo "Error: -s must be set to Male or Female (case sensitive)! Omit -s from command to run both" && print_usage && exit 1; fi
		   sex_kmers="${OPTARG}" ;;
		t) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -t requires an argument" && print_usage && exit 1; fi
			t_num="${OPTARG}" ;;
		k) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -k requries an argument" && print_usage && exit 1; fi
		   if ! [ "$OPTARG" -eq "$OPTARG" ] 2>/dev/null; then echo "Error: -k must be an integer" && print_usage && exit 1; fi
		   k="${OPTARG}" ;;
		H) histogram_plots=1 ;;
		b) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -b requries an argument" && print_usage && exit 1; fi
		   if ! [ "$OPTARG" -eq "$OPTARG" ] 2>/dev/null; then echo "Error: -b must be an integer" && print_usage && exit 1; fi
		   bin_size="${OPTARG}" ;;
		m) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -m requries an argument" && print_usage && exit 1; fi
		   if ! [ "$OPTARG" -eq "$OPTARG" ] 2>/dev/null; then echo "Error: -m must be an integer" && print_usage && exit 1; fi
		   min_size="${OPTARG}" ;;
		p) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -p requries an argument" && print_usage && exit 1; fi
		   plot_prefix="${OPTARG}" ;;
		o) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -o requires an argument" && print_usage && exit 1; fi
			output="${OPTARG}" ;;
		r) redo=1 ;;
		\?) echo "Error: Invalid option -$OPTARG"
			print_usage
			exit 1 ;;
		:)  echo "option -$OPTARG requires an argument"
			print_usage
			exit 1 ;;
		*) print_usage
			exit 1;;
    esac
done

if [ ! -d "$cbs_meryl_out" ]; then echo "Error: directory ${cbs_meryl_out} doesn't exist!" && print_usage && exit 1; fi
if [ "${#ref_genomes[@]}" -eq 0 ]; then echo "Error: No input for reference genome fastas!" && print_usage && exit 1; fi

# Setting path for script locations
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

################################################################################
#                          1. Mapping loops                                    #
################################################################################

date "+%c --- Creating working dir: ${output}"
mkdir -p "${output}"

date "+%c --- Optional parameters input:"
echo 
echo "threads = $t_num"
echo "k-mer size = $k"
echo "memory = $mem_all"
echo "subset combos = $combos"
echo "redo = $redo"
echo "sex kmers = $sex_kmers (empty will run both)"
echo "histogram plot = $histogram_plots"
echo "plot prefix = $plot_prefix"
echo "bin size = $bin_size"
echo "min chr length = $min_size"
echo

date "+%c --- Working from CBS Discovery output dir: ${cbs_meryl_out}"
cd "${cbs_meryl_out}"

for assem in "${ref_genomes[@]}"; 
do
	# index assembly, faidx for quick visualization
	samtools faidx ${assem}

	if [ -s "${assem}.amb" ] && [ -s "${assem}.ann" ] && [ -s "${assem}.bwt" ] && [ -s "${assem}.pac" ] && [ -s "${assem}.sa" ] && [ "$redo" = 0 ]; then
		echo "Existing BWA index files found for ${assem}! Skipping to BWA MEM"
		echo "Use -r to overwrite BWA index files"
	elif [ -s "${assem}.amb" ] && [ -s "${assem}.ann" ] && [ -s "${assem}.bwt" ] && [ -s "${assem}.pac" ] && [ -s "${assem}.sa" ] && [ "$redo" = 1 ]; then
		echo "Existing BWA index files found for ${assem}! WARNING: redo (-r) set - CBS Phase will overwrite these files!"
		bwa index ${assem}
	else
		bwa index ${assem}
	fi	
	
	# will name directories for each input genome assembly within the output directory
	hap=$(basename $assem)
	mkdir -p "${output}"/"${hap}"

	date "+%c --- Making windows out of ${hap}"

	awk -F '\t' -v OFS='\t' '{print $1, $2}' ${assem}.fai > ${output}/${hap}.chromSizes
	slide=$(expr $bin_size / 2) 
	bedtools makewindows -g ${output}/${hap}.chromSizes -w ${bin_size} -s ${slide} > ${output}/${hap}.${bin_size}.s${slide}.windows.bed
	bedtools makewindows -g ${output}/${hap}.chromSizes -n 100 > ${output}/${hap}.n100.windows.bed

	date "+%c --- Mapping to ${hap} with:"

	# only grepping the male/female only kmer databases
	for cf in $(basename -a kmer_lists/*${sex_kmers}*list.fasta | sed 's/.fasta$//')
	do
		if [ -s "kmer_lists/${cf}.fasta" ]; then
			date "+%c ----- $cf"

			if [ -s "${output}/${hap}/${cf}_scaffolds.bam" ] && [ "$redo" = 0 ]; then
				echo "Existing k-mer mapping depth file found! Looks like Phase successfully finished mapping this k-mer list, moving on to next input"
				echo "Use -r to redo k-mer mapping"
			elif [ -s "${output}/${hap}/${cf}_scaffolds.bam" ] && [ "$redo" = 1 ]; then
				echo "Existing k-mer mapping depth file found! WARNING: redo (-r) set - CBS Phase will overwrite mapping output for this k-mer list!"
				bwa mem -t $t_num -k $k -T $k -a -c 10 $assem kmer_lists/$cf.fasta > ${output}/${hap}/${cf}_scaffolds.sam
				samtools sort -@ $t_num ${output}/${hap}/${cf}_scaffolds.sam -o ${output}/${hap}/${cf}_scaffolds.bam
				rm ${output}/${hap}/${cf}_scaffolds.sam

			else
				bwa mem -t $t_num -k $k -T $k -a -c 10 $assem kmer_lists/$cf.fasta > ${output}/${hap}/${cf}_scaffolds.sam
				samtools sort -@ $t_num ${output}/${hap}/${cf}_scaffolds.sam -o ${output}/${hap}/${cf}_scaffolds.bam
				rm ${output}/${hap}/${cf}_scaffolds.sam
			fi

			samtools index ${output}/${hap}/${cf}_scaffolds.bam
			samtools flagstat ${output}/${hap}/${cf}_scaffolds.bam > ${output}/${hap}/${cf}_scaffolds.bam.flagstat

			samtools coverage ${output}/${hap}/${cf}_scaffolds.bam -o ${output}/${hap}/${cf}_scaffolds.${bin_size}.s${slide}.100bp_windows.cov -m -w 100
			bedtools coverage -g ${output}/${hap}.chromSizes -sorted -b ${output}/${hap}/${cf}_scaffolds.bam -a ${output}/${hap}.${bin_size}.s${slide}.windows.bed > ${output}/${hap}/${cf}_scaffolds.${bin_size}.s${slide}.cov
			bedtools coverage -g ${output}/${hap}.chromSizes -sorted -mean -b ${output}/${hap}/${cf}_scaffolds.bam -a ${output}/${hap}.${bin_size}.s${slide}.windows.bed > ${output}/${hap}/${cf}_scaffolds.${bin_size}.s${slide}.depth.txt
			bedtools coverage -g ${output}/${hap}.chromSizes -sorted -mean -b ${output}/${hap}/${cf}_scaffolds.bam -a ${output}/${hap}.n100.windows.bed > ${output}/${hap}/${cf}_scaffolds.n100.windows.depth.txt

		else
			echo "WARNING: ${cf}.fasta is empty! Skipping k-mer mapping." 
			echo "If this is an error, check that your file pathing is correct."
		fi
	done
	# quick histogram plot to help make calls on contig phasing before final Frontier viewing/exploration
	# note, this will require an R environmnet with required packages, see README
	if [ "$histogram_plots" = 1 ]; then
		date "+%c --- Creating quick plot to view k-mer peaks on ${hap} ---"
		echo
		Rscript "${SCRIPT_DIR}"/src/Phase_kmer_plotting.R -r "${output}/${hap}.chromSizes" -p "${output}/${hap}" -o "${output}/${plot_prefix}.${hap}" -b "${bin_size}" -m ${min_size}
	fi
done

# END
