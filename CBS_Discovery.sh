#!/bin/bash
#
#Dev notes: is looking for zipped read files currently, does meryl accept non-zipped inputs? 
#meryl should be exported or conda called in batch script 
################################################################################
#                          CBS tools: CBS Discovery                            #
#                      Identifying sex-specific k-mers                         #
#                        In-progress as of 11/05/2025                          #
################################################################################
#                       0. Reading command-line input                          #
################################################################################

# Optional
output='./CBS_Discovery_out'
t_num=1
k=21
mem_all=30
combos=''
redo=0
males=()
females=()

print_usage() {
    echo "Usage: $0 [Options...] -m {MALE FILENAMES} -f {FEMALE FILENAMES}"
    echo
    echo "To display this screen use -h for help"
    echo
    echo "Required:"
    echo "	-m	List of male individual sequence files separated by spaces"
	echo "			Prefix of file may be used to shorten inputs (/path/ID123_reads.fq.gz can be input as /path/ID123)"
	echo "			If using Illumina reads, only give the prefix before R1/R2, Discovery will find both"
	echo
    echo "	-f	List of female individual sequence files separated by spaces"
	echo "			Prefix of file may be used to shorten inputs (/path/ID123_reads.fq.gz can be input as /path/ID123)"
	echo "			If using Illumina reads, only give the prefix before R1/R2, Discovery will find both"
	echo
	echo "Options:"
    echo "	-c	High diversity option: List of pool sizes to use for male vs female comparisons"
	echo "			Example: -c 3,5,6"
	echo "			Useful for rerunning analysis when all male to all female comparisons have low signal;"
	echo "			this will test every unqiue combination of individuals for the input pool sizes."
	echo "			You can also include a file with the count db listed, but generally not recommended"
	echo "			unless you have specifc combinations to test - see README for how to format combination file."
    echo
    echo "	-t	Number of threads to use, DEFAULT: 1"
	echo 
	echo "	-g	Memory limit in GB, DEFAULT: 30"
	echo
	echo "	-k	K-mer size, DEFAULT: 21"
	echo
    echo "	-o  Output directory name, DEFAULT: ./CBS_out"
	echo
	echo "	-r	Redo/overwrite steps in output directory"
	echo "			CBS Discovery's default behavior won't overwrite k-mer count dbs of the same name in output directory."
	echo " 			This will save resources if rerunning jobs with the high diversity option, -c,"
	echo "			after viewing all male to all female comparisons, or restarting interrupted jobs."
	echo "			Recommended to use this option if chaning input sequencing files."
}

date "+%c Starting CBS Discovery"

while getopts ':m:f:c:t:g:k:o:rh' flag
do
    case "${flag}" in
		h) print_usage && exit 1 ;;
        f) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -f requires an argument" && print_usage && exit 1; fi
           females+=("${OPTARG}")
           while [[ ${!OPTIND} && ${!OPTIND} != -* ]]; do females+=("${!OPTIND}") && ((OPTIND++)); done ;;
        m) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -m requires input argument" && print_usage && exit 1; fi
		   if [ -z "${OPTARG//[0-9]/}" ]; then echo "Error: -m must be an integer" && print_usage && exit 1; fi
           males+=("${OPTARG}")
           while [[ ${!OPTIND} && ${!OPTIND} != -* ]]; do males+=("${!OPTIND}") && ((OPTIND++)); done ;;
        c) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -c requires an argument" && print_usage && exit 1; fi
	   	   combos="${OPTARG}" ;;
        t) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -t requires an argument" && print_usage && exit 1; fi
		   if [ -z "${OPTARG//[0-9]/}" ]; then echo "Error: -t must be an integer" && print_usage && exit 1; fi
           t_num="${OPTARG}" ;;
        o) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -o requires an argument" && print_usage && exit 1; fi
           output="${OPTARG}" ;;
		k) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -k requries an argument" && print_usage && exit 1; fi
		   if [ -z "${OPTARG//[0-9]/}" ]; then echo "Error: -k must be an integer" && print_usage && exit 1; fi
		   k="${OPTARG}" ;;
        g) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -g requires input argument" && print_usage && exit 1; fi
           mem_all="${OPTARG}" ;;
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

if [ -z "$females" ]; then echo "Error: No files input for female samples!" && print_usage && exit 1; fi
if [ -z "$males" ]; then echo "Error: No files input for male samples!" && print_usage && exit 1; fi

# Setting path for script locations
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

date "+%c --- Creating output dir: ${output}"
mkdir -p "${output}"
cd "${output}"

date "+%c --- Optional parameters input:"
echo 
echo "threads = $t_num"
echo "k-mer size = $k"
echo "memory = $mem_all"
echo "subset combos = $combos"
echo "redo = $redo"

date "+%c --- Individuals listed as:"
echo "- Males -"
printf '%s\n' "${males[@]}"
echo
echo "- Females -"
printf '%s\n' "${females[@]}"
echo

################################################################################
#                            1. Meryl count                                    #
################################################################################

mkdir male_count_dbs
mkdir female_count_dbs

date "+%c --- Building kmer dbs:"

# males
for input_path in "${males[@]}"
do
	ID=$(basename $input_path)
	date "+%c -- $ID"
	# check to see if count db already exists and how to proceed
	if [ -d "male_count_dbs/${ID}.${k}mers.count" ] && [ "$redo" = 0 ]; then
		echo "${ID}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
		echo "Use -r to overwrite k-mer count dbs"
	elif [ -d "male_count_dbs/${ID}.${k}mers.count" ] && [ "$redo" = 1 ]; then
		echo "${ID}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"

		fasta_path=$(dirname $input_path)
		fasta_path=${fasta_path//$'\r'/}
		fastas=$(find ${fasta_path} -maxdepth 1 -type f \( -name "*${ID}*.fastq.gz" -o -name "*${ID}*.fq.gz" \))
		echo -e "\nFound files:"
		printf '%s\n' "${fastas[@]}"
		echo
        meryl k=$k count threads=$t_num memory=$mem_all output temp_db $fastas
		mv temp_db male_count_dbs/${ID}.${k}mers.count
	else
		fasta_path=$(dirname $input_path)
		fasta_path=${fasta_path//$'\r'/}
		fastas=$(find ${fasta_path} -maxdepth 1 -type f \( -name "*${ID}*.fastq.gz" -o -name "*${ID}*.fq.gz" \))
		echo -e "\nFound files:"
		printf '%s\n' "${fastas[@]}"
		echo
        meryl k=$k count threads=$t_num memory=$mem_all output temp_db $fastas
		mv temp_db male_count_dbs/${ID}.${k}mers.count
	fi
done

# females
for input_path in "${females[@]}"
do
    ID=$(basename $input_path)
	date "+%c --- $ID"
    # check to see if count db already exists
	if [ -d "female_count_dbs/${ID}.${k}mers.count" ] && [ "$redo" = 0 ]; then
		echo "${ID}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
		echo "Use -r to overwrite k-mer count dbs"
	elif [ -d "female_count_dbs/${ID}.${k}mers.count" ] && [ "$redo" = 1 ]; then
		echo "${ID}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"

        fasta_path=$(dirname $input_path)
	    fasta_path=${fasta_path//$'\r'/}
        fastas=$(find "${fasta_path}" -maxdepth 1 -type f \( -name "*${ID}*.fastq.gz" -o -name "*${ID}*.fq.gz" \))
	    echo -e "\nFound files:"
	    printf '%s\n' "${fastas[@]}"
	    echo
        meryl k=$k count threads=$t_num memory=$mem_all output temp_db $fastas
		mv temp_db female_count_dbs/${ID}.${k}mers.count
    else
        fasta_path=$(dirname $input_path)
	    fasta_path=${fasta_path//$'\r'/}
        fastas=$(find "${fasta_path}" -maxdepth 1 -type f \( -name "*${ID}*.fastq.gz" -o -name "*${ID}*.fq.gz" \))
	    echo -e "\nFound files:"
	    printf '%s\n' "${fastas[@]}"
	    echo
        meryl k=$k count threads=$t_num memory=$mem_all output temp_db $fastas
		mv temp_db female_count_dbs/${ID}.${k}mers.count
    fi
done

# overall for meryl difference at the end
# check to see if count dbs already exists and how to proceed
date "+%c -- All males"
if [ -d "allMales.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allMales.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allMales.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allMales.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"

	meryl k=$k union-sum threads=$t_num memory=$mem_all male_count_dbs/* output temp_db
	mv temp_db allMales.${k}mers.count
else
	meryl k=$k union-sum threads=$t_num memory=$mem_all male_count_dbs/* output temp_db
	mv temp_db allMales.${k}mers.count
fi

# check to see if count dbs already exists and how to proceed
date "+%c -- All females"
if [ -d "allFemales.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allFemales.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allFemales.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allFemales.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"

	meryl union-sum threads=$t_num memory=$mem_all female_count_dbs/* output temp_db
	mv temp_db allFemales.${k}mers.count
else
	meryl union-sum threads=$t_num memory=$mem_all female_count_dbs/* output temp_db
	mv temp_db allFemales.${k}mers.count
fi

################################################################################
#                 2. Generate n choose k combinations                          #
################################################################################

# was the combos option used?
if [ -n "$combos" ]; then
	male_combos=()
	female_combos=()
	# check if combo was put in as integer(s) or a file, if not a file then make combos
	if [ -f "$combos" ]; then
		combo_file_pref=$(echo "${combos}" | sed 's/.txt$//')
		cat "${combos}" | grep "^male" | awk -F'\t' '{print $2}' > "${combo_file_pref}_males.txt"
		cat "${combos}" | grep "^female" | awk -F'\t' '{print $2}' > "${combo_file_pref}_females.txt"
		male_combos+=("${combo_file_pref}_males.txt")
		female_combos+=("${combo_file_pref}_females.txt")
	else
		male_counts=( $(find male_count_dbs -mindepth 1 -maxdepth 1 -type d) )
		female_counts=( $(find female_count_dbs -mindepth 1 -maxdepth 1 -type d) )
		for combo in $(echo ${combos} | sed 's/,/ /g')
		do
			date "+%c --- Generating all choose $combo combinations:"
	
			# males
			if [ ${#male_counts[@]} -le $combo ]; then
				echo "Not enough male individuals input for choose $combo combinations"
			else
            	"${SCRIPT_DIR}"/src/generate_nCk_combos.sh $combo "${male_counts[@]}" > combo_list_${combo}_males.txt
            	male_combos+=("combo_list_${combo}_males.txt")

				#echo "Males:"
				#cat combo_list_${combo}_males.txt
			fi

			# females
			if [ ${#female_counts[@]} -le $combo ]; then
				echo "Not enough female individuals input for choose $combo combinations"
			else	
				"${SCRIPT_DIR}"/src/generate_nCk_combos.sh $combo "${female_counts[@]}" > combo_list_${combo}_females.txt
				female_combos+=("combo_list_${combo}_females.txt")

				#echo "Females:"
				#cat combo_list_${combo}_females.txt
			fi
		done
	fi
fi
################################################################################
#                       3. Meryl intersects                                    #
################################################################################
date "+%c --- Compiling intersects:"

if [ -n "$combos" ]; then

	# males
	for list in "${male_combos[@]}"
	do
		# check is really only for file inputs for combos since these can be empty for one sex or the other (or an error when grepping input) 
		if [ -s "$list" ]; then
			mkdir male_intersect_dbs
			if [ -f "$combos" ]; then
				combo="custom"
			else
				# this grep is for combo files made by the script, if a custom file was inserted then the outputs get labeled as such
				combo=$(echo $list | grep -o "[0-9]*")
			fi
			date "+%c -- Working on male choose-${combo} combinations:"
			iter=1
			while read line 	
			do
				date "+%c - Males.intersect.${combo}_${iter}.${k}mers.count"
				if [ -d "male_intersect_dbs/Males.intersect.${combo}_${iter}.${k}mers.count" ] && [ "$redo" = 0 ]; then
					echo "Males.intersect.${combo}_${iter}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
					echo "Use -r to overwrite k-mer count dbs"
				elif [ -d "male_intersect_dbs/Males.intersect.${combo}_${iter}.${k}mers.count" ] && [ "$redo" = 1 ]; then
					echo "Males.intersect.${combo}_${iter}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
					meryl intersect threads=$t_num memory=$mem_all output temp_db ${line}
					mv temp_db male_intersect_dbs/Males.intersect.${combo}_${iter}.${k}mers.count
				else
					meryl intersect threads=$t_num memory=$mem_all output temp_db ${line}
					mv temp_db male_intersect_dbs/Males.intersect.${combo}_${iter}.${k}mers.count
				fi
				((iter++))
			done < ${list}
		else
			echo "WARN: No male combinations grepped from input file - if this is an error please check the combinations input file format!"
		fi
	done

	# females
	for list in "${female_combos[@]}"
	do
		# check is really only for file inputs for combos since these can be empty for one sex or the other (or an error when grepping input) 
		if [ -s "$list" ]; then
			mkdir female_intersect_dbs
			if [ -f "$combos" ]; then
				combo="custom"
			else
				# this grep is for combo files made by the script, if a custom file was inserted then the outputs get labeled as such
				combo=$(echo $list | grep -o "[0-9]*")
			fi
			date "+%c -- Working on female choose-${combo} combinations:"
			iter=1
			while read line
			do
				date "+%c - Females.intersect.${combo}_${iter}.${k}mers.count"
				if [ -d "female_intersect_dbs/Females.intersect.${combo}_${iter}.${k}mers.count" ] && [ "$redo" = 0 ]; then
					echo "Females.intersect.${combo}_${iter}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
					echo "Use -r to overwrite k-mer count dbs"
				elif [ -d "female_intersect_dbs/Females.intersect.${combo}_${iter}.${k}mers.count" ] && [ "$redo" = 1 ]; then
					echo "Females.intersect.${combo}_${iter}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
					meryl intersect threads=$t_num memory=$mem_all output temp_db ${line}
					mv temp_db female_intersect_dbs/Females.intersect.${combo}_${iter}.${k}mers.count
				else
					meryl intersect threads=$t_num memory=$mem_all output temp_db ${line}
					mv temp_db female_intersect_dbs/Females.intersect.${combo}_${iter}.${k}mers.count
				fi
				((iter++))
			done < ${list}
		else
			echo "WARN: No female combinations grepped from input file - if this is an error please check the combinations input file format!"
		fi
	done
fi

# testing overall intersects
# check to see if count dbs already exists and how to proceed
if [ -d "allMales.intersect.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allMales.intersect.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allMales.intersect.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allMales.intersect.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
	meryl intersect threads=$t_num memory=$mem_all output temp_db male_counts_db/*
	mv temp_db allMales.intersect.${k}mers.count
else
	meryl intersect threads=$t_num memory=$mem_all output temp_db male_count_dbs/*
	mv temp_db allMales.intersect.${k}mers.count
fi

# check to see if count dbs already exists and how to proceed
if [ -d "allFemales.intersect.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allFemales.intersect.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allFemales.intersect.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allFemales.intersect.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
	meryl intersect threads=$t_num memory=$mem_all output temp_db female_count_dbs/*
	mv temp_db allFemales.intersect.${k}mers.count
else
	meryl intersect threads=$t_num memory=$mem_all output temp_db female_count_dbs/*
	mv temp_db allFemales.intersect.${k}mers.count
fi

################################################################################
#                        4. Meryl difference                                   #
################################################################################

date "+%c --- Retrieving male/female only kmers:"

# final result directory for all kmer lists
mkdir kmer_lists

if [ -n "$combos" ]; then

	# males
	# check if empty since custom combinations from input files may not have combos for both sexes
	if [ -d male_intersect_dbs ]; then
		mkdir male_diff_dbs
		for intersect in $(ls -d male_intersect_dbs/*)	
		do
			if [ -f "$combos" ]; then
				x=$(basename $intersect | grep -o "custom_[0-9]*")
			else
				x=$(basename $intersect | grep -o "[0-9]*_[0-9]*")
			fi
			date "+%c - MaleOnly.${x}.${k}mers.count"
			if [ -d "male_diff_dbs/MaleOnly.${x}.${k}mers.count" ] && [ "$redo" = 0 ]; then
				echo "MaleOnly.${x}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
				echo "Use -r to overwrite k-mer count dbs"
			elif [ -d "male_diff_dbs/MaleOnly.${x}.${k}mers.count" ] && [ "$redo" = 1 ]; then
				echo "MaleOnly.${x}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
				meryl difference threads=$t_num memory=$mem_all output temp_db ${intersect} allFemales.${k}mers.count
				meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/MaleOnly.${x}.${k}mers.list.txt
				awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/MaleOnly.${x}.${k}mers.list.txt > kmer_lists/MaleOnly.${x}.${k}mers.list.fasta
				mv temp_db male_diff_dbs/MaleOnly.${x}.${k}mers.count
			else
				meryl difference threads=$t_num memory=$mem_all output temp_db ${intersect} allFemales.${k}mers.count
				meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/MaleOnly.${x}.${k}mers.list.txt
				awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/MaleOnly.${x}.${k}mers.list.txt > kmer_lists/MaleOnly.${x}.${k}mers.list.fasta
				mv temp_db male_diff_dbs/MaleOnly.${x}.${k}mers.count
			fi
		done
	fi
	
	# females
	# check if empty since custom combinations from input files may not have combos for both sexes
	if [ -d female_intersect_dbs ]; then
		mkdir female_diff_dbs
		for intersect in $(ls -d female_intersect_dbs/*)
		do
			if [ -f "$combos" ]; then
				x=$(basename $intersect | grep -o "custom_[0-9]*")
			else
				x=$(basename $intersect | grep -o "[0-9]*_[0-9]*")
			fi
			date "+%c - FemaleOnly.${x}.${k}mers.count"
			if [ -d "female_diff_dbs/FemaleOnly.${x}.${k}mers.count" ] && [ "$redo" = 0 ]; then
				echo "FemaleOnly.${x}.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
				echo "Use -r to overwrite k-mer count dbs"
			elif [ -d "female_diff_dbs/FemaleOnly.${x}.${k}mers.count" ] && [ "$redo" = 1 ]; then
				echo "FemaleOnly.${x}.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
				meryl difference threads=$t_num memory=$mem_all output temp_db "${intersect}" allMales.${k}mers.count
				meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/FemaleOnly.${x}.${k}mers.list.txt
				awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/FemaleOnly.${x}.${k}mers.list.txt > kmer_lists/FemaleOnly.${x}.${k}mers.list.fasta
				mv temp_db female_diff_dbs/FemaleOnly.${x}.${k}mers.count
			else
				meryl difference threads=$t_num memory=$mem_all output temp_db "${intersect}" allMales.${k}mers.count
				meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/FemaleOnly.${x}.${k}mers.list.txt
				awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/FemaleOnly.${x}.${k}mers.list.txt > kmer_lists/FemaleOnly.${x}.${k}mers.list.fasta
				mv temp_db female_diff_dbs/FemaleOnly.${x}.${k}mers.count
			fi
		done
	fi
fi

# difference and print of overall intersects
# check to see if count dbs already exists and how to proceed
if [ -d "allMaleOnly.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allMaleOnly.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allMaleOnly.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allMaleOnly.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
	meryl difference threads=$t_num memory=$mem_all output temp_db allMales.intersect.${k}mers.count allFemales.${k}mers.count
	meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/allMaleOnly.${k}mers.list.txt
	awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/allMaleOnly.${k}mers.list.txt > kmer_lists/allMaleOnly.${k}mers.list.fasta
	mv temp_db allMaleOnly.${k}mers.count
else
	meryl difference threads=$t_num memory=$mem_all output temp_db allMales.intersect.${k}mers.count allFemales.${k}mers.count
	meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/allMaleOnly.${k}mers.list.txt
	awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/allMaleOnly.${k}mers.list.txt > kmer_lists/allMaleOnly.${k}mers.list.fasta
	mv temp_db allMaleOnly.${k}mers.count
fi

# check to see if count dbs already exists and how to proceed
if [ -d "allFemaleOnly.${k}mers.count" ] && [ "$redo" = 0 ]; then
	echo "allFemaleOnly.${k}mers.count found in output directory! CBS Discovery will continue on to the next individual."
	echo "Use -r to overwrite k-mer count dbs"
elif [ -d "allFemaleOnly.${k}mers.count" ] && [ "$redo" = 1 ]; then
	echo "allFemaleOnly.${k}mers.count found in output directory! WARNING: redo (-r) set - CBS Discovery will overwrite this count db!"
	
	meryl difference threads=$t_num memory=$mem_all output temp_db allFemales.intersect.${k}mers.count allMales.${k}mers.count
	meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/allFemaleOnly.${k}mers.list.txt
	awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/allFemaleOnly.${k}mers.list.txt > kmer_lists/allFemaleOnly.${k}mers.list.fasta
	mv temp_db allFemaleOnly.${k}mers.count
else
	meryl difference threads=$t_num memory=$mem_all output temp_db allFemales.intersect.${k}mers.count allMales.${k}mers.count
	meryl print temp_db threads=$t_num memory=$mem_all > kmer_lists/allFemaleOnly.${k}mers.list.txt
	awk 'NF {print ">kmer"NR"\n"$1}' kmer_lists/allFemaleOnly.${k}mers.list.txt > kmer_lists/allFemaleOnly.${k}mers.list.fasta
	mv temp_db allFemaleOnly.${k}mers.count
fi

date  "+%c --- CBS-Discovery complete!"
echo -e "\nUnique k-mer counts:\n"
wc -l kmer_lists/*.txt | head -n -1 | sed 's/kmer_lists\///' | sed 's/.21mers.list.txt//'; echo

################################################################################
#                               End!                                           #
################################################################################

