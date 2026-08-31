# CBS-tools
A collection of k-mer based tools for sex chromosome identification, assembly, and exploration: the Final Frontier!

The following tools are documented in our manuscript [link] and outline the steps we usually take in exploring dioecious assemblies. However they do not have to be used sequentially. Users can use just CBS-Discovery. Users can skip CBS-Phase and pipe the CBS-Discovery k-mer lists straight into CBS-Frontier. Optionally, if users already have k-mer lists but would like to explore the linked genes and assembly distributions, those non-CBS lists can be imported to CBS-Frontier.

**CBS-Discovery**: Create sex-specific k-mers lists from WGS of multiple isolates
- We recommend starting with 6 isolates of each sex, but not required
- Call ZW or XY system

**CBS-Phase**: Helper tools for quickly mapping CBS-Discovery k-mer lists to genome assembly
- Identify sex-specific contigs to assist with manual phasing
- Optional R plots can be generated to help identify sex-linked contigs

**CBS-Frontier**: Interactive visualization of k-mers mapping to assembly
- Option to estimate the SDR/PAR boundary with changepoint analysis
- Option to import gene/repeat gffs to view k-mer linked annotations
- Please see the CBS-Frontier GitHub for instructions! ([link])
## Installation

### CBS-Discovery
Follow meryl installation instructions: https://github.com/marbl/meryl

#### Example for Linux-amd64 from the meryl GitHub
```
wget https://github.com/marbl/meryl/releases/download/v1.4.2/meryl-1.4.2.Linux-amd64.tar.xz
tar -xJf meryl-1.4.2.Linux-amd64.tar.xz
export PATH=/path/to/meryl-1.4.2/bin:$PATH
```

### CBS-Phase
`conda env create --file CBS-Phase.yaml`

OR

`conda create -n CBS-Phase -c bioconda -c conda-forge samtools bedtools bwa r-base r-optparse r-ggplot2 r-ggpubr r-changepoint r-svglite`

>[!NOTE]
>The R k-mer peak plotting is merely a quick visualization to speed up contig identification. Users can skip install the r-packages if they choose. CBS-Phase output will also include samtools k-mer histogram plots and mapping statistics to help users select contigs for manually phasing that do not require the r-packages.

## Usage

### CBS-Discovery

#### WGS prep
Both long-read PacBio HiFi or short-read Illumina DNAseq files will work with CBS-tools. While we don't notice a huge difference, we recommend trimming up files before running CBS-Discovery with your favorite QC tool like [trimmomatic](https://github.com/usadellab/trimmomatic) or [fastp](https://github.com/opengene/fastp).

```
# example trimmomatic command
trimmomatic PE -threads 8 SRA12345_R1*.gz SRA12345_R2*.gz SRA12345_R1_Q30.fq.gz unpaired/SRA12345_R1_Q30_unpaired.fq.gz SRA12345_R2_Q30.fq.gz unpaired/SRA12345_R2_Q30_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:10:30 MINLEN:40
```

#### Required inputs
| Flag   | Description |
| ------------- | ------------- |
| -m | List of male individual sequence files separated by spaces. Prefix of file may be used to shorten inputs (/path/ID123_reads.fq.gz can be input as /path/ID123). If using Illumina reads, only give the prefix before R1/R2, Discovery will find both!  |
| -f | List of female individual sequence files separated by spaces. Prefix of file may be used to shorten inputs (/path/ID123_reads.fq.gz can be input as /path/ID123). If using Illumina reads, only give the prefix before R1/R2, Discovery will find both!  |

#### Options
| Flag | Description |
| ------------- | ------------- |
| -c | High diversity option. List of pool subset sizes to use for male vs female comparisons (e.g.: -c 3,5,6). Useful for rerunning analysis when all male to all female comparisons have low signal. This will test every unqiue combination of individuals for the input pool sizes. You can also include a file with the count db listed, but generally not recommended unless you have specifc combinations to test. |
| -t | Number of threads to use, DEFAULT: 1 |
|	-g | Memory limit in GB, DEFAULT: 30 |
| -k | K-mer size, DEFAULT: 21 |
| -o | Output directory name, DEFAULT: ./CBS_out |
| -r | Redo/overwrite steps in output directory. CBS Discovery's default behavior won't overwrite k-mer count dbs of the same name in output directory. This will save resources if rerunning jobs with the high diversity option, -c, after viewing all male to all female comparisons, or restarting interrupted jobs. Recommended to use this option if changing input sequencing files. |
>[!NOTE]
>Users who would like to skip testing every n-choose-k combination when using the high diversity option (**-c**) or have specific sub-combinations they'd like to test can choose to upload each combo as a file, using `-c combos_file.txt`. Files should be formatted in two columns, the first being the sex, the second being a list of the meryl k-mer count databases in the Discovery output separated by spaces. (Usually this is an option run after the initial all M vs all F comparison is inconclusive, but if sub-combinations are known ahead of time, these can be formatted as [InputID].[k]mers.count, like SRA12345.21mers.count). For example:
>```
>male  SRR35506473.21mers.count SRR35506472.21mers.count SRR35506471.21mers.count
>male  SRR35506472.21mers.count SRR35506471.21mers.count SRR35506470.21mers.count
>male  SRR35506471.21mers.count SRR35506470.21mers.count SRR35506468.21mers.count
>female  SRR35506469.21mers.count SRR35506467.21mers.count SRR35506466.21mers.count
>female  SRR35506467.21mers.count SRR35506466.21mers.count SRR35506462.21mers.count
>female  SRR35506466.21mers.count SRR35506462.21mers.count SRR35506459.21mers.count
>```

#### Example commands - lots of options!
>[!NOTE]
>Make sure the meryl bin is in the path if not already done when installing meryl!
>`export PATH=/path/to/meryl-1.4.2/bin:$PATH`

```
# export the CBS-tools directory to path
export PATH=$PATH:path/to/CBS_tools/
```

Running the command out of the SRA file locations could look like this:

```
MALE_PREFIXES=(SRR35506473 SRR35506472 SRR35506471 SRR35506470 SRR35506468)
FEMALE_PREFIXES=(SRR35506469 SRR35506467 SRR35506466 SRR35506462 SRR35506459)

CBS_Discovery.sh -k 21 -m ${MALE_PREFIXES[@]} -f ${FEMALE_PREFIXES[@]} -t 12 -g 36 -o /path/to/Discovery_out
```

Storing WGS in different directories could look like this:

```
# set FASTA input directories (leave slashes at the end!)
MALE_DIR=/home/male_sras/
FEMALE_DIR=/outside_collaborator/share/female_sras/

# append the paths to the IDs
MALES=("${MALE_PREFIXES[@]/#/${MALE_DIR}}")
FEMALES=("${FEMALE_PREFIXES[@]/#/${FEMALE_DIR}}")

CBS_Discovery.sh -k 21 -m ${MALES[@]} -f ${FEMALES[@]} -t 12 -g 36 -o /path/to/Discovery_out
```

Even if one file is stored elsewhere this should work:

```
CBS_Discovery.sh -k 21 -m ${MALES[@]} /path/to/other/sra/to/test -f ${FEMALES[@]} -t 12 -g 36 -o /path/to/Discovery_out -r
```

### CBS-Phase

#### Required inputs
| Flag   | Description |
| ------------- | ------------- |
| -d | Directory name of CBS-Discovery results, provide the full path |
| -g | Reference genome assembly fasta(s) to align k-mers against. Multiple haplotypes or genomes can be listed separated by spaces. |
| -k | K-mer size |

#### Mapping options
| Flag | Description |
| ------------- | ------------- |
| -s | [Male|Female] Sex-specific k-mers to map (case sensitive). Can save time by mapping k-mer lists of one sex, for example if heterogametic sex is already known. If results were inconclusive, Phase can map both sets of sex-specific k-mers by omitting this option (DEFAULT). |
|	-t | Number of threads to use, DEFAULT: 1 |

#### Plotting options (if using R-packages)
| Flag | Description |
| ------------- | ------------- |
|	-H | Quick histogram plots to check k-mer peaks for phasing sex-linked contigs between haplotypes. |
| -b | Bin-size used in ggcoverage histogram plots, DEFAULT: 1000000. |
|	-m | Minimum scaffold or contig size to map to in reference assembly, DEFAULT: 2000000. |
| -p | Prefix to add to plot file names, DEFAULT: kmer_hist. |

#### General options
| Flag | Description |
| ------------- | ------------- |
|	-o| Output directory path and name, DEFAULT: ./CBS_Phaser_out. |
| -r | Redo/overwrite results in output directory. CBS Phaser's default behavior won't overwrite files of the same name in output directory. If additional k-mer lists are added to an existing Phaser run, it will skip existing outputs and only run on new inputs. If an error occurred during the job, it would be best to use -r or delete the suspect files manually and rerun without -r. |

#### Example commands

