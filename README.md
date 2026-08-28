# Installation

## CBS-Discovery
Follow meryl installation instructions: https://github.com/marbl/meryl

#### Example for Linux-amd64
wget https://github.com/marbl/meryl/releases/download/v1.4.2/meryl-1.4.2.Linux-amd64.tar.xz
tar -xJf meryl-1.4.2.Linux-amd64.tar.xz
export PATH=/path/to/meryl-1.4.2/bin:$PATH

## CBS-Phase:
conda env create --file CBS-Phase.yaml

OR

conda create -n CBS-Phase -c bioconda -c conda-forge samtools bedtools bwa r-base r-optparse r-ggplot2 r-ggpubr r-changepoint r-svglite
