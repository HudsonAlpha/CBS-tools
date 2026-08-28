#!/usr/bin/env Rscript
### testing ggcoverage

#install.packages("remotes")
#remotes::install_github("showteeth/ggcoverage")

rm(list=ls())

library(optparse)
library(ggplot2)
library(gridExtra)
library(ggcoverage)

# command line options for use on the server

options_list <- list(
  make_option(c("-r", "--ref_file"), type = "character", default = NULL,
              help = "Path to FAI samtools file for assembly reference", metavar = "FILE"),
  make_option(c("-p", "--phaser_path"), type = "character", default = NULL,
              help = "Path to Phaser output folder for assembly reference", metavar = "PATH"),
  make_option(c("-b", "--bin_size"), type="integer", default=100000, 
              help="Histogram bin size [default %default]", metavar="INTEGER"),
  make_option(c("-m", "--min_chr_size"), type="integer", default=2000000, 
              help="Minimum contig/scaffold size to show in plot [default %default]", metavar="INTEGER"),
  make_option(c("-o", "--outpref"), type = "character", default = 'kmer_hist',
              help = "Output file name prefix [default %default]", metavar = "STRING")
)

# Create the parser object
input_opts <- OptionParser(option_list = options_list)
# Parse the arguments
opts_parsed <- parse_args(input_opts)
# print arguments
print("Options input:")
str(opts_parsed)

# Check for required options
if (is.null(opts_parsed$ref_file)) {
  print_help(input_opts)
  stop("--ref_file is required.", call.=FALSE)
}
if (is.null(opts_parsed$phaser_path)) {
  print_help(input_opts)
  stop("--phaser_path is required.", call.=FALSE)
}

################################################################################
# Inputs
################################################################################

# prefix name for file output
#outpref <- "Ursus_only_chr"

# CBS-Phaser output folder name with bam/bam.bai files
#track_folder <- "GCA_023065955.2_UrsArc2.0_genomic.fna"

# This is our recommended starting size but can be tweaked depending on genome size
#bin_size <- 100000

# Minimum chr size for plotting, if you want to show mappings for all contigs or
# fragments (like if comparing to unscaffolded assembly) then set to 1
#min_chr_size <- 2000000

################################################################################
# Run scripts below for output, shouldn't need to change anything below 
################################################################################

# this index file should be coming straight from Phaser as well and should already match the folder name
#ref_file <- list.files("./",pattern = paste0(track_folder,".fai"), full.names = T)

# Pulls in samtools faidx index file file for scaffold names and lengths

ref_index <- read.delim(file = opts_parsed$ref_file, header = F, sep = "\t")

# optional size subset
ref_index <- ref_index[ref_index$V2 > opts_parsed$min_chr_size,]
print(head(ref_index))

meta_file <- data.frame(SampleName=gsub(".bam$", "", list.files(path = opts_parsed$phaser_path, pattern = "*.bam$")), 
                        Type=gsub(".*(M|F).*aleOnly.([0-9]_[0-9]*|custom_[0-9]*).*", "\\1 \\2", list.files(path = opts_parsed$phaser_path, pattern = "*.bam$")), 
                        Group=gsub(".*(M|F).*aleOnly.*", "\\1", list.files(path = opts_parsed$phaser_path, pattern = "*.bam$")))

# removing custom for plotting facet labels
meta_file$Type <- gsub("custom_", "", meta_file$Type)
meta_file$Type[grep("all", meta_file$SampleName)] <- paste("all", meta_file$Group[grep("all", meta_file$SampleName)])
print(meta_file)

PlotList <- apply(ref_index, 1, function(x){
  print(paste0(x[[1]],":1-",x[[2]]))
  track_df <- LoadTrackFile(
    track.folder = opts_parsed$phaser_path,
    meta.info = meta_file,
    format = "bam",
    region = paste0(x[[1]],":1-",x[[2]]),
    norm.method = "None",
    bin.size = opts_parsed$bin_size
  )
  
  #might be good for full run:
  bincoverage <- ggcoverage(
    data = track_df,
    color = c("#356920","#8A4F21"),
    mark.region = NULL,
    range.position = "out",
    plot.type = "facet",
    facet.key = "Type",
    group.key = "Group",
    facet.y.scale = "fixed"
    #facet.order = 
    #joint.avg = TRUE
    #plot.space = ""
  ) + 
    scale_y_continuous(limits = c(0,10), breaks = c(5,10)) +
    scale_x_continuous(limits = c(1,as.integer(x[[2]]))) +
    labs(title = x[[1]], subtitle = paste0("bin size = ", opts_parsed$bin_size, " bp")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(bincoverage)
})

allPlots <- marrangeGrob(grobs = PlotList, ncol = 1, nrow = 1)

# attempt to make plot height automatic for printing, but this may need to be adjusted
plot_h <- nrow(meta_file)/2 + 1

# is there a way to calculate best height so it can be automatic?
ggsave(paste0(opts_parsed$outpref, ".pdf"), allPlots, width = 8, height = plot_h, units = "in")