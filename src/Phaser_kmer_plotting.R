### CBS-Frontier draft

# this script is still in development and not quite "pipelinized" yet
# the items in the hash box below should be the only things you need to change if
# using default output from Discovery and Phaser
# However you can also customize the figures in the PlotList loop below

rm(list = ls())

library(optparse)
library(ggplot2)
library(ggpubr)

################################################################################
# Command line options
################################################################################

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
# Edit the inputs for local use
################################################################################

# working directory
setwd(paste(opts_parsed$phaser_path))

# prefix name for file output
outpref <- opts_parsed$outpref

# CBS-Phaser output folder name with bam/bam.bai files
ref_fai <- paste(opts_parsed$ref_file)

bin_size <- as.integer(opts_parsed$bin_size)

# Minimum chr size for plotting, if you want to show mappings for all contigs or
# fragments (like if comparing to unscaffolded assembly) then set to 1
min_chr_size <- as.integer(opts_parsed$min_chr_size)

################################################################################
# Run scripts below for output, shouldn't need to change anything below 
################################################################################

# chromSizes for y-lim
ref_chromSizes <- read.delim(file = ref_fai, header = F, sep = "\t")
ref_chromSizes <- ref_chromSizes[ref_chromSizes$V2 > min_chr_size,]
ref_chromSizes <- ref_chromSizes[ref_chromSizes$V2 > bin_size,]

# covs <- foreach(i=list.files(path = ".", pattern = "*.depth.txt"), .combine = rbind) %do% {
#   f <- read.csv(i, header = F, sep = "\t")
#   #colnames(f) <- c("chr", "start", "end", "depth")
#   colnames(f) <- c("chr", "start", "end", "depth")
#   f$sample <- gsub(".[0-9][0-9]mers.list_scaffolds.depth.txt","", i)
#   f$sex <- gsub(".*(Male|Female).*Only.*", "\\1", f$sample)
#   return(f)
# }
depth_files <- list.files(path = ".", pattern = paste0(bin_size, ".s[0-9]*.depth.txt"))
covs_list <- lapply(depth_files, function(i) {
  f <- read.csv(i, header = FALSE, sep = "\t")
  colnames(f) <- c("chr", "start", "end", "depth")
  f$sample <- gsub(".[0-9][0-9]mers.list_scaffolds.[0-9]*.s[0-9]*.depth.txt", "", i)
  f$sex <- gsub(".*(Male|Female).*Only.*", "\\1", f$sample)
  f
})
covs <- do.call(rbind, covs_list)
covs <- covs[covs$chr %in% ref_chromSizes$V1,]

# 100 windows to make comparable with samtools histogram
n100_depth_files <- list.files(path = ".", pattern = ".n100.windows.depth.txt")
n100_covs_list <- lapply(n100_depth_files, function(i) {
  f <- read.csv(i, header = FALSE, sep = "\t")
  colnames(f) <- c("chr", "start", "end", "depth")
  f$sample <- gsub(".[0-9][0-9]mers.list_scaffolds.n100.windows.depth.txt", "", i)
  f$sex <- gsub(".*(Male|Female).*Only.*", "\\1", f$sample)
  f
})
n100_covs <- do.call(rbind, n100_covs_list)
n100_covs <- n100_covs[n100_covs$chr %in% ref_chromSizes$V1,]

# setting color palettes and linetypes

color_frame <- covs[!duplicated(covs[,c("sample", "sex")]), c("sample", "sex")]
color_frame$color_vals <- "black"

# Colors come from RColorBrewer palettes Green and Purple
# rev(brewer.pal(n=9, name = "Greens"))
# rev(brewer.pal(n=9, name = "Purples"))

female_colors <- adjustcolor(c("#00441B", "#006D2C", "#238B45", "#41AB5D", "#74C476", "#A1D99B"), alpha.f = 0.75)
male_colors <- adjustcolor(c("#3F007D", "#54278F", "#6A51A3", "#807DBA", "#9E9AC8", "#BCBDDC"), alpha.f = 0.75)

color_frame$color_vals[color_frame$sex=="Female"] <- 
  female_colors[c(1:length(which(color_frame$sex=="Female")))]
color_frame$color_vals[color_frame$sex=="Male"] <- 
  male_colors[c(1:length(which(color_frame$sex=="Male")))]

covs <- merge(covs, color_frame, by = c("sample", "sex"))
n100_covs <- merge(n100_covs, color_frame, by = c("sample", "sex"))

chr.split <- split(covs, covs$chr)
n100_chr.split <- split(n100_covs, n100_covs$chr)

# user chosen bins (or default 1+06)
cov_plots <- lapply(names(chr.split), function(chr) {
  c <- chr.split[[chr]]
  ggplot(c, aes(x = start, y = depth, group = sample, color = color_vals, linetype = sex)) +
    geom_line(show.legend = T) +
    theme_bw(base_size = 8) +
    labs(title = chr, subtitle = paste0("bin size = ", bin_size, " bp"), x = "position (bp)", y = "coverage depth") +
    scale_color_identity(guide = "legend", breaks = color_frame$color_vals, labels = color_frame$sample, 
                         limits = color_frame$color_vals, name = "k-mer intersect", drop = FALSE) +
    scale_linetype_manual(guide = "legend", values = c("Female" = "solid", "Male" = "dashed"), 
                          limits = names(c("Female" = "solid", "Male" = "dashed")), drop = FALSE) +
    scale_x_continuous(breaks = c(1, ref_chromSizes$V2[ref_chromSizes$V1== chr]),
                       labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    xlim(c(1, ref_chromSizes$V2[ref_chromSizes$V1 == chr])) +
    theme(plot.title = element_text(hjust = 0.5))
})

# scaled to y plot
cov_plots_scaled <- lapply(names(chr.split), function(chr) {
  c <- chr.split[[chr]]
  ggplot(c, aes(x = start, y = depth, group = sample, color = color_vals, linetype = sex)) +
    geom_line(show.legend = T) + 
    theme_bw(base_size = 8) +
    labs(title = chr, subtitle = paste0("bin size = ", bin_size, " bp"), x = "position (bp)", y = "coverage depth") +
    scale_color_identity(guide = "legend", breaks = color_frame$color_vals, labels = color_frame$sample, 
                         limits = color_frame$color_vals, name = "k-mer intersect", drop = FALSE) +
    scale_linetype_manual(guide = "legend", values = c("Female" = "solid", "Male" = "dashed"), 
                          limits = names(c("Female" = "solid", "Male" = "dashed")), drop = FALSE) +
    scale_x_continuous(breaks = c(1, ref_chromSizes$V2[ref_chromSizes$V1==c$chr[1]]),
                       labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    xlim(c(1, ref_chromSizes$V2[ref_chromSizes$V1== chr])) +
    ylim(c(0, max(covs$depth))) +
    theme(plot.title = element_text(hjust = 0.5))
  
}) 

allPlots <- ggarrange(plotlist = cov_plots,
                      ncol = 2, nrow = 3,
                      common.legend = TRUE,
                      legend = "bottom")

allPlots_scaled <- ggarrange(plotlist = cov_plots_scaled, 
                             ncol = 2, nrow = 3, 
                             common.legend = TRUE, 
                             legend = "bottom")

# save outputs
ggexport(allPlots, filename = paste0(outpref, ".", bin_size, "bp.pdf"))
ggexport(allPlots_scaled, filename = paste0(outpref, ".", bin_size, "bp_scaled.pdf"))

# n100 windows
n100_cov_plots <- lapply(names(n100_chr.split), function(chr) {
  c <- n100_chr.split[[chr]]
  ggplot(c, aes(x = start, y = depth, group = sample, color = color_vals, linetype = sex)) +
    geom_line(show.legend = T) +
    theme_bw(base_size = 8) +
    labs(title = chr, subtitle = paste0("100 bins"), x = "position (bp)", y = "coverage depth") +
    scale_color_identity(guide = "legend", breaks = color_frame$color_vals, labels = color_frame$sample, 
                         limits = color_frame$color_vals, name = "k-mer intersect", drop = FALSE) +
    scale_linetype_manual(guide = "legend", values = c("Female" = "solid", "Male" = "dashed"), 
                          limits = names(c("Female" = "solid", "Male" = "dashed")), drop = FALSE) +
    scale_x_continuous(breaks = c(1, ref_chromSizes$V2[ref_chromSizes$V1== chr]),
                       labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    xlim(c(1, ref_chromSizes$V2[ref_chromSizes$V1 == chr])) +
    theme(plot.title = element_text(hjust = 0.5))
})

# scaled to y plot
n100_cov_plots_scaled <- lapply(names(n100_chr.split), function(chr) {
  c <- n100_chr.split[[chr]]
  ggplot(c, aes(x = start, y = depth, group = sample, color = color_vals, linetype = sex)) +
    geom_line(show.legend = T) + 
    theme_bw(base_size = 8) +
    labs(title = chr, subtitle = paste0("100 bins"), x = "position (bp)", y = "coverage depth") +
    scale_color_identity(guide = "legend", breaks = color_frame$color_vals, labels = color_frame$sample, 
                         limits = color_frame$color_vals, name = "k-mer intersect", drop = FALSE) +
    scale_linetype_manual(guide = "legend", values = c("Female" = "solid", "Male" = "dashed"), 
                          limits = names(c("Female" = "solid", "Male" = "dashed")), drop = FALSE) +
    scale_x_continuous(breaks = c(1, ref_chromSizes$V2[ref_chromSizes$V1==c$chr[1]]),
                       labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    xlim(c(1, ref_chromSizes$V2[ref_chromSizes$V1== chr])) +
    ylim(c(0, max(n100_covs$depth))) +
    theme(plot.title = element_text(hjust = 0.5))
  
}) 

n100_allPlots <- ggarrange(plotlist = n100_cov_plots,
                      ncol = 2, nrow = 3,
                      common.legend = TRUE,
                      legend = "bottom")

n100_allPlots_scaled <- ggarrange(plotlist = n100_cov_plots_scaled, 
                             ncol = 2, nrow = 3, 
                             common.legend = TRUE, 
                             legend = "bottom")

# save outputs
ggexport(n100_allPlots, filename = paste0(outpref, "100_bins.pdf"))
ggexport(n100_allPlots_scaled, filename = paste0(outpref, "100_bins_scaled.pdf"))

# saving un-scaled cov plots as RData object for manual edits
save(covs, n100_covs, cov_plots, n100_cov_plots, file = paste0(outpref, "_plotList.RData"))