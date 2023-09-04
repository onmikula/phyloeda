### TOOLS
library(ape)
library(dplyr)
library(vioplot)
library(maps)

source("func/read_loci.R")
source("func/write_loci.R")
source("func/remove_repeats.R")
source("func/ddrad_ascertbias.R")
source("func/loci_utils.R")
source("func/snps_utils.R")
source("func/even_subsampling.R")
source("func/plot_gencoord.R")

source("func/collapse_diploid_haplotypes.R")
source("func/concatenate.R")
source("func/branchcutting.R")
source("func/read.fasta.R")
source("func/write.delim.R")

source("func/binary_snps.R")
source("func/pcomp.R")
source("func/elbow.R")

source("func/plot_str_output.R")
source("func/writeSTRUCTURE.R")

source("func/plot_phylomap.R")


### SEQUENCES
# This pipeline assumes an input in the form of a list of sequences alignments, each corresponding to a discrete locus. The parental copies are distinguished (heterozygous positions phased) within loci, but not across them. In a diploid organism, therefore, every individual is represented by two sequences at every locus. This is included, e.g., in '.alleles.fas' output of Stacks or '. alleles' output of ipyrad. The functions 'read_ipyrad' and 'read_stacks' read in data from the respective assemblers and produce these lists of locus-specific alignments. The loci can be the filtered to retain, e.g., just those with at least 80% occupancy (= 80% sequencing success across individuals), loci of some minimum length or just variable loci (containing at least one SNP). Also, one can retain just a randomly selected subset of loci or - in case of assembly by mapping to reference - a specified number of loci evenly spaced along the reference genome.


# LOADING & FILTERING DATA
# loading of a small sample of 1500 loci (the real numbers can be as high as tens of thousands)
loci <- read_ipyrad("data/georychus1500.alleles")

# possibly you can get rid off sequence repeats using the function 'remove_repeats' (but it is time-consuming with large data sets)
# loci <- remove_repeats(loci)

# The choice of the occupancy threshold is assisted by the assessment of possible ascertainment bias in increasingly filtered sets of loci. The diagnostic plots show distributions of mean pairwise genetic distance, no. of alleles and no. of SNPs.
# with a large data set it is wise to apply sample thinning a analyze, e.g., every 10th locus (thinning=10)

ascert <- ddrad_ascertbias(loci, thinning=1)
if (.Platform$OS.type == "unix") device <- match.fun("quartz") else device <- match.fun("x11")
device(width=12, height=4); par(mfrow=c(1,3), mai=c(1.02,1.02,0.22,0.82))
plotAscertBias(x="occupancy", y="mdist", ascert=ascert, xlab="min. occupancy", ylab="mean pairwise distance", cex.lab=0.9)
plotAscertBias(x="occupancy", y="nalleles", ascert=ascert, xlab="min. occupancy", ylab="no. of alleles per locus", cex.lab=0.9)
plotAscertBias(x="occupancy", y="nsnps", ascert=ascert, xlab="min. occupancy", ylab="no. of SNPs per locus", cex.lab=0.9)


# at least 80% occupancy
loci <- subset_loci(loci, occupancy(loci) >= 0.80)

# at least 100bp-long loci (when rm.gaps=TRUE gaps are excluded from the calculation - often meaningful, but can take long)
loci <- subset_loci(loci, locilengths(loci, rm.gaps=TRUE) >= 100)

# subset of loci for the purpose of demonstration
# 500 loci evenly spaced along the references genome (+ plot of their position along 20 longest scaffolds)
bedtable <- attr(loci, "bedtable")
chrsizes <- read.delim("data/Fukomys_damarensis_GCA_000743615.1_chromosizes.txt", header=FALSE)
even <- even_subsampling(bed=bedtable, nloci=500, chromsize=chrsizes)
plot_gencoord(bedtable, chrsizes, subset=even[,4], no=20, order=TRUE)
loci <- subset_loci(loci, subset=even[,4])
write_ipyrad(loci, "data/georychus500.alleles")

# alternatively - 500 randomly selected loci
#loci <- subset_loci(loci, n=500, seed=123)
#write_ipyrad(loci, "data/georychus500.alleles")



# EXTRACTING SNPS
# read in pre-selected loci
loci <- read_ipyrad("data/georychus500.alleles")

# extraction of biallelic (nalleles=2), parsimony-informaive (pis=TRUE) SNPs
snps <- get_snps(loci, pis=TRUE, nalleles=2)
# retaining just a single SNP per locus, while minimizing amount of missing data
# a seed for stochastic component of the selection is specified to obtain always the same set of SNPs
ssnps <- single_snps(snps, minmiss=TRUE, seed=123)



### CONCATENATED ML TREE
# assumes 'iqtree' binary in 'bin' folder
# ML tree is estimated and partitioned into OTUs using branchcutting method
# ML tree is displayed with tips colored according to pre-defined OTUs

# file names
prefix <- "georychus_ssnps"
iqt_input <- paste0("data/", prefix, ".fasta")
iqt_output <- paste0("results/", prefix)
iqtreefile <- paste0(iqt_output, ".iqtree")
mltreefile <- paste0(iqt_output, ".treefile")
otus_bcut <- paste0(iqt_output, "_bcut.txt")

# input data
collapsed <- concatenate(collapse_diploid_haplotypes(ssnps, allele="_[01]$"))
write.fasta(collapsed, iqt_input)

# analysis
iqtree <- ifelse(.Platform$OS.type == "unix", "./bin/iqtree2", "bin/iqtree2.exe")
system2(iqtree, args=c("-s", iqt_input, "-m", "MFP", "-B", "1000", "-T", "AUTO", "--prefix", iqt_output))

# results
path <- sub("\\/$", "", regmatches(iqt_output, regexpr("^.+\\/", iqt_output)))
outputs <- list.files(path, pattern=prefix, full.names=TRUE)
invisible(file.remove(setdiff(outputs, c(iqtreefile, mltreefile))))
mltree <- ape::read.tree(mltreefile)

# delimitation of OTUs by the branchcuttting algorithm
bcut <- branchcutting(mltree)
otus <- write.bcut(bcut, otus_bcut, return=TRUE)

# tree with tips colored by estimated OTUs
color <- grDevices::palette()[otus$OTU[match(mltree$tip.label, otus$ID)] + 1]
plot(mltree, type="unrooted", lab4ut="axial", no.margin=TRUE, label.offset=0.01, cex=0.7)
ape::tiplabels(pch=21, col=1, bg=color, cex=1.5)


### Principal Component Analysis (PCA)
# performed on counts of binary-coded alleles of biallelic SNPs 
# as described in Patterson et al. 2006 (doi:10.1371/journal.pgen.0020190 )

# converting biallelic SNPs into binary (0, 1) data 
# (oneperind==TRUE means rows contain counts of '1' alleles, i.e. 0 & 2 ~ homozygotes, 1 ~ heterozygote)
binary <- make_binary_snps(snps, onerowperind=TRUE, allele="_[01]$")

# the centring (center==TRUE) in 'make_binary_snps; converts (0, 1, 2) into (-1, 0, 1)
# scaling by sqrt(p * (1 - p)), p being the frequency of '1 allele'
# & centering by substraction of sample means
binary <- scale_binary_snps(binary, center=TRUE)

# PCA by spectral decomposition of a matrix of individual genotype similarities
# (cross-products of individual genotype vectors, after pairwise deletion of missing data)
# the elbow method of Salvador & Chan (2004) is used to determine the number of PCs to be retained
pc <- pcomp(binary)
k <- elbow(pc$values)$best

# PC scores are displayed for axes 1-2 and for all pairs of axes 1-k 
# colors indicate classification of individuals into branchcutting units
prefix <- "georychus_ssnps"
otus <- read.delim(paste0("results/", prefix, "_bcut.txt"))
palette <- setNames(grDevices::palette()[1:max(otus$OTU)+1], seq(max(otus$OTU)))
color <- palette[otus$OTU[match(rownames(pc$vectors), otus$ID)]]
plot(pc$vectors, pch=21, cex=2, bg=color, cex.lab=1.3, cex.axis=1.25)
legend("topleft", legend=paste("BCUT", names(palette)), pch=21, pt.bg=palette, pt.cex=2, bty="n")
pairs(pc$vectors[,1:k], pch=21, cex=2, bg=color, lwd=0.5)


### STRUCTURE
# assumes 'STRUCTURE' binary in 'bin' folder
# parses STRUCTURE control files
# estimates admixture coefficients for the specified K, which can be:
# - the number of pre-defined groups (OTUs, lineages, species, populations...)
# - the number of clusters in the UPGMA dendrogram of PC scores chosen by the maximum cross-section branch length criterion
# - a series of numbers (K is indicated in output names)
# loads the results and displays them in the form of barplot with color opaqueness indicating confidence about the estimates (argument 'INTERVAL=TRUE')

# setting K
prefix <- "georychus_ssnps"
otus <- read.delim(paste0("results/", prefix, "_bcut.txt"))
K <- nlevels(factor(otus$OTU))

# file names
prefix <- "georychus_unsps"
str_input <- paste0("data/", prefix, ".str")
str_params <- paste0("results/", prefix, paste0("_K",K), c("_mainpar", "_extrapar"), ".txt")
str_output <- paste0("results/", prefix, paste0("_K",K), ".res")
str_output_f <- paste0(str_output, "_f")

# parsing files
strdat <- write_str_data(loci=ssnps, file=str_input, allele="_[01]$", return=TRUE)
strpar <- data.frame(
Arg=c("infile", "outfile", "MAXPOPS", "NUMINDS", "NUMLOCI", "ONEROWPERIND", "LABEL", "POPDATA", "POPFLAG", "LOCDATA", "MAPDISTANCES", "NOADMIX", "LINKAGE", "USEPOPINFO", "LOCPRIOR", "LOCISPOP", "ANCESTDIST", "ANCESTPINT","BURNIN","NUMREPS"),
Value=c(str_input, str_output, K, nrow(strdat), ncol(strdat)/2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.90, "3000", "10000"))
parse_str_params(df=strpar, templates=c("data/mainparams", "data/extraparams"), outputs=str_params)

# analysis & results
structure <- ifelse(.Platform$OS.type == "unix", "./bin/structure", "bin/structure.exe")
system2(structure, args=c("-m", str_params[1], "-e", str_params[2]))

str_res <- read_str_output(str_output_f)

palette <- setNames(grDevices::palette()[1:max(otus$OTU)+1], seq(max(otus$OTU)))
str_palette <- custom_str_palette(str_res, otus, palette, ind="ID", pop="OTU")
plot_str_output(str_res, interval=TRUE, palette=str_palette, show.ind.label=TRUE)



### PHYLOGEOGRAPHIC MAP
# information about the individuals (incl. geographical coordinates)
info <- read.delim("data/georychus_info.txt", stringsAsFactors=FALSE)

# branchcutting classification of the individuals
prefix <- "georychus_ssnps"
otus <- read.delim(paste0("results/", prefix, "_bcut.txt"))
palette <- setNames(grDevices::palette()[1:max(otus$OTU)+1], seq(max(otus$OTU)))

# map
plot_phylomap(info, otus, color=palette)
legend("bottomright", legend=paste("BCUT", names(palette)), pch=21, pt.bg=palette, pt.cex=2, bty="n")
