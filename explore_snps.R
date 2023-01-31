### TOOLS
library(ape)
library(dplyr)
library(vioplot)
source("func/readIpyrad.R")
source("func/readStacks.R")
source("func/ddrad_ascertbias.R")
source("func/snps_utils.R")
source("func/binary_snps.R")
source("func/concatenate.R")
source("func/readNexus.R")
source("func/writeNexus.R")
source("func/read.fasta.R")
source("func/oneperind.R")
source("func/pcoord.R")
source("func/elbow.R")
source("func/writeSTRUCTURE.R")
source("func/write.delim.R")
source("func/max_crosssection.R")
source("func/branchcutting.R")



### SEQUENCES
# The below pipeline  assumes an input in the form of list of alignments corresponding to discrete loci with phased alleles (i.e., in diploid organism, every individual is represented at every locus by two sequences). This can be obtained from the appropriate output of ddRAD-loci assembler like ipyrad or Stacks. The functions 'readIpyrad' and 'readStacks' read in such data and produce the list of locus-specific alignments. The loci can be the filtered to retain, e.g., just those with at least 80% occupancy (= 80% sequencing success across individuals), loci of some minimum length or just variable loci (containing at least one SNP). Also, one can retain just a randomly selected subset of loci.

# LOADING & FILTERING DATA
loci <- readStacks("data/georychus1000.alleles.fas")

# The choice of the occupancy threshold is assisted by the assessment of possible ascertainment bias in increasingly filtered sets of loci. The diagnostic plots show distributions of mean pairwise genetic distance and no. of alleles.
ascert <- ddrad_ascertbias(loci, thinning=1)
if (.Platform$OS.type == "unix") device <- match.fun("quartz") else device <- match.fun("x11")
device(width=10, height=5); par(mfrow=c(1,2))
plotAscertBias(x="occupancy", y="mdist", ascert=ascert, xlab="min. occupancy", ylab="mean pairwise distance", cex.lab=1.1, cex.axis=0.9)
plotAscertBias(x="occupancy", y="nalleles", ascert=ascert, xlab="min. occupancy", ylab="no. of alleles per locus", cex.lab=1.1, cex.axis=0.9)

# at least 80% occupancy
loci <- loci[occupancy(loci) >= 0.85]
# at least 100bp-long loci
loci <- loci[locilengths(loci) >= 100]
# subset of randomly selected loci
loci <- subset_loci(loci, n=100, seed=123)

# save the final selection of loci
save(loci, file="data/georychus_subset_loci.R")
# or, alternatively:
# saveRDS(loci, file="data/loci_subset.rds")


# EXTRACTING SNPS
# loading of pre-selected loci
load("data/georychus_subset_loci.R")

# extraction of biallelic (nalleles=2), parsimony-informaive (pis=TRUE) SNPs
snps <- get_snps(loci, pis=TRUE, nalleles=2)
# retaining just a single SNP per locus, while minimizing amount of missing data, a seed for stochastic component of the selection is specified to obtain always the same set of SNPs
usnps <- single_snps(snps, minmiss=TRUE, seed=123)



### INFO
# loading of information about classification of individuals
# to assisst interpretation of the observed patterns
# not necessary for the analyses themselves!
# setting up a color palette for plotting

info <- read.delim("data/georychus_info.txt", stringsAsFactors=FALSE)
palette <- c(Cape="#e6194b", Struisbaai="#f58231", Oudshoorn="#3cb44b", KwaZuluNatal="#4363d8", Mpumalanga="#46f0f0")



### PCA & PCoA
# Principal Component Analysis (PCA) on counts of binary-coded alleles of biallelic SNPs
# or Principal Coordinates Analysis (PCoA) on a matrix of pairwise genotype dissimilarities
# The dissimilarity is an arithmetic mean of differences in counts of alleles coded by '1s' after pairwise deletion of missing data. The missing data can be then imputed as an arithmetic mean of counts of '1' alleles observed in k nearest-neighbors at given SNP, with the pairwise dissimilarities being used to find the neighbors.
# PCA is performed without scaling as all variables are counts and hence mutually comparable
# and without centering as they are pre-centered so the heterozygote (natural centers) ~ 0.
# The centering by sample mean is implicit in PCoA and its results can differ from those of PCA, therefore.
# The elbow method of Salvador & Chan (2004) is used to determine the number of PCs/PCos to be retained (k)
# PC scores are displayed for axes 1-2 and for all pairs of axes 1-k 


# PREPARING BINARY SNP MATRIX
# converting biallelic SNPs into binary (0, 1) data 
# if oneperind==TRUE, then rows contain counts of '1' alleles (0 & 2 ~ homozygotes, 1 ~ heterozygote)
# if center==TRUE, then (0, 1, 2) is converted into (-1, 0, 1)
binary <- make_binary_snps(concatenate(usnps), center=TRUE, format="snps", onerowperind=TRUE, allele="_A[[:digit:]]$")
# calculating dissimilarity matrix
dissmat <- dissim.snps(binary)
# missing data imputation 
binary <- impute_binary_snps(binary, d=dissmat, k=3)
# deleting columns with any missing data
binary <- clean_data(binary, remove="columns")
# export 
write.delim(binary, "data/georychus_subset_binary.txt", row.names=TRUE)


# PRINCIPAL COORDINATES ANALYSIS (PCoA)
pc <- pcoord(dissmat)
k <- elbow(pc$values)$best
pcs <- pc$vectors
bg <- palette[info$OTU[match(rownames(pcs), info$ID)]]

par(mai=c(1.02,1.02,0.42,0.42))
plot(pcs, pch=21, cex=2, bg=bg, cex.lab=1.5, cex.axis=1.25)
legend("topright", legend=names(palette), pch=21, pt.bg=palette, pt.cex=2, bty="n")

pairs(pcs[,1:k], pch=21, cex=2, bg=bg, lwd=0.5)


# PRINCIPAL COMPONENT ANALYSIS (PCA)
# note the PCs may be correlated when no further centering is applied, explore the effect of setting center=TRUE 
pc <- prcomp(binary, center=FALSE)
k <- elbow(pc$sdev^2)$best
pcs <- pc$x
bg <- palette[info$OTU[match(rownames(pcs), info$ID)]]

par(mai=c(1.02,1.02,0.42,0.42))
plot(pcs, pch=21, cex=2, bg=bg, cex.lab=1.5, cex.axis=1.25)
legend("topleft", legend=names(palette), pch=21, pt.bg=palette, pt.cex=2, bty="n")

pairs(pcs[,1:k], pch=21, cex=2, bg=bg, lwd=0.5)



### STRUCTURE
# assumes 'STRUCTURE' binary in 'bin' folder
# parses STRUCTURE control files
# estimates admixture coefficients for the specified K, which can be:
# - the number of pre-defined groups (OTUs, lineages, species, populations...)
# - the number of clusters in the UPGMA dendrogram of PC scores chosen by the maximum cross-section branch length criterion
# - a series of numbers (K is indicated in output names)
# loads the results and displays them in the form of barplot with color opaqueness indicating confidence about the estimates (argument 'INTERVAL=TRUE')

# setting K
K <- nlevels(factor(info$OTU))
#upgma <- ape::as.phylo(stats::hclust(dist(pcs[,1:k]), method="average"))
#K <- max(max_crosssection(upgma)[,2])

# file names
str_prefix <- "georychus_unsps"
str_input <- paste0("data/", str_prefix, ".str")
str_params <- paste0("results/", str_prefix, paste0("_K",K), c("_mainpar", "_extrapar"), ".txt")
str_output <- paste0("results/", str_prefix, paste0("_K",K), ".res")
str_output_f <- paste0(str_output, "_f")

# parsing files
strdat <- write_str_data(loci=usnps, file=str_input, allele="_A[[:digit:]]$", return=TRUE)
strpar <- data.frame(
Arg=c("infile", "outfile", "MAXPOPS", "NUMINDS", "NUMLOCI", "ONEROWPERIND", "LABEL", "POPDATA", "POPFLAG", "LOCDATA", "MAPDISTANCES", "NOADMIX", "LINKAGE", "USEPOPINFO", "LOCPRIOR", "LOCISPOP", "ANCESTDIST", "ANCESTPINT","BURNIN","NUMREPS"),
Value=c(str_input, str_output, K, nrow(strdat), ncol(strdat)/2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.90, "3000", "10000"))
parse_str_params(df=strpar, templates=c("data/mainparams", "data/extraparams"), outputs=str_params)

# analysis & results
structure <- ifelse(.Platform$OS.type == "unix", "./bin/structure", "bin/structure.exe")
system2(structure, args=c("-m", str_params[1], "-e", str_params[2]))

str_res <- read_str_output(str_output_f)
str_palette <- custom_str_palette(str_res, info, palette, ind="ID", pop="OTU")
plot_str_output(str_res, interval=TRUE, palette=str_palette, show.ind.label=TRUE)



### CONCATENATED ML TREE
# assumes 'iqtree' binary in 'bin' folder
# ML tree is estimated and partitioned into OTUs using branchcutting method
# ML tree is displayed with tips colored according to pre-defined OTUs

# file names
iqt_prefix <- "georychus_unsps_concat"
iqt_input <- paste0("data/", iqt_prefix, ".fasta")
iqt_output <- paste0("results/", iqt_prefix)
iqtreefile <- paste0(iqt_output, ".iqtree")
mltreefile <- paste0(iqt_output, ".treefile")
otus_bcut <- paste0(iqt_output, "_otus_bcut.txt")

# input data
concat <- concatenate(oneperind(usnps, allele="_A[[:digit:]]$", remove=TRUE))
write.fasta(concat, iqt_input)

# analysis & results
iqtree <- ifelse(.Platform$OS.type == "unix", "./bin/iqtree2", "bin/iqtree2.exe")
system2(iqtree, args=c("-s", iqt_input, "-m", "MFP", "-B", "1000", "-T", "AUTO", "--prefix", iqt_output))

path <- sub("\\/$", "", regmatches(iqt_output, regexpr("^.+\\/", iqt_output)))
outputs <- list.files(path, pattern=iqt_prefix, full.names=TRUE)
invisible(file.remove(setdiff(outputs, c(iqtreefile, mltreefile))))

mltree <- ape::read.tree(mltreefile)
bg <- palette[info$OTU[match(mltree$tip.label, info$ID)]]
plot(mltree, type="unrooted", lab4ut="axial", no.margin=TRUE, label.offset=0.01, cex=0.7)
ape::tiplabels(pch=21, col=1, bg=bg, cex=1.5)

# delimitation of OTUs by the branchcuttting algorithm
bcut <- branchcutting(mltree)
otus <- write.bcut(bcut, otus_bcut, return=TRUE)


