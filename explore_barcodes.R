### NOTES & COMMENTS

# The pipeline:
# - selects a number of the most divergent unique haplotypes (longer than a specified length)
# - infers maximum likelihood (ML) tree from the selected sequences (substitution model selection included)
# - delimits operational taxonomic units (OTUs)
# - estimates ML phylogenetic placement of the remaining (=query) sequences
# - classifies query sequences into the OTUs
# - displays OTU-colored ML tree
# - optionally displays OTU-colored map of sampling locations

# In its basic form it requires specifying just the name of input alignment and outgroup indicator,
# for the mapping also an info file with geographical coordinates has to be supplied.
# Every section (introduced by '###' header) can be run separately or in a series.
# Every section starts with automated specification of necessary file names,
# then it calls either R function or third party software (IQTREE, mptp, RAxML)
# and also inputs and outputs are managed along the way
# (e.g., partition file is created, tab-delimited file is exported etc.)

# The worked out example is based on data from Uhrová et al. 2022 (Mol. Phyl Evol. 167: 107337).



### TOOLS
library(ape)
library(maps)
source("func/haplotypes.R")
source("func/raxml_partitions.R")
source("func/epatools.R")
source("func/read_ptp.R")
source("func/rootanddrop.R")
source("func/branchcutting.R")
source("func/plot_clades.R")
source("func/plot_phylomap.R")
source("func/read.fasta.R")
source("func/write.delim.R")



### SEQUENCE FILE NAME & OTHER PARAMETERS
# seqfile: the name of input alignment
# infofile: tab-delimited file with geographical coordinates, necessary only for 'plot_phylomap'
# outgroup: names of outgroup sequences or their unique identifier(s), e.g., "outgroup"
# minlength: minimum length for sequences to be included in the tree, zero indicates including everything
# no: number of the most divergent haplotypes to be retained. 'Inf' means "retain all"
# platform: either 'unix' or 'windows'
# coding: whether to consider partitioning by codon position
# method: OTU-picking method, either "bcut" (in R, all platforms) or "mptp" (calling 'mptp', unix only)
# results: where to store results, "." means the working directory 

seqfile <- "data/heliophobius_cytb.fasta"
infofile <- "data/heliophobius_info.txt"
outgroup <- "outgroup"
minlength <- 700
no <- 150
platform <- .Platform$OS.type
coding <- TRUE
method <- "bcut"
results <- "results"


### HAPLOTYPES
# data
cytb <- read.fasta(seqfile)
outgroups <- rownames(cytb)[grepl(outgroup, rownames(cytb))]

# file names
hapmap <- sub("\\..+$", "_haps.txt", seqfile)
hapseq <- sub("txt$", sub("^.+\\.", "", seqfile), hapmap)

# unique haplotypes
haplist <- haplotypes(seqfile, outgroup=outgroups)
writeLines(sapply(haplist$assign, paste, collapse=" "), hapmap)

# most divergent sequences
haps <- haplist$haplotypes[seqlen(haplist$haplotypes) >= minlength,,drop=FALSE]
if (no < (nrow(haps) - length(outgroups))) {
	haps <- divergentseq(input=haps, n=no, outgroup=outgroups, model="raw")
}
write.fasta(haps, hapseq)



### MAXIMUM LIKELIHOOD TREE (IQ-TREE)
# file names & parsing partition file
extension <- regmatches(seqfile, regexpr("\\.[[:alpha:]]+$", seqfile))
filename <- sub(extension, "", seqfile, fixed=TRUE)
hapseq <- paste0(filename, "_haps", extension)
partfile <- paste0(filename, "_part", ".txt")
nbp <- ncol(read.fasta(hapseq))
part <- ifelse(isTRUE(coding), "1-2-3", "123")
write_raxml_partitions(data.frame(start=1, end=nbp, part=part, row.names="gene"), partfile)
prefix <- paste0(results, "/", sub("^.*\\/", "", filename), "_iqt")

# iq-tree
iqtree <- ifelse(platform == "unix", "./bin/iqtree2", "bin/iqtree2.exe")
system2(iqtree, args=c("-s", hapseq, "-p", partfile, "-m", "TESTMERGE", "-mset", "HKY", "-mfreq", "F", "-mrate", "E,I,G,I+G", "-B", "1000", "-T", "AUTO", "--prefix", prefix))

# post-processing of outputs
outputs <- list.files(results, pattern=sub("^.*\\/", "", prefix), full.names=TRUE)
best_scheme <- outputs[grep("\\.best_scheme$", outputs)]
best_model <- outputs[grep("\\.best_model\\.nex$", outputs)]
treefile <- outputs[grep("\\.treefile$", outputs)]
writeLines(gsub("DNAF", "DNA", readLines(best_scheme)), best_scheme)
retained <- c(best_scheme, best_model, treefile)
invisible(file.remove(setdiff(outputs, retained)))
invisible(file.rename(retained, sub("_iqt", "", retained)))
invisible(file.remove(partfile))



### OTU PICKING (branchcutting / mPTP)
# method
if (platform != "unix" & method == "mptp") {
	method <- "bcut"
	warning("on a unix-type platform 'method' was switched to 'bcut'")	
} 

# file names
treefile <- paste0(results, "/", gsub("^.*\\/|\\..+$", "", seqfile), ".treefile")
rtreefile <- sub("\\.treefile$", "_rooted.treefile", treefile)
prefix <- paste0(results, "/", gsub("^.*\\/|\\..+$", "", seqfile), "_", method)
otufile <- sub(method, "otus.txt", prefix)

# input tree
tree <- ape::read.tree(treefile)
rtree <- rootanddrop(tree, outgroup=outgroup)
ape::write.tree(rtree, rtreefile)

# otu picking
if (method == "mptp") {
	system2("./bin/mptp", args=c("--ml", "--multi", "--tree_file", rtreefile, "--output_file", prefix))
	otus <- read_ptp(paste0(prefix, ".txt"))
}
if (method == "bcut") {
	bcut <- branchcutting(rtree)
	otus <- write.bcut(bcut)
}
names(otus) <- c("ID", method)
write.delim(otus, otufile)	



### PHYLOGENETIC PLACEMENTS (EPA)
# file names and copying
prefix <- gsub("^.*\\/|\\..+$", "", seqfile)
treefile <- paste0(results, "/", prefix, ".treefile")
best_scheme <- paste0(results, "/", prefix, ".best_scheme")
inputs <- c(seq=seqfile, tree=treefile, part=best_scheme)
copies <- sub("^.*\\/", "", inputs)
invisible(file.copy(inputs, copies))
jplacefile <- paste0(results, "/", prefix, ".jplace")

# EPA (in RAxML)
raxml <- ifelse(platform == "unix", "./bin/raxmlHPC-SSE3", "bin/raxmlHPC.exe")
system2(raxml, args=c("-f", "v", "-m", "GTRGAMMA", "--HKY85", "--epa-keep-placements=100", "--epa-prob-threshold=0.01", "-s", copies["seq"], "-t", copies["tree"], "-q", copies["part"], "-n", prefix))

# post-processing of outputs
outputs <- list.files(pattern=paste0("RAxML.*\\.", prefix))
invisible(file.rename(outputs[grep("\\.jplace$", outputs)], jplacefile))
invisible(suppressWarnings(file.remove(outputs)))
invisible(suppressWarnings(file.remove(c(copies, paste0(copies, ".reduced")))))



### CLASSIFICATION
# file names
prefix <- gsub("^.*\\/|\\..+$", "", seqfile)
jplacefile <- paste0(results, "/", prefix, ".jplace")
otufile <- paste0(results, "/", prefix, "_otus.txt")

# read in .jplace file, root the tree & classify its branches & query sequences to OTUs
jplace <- read_jplace(jplacefile)
jplace <- root_jplace(jplace, outgroup=outgroup)
otus <- read.delim(otufile)
jplace_otus <- classify_jplace(jplace, otus, ancestral=TRUE)
class_otus <- classify_sequences(jplace_otus)
class_otus_ml <- ml_classification(class_otus)

# post-processing outputs
otus <- data.frame(otus, position="crown", probability=1, totalprob=1, query=FALSE, stringsAsFactors=FALSE)
query <- data.frame(class_otus_ml[,c("id","clade","position","probability","totalprob")], query=TRUE, stringsAsFactors=FALSE)
names(query)[1:5] <- c("ID", names(otus)[2], "position", "probability","totalprob")
otus <- rbind(otus, query)
write.delim(otus, otufile)



### PLOT TREE
# file names
prefix <- gsub("^.*\\/|\\..+$", "", seqfile)
rtreefile <- paste0(results, "/", prefix, "_rooted.treefile")
otufile <- paste0(results, "/", prefix, "_otus.txt")

# data & plot
mltree <- ape::read.tree(rtreefile)
otus <- read.delim(otufile)
plot_clades(mltree, otus)



### PLOT MAP
# file names
otufile <- paste0(results, "/", gsub("^.*\\/|\\..+$", "", seqfile), "_otus.txt")

# data & map
otus <- read.delim(otufile)
info <- read.delim(infofile)
plot_phylomap(info, otus)


