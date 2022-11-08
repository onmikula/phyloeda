### TOOLS
library(ape)
library(dplyr)
source("func/readIpyrad.R")
source("func/readStacks.R")
source("func/snps_utils.R")
source("func/concatenate.R")
source("func/oneperind.R")
source("func/elbow.R")
source("func/writeSTRUCTURE.R")
source("func/branchcutting.R")
read.fasta <- function(file) toupper(ape::read.dna(file, format="fasta", as.matrix=TRUE, as.character=TRUE))
write.fasta <- function(seq, file) ape::write.dna(toupper(as.matrix(seq)), file, format="fasta")
write.delim <- function(x, file, row.names=FALSE) write.table(x, file, row.names=row.names, quote=FALSE, sep="\t")



### SEQUENCES
# The below pipeline  assumes an input in the form of list of alignments corresponding to discrete loci with phased alleles (i.e., in diploid organism, every individual is represented at every locus by two sequences). This can be obtained from the appropriate output of ddRAD-loci assembler like ipyrad or Stacks. The functions 'readIpyrad' and 'readStacks' read in such data and produce the list of locus-specific alignments. The loci can be the filtered to retain, e.g., just those with at least 80% occupancy (=sequencing success across individuals), loci of some minimum length or just variable loci (containing at least one SNP). Also, one can retain just a randomly selected subset of loci.

#loci <- readIpyrad("apodemus_agrarius_80occup_4del.alleles")
loci <- readStacks("data/georychus1000.alleles.fas")

# at least 80% occupancy
loci <- loci[occupancy(loci) >= 0.80]
# at least 100bp-long loci
loci <- loci[locilengths(loci) >= 100]
# only variable loci (especially this step must be properly justified!)
# loci <- loci[count_snps(find_snps(loci), type="var") > 0]

# subset of randomly selected loci
nloci <- 100
set.seed(123)
subset <- sort(sample(length(loci), nloci))
loci <- loci[subset]

# save the final selection of loci
save(loci, file="data/loci_subset.R")
# or, alternatively:
# saveRDS(loci, file="data/loci_subset.rds")


### SNPS
# loading of pre-selected loci, calling of biallelic SNPs (one per locus)

load("data/loci_subset.R")

usnps <- get_snps(loci, nalleles=2, type="pis", single=TRUE, as.matrix=FALSE)


### INFO
# loading of information about classification of individuals (not necessary for purely exploratory analyses)
# setting up a color palette for subsequent plotting

info <- read.delim("data/georychus_info.txt", stringsAsFactors=FALSE)
palette <- c(Cape="#e6194b", Struisbaai="#f58231", Oudshoorn="#3cb44b", KwaZuluNatal="#4363d8", Mpumalanga="#46f0f0")


### PCA
# PCA on counts of biallelic SNPs at different loci as input variables
# it is performed without scaling as all variables are counts (and hence mutually comparable)
# and without centering as they are pre-centered so the heterozygote (natural centers) ~ 0
# elbow method of Salvador & Chan (2004) is used to determine number of PCs to be retained (k)
# PS scores are displayed for PC1-2, for all pairs of PC1-k and then used for UPGMA clustering

binary <- make_binary_snps(concatenate(usnps), center=FALSE, format="snps", onerowperind=TRUE, allele="_A[[:digit:]]$")
dst <- dissim.snps(binary)
pcoa <- cmdscale(dst, eig=TRUE)
k <- elbow(pcoa$eig[pcoa$eig > 0])$best
pcs <- cmdscale(dst, k=k)

bg <- palette[info$Lineage[match(rownames(pcs), info$ddRAD)]]

par(mai=c(1.02,1.02,0.42,0.42))
plot(pcs, pch=21, cex=2, bg=bg, cex.lab=1.5, cex.axis=1.25)
legend("topright", legend=names(palette), pch=21, pt.bg=palette, pt.cex=2, bty="n")

pairs(pcs[,1:k], pch=21, cex=2, bg=bg, lwd=0.5)



### STRUCTURE
# parses STRUCTURE control files
# estimates admixture coefficients for the specified K (here 5)
# loads the results and displays them in the form of barplot with color opaqueness indicating confidence about the estimates (argument 'INTERVAL=TRUE')

strdat <- write_str_data(loci=usnps, file="data/georychus_unsps.str", allele="_A[[:digit:]]$", return=TRUE)

K <- 5
strpar <- 
data.frame(
Arg=c("infile", "outfile", "MAXPOPS", "NUMINDS", "NUMLOCI", "ONEROWPERIND", "LABEL", "POPDATA", "POPFLAG", "LOCDATA", "MAPDISTANCES", "NOADMIX", "LINKAGE", "USEPOPINFO", "LOCPRIOR", "LOCISPOP", "ANCESTDIST", "ANCESTPINT","BURNIN","NUMREPS"),
Value=c("data/georychus_unsps.str", "results/georychus_unsps.res", K, nrow(strdat), ncol(strdat)/2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.90, "30000", "100000"))

parse_str_params(df=strpar, templates=c("data/mainparams", "data/extraparams"), outputs=c("data/georychus_mainpar.txt", "data/georychus_extrapar.txt"))

structure <- ifelse(.Platform$OS.type == "unix", "./bin/structure", "bin/structure.exe")
system2(structure, args=c("-m", "data/georychus_mainpar.txt", "-e", "data/georychus_extrapar.txt"))

str_res <- read_str_output("results/georychus_unsps.res_f")
cluster <- apply(str_res[,seq(attr(str_res,"K"))], 1, which.max)
mscspec <- info$Lineage[match(names(cluster), info$ddRAD)]
matching <- sort(apply(table(mscspec, cluster), 1, which.max))
plot_str_output(str_res, interval=TRUE, palette=palette[names(matching)])





### CONCATENATED ML TREE
# assumes 'iqtree' binary in 'bin' folder
# ML tree is estimated and partitioned into OTUs using branchcutting method
# ML tree is displayed with tips colored according to OTUs


concat <- concatenate(oneperind(usnps, allele="_A[[:digit:]]$", remove=TRUE))
write.fasta(concat, "data/georychus_unsps_concat.fasta")

iqtree <- ifelse(.Platform$OS.type == "unix", "./bin/iqtree2", "bin/iqtree.exe")
system2(iqtree, args=c("-s", "data/georychus_unsps_concat.fasta", "-m", "MFP", "-B", "1000", "-T", "AUTO", "--prefix", "results/georychus_unsps_concat"))

mltree <- ape::read.tree("results/georychus_unsps_concat.treefile")
bcut <- branchcutting(mltree)
otus <- write.bcut(bcut, "results/georychus_bcut_otus.txt", return=TRUE)

bg <- palette[info$Lineage[match(mltree$tip.label, info$ddRAD)]]
#bg <- otus$OTU[match(mltree$tip.label, otus$ID)]
plot(mltree, type="unrooted", lab4ut="axial", no.margin=TRUE, label.offset=0.01, cex=0.7)
ape::tiplabels(pch=21, col=1, bg=bg, cex=1.5)



