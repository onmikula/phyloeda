### TOOLS
library(phyloeda)
library(ape)
library(branchcutting)
library(epaclades)
library(maps)


### PIPELINE
# input specification
seqfile <- "Aethomys_2024-08-15.fasta"
infofile <- "Aethomys_2024-08-15.txt"
prefix <- "Aethomys_cytb"
outdir <- "results"
outgroup <- "outgroup"
minlength <- 700
no <- 75

# file names
refseq <- paste0(outdir, "/", prefix, ".fasta")
treefile <- paste0(outdir, "/", prefix, ".treefile")
bcutfile <- paste0(outdir, "/", prefix, "_bcut.txt")
epafile <- paste0(outdir, "/", prefix, "_epa.fasta")
jplacefile <- paste0(outdir, "/", prefix, ".jplace")
finalclass <- paste0(outdir, "/", prefix, "_otus.txt")

# haplotypes
cytb <- read.fasta(seqfile)
haps <- haplotypes(cytb, outgroup=outgroup)
divseq <- divergentseq(input=haps, n=no, minlen=minlength, outgroup=outgroup, fasta=refseq)

# maximum likelihood tree (IQ-TREE)
args <- data.frame(args=c("-s", "-m", "-mset", "-mfreq", "-mrate", "-B", "-T"),
values=c(refseq, "TESTMERGE", "mrbayes", "FO,FQ", "E,I,G,I+G", "1000", "AUTO"))
run_iqtree(command="./bin/iqtree2", args=args, prefix=prefix, outdir=outdir, coding=TRUE, retain=c("\\.best_model", "\\.treefile$"))

# species delimitation (branchcutting)
bc <- branchcutting(treefile, outgroup=outgroup)
taxa <- export(bc, bcutfile)

# phylogenetic placement (EPA in RAxML) + haplotype matches 
query <- findqueryseq(ref=treefile, haps)
if (length(query) > 0) {
	epainput(ref=treefile, query=query, seq=cytb, file=epafile)
	run_epa(command="./bin/raxmlHPC-SSE3", args=c("-s", epafile, "-t", treefile), prefix=prefix, outdir=outdir)
	jp <- read_jplace(jplacefile)
	jp <- root_jplace(jp, outgroup=outgroup)
	jp <- classify_jplace(jp, clades=taxa)
	otus <- epaformat(taxa=taxa, haps=haps, epa=jp, file=finalclass, return=TRUE, na.rm=FALSE, method="bcut")
} else {
	otus <- epaformat(taxa=taxa, haps=haps, file=finalclass, return=TRUE, na.rm=FALSE, method="bcut")
}

# tree
tree <- rootanddrop(read_tree(treefile), outgroup)
plot_clades(tree, taxa)

# map
info <- read.delim(infofile)
info$bcut <- otus$clade[match(info$ID, otus$id)]
plot_phylomap(info, lin="bcut")

