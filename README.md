# phyloeda
Exploratory analyses of phylogenetic / phylogeographic data.

Currently, it includes two pipelines:
1.  explore_barcodes starts with alignment of barcode genes (typically CYTB, COI, dloop, 16S, ...) and goes through haplotype selection, ML tree inference, OTU picking and phylogenetic placement of sequences to the visualization of OTUs (colored tree) and their distributional ranges (simple map with colored points)
2. explore_snps with exploratory analyses of biallelic SNPs: PCA, STRUCTURE and ML inference of concatenated tree
