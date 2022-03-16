# Test whether abundance-weighted FracMinHash sketches can be used to regenerate findings from greenelab/core-accessory-interactome

The [greenelab/core-accessory-interactome](https://github.com/greenelab/core-accessory-interactome/) repository (core-acc) uses compendia of *P. aeruginosa* gene expression to identify relationships between the expression of core and accessory genes.
The compendia were created using publicly available RNA-seq data (see [preprint](https://doi.org/10.1101/2022.01.24.477642), [data repository](https://osf.io/s9gyu/), and [GitHub repository](https://github.com/georgiadoing/pa-seq-compendia)) and represent gene counts for two *P. aeruginosa* strains, PAO1 and PA14.
The gene counts were derived from quasi-mapping SRA fastq files against both reference transcriptomes.
Using these compendia, the core-acc project identified 545 samples that best mapped to the PA14 strain, and 861 samples that best mapped to the PAO1 strain. 
The core-acc project performed first- and second-order correlation analyses of gene expression between the two strain types to identify: 
1) core genes that had consistent co-expression with consistent sets of other genes across RNA seq libraries and across strains 
2) modules of accessory genes that were co-expressed.

The goal of this repository is to determine whether the same relationships can be recovered using abundance-weighted FracMinHash sketches instead of gene counts.
If it can, then that may open the door to FracMinHash sketches being used to identify similar relationships in other data sets.

FracMinHash sketches capture a fraction of the k-mers in a given data set (see [here](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract)). 

The benefits of FracMinHash sketches include:
+ dramatic data reduction 
+ streaming-compatible, so data (FASTQ/FASTA) can be piped directly to sketching without writing to hard disk first
+ faster processing from raw data to abundances than quasi-mapping- or mapping-based quantification pipelines
+ reference- and annotation-free, so k-mers that are in a sample but are not yet in a database can be included in analysis

However, FracMinHash sketches may not capture the same information as gene count matrices because:
+ nucleotide k-mers may be too specific to capture information that is represented by gene counts; nucleotide k-mers are more brittle to evolutionary distance than (quasi-)mapping.
In general, a k-mer size of 21 is fairly specific to a genus of organisms, a k-mer of size 31 is specific to a species of organisms, and a k-mer of size 51 is specific to a strain/genome of an organism (see [here](https://bluegenes.github.io/2022-paper-protein-kmers/#amino-acid-k-mer-length-selection) and [here](https://journals.asm.org/doi/full/10.1128/mSystems.00020-16)).
This means that even for fairly closely related organisms, there may be relatively few nucleotide k-mers in common between genomes.
Anecdotally, we expect a Jaccard similarity of approximately 0.4 for genomes of the same species, and we hope that that is enough k-mer overlap to represent the same information contained in gene counts.
However, this could create problems if not enough k-mers are shared between related genomes.
+ FracMinHashes are randomly (but consistently) *fractional*. 
We use a scaled value of 1000 to calculate sketches, meaning approximately 1/1000th of all k-mers are included in the final sketch. 
These are approximately distributed across the genome, but the hash function is not designed to guaruntee sampling from all windows of the genome, so some sections may be missed.
The probability that a k-mer will be observed in each 1000 base pair window of a genome sequence approximately follows a Poisson distribution (I think).
The hash function also does not guaruntee biological relevance of a sequence, so a k-mer may represent a gene, a promoter, intergenic sequence, etc.
+ FracMinHash abundances are abundances for a single k-mer instead of the number of reads that map to a gene. 
Given the small number of base pairs that these abundances are constructed from, it is possible that the abundances may not reflect the abundances of surrounding k-mers and therefore may not be useful for correlations.
+ There are many more hashes than genes, so using a scaled value of 1000 may create a data set that is too high dimensional to analyze as was done in core-acc.
