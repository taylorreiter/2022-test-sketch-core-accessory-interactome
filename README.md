# Test whether abundance-weighted FracMinHash sketches can be used to regenerate findings from greenelab/core-accessory-interactome

The goal of this repository is to determine whether gene co-expression correlations can be recovered using abundance-weighted FracMinHash sketches in place of gene counts.
If FracMinHash sketches can be used for this purpose, this may open the door for these sketches to be used to identify similar relationships in other data sets.

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
These are approximately distributed across the genome, but the hash function is not designed to guarantee sampling from all windows of the genome, so some sections may be missed.
The probability that a k-mer will be observed in each 1000 base pair window of a genome sequence approximately follows a Poisson distribution (I think).
The hash function also does not guarantee biological relevance of a sequence, so a k-mer may represent a gene, a promoter, intergenic sequence, etc.
+ FracMinHash abundances are abundances for a single k-mer instead of the number of reads that map to a gene. 
Given the small number of base pairs that these abundances are constructed from, it is possible that the abundances may not reflect the abundances of surrounding k-mers and therefore may not be useful for correlations.
+ There are many more hashes than genes, so using a scaled value of 1000 may create a data set that is too high dimensional to analyze as was done in core-acc.

## Test data set

To test whether sub sampled k-mer abundances can be used to recover gene co-expression correlations, this repository uses *P. aeruginosa* gene expression compendia.
The compendia were created using publicly available RNA-seq data (see [preprint](https://doi.org/10.1101/2022.01.24.477642), [data repository](https://osf.io/s9gyu/), and [GitHub repository](https://github.com/hoganlab-dartmouth/pa-seq-compendia)) and represent gene counts for two *P. aeruginosa* strains, PAO1 and PA14.
The gene counts were derived from quasi-mapping thousands of SRA fastq files against both reference transcriptomes.
While each compendia contains gene counts for thousands of samples, 545 samples were identified as PA14 strain while 861 samples were identified as PAO1 strain (documented [here](https://github.com/greenelab/core-accessory-interactome/blob/master/data/metadata/SRA_annotations.tsv)). 

During compendia construction, Doing et al. noted that known gene-gene relationships could be recovered using correlation analysis of gene counts.
As depicted in the image below, genes that are part of the same operon and thus presumably co-transcribed on polycistronic transcripts have gene counts with high correlation coefficients after normalization.

![](https://i.imgur.com/jHkwLn6.png)

This repository tests whether a similar pattern can be recovered using FracMinHash abundances.
We initially focus on the PA14 compendium as it is smaller. 

## Things tested

1. Sub sampled k-mer abundance distributions. Do these distributions look similar to gene count distributions (negative binomial)? 
    + [20220322_compare_kmer_and_gene_count_dists.ipynb](notebooks/20220322_compare_kmer_and_gene_count_dists.ipynb)
2. K-mer coverage for genes and operons. Do most genes and most operons have at least one hash (k-mer) in the FracMinHash sketch?
    + [20220325_hashes_in_transcripts_and_operons.ipynb](notebooks/20220325_hashes_in_transcripts_and_operons.ipynb)
3. Normalization techniques and correlation analysis. Can we recover co-expression patterns using correlation analysis of normalized k-mer abundances?
    + DESeq2 median of ratios normalization: [20220328_correlations_mrnorm.ipynb](notebooks/20220328_correlations_mrnorm.ipynb)
    + DESeq2 variance stabilized transformation: [20220330_correlations_vstnorm.ipynb](notebooks/20220330_correlations_vstnorm.ipynb)
    + Seurat log normalization: [20220330_correlations_lognorm.ipynb](notebooks/20220330_correlations_lognorm.ipynb)
    + Seurat centered log ratio normalization: [20220330_correlations_clrnorm.ipynb](notebooks/20220330_correlations_clrnorm.ipynb)
    + Seurat counts per million normalization: [20220330_correlations_rcnorm.ipynb](notebooks/20220330_correlations_rcnorm.ipynb) 

## Results and conclusions

**Results**

See notebooks for more details.

+ K-mer abundance distributions are weird
+ Many genes don't have a k-mer, but most operons do
+ Many spurious correlations (e.g. high correlation coefficients for k-mers that are not in the same operons) remain after normalization, obscuring signal from intra-operon correlations. Median of ratios and centered log ratio normalization appear to have performed the best, however many spurious correlations remain after normalization. In situations where we won't have gene annotations available, this would make it difficult to determine which k-mers are actually part of the same operon and which are not. 

**Conclusions**

I don't think k-mer abundances can be used to recover co-expression patterns; too many spurious correlations remain after normalization.
Brainstorming possible reasons why this might be:
1. The RNAseq samples in the PA14 compendium may be too dissimilar, meaning that even at a k-mer size of 21 there are too few shared k-mers between many of samples to capture co-expression relationships. There is some evidence that this may be the case; the compendium was constructed using a k-mer size of 15 to seed the mappings. When the default k-mer size of 31 was used, the mapping rates were very low for most RNAseq samples. Dropping the k-mer size to 15 increased the mapping rate by 20%. This indicates that there are many SNPs in the RNAseq samples. Because k-mers have to exactly match to enable comparisons this may be artificially reducing gene-gene correlations.
2. K-mer count distributions are weird, so normalizing them doesn't lead to intended outcomes. 
3. K-mer abundances may not be reflective of gene abundances, since gene abundances are inferred from number of reads mapped and k-mer abundances are an exact count of that k-mer. 

**Next steps**

It is possible that there is some form of normalization that I haven't applied that would work to remove the spurious correlations.
If I become aware of any other methods, I'll give them a try.

## Running the code in this repository

This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
To setup the environment required for this repository, run the following commands:
```
conda env create --name core_acc --file environment.yml
conda activate core_acc
```


Snakemake can parallelize job submission and modulate resource usage (RAM, CPUs). 
We used the command below in a slurm cluster, but other cluster engines [are also supported](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

```
snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=60000 --cluster "sbatch -t 720 -J cacc -p shas -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```

Alternatively, snakemake can be executed without a cluster:
```
snakemake -j 2 --use-conda --rerun-incomplete -k -n
```
These parameters are described below:

+`-j 2` parallelizes the snakefile over two cores (drop to `-j 1` to only run one process at a time)
+ `--use-conda` tells snakemake to use conda to manage software environments for each rule
+ `--rerun-incomplete` tells snakemake to rerun rules when it thinks a file is incomplete, e.g. as may occur if a file is half-finished when a snakemake process is terminated.
+ `-k` indicates for the snakefile to keep running even if a rule fails. Snakemake will attempt to run all rules that don't depend on the output of the failed rule.
+ `-n` specifies a dry run. Remove this to actually execute the snakefile. The dry run is useful to make sure snakemake is running the desired rules for the desired number of times.

After running the code in the Snakefile, the notebooks in the `notebooks` directory analyze the output of the snakemake pipeline. 
The notebooks are sorted by date of analysis, and `20220324_make_wide_kmer_count_df.ipynb` needs to be run prior to the `*norm.ipynb` notebooks.  
