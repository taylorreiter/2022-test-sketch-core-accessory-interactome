import pandas as pd

m = pd.read_csv("inputs/metadata.tsv", sep = "\t")
SRX = list(m['experiment'])
KSIZE = [21]
STRAIN = ['pao1', 'pa14']

rule all:
    input:
        expand("outputs/srx_sourmash_sketch_filtered_acc/{srx}_k{ksize}.sig", srx = SRX, ksize = KSIZE),
        expand("outputs/srx_sourmash_sketch_filtered_core/{srx}_k{ksize}.sig", srx = SRX, ksize = KSIZE)


rule convert_srx_to_srr:
    output: "outputs/srx_to_srr/{srx}.tsv"
    conda: "envs/pysradb.yml"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    pysradb srx-to-srr --detailed --save-to {output} {wildcards.srx}
    '''
    
rule srx_sketch:
    input: "outputs/srx_to_srr/{srx}.tsv"
    output: "outputs/srx_sourmash_sketch/{srx}.sig"
    threads: 1
    resources: mem_mb = 1000
    run:
        srx_to_srr = pd.read_csv(input, sep = "\t")
        row = srx_to_srr.loc[srx_to_srr['experiment_accession'] == wildcards.srx]
        srr = row['run_accession'].values[0]
        shell("fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {srr} |"
              "sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.srx} -o {output} -")

rule srx_filter_sketches:
    input: "outputs/srx_sourmash_sketch/{srx}.sig"
    output: "outputs/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig filter --name {wildcards.srx} -o {output} --min-abundance 2 -k {wildcards.ksize} --dna {input}
    '''

##########################################################
## Define core transcriptome
##########################################################

rule txome_download_pa14:
    output: "inputs/txomes/pa14.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_16_collection/pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625/cdna/Pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625.ASM1462v1.cdna.all.fa.gz
    '''

rule txome_download_pao1:
    output: "inputs/txomes/pao1.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/cdna/Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.cdna.all.fa.gz
    '''

rule txome_sketch:
    input: "inputs/txomes/{strain}.cdna.all.fa.gz"
    output: "outputs/txomes_sourmash_sketch/{strain}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.strain} -o {output} {input}
    '''

rule txome_intersect_strain:
    input: expand("outputs/txomes_sourmash_sketch/{strain}.sig", strain = STRAIN)
    output: "outputs/txomes_sourmash_sketch_core/core_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig intersect --dna -o {output} -k {wildcards.ksize} {input}
    '''

#############################################################
## Use core transcriptome to identify "core" and "acc" hashes
## in signatures
#############################################################

rule srx_identify_core:
    input: 
        core = "outputs/txomes_sourmash_sketch_core/core_k{ksize}.sig",
        srx =  "outputs/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    output: "outputs/srx_sourmash_sketch_filtered_core/{srx}_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig intersect -A {input.srx} --dna -o {output} -k {wildcards.ksize} {input}
    '''

rule srx_identify_acc:
    input: 
        core = "outputs/txomes_sourmash_sketch_core/core_k{ksize}.sig",
        srx =  "outputs/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    output: "outputs/srx_sourmash_sketch_filtered_acc/{srx}_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig subtract --dna -o {output} -k {wildcards.ksize} {input.srx} {input.core}
    '''

############################################################
## create hash:annotation map
############################################################

rule txome_sketch_singleton:
    input: "inputs/txomes/{strain}.cdna.all.fa.gz"
    output: "outputs/txomes_sourmash_sketch_singleton/{strain}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name-from-first --singleton -o {output} {input}
    '''
