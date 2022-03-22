import pandas as pd

outdir = "/projects/treiter@xsede.org/2022-test-sketch-core-accessory-interactome/outputs"
indir = "/projects/treiter@xsede.org/2022-test-sketch-core-accessory-interactome/inputs"

m = pd.read_csv("inputs/metadata.tsv", sep = "\t")
SRX = list(m['experiment'])
KSIZE = [21]
STRAIN = ['pao1', 'pa14']

rule all:
    input:
        expand(f"{outdir}/srx_sourmash_sketch_filtered_acc/{{srx}}_k{{ksize}}.sig", srx = SRX, ksize = KSIZE),
        expand(f"{outdir}/srx_sourmash_sketch_filtered_core/{{srx}}_k{{ksize}}.sig", srx = SRX, ksize = KSIZE),
        expand(f"{outdir}/srx_sourmash_sketch_filtered_csv/{{srx}}_k{{ksize}}.csv", srx = SRX, ksize = KSIZE)


rule convert_srx_to_srr:
    output: outdir + "/srx_to_srr/{srx}.tsv"
    conda: "envs/pysradb.yml"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    pysradb srx-to-srr --detailed --saveto {output} {wildcards.srx}
    '''
    
rule srx_sketch:
    input: outdir + "/srx_to_srr/{srx}.tsv"
    output: outdir + "/srx_sourmash_sketch/{srx}.sig"
    threads: 1
    resources: mem_mb = 1000
    run:
        srx_to_srr = pd.read_csv(str(input), sep = "\t")
        row = srx_to_srr.loc[srx_to_srr['experiment_accession'] == wildcards.srx]
        srr = row['run_accession'].values[0]
        shell("fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {srr} |"
              "sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.srx} -o {output} -")

rule srx_filter_sketches:
    input: outdir + "/srx_sourmash_sketch/{srx}.sig"
    output: outdir + "/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig filter --name {wildcards.srx} -o {output} --min-abundance 2 -k {wildcards.ksize} --dna {input}
    '''

rule convert_srx_filtered_sketches_to_csv:
    input:outdir + "/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    output: outdir + "/srx_sourmash_sketch_filtered_csv/{srx}_k{ksize}.csv"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

##########################################################
## Define core transcriptome
##########################################################

rule txome_download_pa14:
    output: indir + "/txomes/pa14.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_16_collection/pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625/cdna/Pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625.ASM1462v1.cdna.all.fa.gz
    '''

rule txome_download_pao1:
    output: indir + "/txomes/pao1.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/cdna/Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.cdna.all.fa.gz
    '''

rule txome_sketch:
    input: indir + "/txomes/{strain}.cdna.all.fa.gz"
    output: outdir + "/txomes_sourmash_sketch/{strain}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.strain} -o {output} {input}
    '''

rule txome_intersect_strain:
    input: expand(f"{outdir}/txomes_sourmash_sketch/{{strain}}.sig", strain = STRAIN)
    output: outdir + "/txomes_sourmash_sketch_core/core_k{ksize}.sig"
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
        core = outdir + "/txomes_sourmash_sketch_core/core_k{ksize}.sig",
        srx =  outdir + "/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    output: outdir + "/srx_sourmash_sketch_filtered_core/{srx}_k{ksize}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig intersect -A {input.srx} --dna -o {output} -k {wildcards.ksize} {input}
    '''

rule srx_identify_acc:
    input: 
        core = outdir + "/txomes_sourmash_sketch_core/core_k{ksize}.sig",
        srx =  outdir + "/srx_sourmash_sketch_filtered/{srx}_k{ksize}.sig"
    output: outdir + "/srx_sourmash_sketch_filtered_acc/{srx}_k{ksize}.sig"
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
    input: indir + "/txomes/{strain}.cdna.all.fa.gz"
    output: outdir + "/txomes_sourmash_sketch_singleton/{strain}.sig"
    threads: 1
    resources: mem_mb = 1000
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name-from-first --singleton -o {output} {input}
    '''
