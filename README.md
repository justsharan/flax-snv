This repository contains a pipeline and some data analysis scripts I've written for a variant calling project on EMS-mutagenized flax. Here are some instructions on how to execute the pipeline.

## Setup

1. Clone the repository

```sh
git clone https://github.com/justsharan/flax-snv.git
cd flax-snv
```

2. Install and setup SnpEff

[SnpEff](https://pcingola.github.io/SnpEff/) is the tool I have used to integrate the variant calling data with genomic information and annotate the vcf file. It can be installed as such. The `-n` flag preserves the config I have included in this repository.

```sh
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip -n snpEff_latest_core.zip
```

3. Provide flax FASTA files for SnpEff

```sh
mkdir -p snpEff/data/flax/
```

In that location, SnpEff expects 4 files with these exact filenames. The content of the files is self-explanatory.
* `cds.fa`
* `genes.gff`
* `protein.fa`
* `sequences.fa`

For this project in particular, the requisite files were obtained from these sources:
* https://phytozome-next.jgi.doe.gov/info/Lusitatissimum_v1_0
* You, F. M., Xiao, J., Li, P., Yao, Z., Jia, G., He, L., Zhu, T., Luo, M. C., Wang, X., Deyholos, M. K., & Cloutier, S. (2018). Chromosome-scale pseudomolecules refined by optical, physical and genetic maps in flax. The *Plant journal : for cell and molecular biology*, *95*(2), 371–384. https://doi.org/10.1111/tpj.13944

My pipeline also uses the `sequences.fa` file as the reference for aligning the WGS reads. This behavior can be changed.

4. Install dependencies

This pipeline requires the following dependencies:
* nextflow
* fastqc
* multiqc
* bwa-mem2
* samtools
* bcftools

If you are using conda, you can do the following:

```sh
conda create -n flaxsnv nextflow fastqc multiqc bwa-mem2 samtools bcftools
conda activate flaxsnv
```

## Running the pipeline

> It is recommended to run this within `tmux`.

The basic command is `nextflow run snv_pipeline.nf`. It takes the following custom arguments:

|Option|Default|Description|
|:-----|:------|:----------|
|`reads`|`$PWD/data/*_{1,2}.fastq`|The location of the WGS reads. The default option looks at all paired-end reads in the data folder. For example, `SRR12345_1.fastq` and `SRR12345_2.fastq` will be grouped together.|
|`outdir`|`$PWD/out`|The output directory|
|`snpeff`|`$PWD/snpEff`|The location of the SnpEff install. Must include the config file and data directory as specified in setup.|
|`genome`|`$snpeff/data/flax/sequences.fa`|The reference genome to align the sequences to. Uses the `sequences.fa` file from the SnpEff setup.|
|`platform`|`DNBSEQ`|Used in the read groups of the alignment files. Necessary for processing of alignment files by samtools. [[Read more]](https://samtools.github.io/hts-specs/SAMv1.pdf)|

These can be passed in as arguments to the nextflow command. To change the location of the reads, for example, you can do:

```sh
nextflow run snv_pipeline.nf --reads "some/other/path/*_{1,2}.fq.gz"
```

Other built-in Nextflow options you may find useful:

* `-preview`: Perform a dry run of the pipeline. Useful to check if your pipeline code has any syntax issues.
* `-with-dag`: Produce a flowchart of the different pipeline steps.
* `-with-report`: Produce a report of memory/CPU usage for each task, etc.
* `-with-timeline`: Produce a timeline of when each command was ran and how long it took.

[[More info]](https://www.nextflow.io/docs/latest/tracing.html)
