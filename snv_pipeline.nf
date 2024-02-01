params.reads = "$projectDir/data/*_{1,2}.fastq"
params.genome = "$projectDir/Phytozome/PhytozomeV10/Lusitatissimum/assembly/Lusitatissimum_200_BGIv1.0.fa"
params.annotation = "$projectDir/Phytozome/PhytozomeV10/Lusitatissimum/annotation/Lusitatissimum_200_v1.0.gene_exons.gff3"
params.outdir = "$projectDir/out"

println """
    Genome:     $params.genome
    WGS Reads:  $params.reads
    Annotation: $params.annotation
    """.stripIndent()

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
    index = BUILD_BWA_INDEX(reference)
    bamfile = BWA_ALIGN_SORT(reads, index.first())
    bamfile_stats = BAMTOOLS_STATS(bamfile)
    index_ref = SAMTOOLS_REF_INDEX(reference)
    vcffile = BCFTOOLS_CALL(index_ref.first(), bamfile)
    vcf_stats = BCFTOOLS_STATS(vcffile)
    MULTIQC(fastqc.mix(bamfile_stats, vcf_stats).collect())
}

process FASTQC {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    fastqc -q $reads
    """
}

process MULTIQC {
    publishDir params.outdir, mode:"copy"

    input:
    path("*")

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}

process BUILD_BWA_INDEX {
    publishDir "$params.outdir/bwa_index", mode:"copy"

    input:
    path(genome)

    output:
    tuple path(genome), path("*")

    script:
    """
    bwa index $genome
    """
}

process BWA_ALIGN_SORT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    tuple path(genome), path("*")

    output:
    path("${sample_id}.bam")

    script:
    """
    bwa mem -t $task.cpus $genome $reads | samtools sort -o ${sample_id}.bam
    """
}

process BAMTOOLS_STATS {
    tag "Getting stats on ${bamfile.getName()}"

    input:
    path(bamfile)

    output:
    path("${bamfile.getBaseName()}.txt")

    script:
    """
    bamtools stats -in $bamfile > ${bamfile.getBaseName()}.txt
    """
}

process SAMTOOLS_REF_INDEX {
    tag "Building genome .fai"

    input:
    path(genome)

    output:
    tuple path(genome), path("${genome.getName()}.fai")

    script:
    """
    samtools faidx $genome
    """
}

process BCFTOOLS_CALL {
    tag "${bamfile.getBaseName()}"
    publishDir "$params.outdir/vcf", mode:"copy"

    input:
    tuple path(genome), path("*")
    path(bamfile)

    output:
    path("${bamfile.getBaseName()}.vcf.gz")

    script:
    """
    bcftools mpileup -Ou -f $genome $bamfile | bcftools call -m -v --ploidy 2 -Ov -o ${bamfile.getBaseName()}.vcf.gz
    """
}

process BCFTOOLS_STATS {
    tag "Getting stats on ${vcffile.getName()}"

    input:
    path(vcffile)

    output:
    path("${vcffile.getBaseName()}.txt")

    script:
    """
    bcftools stats $vcffile > ${vcffile.getBaseName()}.txt
    """
}
