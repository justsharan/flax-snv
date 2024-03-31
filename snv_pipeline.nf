params.reads = "$projectDir/data/*_{1,2}.fastq"
params.outdir = "$projectDir/out"
params.snpeff = "$projectDir/snpEff"
params.genome = "$params.snpeff/data/flax/sequences.fa"
params.platform = "DNBSEQ"

println """
    F L A X    V A R I A N T S    P I P E L I N E
    =============================================
    Genome:     $params.genome
    WGS Reads:  $params.reads
    SnpEff:     $params.snpeff
    """.stripIndent()

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)
snpeff = Channel.fromPath(params.snpeff, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
    index = BUILD_INDEX(reference)
    alignments = BWA_ALIGN(reads, index.first())
    dedup = SAMTOOLS_PROCESS(alignments)
    variants = BCFTOOLS_CALL(dedup.collect(), index.first())
    snpeff = SNPEFF_ANNOTATE(variants, BUILD_SNPEFF_DB(snpeff).first())
    MULTIQC(fastqc.mix(BCFTOOLS_STATS(variants), snpeff).collect())
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

process BUILD_INDEX {
    input:
    path(genome)

    output:
    tuple path(genome), path("*")

    script:
    """
    bwa-mem2 index $genome
    samtools faidx $genome
    """
}

process BWA_ALIGN {
    tag "$sample_id"
    publishDir "$params.outdir/bam"

    input:
    tuple val(sample_id), path(reads)
    tuple path(genome), path("*")

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa-mem2 mem -R '@RG\\tID:$sample_id\\tPL:$params.platform\\tLB:$sample_id' $genome $reads | samtools sort > ${sample_id}.bam
    """
}

process SAMTOOLS_PROCESS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bamfile)

    output:
    path("${sample_id}.dedup.bam")

    script:
    """
    samtools collate -Ou $bamfile \
        | samtools fixmate -m - - \
        | samtools sort - - \
        | samtools markdup -d 100 - ${sample_id}.dedup.bam
    """
}

process BCFTOOLS_CALL {
    publishDir params.outdir

    input:
    path("*")
    tuple path(genome), path("*")

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -f $genome *.dedup.bam \
        | bcftools call -f GQ,GP --ploidy 2 -m -Ob -o variants.vcf.gz
    """
}

process BCFTOOLS_STATS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcffile)

    output:
    path("${sample_id}.txt")

    script:
    """
    bcftools stats $vcffile > ${sample_id}.txt
    """
}

process BUILD_SNPEFF_DB {
    input:
    path(snpEff)

    output:
    path(snpEff)

    script:
    """
    java -jar $snpEff/snpEff.jar build -gff3 -v flax -noCheckCds -noCheckProtein
    """
}

process SNPEFF_ANNOTATE {
    tag "$sample_id"
    publishDir params.outdir, mode:"copy"

    input:
    path(vcffile)

    output:
    path("annotated.vcf.gz")
    path("stats.csv")

    script:
    """
    java -jar $params.snpeff/snpEff.jar flax $vcffile -csvStats stats.csv > annotated.vcf.gz
    """
}
