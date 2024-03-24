params.reads = "$projectDir/data/*_{1,2}.fastq"
params.genome = "$projectDir/data/*.fasta"
params.outdir = "$projectDir/out"
params.snpeff = "$projectDir/snpEff"
params.platform = "illumina"

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
    variants = BCFTOOLS_CALL(alignments, reference.first())
    variant_indexes = BCFTOOLS_INDEX(variants)
    merged = variants.map { v -> v[1] }.mix(variant_indexes).collect() | BCFTOOLS_MERGE
    snpeff_db = BUILD_SNPEFF_DB(snpeff)
    snpeff = SNPEFF_ANNOTATE(merged, snpeff_db.first())
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
    bwa-mem2 mem $genome $reads | samtools sort > ${sample_id}.bam
    """
}

process BCFTOOLS_CALL {
    tag "$sample_id"
    publishDir "$params.outdir/vcf"

    input:
    tuple val(sample_id), path(bamfile)
    path(genome)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    bcftools mpileup -Ou -f $genome $bamfile \
        | bcftools call -m -v --ploidy 2 -Ov -o ${sample_id}.vcf.gz
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

process BCFTOOLS_INDEX {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcffile)

    output:
    path("*.tbi")

    script:
    """
    bcftools index --threads $task.cpus -t $vcffile
    """
}

process BCFTOOLS_MERGE {
    publishDir params.outdir
    input:
    path("*")

    output:
    path("variants.vcf.gz")

    script:
    """
    bcftools merge *.vcf.gz -o variants.vcf.gz
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
    path(snpEff)

    output:
    path("annotated.vcf.gz")
    path("stats.csv")

    script:
    """
    java -jar $snpEff/snpEff.jar flax $vcffile -csvStats stats.csv > annotated.vcf.gz
    """
}
