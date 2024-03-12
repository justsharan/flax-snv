params.reads = "$projectDir/data/*_{1,2}.fastq"
params.genome = "$projectDir/data/*.fasta"
params.annotation = "$projectDir/data/*.gff3"
params.outdir = "$projectDir/out"
params.snpeff = "$projectDir/snpEff"
params.platform = "illumina"

println """
    F L A X    V A R I A N T S    P I P E L I N E
    =============================================
    Genome:     $params.genome
    WGS Reads:  $params.reads
    Annotation: $params.annotation
    SnpEff:     $params.snpeff
    """.stripIndent()

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)
snpeff = Channel.fromPath(params.snpeff, checkIfExists: true)

workflow {
    // Quality control on sequencing data
    fastqc = FASTQC(reads)
    // Build bwa index for reference genome
    index = BUILD_INDEX(reference)
    // Align samples to reference
    alignments = BWA_ALIGN_POSTPROCESS(reads, index.first())
    // Call variants and merge
    variants = BCFTOOLS_CALL(alignments, reference.first())
    merged = BCFTOOLS_MERGE(variants.collect())
    variant_stats = BCFTOOLS_STATS(variants)
    // Annotate variants
    snpeff_db = BUILD_SNPEFF_DB(snpeff)
    snpeff = SNPEFF_ANNOTATE(merged, snpeff_db.first())
    // Generate report
    MULTIQC(fastqc.mix(variant_stats, snpeff).collect())
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
    bwa index $genome
    samtools faidx $genome
    """
}

process BWA_ALIGN_POSTPROCESS {
    tag "$sample_id"
    publishDir "$params.outdir/bam", mode:"copy"

    input:
    tuple val(sample_id), path(reads)
    tuple path(genome), path("*")

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem -t $task.cpus $genome $reads -r "@RG\\tID:$sample_id\\tPL:$params.platform\\tLB:$sample_id" \
        | samblaster \
        | samclip --ref $genome \
        | samtools sort > ${sample_id}.bam
    """
}

process BCFTOOLS_CALL {
    tag "$sample_id"
    publishDir "$params.outdir/vcf", mode:"copy"

    input:
    tuple val(sample_id), path(bamfile)
    path(genome)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    bcftools mpileup -Ou -f $genome $bamfile | bcftools call -m -v --ploidy 2 -Ov -o ${bamfile.getBaseName()}.vcf.gz
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

process BCFTOOLS_MERGE {
    publishDir params.outdir, mode:"copy"
    input:
    path("*")

    output:
    path("variants.vcf.gz")

    script:
    """
    bcftools merge **/*.vcf.gz -o variants.vcf.gz
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
