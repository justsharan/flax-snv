params.reads = "$projectDir/data/*_{1,2}.fastq"
params.genome = "$projectDir/data/*.fasta"
params.annotation = "$projectDir/data/*.gff3"
params.outdir = "$projectDir/out"
params.snpeff = "$projectDir/snpEff"

println """
    Genome:     $params.genome
    WGS Reads:  $params.reads
    Annotation: $params.annotation
    SnpEff:     $params.snpeff
    """.stripIndent()

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)
snpeff = Channel.fromPath(params.snpeff, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
    index = BUILD_BWA_INDEX(reference)
    bamfile = BWA_ALIGN_SORT(reads, index.first())
    bamfile_stats = BAMTOOLS_STATS(bamfile)
    index_ref = SAMTOOLS_REF_INDEX(reference)
    vcffile = BCFTOOLS_CALL(index_ref.first(), bamfile)
    vcf_stats = BCFTOOLS_STATS(vcffile)
    snpeff_db = BUILD_SNPEFF_DB(snpeff)
    snpeff = SNPEFF_ANNOTATE(vcffile, snpeff_db.first())
    MULTIQC(fastqc.mix(bamfile_stats, vcf_stats, snpeff).collect())
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
    tag "${vcffile.getSimpleName()}"
    publishDir "$params.outdir/annotated", mode:"copy"

    input:
    path(vcffile)
    path(snpEff)

    output:
    path("${vcffile.getSimpleName()}.ann.vcf.gz")
    path("${vcffile.getSimpleName()}.csv")

    script:
    """
    java -jar $snpEff/snpEff.jar flax $vcffile -csvStats ${vcffile.getSimpleName()}.csv > ${vcffile.getSimpleName()}.ann.vcf.gz
    """
}
