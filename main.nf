println """
    Genome:    ${params.genome}
    WGS Reads: ${params.reads}
    """.stripIndent()

reads     = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
    BWA_INDEX(reference)
    BWA_ALIGNMENT(reads, BWA_INDEX.out.bwa_index.first())
    MULTIQC(fastqc.collect())
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q $reads
    """
}

process BWA_INDEX {
    publishDir "${params.outdir}/bwa_index", mode:"copy"

    input:
    file(genome)

    output:
    tuple path(genome), path("*"), emit: bwa_index

    script:
    """
    bwa index $genome
    """
}

process BWA_ALIGNMENT {
    tag "$sample_id to reference"

    input:
    tuple val(sample_id), path(reads)
    tuple path(genome), path("*")

    output:
    "*.bam"

    script:
    """
    bwa mem $genome $reads \
        | samtools view -b \
        > ${sample_id}.bam
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
