println """
    Genome:    ${params.genome}
    WGS Reads: ${params.reads}
    """.stripIndent()

reads     = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
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
