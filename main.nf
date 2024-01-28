println """
    Genome:    ${params.genome}
    WGS Reads: ${params.reads}
    """.stripIndent()

reads     = Channel.fromFilePairs(params.reads, checkIfExists: true)
reference = Channel.fromPath(params.genome, checkIfExists: true)

workflow {
    fastqc = FASTQC(reads)
    BWA_INDEX(reference)
    bamfile = BWA_ALIGNMENT(reads, BWA_INDEX.out.bwa_index.first())
    VARIANT_CALLING(reference, bamfile.collect())
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
    path("${sample_id}.bam"), emit: bamfile

    script:
    """
    bwa mem $genome $reads \
        | samtools view -b \
        > ${sample_id}.bam
    """
}

process VARIANT_CALLING {
  tag "${bamfile.getBaseName()}"
  publishDir "${params.outdir}/bcf", mode:"copy"

  input:
  file(genome)
  file(bamfile)

  output:
  path("${bamfile.getBaseName()}.bcf")

  script:
  """
    bcftools mpileup -Ou -f $genome $bamfile | bcftools call -mv --ploidy 2 -Ob -o ${bamfile.getBaseName()}.bcf
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
