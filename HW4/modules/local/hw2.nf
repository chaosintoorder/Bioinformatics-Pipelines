nextflow.enable.dsl = 2

process DOWNLOAD_SRA {
    tag "$meta.id"
    publishDir "${params.outdir}/reads/raw", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}_*.fastq.gz")

    script:
    """
    fasterq-dump --outdir . ${meta.id}
    gzip -f ${meta.id}_1.fastq ${meta.id}_2.fastq
    """
}

process RUN_QC {
    tag "${reads_type}:${meta.id}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(reads_type), val(meta), path(reads)

    output:
    path "${reads_type}_${meta.id}_fastqc"

    script:
    """
    mkdir -p ${reads_type}_${meta.id}_fastqc
    fastqc -o ${reads_type}_${meta.id}_fastqc ${reads}
    """
}

process TRIMM {
    tag "$meta.id"
    publishDir "${params.outdir}/reads/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R*_trimmed.fq.gz")

    script:
    """
    ADAPTERS=\$(find \${CONDA_PREFIX:-/opt/conda} /usr/local -name TruSeq3-PE.fa -print -quit 2>/dev/null)
    trimmomatic PE ${reads[0]} ${reads[1]} \\
        ${meta.id}_R1_trimmed.fq.gz ${meta.id}_R1_unpaired.fq.gz \\
        ${meta.id}_R2_trimmed.fq.gz ${meta.id}_R2_unpaired.fq.gz \\
        ILLUMINACLIP:\$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process BWA_ALIGN {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("${meta.id}.sam")

    script:
    """
    bwa index ${reference}
    bwa mem ${reference} ${reads[0]} ${reads[1]} > ${meta.id}.sam
    """
}

process SORT_BAM {
    tag "$meta.id"
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), path("${meta.id}.flagstat.txt")

    script:
    """
    samtools sort -o ${meta.id}.sorted.bam ${sam}
    samtools index ${meta.id}.sorted.bam
    samtools flagstat ${meta.id}.sorted.bam > ${meta.id}.flagstat.txt
    """
}

process PLOT_COVERAGE {
    tag "$meta.id"
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat)

    output:
    path "${meta.id}.coverage.tsv"
    path "${meta.id}.coverage_summary.tsv"
    path "${meta.id}.coverage.png"

    script:
    """
    samtools depth -aa ${bam} > ${meta.id}.coverage.tsv
    Rscript ${projectDir}/bin/plot_coverage.R ${meta.id}
    """
}

workflow HW2_PIPELINE {
    main:
    if( params.samplesheet ) {
        input_reads = Channel.fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple([id: row.sample, group: row.group], [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)])
            }
    } else if( params.accession ) {
        input_reads = DOWNLOAD_SRA(Channel.of([id: params.accession, group: params.accession]))
    } else {
        input_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
            .map { reads_label, reads -> tuple([id: reads_label.toString(), group: reads_label.toString()], reads) }
    }

    raw_qc_input = input_reads.map { meta, reads -> tuple('raw', meta, reads) }
    trimmed_reads = TRIMM(input_reads)
    trimmed_qc_input = trimmed_reads.map { meta, reads -> tuple('trimmed', meta, reads) }
    RUN_QC(raw_qc_input.mix(trimmed_qc_input))

    if( params.reference ) {
        reference = Channel.value(file(params.reference, checkIfExists: true))
        reference_for_mapping = trimmed_reads.combine(reference).map { meta, reads, ref ->
            tuple(meta, reads, ref)
        }
        reference_for_variants = trimmed_reads.combine(reference).map { meta, reads, ref ->
            tuple(meta, ref)
        }
    } else {
        error "Reference genome is required"
    }

    aligned_reads = BWA_ALIGN(reference_for_mapping)
    sorted_bam = SORT_BAM(aligned_reads)
    PLOT_COVERAGE(sorted_bam)

    emit:
    bam = sorted_bam
    reference = reference_for_variants
}