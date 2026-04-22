params.reads = null
params.sra_id = null
params.reference = null

// =======================
// PROCESS: FETCH READS
// =======================
process fetch_reads {
    input:
        val sra_id

    output:
        tuple val(sra_id), path("${sra_id}_1.fastq"), path("${sra_id}_2.fastq")

    script:
    """
    fasterq-dump ${sra_id}
    """
}

// =======================
// PROCESS: QC
// =======================
process run_qc_raw {
    input:
        tuple val(id), path(r1), path(r2)

    output:
        path "raw_qc/", type: 'dir'

    script:
    """
    mkdir raw_qc
    fastqc -o raw_qc $r1 $r2
    """
}

process run_qc_trimmed {
    input:
        tuple val(id), path(r1), path(r2)

    output:
        path "trimmed_qc/", type: 'dir'

    script:
    """
    mkdir trimmed_qc
    fastqc -o trimmed_qc $r1 $r2
    """
}
// =======================
// PROCESS: TRIMMING
// =======================
process trimm {
    input:
        tuple val(id), path(r1), path(r2)

    output:
        tuple val(id), path("R1_p.fq"), path("R2_p.fq")

    script:
    """
    trimmomatic PE \
        $r1 $r2 \
        R1_p.fq R1_u.fq \
        R2_p.fq R2_u.fq \
        LEADING:3 TRAILING:3 MINLEN:36
    """
}

// =======================
// PROCESS: ASSEMBLY
// =======================
process assembly {
    input:
        tuple val(id), path(r1), path(r2)

    output:
        path "assembly/contigs.fasta"

    script:
    """
    spades.py -1 $r1 -2 $r2 -o assembly
    """
}

// =======================
// PROCESS: MAPPING
// =======================
process mapping {
    input:
        tuple val(id), path(r1), path(r2)
        path ref

    output:
        path "aligned.bam"

    script:
    """
    bwa index $ref
    bwa mem $ref $r1 $r2 | samtools view -bS - | samtools sort -o aligned.bam
    """
}

// =======================
// PROCESS: COVERAGE
// =======================
process coverage {
    input:
        path bam

    output:
        path "coverage.txt"

    script:
    """
    samtools depth $bam > coverage.txt
    """
}

// =======================
// WORKFLOW
// =======================
workflow pipeline {

    take:
        reads_ch

    main:
        qc1 = run_qc_raw(reads_ch)

        trimmed = trimm(reads_ch)

        qc2 = run_qc_trimmed(trimmed)

        if (params.reference) {
            ref = channel.of(file(params.reference))
        } else {
            ref = assembly(trimmed).map{ it }
        }

        mapped = mapping(trimmed, ref)

        cov = coverage(mapped)

    emit:
        qc1 = qc1
        qc2 = qc2
        coverage = cov
}

// =======================
// ENTRY WORKFLOW
// =======================
workflow {

    main:

    if (params.sra_id) {
        reads = fetch_reads(params.sra_id)
    } else {
        reads = channel.fromFilePairs(params.reads)
    }

    pipeline(reads)

    publish:
        qc1 = pipeline.out.qc1
        qc2 = pipeline.out.qc2
        coverage = pipeline.out.coverage
}

// =======================
// OUTPUT
// =======================
output {

    qc1 {
        path 'qc_initial'
    }

    qc2 {
        path 'qc_trimmed'
    }

    coverage {
        path 'coverage'
    }
}
