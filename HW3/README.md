deadline: 6.05
# Assignment

You have to continue the development of a pipeline: 
1. Take all processes and workflows from HW2 and import them to new pipeline 
2. Import process(es) for variant calling from nf-core and call variants on resulting .bam file (use bcftools or smth similar) 
3. Create configuration file that defines: 
   3.1. All required and optional input for pipeline. 
   3.2. 3 configuration profiles: 
   a) local – for execution on current machine, use path to local conda environment for dependencies handling (pass as parameter) 
   b) cluster – for execution on remote cluster, use conda .yml file or explicit conda package list for dependencies handling 
   c) container – for execution of all processes within containers, Docker or Singularity.
   3.3. Each configuration profile should also define appropriate hardware requirements 
5. Create and upload to docker hub all required for your pipeline docker images. 
6. Upload your pipeline code and conda .yml files to GitHub in separate branch – HW3.


# 1. Importing processes from HW2 with correction

Check current versions of tools for reproducibility:
```
$ conda --version
conda 26.1.1
$ docker --version
Docker version 29.4.0, build 9d7ad9f
$ docker ps
CONTAINER ID   IMAGE     COMMAND   CREATED   STATUS    PORTS     NAMES
$ nextflow -version

      N E X T F L O W
      version 23.10.0 build 5891
      created 15-10-2023 15:14 UTC (18:14 MSD)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

Create the HW3 structure and copy files from HW2:
```
$ cd ~
$ mkdir -p hw3/bin
$ cp ~/hw2_pipeline/main.nf hw3/
$ cp ~/hw2_pipeline/nextflow.config hw3/
$ cp ~/hw2_pipeline/plot_coverage.R hw3/bin/
$ ls -la hw3/
total 24
drwxr-xr-x  3 student student 4096 May  6 03:51 .
drwxr-x--- 24 student student 4096 May  6 03:51 ..
drwxr-xr-x  2 student student 4096 May  6 03:51 bin
-rw-r--r--  1 student student 4122 May  6 03:51 main.nf
-rw-r--r--  1 student student  160 May  6 03:51 nextflow.config

```

Let's split the .nf files into three parts: main.nf, hw2.nf, and main.nf in the "mpileup" folder (see step 2). Let's adjust the content of main.nf in the "hw3" folder:
```
$ nano ~/hw3/main.nf
```

```
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hw2_pipeline } from './modules/local/hw2'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup'

params.reads = null
params.sra_id = null
params.reference = null
params.outdir = "results"

if (params.reads && params.sra_id) {
    error "Please provide either '--reads' or '--sra_id', not both"
}
if (!params.reads && !params.sra_id) {
    error "Please provide either '--reads' or '--sra_id'"
}
if (!params.reference) {
    error "Please provide '--reference' to a genome fasta file"
}

workflow {
    // 1. Run the main pipeline (from HW2)
    hw2_pipeline()
    
    // 2. Prepare input data for bcftools mpileup
    //    Format: tuple(meta, bam, intervals_mpileup, intervals_call)
    bam_input = hw2_pipeline.out.bam.map { sample, bam, bai ->
        tuple([id: sample], bam, [], [])
    }
    
    // 3. Reference with index
    //    Create .fai if it doesn't exist
    ref_fai = file(params.reference + ".fai")
    if (!ref_fai.exists()) {
        "samtools faidx ${params.reference}".execute()
    }
    
    reference_input = hw2_pipeline.out.reference.map { ref ->
        tuple([id: ref.baseName], ref, ref_fai)
    }
    
    // 4. Call variants using nf-core module
    BCFTOOLS_MPILEUP(bam_input, reference_input, false)
}
```

hw2.nf (hw3\modules\local):
```
process fetch_reads {
    tag "$sra_id"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    
    input:
        val sra_id

    output:
        tuple val(sra_id), path("${sra_id}_1.fastq"), path("${sra_id}_2.fastq")

    script:
    """
    fasterq-dump --split-files $sra_id
    """
}

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/qc_${reads_type}", mode: 'copy'
    
    input:
        val reads_type
        tuple val(sample_id), path(r1), path(r2)

    output:
        path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}

process trimm {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq"
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_R1_p.fastq"), path("${sample_id}_R2_p.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        $r1 $r2 \
        ${sample_id}_R1_p.fastq ${sample_id}_R1_u.fastq \
        ${sample_id}_R2_p.fastq ${sample_id}_R2_u.fastq \
        LEADING:3 TRAILING:3 MINLEN:36
    """
}

workflow trimmed_qc_wf {
    take:
        trimmed_reads
    main:
        fastqc('trimmed', trimmed_reads)
    emit:
        fastqc.out
}

process alignment {
    tag "$sample_id"
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2), path(ref)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa index $ref
    bwa mem -t ${task.cpus} $ref $r1 $r2 | \
        samtools view -@ ${task.cpus} -bS - | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process plot_coverage {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam), path(bai)
        path r_script
    
    output:
        tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")
    
    script:
    """
    samtools depth $bam > ${sample_id}_depth.txt
    Rscript $r_script ${sample_id}
    """
}

workflow hw2_pipeline {
    main:
    if (params.sra_id) {
        raw_reads_ch = fetch_reads(params.sra_id)
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Please provide either '--sra_id' or '--reads'"
    }
    
    // raw QC
    fastqc('raw', raw_reads_ch)
    
    // trimming
    trimmed_reads_ch = trimm(raw_reads_ch)
    
    // trimmed QC (через workflow обертку)
    trimmed_qc_wf(trimmed_reads_ch)
    
    // reference
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    // mapping
    mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    bam_ch = alignment(mapping_in_ch)
    
    // coverage
    plot_coverage(bam_ch, file("${projectDir}/bin/plot_coverage.R"))
    
    emit:
    bam = bam_ch
    reference = ref_ch
}
```

Create environment.yml (the "hw3" folder):
```
$ nano ~/hw3/environment.yml
```

```
name: hw3-pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc=0.12.1
  - trimmomatic=0.39
  - sra-tools=3.0.0
  - spades=3.15.5
  - bwa=0.7.17
  - samtools=1.16.1
  - bcftools=1.17
  - r-base=4.2
  - r-ggplot2
```

Create Dockerfile:
```
$ nano ~/hw3/Dockerfile
```

```
FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a

COPY bin/plot_coverage.R /usr/local/bin/plot_coverage.R
RUN chmod +x /usr/local/bin/plot_coverage.R

ENV PATH /opt/conda/envs/hw3-pipeline/bin:$PATH

WORKDIR /data
```

Сreate README.md:
```
$ nano ~/hw3/README.md
```

```
# HW3: NGS Pipeline with variant calling

## Description
Pipeline for paired-end NGS data processing with variant calling.

## Quick Start

### With SRA accession:
```bash
nextflow run main.nf -profile local --sra_id SRR1175163 --reference dengue.fasta --outdir results
and so on.............
```
---

# 2. Importing nf-core processes for variant calling

```
$ pip install nf-core
$ nf-core modules install bcftools/mpileup
```

Write the content of main.nf in the "mpileup" folder:
```
process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(bam), path(intervals_mpileup), path(intervals_call)
    tuple val(meta2), path(fasta), path(fai)
    val save_mpileup
    
    output:
    path "${meta.id}.vcf.gz", emit: vcf
    path "${meta.id}.vcf.gz.tbi", emit: tbi
    path "${meta.id}.stats.txt", emit: stats

    script:
    """
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --output-type u \\
        ${bam} \\
        | bcftools call --multiallelic-caller --variants-only \\
        | bcftools view --output-file ${meta.id}.vcf.gz --output-type z
    
    tabix -p vcf ${meta.id}.vcf.gz
    bcftools stats ${meta.id}.vcf.gz > ${meta.id}.stats.txt
    """
}
```

---

# 3. Config with 3 profiles

Create a new nextflow.config with 3 profiles - local, cluster and container:
```
// =======================
// DEFAULT PARAMETERS
// =======================
params {
    reads = null
    sra_id = null
    reference = null
    outdir = "results"
    save_mpileup = false
    conda_env = null
}

// Default process settings
process {
    cpus = 4
    memory = '8 GB'
    time = '4h'
    errorStrategy = 'retry'
    maxRetries = 3
}

// =======================
// PROFILE 1: LOCAL (uses conda environment)
// =======================
profiles {
    local {
        process {
            cpus = 4
            memory = '2 GB'
            time = '8h'
            conda = params.conda_env ?: '/home/student/miniconda3/envs/hw3'
        }
        conda.enabled = true
        docker.enabled = false
    }
    
    // =======================
    // PROFILE 2: CLUSTER (for SLURM cluster)
    // =======================
    cluster {
        process {
            executor = 'slurm'
            cpus = 16
            memory = '64 GB'
            time = '48h'
            queue = 'normal'
            conda = "${projectDir}/environment.yml"
        }
        conda.enabled = true
        docker.enabled = false
    }
    
    // =======================
    // PROFILE 3: CONTAINER (Docker)
    // =======================
    container {
        process {
            cpus = 4
            memory = '16 GB'
            time = '8h'
            container = 'chaosintoorder/hw3-pipeline:latest'
        }
        conda.enabled = false
        docker.enabled = true
        docker.runOptions = '--rm --platform linux/amd64'
    }
}

// Results publishing settings for specific processes
process {
    withName: 'BCFTOOLS_MPILEUP' {
        publishDir = [ 
            path: { "${params.outdir}/variants" }, 
            mode: 'copy'
        ]
    }
    
    withName: 'plot_coverage' {
        publishDir = [ 
            path: { "${params.outdir}/coverage" }, 
            mode: 'copy' 
        ]
    }
    
    withName: 'fastqc' {
        publishDir = [ 
            path: { "${params.outdir}/qc" }, 
            mode: 'copy' 
        ]
    }
    
    withName: 'fetch_reads' {
        publishDir = [ 
            path: { "${params.outdir}/raw_reads" }, 
            mode: 'copy' 
        ]
    }
    
    withName: 'trimm' {
        publishDir = [ 
            path: { "${params.outdir}/trimmed_reads" }, 
            mode: 'copy',
            pattern: "*_p.fastq"
        ]
    }
}
```

Let's do a test run with Conda. But first, let's download the reference genome for the dengue virus. I downloaded the complete virus genome from the link: https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.

Create a conda environment:
```
$ conda env create -f environment.yml -n hw3
$ conda activate hw3
```

Check that everything is installed:
```
$ which fastqc trimmomatic bwa samtools bcftools spades
/home/student/miniconda3/envs/hw3/bin/fastqc
/home/student/miniconda3/envs/hw3/bin/trimmomatic
/home/student/miniconda3/envs/hw3/bin/bwa
/home/student/miniconda3/envs/hw3/bin/samtools
/home/student/miniconda3/envs/hw3/bin/bcftools
```

Run:
```
$ nextflow run main.nf -profile local --sra_id SRR1175163 --reference dengue_ref.fasta --outdir results_dengue
```

В процессе запуска я столкнулся с ошибкой: fetch_reads требует 16 GB RAM, но у меня доступно только 3.7 GB. Поэтому в файле nextflow.config я исправил настройки памяти.

During the launch process, I encountered an error: fetch_reads requires 16 GB RAM, but I only have 3.7 GB available. Therefore, in the nextflow.config file, I fixed the memory settings.

Clear the Nextflow cache, check each profile:
```
$ nextflow clean -f
$ nextflow config -profile local
$ nextflow config -profile cluster
$ nextflow config -profile container
```

And run again.
```
$ nextflow run main.nf -profile local --sra_id SRR1175163 --reference dengue_ref.fasta --outdir results_dengue
N E X T F L O W  ~  version 23.10.0
Launching `main.nf` [thirsty_einstein] DSL2 - revision: 01943ac036
[42/fb28af] Submitted process > hw2_pipeline:fetch_reads (SRR1175163)
[2e/d46550] Submitted process > hw2_pipeline:trimm (SRR1175163)
[a8/f9df7d] Submitted process > hw2_pipeline:fastqc (SRR1175163)
[a6/3a266b] Submitted process > hw2_pipeline:trimmed_qc_wf:fastqc (SRR1175163)
[e8/e28fdd] Submitted process > hw2_pipeline:alignment (SRR1175163)
[e2/2773c5] Submitted process > hw2_pipeline:plot_coverage (SRR1175163)
[6e/0e1dd9] Submitted process > BCFTOOLS_MPILEUP (SRR1175163)
```

---

# 4. Docker image on Docker Hub

Log in to Docker:
```
$ docker login
Authenticating with existing credentials... [Username: chaosintoorder]

i Info → To login with a different account, run 'docker logout' followed by 'docker login'


Login Succeeded
```

Build the image:
```
$ docker build -t chaosintoorder/hw3-pipeline:latest .
```

Check:
```
$ docker images | grep hw3-pipeline
WARNING: This output is designed for human readability. For machine-readable output, please use --format.
chaosintoorder/hw3-pipeline:latest   90d0ebc3c2bc       4.95GB         1.38GB
```

Push:
```
$ docker push chaosintoorder/hw3-pipeline:latest
The push refers to repository [docker.io/chaosintoorder/hw3-pipeline]
51521d07e119: Pushed
7708bc2ffc3e: Pushed
447fd0a1fb70: Pushed
3531af2bc2a9: Pushed
4f4fb700ef54: Pushed
3b4f41cf09c7: Pushed
18f3fce960ec: Pushed
fa93fe8f5025: Pushed
3bb161e09909: Pushed
3fd3d29769a2: Pushed
0bbe34ac03d8: Pushed
latest: digest: sha256:90d0ebc3c2bc5b63dd7a39a0334d53b894fcb0ceb78ff15d887162fb81588a5a size: 856
```
---

# 5. Pushing to GitHub

Done.