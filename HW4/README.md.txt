# HW4: Multi-sample NGS Pipeline with CSV Input and Variant Filtering

### 1. CSV Input (main.nf)
```nextflow
Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        tuple([id: row.sample, group: row.group], 
              [file(row.fastq_1), file(row.fastq_2)])
    }
```

### 2. Group-based Processing (main.nf)
```nextflow
bam_for_variants = results.map { group, sample_id, bam_ch ->
    bam_ch.map { sid, bam, bai ->
        tuple([id: sample_id, group: group], bam, [], [])
    }
}.flatten()
```

### 3. Variant Calling (bcftools/mpileup/main.nf)
```nextflow
bcftools mpileup --fasta-ref ${fasta} ${bam} |
    bcftools call --multiallelic-caller --variants-only |
    bcftools view --output-file ${prefix}.vcf.gz --output-type z
```

### 4. Join Back by Group (main.nf)
```nextflow
variants_by_group = BCFTOOLS_MPILEUP.out.vcf
    .map { meta, vcf -> tuple([group: meta.group], vcf) }
    .groupTuple()
```

### 5. Variant Filtering (main.nf)
```nextflow
process FILTER_VARIANTS {
    script:
    """
    bcftools view -i 'QUAL>=20' -Oz -o filtered.vcf.gz input.vcf
    """
}
```

### 6. Stub Mode (hw2.nf)
```nextflow
process DOWNLOAD_SRA {
    stub:
    """
    touch ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz
    """
}
```

## CSV Format (samplesheet.csv)
```csv
sample,group,fastq_1,fastq_2
dengue1,virus_A,SRR1175163,
dengue2,virus_A,SRR1175164,
dengue3,virus_B,SRR1175165,
```

## Usage

### Single sample (backward compatible)
```bash
nextflow run main.nf --accession SRR1175163 --reference dengue.fasta
```

### Multiple samples from CSV
```bash
nextflow run main.nf --samplesheet samplesheet.csv --reference dengue.fasta
```

### Stub mode for development
```bash
nextflow run main.nf --samplesheet samplesheet.csv --reference dengue.fasta -profile stub
```

## Output Structure
```
results/
├── reads/raw/           # Downloaded FASTQ
├── reads/trimmed/       # Trimmed reads
├── qc/                  # FastQC reports
├── mapping/             # BAM files
├── coverage/            # Coverage plots (PNG + TSV)
├── variants/            # Raw VCF files
└── filtered_variants/   # QUAL-filtered VCFs by group
    ├── virus_A_filtered/
    └── virus_B_filtered/
```

## Known Issues & Solutions

### Issue: Stub mode not recognized
**Error**: `Unknown configuration profile: 'stub'`

**Solution**: Update Nextflow to >=24.10.0
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Issue: Statements outside workflow
**Error**: `Statements cannot be mixed with script declarations`

**Solution**: Move `def input_count` inside workflow block

## Dependencies (envs/hw3.yml)
```yaml
name: hw3
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.2
  - r-ggplot2
  - fastqc
  - trimmomatic
  - sra-tools
  - bwa
  - samtools
  - bcftools
  - htslib
```