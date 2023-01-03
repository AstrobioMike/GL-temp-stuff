# Bioinformatics pipeline for Methylation Sequencing (Methyl-Seq) data

> **This document holds an overview and some example commands of how GeneLab processes bisulfite sequencing (methylseq) datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**

---

**Date:**   
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee (GeneLab Analysis Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

---

# Table of contents

- [Bioinformatics pipeline for Methylation Sequencing (Methyl-Seq) data](#bioinformatics-pipeline-for-methylation-sequencing-methyl-seq-data)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [2. Adapter trimming/quality filtering](#2-adapter-trimmingquality-filtering)
    - [If not RRBS or if RRBS using MseI digestion](#if-not-rrbs-or-if-rrbs-using-msei-digestion)
    - [If RRBS with MspI digestion](#if-rrbs-with-mspi-digestion)
    - [If RRBS with NuGEN ovation kit](#if-rrbs-with-nugen-ovation-kit)
      - [First adapter-trimming/quality-filtering with trimgalore](#first-adapter-trimmingquality-filtering-with-trimgalore)
      - [Now running NuGEN-specific script](#now-running-nugen-specific-script)
  - [3. Filtered/Trimmed Data QC](#3-filteredtrimmed-data-qc)
    - [3a. Trimmed Data QC](#3a-trimmed-data-qc) 
    - [3b. Compile Trimmed Data QC](#3b-compile-trimmed-data-qc)
  - [4. Alignment](#4-alignment)
    - [4a. Generate reference](#4a-generate-reference)
    - [4b. Align](#4b-align)
    - [4c. Sort Alignment Files](#4c-sort-alignment-files)
  - [5. Alignment QC](#5-alignment-qc)
  - [6. Deduplicate (skip if data are RRBS)](#6-deduplicate-skip-if-data-are-rrbs)
  - [7. Extract methylation calls](#7-extract-methylation-calls)
  - [8. Generate individual sample report](#8-generate-individual-sample-report)
  - [9. Generate combined summary report](#9-generate-combined-summary-report)
  - [10. Generate reference genome annotation information](#10-generate-reference-genome-annotation-information)
    - [10a. GTF to BED conversion](#10a-gtf-to-bed-conversion)
    - [10b. Making a mapping file of genes to transcripts](#10b-making-a-mapping-file-of-genes-to-transcripts)
  - [11. Differential methylation analysis](#11-differential-methylation-analysis)

---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|FastQC| 0.11.9 |[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC| 1.12 |[https://multiqc.info/](https://multiqc.info/)|
|TrimGalore!| 0.6.7 |[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|Bismark| 0.23.1 |[https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)|
|bowtie2| 2.4.4 |[https://github.com/BenLangmead/bowtie2#overview](https://github.com/BenLangmead/bowtie2#overview)|
|samtools| 1.16.1 |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
|qualimap| 2.2.2d |[http://qualimap.conesalab.org/](http://qualimap.conesalab.org/)|
|methylKit|1.20.0|[https://bioconductor.org/packages/release/bioc/html/methylKit.html](https://bioconductor.org/packages/release/bioc/html/methylKit.html)|

---

# General processing overview with example commands

> Exact processing commands for specific datasets that have been released are provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 
>
> Output files listed in **bold** below are included with the processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC

```bash
fastqc -o raw_fastqc_output *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input data:**

* *raw.fastq.gz (raw reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)

<br>

### 1b. Compile Raw Data QC

```bash
multiqc --interactive -o raw_multiqc_data -n raw_multiqc -z raw_fastqc_output/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* raw_fastqc_output/*fastqc.zip (FastQC output data from [Step 1a](#1a-raw-data-qc))

**Output data:**

* **raw_multiqc.html** (multiqc output html summary)
* **raw_multiqc_data.zip** (zipped directory containing multiqc output data)

<br>  

---

## 2. Adapter trimming/quality filtering
See `trim_galore --help` [menu](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3035) for more info on any of the below.

<br>

### If not RRBS or if RRBS using MseI digestion
Note that the `--rrbs` option is **not** appropriate when RRBS (reduced representation bisulfite sequencing) libraries were prepared with MseI digestion (see `trim_galore --help` menu [(starting at this line)](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3337).

**Single-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  --output_dir /path/to/TrimGalore/output/directory \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

### If RRBS with MspI digestion
Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end; see `trim_galore --help` menu [e.g., here](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3315)). 

**Single-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  --rrbs \
  --output_dir /path/to/TrimGalore/output/directory \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  --rrbs \
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

### If RRBS with NuGEN ovation kit
Libraries prepared with the NuGEN ovation kit need to be procesed with an additional script provided by the company's [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq). 

Following their instructions, we first run an adapter-trimming/quality-filtering step with trimgalore. Note that the `--rrbs` option is not appropriate to pass to trimgalore when this kit is used (see `trim_galore --help` menu [(starting at this line)](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3329). Then we utilize the company's script to remove the random diversity sequences added by the kit. 

#### First adapter-trimming/quality-filtering with trimgalore

**Single-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  -a AGATCGGAAGAGC \
  --output_dir /path/to/TrimGalore/output/directory \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores 4 \
  --phred33 \
  --paired \
  -a AGATCGGAAGAGC \
  -a2 AAATCAAAAAAAC \
  --output_dir /path/to/TrimGalore/output/directory \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

#### Now running NuGEN-specific script

The NuGEN-specific script can be downloaded from their [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq) with the following:

```bash
curl -LO https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/trimRRBSdiversityAdaptCustomers.py
```

**Single-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py \
  -1 sample-1_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_trimmed.fastq_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py \
  -1 sample-1_R1_trimmed.fastq.gz \
  -2 sample-1_R2_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_trimmed.fastq_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_trimmed.fastq_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

**Parameter Definitions for `trim_galore`:**  

* `--gzip` - gzip compress the output(s)
* `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
* `--cores` - number of cores to use (this value is dependent on the number of threads available on the system running trim galore)
* `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
* `--rrbs` - specific trimming suitable for RRBS data generated with MspI digestion only
* `-a` - specific adapter sequence to be trimmed off of forward or single reads (applicable for libraries prepared with the NuGEN ovation kit)
* `-a2` - specific adapter sequence to be trimmed off of reverse reads (applicable for libraries prepared with the NuGEN ovation kit)
* `--paired` - specifies data are paired-end
* `--output_dir` - the output directory to store results
* positional arguments represent the input read files, 2 of them if paired-end data


**Parameter Definitions for `trimRRBSdiversityAdaptCustomers.py `:**  

- `-1` - forward or single input reads
- `-2` - reverse reads if paired-end data

**Input Data:**

* \*fastq.gz (raw reads)

**Output Data:**

* **gzip-compressed fastq files** (adapter-trimmed/quality-filtered reads)
* **\*trimming_report.txt** (trimming report)

<br>

---

## 3. Filtered/Trimmed Data QC

<br>

### 3a. Trimmed Data QC

```bash
fastqc -o trimmed_fastqc_output/ *trimmed.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*trimmed.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input data:**

* *trimmed.fastq.gz (filtered/trimmed reads, output from [Step 2](#2-adapter-trimmingquality-filtering))

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)

<br>

### 3b. Compile Trimmed Data QC

```bash
multiqc --interactive -o trimmed_multiqc_data -n trimmed_multiqc -z trimmed_fastqc_output/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`trimmed_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* trimmed_fastqc_output/*fastqc.zip (FastQC output data from [Step 3a](#3a-trimmed-data-qc))

**Output data:**

* **trimmed_multiqc.html** (multiqc output html summary)
* **trimmed_multiqc_data.zip** (zipped directory containing multiqc output data)

<br>

---

## 4. Alignment

<br>

### 4a. Generate reference
The reference will need to be specific to the organism that was sequenced. Bismark operates on a directory holding the target reference genome in fasta format.

```bash
# creating directory to hold reference and moving it into there
mkdir bismark_reference_genome
mv ref-genome.fasta bismark_reference_genome/

bismark_genome_preparation --path_to_aligner /path/to/bowtie2 \
  --bowtie2 \
  --parallel 4 \
  bismark_reference_genome/
```

**Parameter Definitions:**

*	`--path_to_aligner` - specify path to bowtie2 software if it is not in your `$PATH`
*	`--bowtie2` - specify bismark to create bisulfite indexes for use with Bowtie2
*	`--parallel` – specifies how many threads to use (note these will be doubled as it operates on both strands simultaneously)
*  positional argument specifing the directory holding the reference genome (should end in ".fa" or ".fasta", can be gzipped and including ".gz")

**Input data:**

* a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file))

**Output data:**

* the reference genome directory that was provided as input will now hold indexes for the bisulfite-converted reference genome (all `*.bt2` files are indexes, all `*.fa` files are converted versions of the reference genome)
* bismark_reference_genome/Bisulfite_Genome/
  * **CT_converion/**
    * **BS_CT.1.bt2**
    * **BS_CT.2.bt2**
    * **BS_CT.3.bt2**
    * **BS_CT.4.bt2**
    * **BS_CT.rev.1.bt2**
    * **BS_CT.rev.2.bt2**
    * **genome_mfa.CT_conversion.fa**
  * **GA_conversion/**
    * **BS_GA.1.bt2**
    * **BS_GA.2.bt2**
    * **BS_GA.3.bt2**
    * **BS_GA.4.bt2**
    * **BS_GA.rev.1.bt2**
    * **BS_GA.rev.2.bt2**
    * **genome_mfa.GA_conversion.fa**
  * **\*.txt** (captured standard output from the command)



> **NOTE**  
> If using RNA, the `--hisat2` flag is added, which specifies to use the splice-aware aligner [HISAT2](https://github.com/DaehwanKimLab/hisat2#hisat2), and the outputs would include 8 `*ht2` files in separate sub-directories along with each reference-genome conversion.

<br>

### 4b. Align

Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end). 

**Single-end example**  

```bash
bismark --path_to_bowtie2 /path/to/bowtie2 \
  --samtools_path /path/to/samtools \
  --bowtie2 \
  --bam \
  --parallel 4 \
  --non_bs_mm \
  --basename sample-1_bismark_bt2 \
  --output_dir /path/to/alignment/output/files/sample-1 \
  --genome_folder bismark_reference_genome/ \
  sample-1_trimmed.fastq.gz

# renaming output report file so works with sorted bam file/auto-detection of bismark2summary later
mv sample-1_bismark_bt2_SE_report.txt sample-1_bismark_bt2_sorted_SE_report.txt
```

**Paired-end example**  

```bash
bismark --path_to_bowtie2 /path/to/bowtie2 \
  --samtools_path /path/to/samtools \
  --bowtie2 \
  --bam \
  --parallel 4 \
  --non_bs_mm \
  --basename sample-1_bismark_bt2 \
  --output_dir /path/to/alignment/output/files/sample-1 \
  --genome_folder bismark_reference_genome/ \
  -1 sample-1_R1_trimmed.fastq.gz \
  -2 sample-1_R2_trimmed.fastq.gz

# renaming output report file so works with sorted bam file/auto-detection of bismark2summary later
mv sample-1_bismark_bt2_PE_report.txt sample-1_bismark_bt2_PE_report.txt
```

**Parameter Definitions:**

* `--path_to_bowtie2` - specifies path to bowtie2 software if it is not in your `$PATH` 
* `--samtools_path` - specifies path to samtools software if it is not in your `$PATH`
* `--bowtie2` - specifies to use bowtie2 for alignment (limited to end-to-end alignments)
* `--bam` - specifies to convert the default output sam format into compressed bam format
* `--parallel` - allows us to specify the number of threads to use (note: will consume 3-5X this value)
* `--non_bs_mm` - outputs an extra column specifying the number of non-bisulfite mismatches each read has
* `--basename` - base file name to use for all output files
* `--genome_folder` - specifies the directory holding the reference genome indexes (the same that was provided in [Step 4a.](#4a-generate-reference) above)
* input trimmed-reads are provided as a positional argument if they are single-end data
* `-1` - where to specify the forward trimmed reads if paired-end data
* `-2` - where to specify the reverse trimmed reads if paired-end data


**Input data:**
* bismark_reference_genome/ (directory holding indexes of reference genome)
* gzip-compressed fastq files (adapter-trimmed/quality-filtered reads)

**Output data:**  

* sample-1_bismark_bt2.bam (mapping file) 
* **\*report.txt** (bismark mapping report file)


> **NOTE**  
> If using RNA, need to add the `--hisat2` and `--path_to_hisat2` flags.

<br>

### 4c. Sort Alignment Files

```bash
samtools sort -@ 4 \
  -o sample-1_bismark_bt2.sorted.bam \
  sample-1_bismark_bt2.bam
```

**Parameter Definitions:**

* `sort` - specifies to use the `sort` sub-program of `samtools`
* `-@` - specifies the number of threads to use
* `-o` - specifies the output file name
* sample-1_bismark_bt2.bam - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2.bam (bismark bowtie2 alignment bam file, output from [Step 4b.](#4b-align) above)

**Output data:**

* **sample-1_bismark_bt2.sorted.bam** (sorted bam file)

<br>

---

## 5. Alignment QC

```bash
qualimap bamqc -bam sample-1_bismark_bt2.sorted.bam \
  -outdir /path/to/alignment/qc/files/sample-1_qualimap \
  --collect-overlap-pairs \
  --java-mem-size=6G \
  -nt 4
```

**Parameter Definitions:**

* `bamqc` - specifies the `bamqc` sub-program of `qualimap`
* `-bam` - specifies the input bam file
* `-outdir` - specifies the directory to print the output files to
* `--collect-overlap-pairs` - includes statistics of overlapping paired-end reads (if data were paired-end, no effect if single-end)
* `--java-mem-size=6G` - specifies the amount of memory to use (here 6G; see [qualimap FAQ here](http://qualimap.conesalab.org/doc_html/faq.html?highlight=java-mem-size))
* `-nt` - specifies the number of threads to use

**Input data:**

* sample-1_bismark_bt2.sorted.bam (sortded bam file, output from [Step 4c.](#4c-sort-alignment-files) above)

**Output data:** 

* **sample-1_bismark_bt2.sorted_qualimap/** - subdirectory of many alignment QC output files and formatting files for presenting in an html file (see [qualimap documentation](http://qualimap.conesalab.org/doc_html/analysis.html#output))


<br>

---

## 6. Deduplicate (skip if data are RRBS)
> **NOTE**  
> This step should **not** be done if the data are RRBS (reduced representation bisulfite sequencing; see e.g., [bismark documentation](https://github.com/FelixKrueger/Bismark/tree/master/Docs#iii-running-deduplicate_bismark)).

```bash
deduplicate_bismark sample-1_bismark_bt2.sorted.bam
```

**Parameter Definitions:**

* positional argument is the sorted bam file produced in step 4 above

**Input data:**

* sample-1_bismark_bt2.sorted.bam - sorted mapping file produced in step 4 above

**Output data:**

* **\*.deduplicated.bam** - a deduplicated mapping file
* **\*.deduplication_report.txt** - report file of deduplication 


<br>

---

## 7. Extract methylation calls


**Single-end example**  

```bash
bismark_methylation_extractor --bedGraph --gzip --comprehensive sample-1_trimmed_sorted_bismark_bt2.bam
    # note, input should be the deduplicated version produced in step 6 above if *not working with RRBS data
```

**Paired-end example**  

```bash
bismark_methylation_extractor --bedGraph --gzip --comprehensive --ignore_r2 2 --ignore_3prime_r2 2 sample-1_trimmed_sorted_bismark_bt2.bam
    # note, input should be the deduplicated version produced in step 6 above if *not working with RRBS data
```


**Parameter Definitions:**

* `--bedGraph` - specifies to generate a bedGraph-formatted file of methylated CpGs (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output))
* `--gzip` - specifies to gzip-compress the larger output files
* `--comprehensive` - specifies to merge all four possible strand-specific methylation outputs into context-dependent output files
* `--ignore_r2` - allows specifying how many bases to ignore from the 5' end of the reverse reads (bismark docs recommend 2, see [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#appendix-iii-bismark-methylation-extractor))
* `--ignore_3prime_r2` - allows specifying how many bases to ignore from the 3' end of the reverse reads (this is utilized in the [nf-core methylseq workflow](https://nf-co.re/methylseq), set at [this line](https://github.com/nf-core/methylseq/blob/03972a686bedeb2920803cd575f4d671e9135af0/main.nf#L643)) 
* the positional argument is an input bam file

**Input data:**

* sample-1_trimmed_sorted_bismark_bt2.bam - sorted bam file produced above (in step 4 if data are RRBS, or step 6 if not and the bam file was deduplicated)

> **NOTE**  
> If data are **not** RRBS, the input bam file should be the deduplicated one produced by step 6 above. 


**Output data:**

* **\*\_context\_\*.txt.gz** - bismark methylation-call files for CpG and non-CpG contexts that were detected; see [bismark documentation](https://github.com/FelixKrueger/Bismark/tree/master/Docs), namely [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#methylation-call) for symbols, and [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#iv-bismark-methylation-extractor) for file format
* **\*.bedGraph.gz** - gzip-compressed bedGraph-formatted file of methylation percentages of each CpG site (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output))
* **\*.bismark.cov.gz** - gzip-compressed bedGraph-formatted file like above "\*.bedGraph.gz", but also including 2 more columns of methylated and unmethylated counts at the specified position (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output))
* **\*.M-bias.txt** - text file with methylation information in the context of the position in reads, helpful for investigating bias as a function of base position in the read (see bismark documentation [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#m-bias-plot))
* **\*report.txt** - text file of general detected methylation information

<br>

---

## 8. Generate individual sample report


```bash
bismark2report --alignment_report sample-1_trimmed_sorted_bismark_bt2_SE_report.txt \
               --splitting_report sample-1_trimmed_sorted_bismark_bt2_splitting_report.txt \
               --mbias_report sample-1_trimmed_sorted_bismark_bt2.M-bias.txt
```

**Parameter Definitions:**

* `--alignment_report` - where to provide the alignment report generated by step 4 above 
* `--splitting_report` - where to provide the splitting report generated by step 7 above 
* `--mbias_report` - where to provide the bias report generated by step 7 above 

**Input data:**

* sample-1_trimmed_sorted_bismark_bt2_SE_report.txt - alignment report generated by step 4 above
* sample-1_trimmed_sorted_bismark_bt2_splitting_report.txt - splitting report generated by step 7 above
* sample-1_trimmed_sorted_bismark_bt2.M-bias.txt - bias report generated by step 7 above

> **NOTE**  
> If data are **not** RRBS, the deduplication report from step 6 above should also be provided to the above command, e.g.: `--dedup_report sample-1_trimmed_sorted_bismark_bt2.deduplication_report.txt` 

**Output data:**

* **\*_report.html** - a bismark summary html file for the given sample

<br>

---

## 9. Generate combined summary report

```bash

bismark2summary sample-1_bismark_bt2.sorted.bam
```

**Input data:**  

* autodetects appropriate files in current working directory intially based on specified bam files generated in step 4 above (here `sample-1_bismark_bt2.sorted.bam` given as a positional argument; other files can't be explicitly provided, but it looks for these:
  * sample-1_trimmed_sorted_bismark_bt2_SE_report.txt generated from step 4 above
  * sample-1_trimmed_sorted_bismark_bt2_splitting_report.txt generated from step 7 above
  * and the deduplication report files if deduplication was done (these are not required)

**Output data:**  

* **bismark_summary_report.txt** - summary table of general information on all samples
* **bismark_summary_report.html** - html summary of general information on all samples

<br>

---

## 10. Generate reference genome annotation information

### 10a. GTF to BED conversion
A bed-formatted annotation file is needed for adding annotation information to the output from the differential methylation analysis. We utilize gtf files from [Ensembl](https://www.ensembl.org/) and convert them as in the following example:

```bash
# downloading mouse reference gtf for this example
curl -LO https://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz

gunzip Mus_musculus.GRCm39.107.gtf.gz

gtfToGenePred Mus_musculus.GRCm39.107.gtf Mus_musculus.GRCm39.107.genePred

genePredToBed Mus_musculus.GRCm39.107.genePred Mus_musculus.GRCm39.107.bed

# removing intermediate files
rm Mus_musculus.GRCm39.107.genePred Mus_musculus.GRCm39.107.gtf
```

**Input data:**

* a reference gtf file ("Mus_musculus.GRCm39.107.gtf.gz" in the above example)

**Output data:**

* **the generated bed file** ("Mus_musculus.GRCm39.107.bed" in the above example)

### 10b. Making a mapping file of genes to transcripts
Making a mapping file of gene names to transcript names, which we need to link functional annotations in a primary output table. 

```bash
awk ' $3 == "transcript" ' Mus_musculus.GRCm39.107.gtf | cut -f 9 | tr -s ";" "\t" | \
    cut -f 1,3 | tr -s " " "\t" | cut -f 2,4 | tr -d '"' \
    > Mus_musculus.GRCm39.107-gene-to-transcript-map.tsv
```

**Input data:**

* a reference gtf file ("Mus_musculus.GRCm39.107.gtf" in the above example)

**Output data:**

* **the generated gene-to-transcript mapping file** with gene IDs in the first column and transcript IDs in the second ("Mus_musculus.GRCm39.107-gene-to-transcript-map.tsv" in the above example)

<br>

---

## 11. Differential methylation analysis

Example data for the R code below can be downloaded and unpacked with the following:

```bash
curl -L -o subset-test-results.tar https://figshare.com/ndownloader/files/38004075
tar -xvf subset-test-results.tar && rm subset-test-results.tar
```

The remainder of this section is performed in R. 

```R
library(tidyverse)
library(methylKit)

## reading in data
file.list <-list("subset-test-results/F-SRR12865062-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/F-SRR12865063-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/F-SRR12865064-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865070-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865071-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865072-sub_trimmed_bismark_bt2.bismark.cov.gz")

sample.list <- list("F-SRR12865062",
                    "F-SRR12865063",
                    "F-SRR12865064",
                    "G-SRR12865070",
                    "G-SRR12865071",
                    "G-SRR12865072")

# reading into memory
myobj <- methRead(location = file.list,
                  sample.id = sample.list,
                  assembly = "Mmus_GRCm39",
                  pipeline = "bismarkCoverage",
                  header = FALSE,
                  treatment = c(1,1,1,0,0,0),
                  mincov = 10)

# example of how to store as tables if memory requirements are too high
# myobj_storage <- methRead(location = file.list,
#                   sample.id = sample.list,
#                   assembly = "Mmus_GRCm39",
#                   dbtype = "tabix",
#                   pipeline = "bismarkCoverage",
#                   header = FALSE,
#                   treatment = c(1,1,1,0,0,0),
#                   dbdir = "methylkit-dbs/",
#                   mincov = 10)

### Individual-base analysis ###
# merging samples
meth <- unite(myobj)

# calculating differential methylation
myDiff <- calculateDiffMeth(meth, mc.cores = 4)

# getting hyper methylated bases
myDiff25p.hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
# making table for writing out
sig_bases_hyper_out_tab <- getData(myDiff25p.hyper) %>% arrange(qvalue)

# getting hypo methylated bases
myDiff25p.hypo <- getMethylDiff(myDiff, difference=25, qvalue = 0.01, type = "hypo")
# making table for writing out
sig_bases_hypo_out_tab <- getData(myDiff25p.hypo) %>% arrange(qvalue)

# getting all differentially methylated bases
myDiff25p <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)
# making table for writing out
sig_bases_all_out_tab <- getData(myDiff25p) %>% arrange(qvalue)

### Tile analysis ###
# tiling
tiles_obj <- tileMethylCounts(myobj, win.size = 1000, step.size = 1000, cov.bases = 10)

# merging tiled samples
tiles_meth <- unite(tiles_obj)

# calculating differential methylation on tiles
tiles_diff <- calculateDiffMeth(tiles_meth, mc.cores = 4)

# getting hyper methlated tiles
# note that with the subset example data, there are no significant tiles
# 25 is the difference cutoff we use, but if wanting to follow this code in full
# with the example data, change that to `difference = 15` in the following
# 3 getMethylDiff() calls, otherwise later functions on the tile objects won't work as those tables will be empty
tiles_myDiff25p.hyper <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01, type = "hyper")
# making table for writing out
tiles_sig_hyper_out_tab <- getData(tiles_myDiff25p.hyper) %>% arrange(qvalue)

# getting hypo methylated tiles
tiles_myDiff25p.hypo <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01, type = "hypo")
# making table for writing out
tiles_sig_hypo_out_tab <- getData(tiles_myDiff25p.hypo) %>% arrange(qvalue)

# getting all differentially methylated tiles
tiles_myDiff25p <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01)
# making table for writing out
sig_tiles_all_out_tab <- getData(tiles_myDiff25p) %>% arrange(qvalue)


### Adding feature information ###
library(genomation)

gene.obj <- readTranscriptFeatures("subset-test-results/Mus_musculus.GRCm39.107.bed", up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, 
                                   unique.prom = TRUE)

## adding features to individual-base objects
diffAnn <- annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
diffAnn.hyper <- annotateWithGeneParts(as(myDiff25p.hyper, "GRanges"), gene.obj)
diffAnn.hypo <- annotateWithGeneParts(as(myDiff25p.hypo, "GRanges"), gene.obj)

# making base-level sig table with features 
sig_all_bases_tab_with_features <- cbind(data.frame(myDiff25p), 
                                         getAssociationWithTSS(diffAnn), 
                                         as.data.frame(getMembers(diffAnn))) %>% .[,-c(8)]

# making base-level sig hyper-methylated table with features
sig_bases_hyper_tab_with_features <- cbind(data.frame(myDiff25p.hyper), 
                                           getAssociationWithTSS(diffAnn.hyper), 
                                           as.data.frame(getMembers(diffAnn.hyper))) %>% .[,-c(8)]

# making base-level sig hypo-methylated table with features
sig_bases_hypo_tab_with_features <- cbind(data.frame(myDiff25p.hypo), 
                                          getAssociationWithTSS(diffAnn.hypo), 
                                          as.data.frame(getMembers(diffAnn.hypo))) %>% .[,-c(8)]


## adding features to tiles objects
tiles_diffAnn <- annotateWithGeneParts(as(tiles_myDiff25p, "GRanges"), gene.obj)
tiles_diffAnn.hyper <- annotateWithGeneParts(as(tiles_myDiff25p.hyper, "GRanges"), gene.obj)
tiles_diffAnn.hypo <- annotateWithGeneParts(as(tiles_myDiff25p.hypo, "GRanges"), gene.obj)


# making tiles sig table with features 
tiles_sig_all_out_tab_with_features <- cbind(data.frame(tiles_myDiff25p), 
                                             getAssociationWithTSS(tiles_diffAnn), 
                                             as.data.frame(getMembers(tiles_diffAnn))) %>% .[,-c(8)]

# making tiles sig hyper-methylated table with features
tiles_sig_hyper_tab_with_features <- cbind(data.frame(tiles_myDiff25p.hyper), 
                                           getAssociationWithTSS(tiles_diffAnn.hyper), 
                                           as.data.frame(getMembers(tiles_diffAnn.hyper))) %>% .[,-c(8)]

# making tiles sig hypo-methylated table with features
tiles_sig_hypo_tab_with_features <- cbind(data.frame(tiles_myDiff25p.hypo), 
                                          getAssociationWithTSS(tiles_diffAnn.hypo), 
                                          as.data.frame(getMembers(tiles_diffAnn.hypo))) %>% .[,-c(8)]


### Adding functional annotations ###
# reading in annotation table appropriate for current organism
functional_annots_tab <- 
    read.table("https://figshare.com/ndownloader/files/35939642", sep = "\t", 
               quote = "", header = TRUE)

# reading in gene to transcript mapping file
gene_transcript_map <- 
    read.table("subset-test-results/Mus_musculus.GRCm39.107-gene-to-transcript-map.tsv", sep = "\t", 
               col.names = c("gene_ID", "feature.name"))


## for individual-base output
# for each transcript ID in the sig_all_bases_tab_with_features table, getting 
# its corresponding gene ID and adding that to the table
sig_all_bases_tab_with_features_and_gene_IDs <- 
    left_join(sig_all_bases_tab_with_features, gene_transcript_map)

sig_hyper_bases_tab_with_features_and_gene_IDs <-
    left_join(sig_bases_hyper_tab_with_features, gene_transcript_map)

sig_hypo_bases_tab_with_features_and_gene_IDs <-
    left_join(sig_bases_hypo_tab_with_features, gene_transcript_map)

# now adding full annotations
sig_all_bases_tab_with_features_and_annots <- 
    left_join(sig_all_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hyper_bases_tab_with_features_and_annots <- 
    left_join(sig_hyper_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hypo_bases_tab_with_features_and_annots <- 
    left_join(sig_hypo_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

# and writing out
write.table(sig_all_bases_tab_with_features_and_annots, 
            "sig-diff-methylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hyper_bases_tab_with_features_and_annots, 
            "sig-diff-hypermethylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hypo_bases_tab_with_features_and_annots, 
            "sig-diff-hypomethylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)


# making and writing out a table of base-level percent methylated
perc.meth <- percMethylation(meth, rowids = TRUE)
perc.meth <- perc.meth %>% data.frame() %>% rownames_to_column("location")

write.table(perc.meth, "base-level-percent-methylated.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)


## for tiles output
# for each transcript ID in the tiles_sig_all_out_tab_with_features table, getting 
# its corresponding gene ID and adding that to the table
sig_all_tiles_tab_with_features_and_gene_IDs <- 
    left_join(tiles_sig_all_out_tab_with_features, gene_transcript_map)

sig_hyper_tiles_tab_with_features_and_gene_IDs <-
    left_join(tiles_sig_hyper_tab_with_features, gene_transcript_map)

sig_hypo_tiles_tab_with_features_and_gene_IDs <-
    left_join(tiles_sig_hypo_tab_with_features, gene_transcript_map)


# now adding full annotations
sig_all_tiles_tab_with_features_and_annots <- 
    left_join(sig_all_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hyper_tiles_tab_with_features_and_annots <- 
    left_join(sig_hyper_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hypo_tiles_tab_with_features_and_annots <- 
    left_join(sig_hypo_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

# and writing out
write.table(sig_all_tiles_tab_with_features_and_annots, 
            "sig-diff-methylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hyper_tiles_tab_with_features_and_annots, 
            "sig-diff-hypermethylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hypo_tiles_tab_with_features_and_annots, 
            "sig-diff-hypomethylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)


# making and writing out a table of tile-level percent methylated
tiles_perc.meth <- percMethylation(tiles_meth, rowids = TRUE)
tiles_perc.meth <- tiles_perc.meth %>% data.frame() %>% rownames_to_column("location")

write.table(tiles_perc.meth, "tile-level-percent-methylated.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)


### Overview figure of percent diff. methylation across features ###
## based on individual bases
pdf("sig-diff-methylated-bases-across-features.pdf")
plotTargetAnnotation(diffAnn, precedence = TRUE, 
                     main = "% of sig. diff. methylated sites across features")
dev.off()

## based on tiles
pdf("sig-diff-methylated-tiles-across-features.pdf")
plotTargetAnnotation(tiles_diffAnn, precedence = TRUE, 
                     main = "% of sig. diff. methylated tiles across features")
dev.off()
```

**Input data:**
* \*.bismark.cov.gz - gzip-compressed bedGraph-formatted files generated in Step 7 above
* Mus_musculus.GRCm38.101.bed - bed file generated in Step 11a above
* Mus_musculus.GRCm38.101-gene-to-transcript-map.tsv - gene-to-transcript mapping file generated in Step 11b above

**Output data:**
* **\*sig-diff-methylated-bases.tsv** - all significantly differentially methylated cytosines in contrast listed in filename (e.g., will be something like "A_vs_B-sig-diff-methylated-bases.tsv")
* **\*sig-diff-hypermethylated-bases.tsv** - cytosines with significantly elevated methylation levels in contrast listed in filename
* **\*sig-diff-hypomethylated-bases.tsv** - cytosines with significantly reduced methylation levels in contrast listed in filename
* **\*sig-diff-methylated-tiles.tsv** - all significantly differentially methylated tiles in contrast listed in filename
* **\*sig-diff-hypermethylated-tiles.tsv** - tiles with significantly elevated methylation levels in contrast listed in filename
* **\*sig-diff-hypomethylated-tiles.tsv** - tiles with significantly reduced methylation levels in contrast listed in filename
* **base-level-percent-methylated.tsv** - table of methylation levels across all cytosines and samples
* **tile-level-percent-methylated.tsv** - table of methylation levels across all tiles and samples
* **\*sig-diff-methylated-bases-across-features.pdf** - overview figure of what percent of identified significantly differentially methylated cytosines are in specific features (promoter, exon, intron)
* **\*sig-diff-methylated-tiles-across-features.pdf** - overview figure of what percent of identified significantly differentially methylated tiles are in specific features (promoter, exon, intron)


> NOTE:   
> All of these files, except "\*percent-methylated.tsv", will be prefixed with contrasted groups, e.g., A_vs_B-\*.

