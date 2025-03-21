# CAR-NGS Frontend

## "CAR-NGS: Buckle Up, We're Taking Your DNA for a Ride!"

CAR-NGS is a toolkit designed to simplify high-throughput genomic data analysis by providing user-friendly R-based frontend functions. These functions serve as an interface to automated data analysis pipelines, enabling seamless execution of genomic workflows without the need for extensive command-line expertise.

## Installation

You can install the tool directly from GitHub using the following commands in R:

```r
install.packages("devtools")
library(devtools)
install_github("https://github.com/Reproducible-Bioinformatics/CAR-NGS", ref="main")
```

Once installed, you can access and use the provided functions to interact with the CAR-NGS pipelines efficiently.

## 6S rRNA Gene Analysis (`sixteenS` function)

The `sixteenS` function is an R-based frontend for executing the 16S rRNA gene sequencing analysis pipeline for microbial community profiling. This function simplifies execution by managing Docker commands directly from R, ensuring reproducibility and ease of use.

### Pipeline Overview

The 16S rRNA analysis pipeline consists of the following steps:

1. **Data Import**: Converts raw paired-end FASTQ files into a QIIME2 artifact format (.qza).
2. **Quality Control**: Summarizes sequencing read quality to assess data integrity.
3. **Denoising (DADA2)**: Filters and corrects sequencing errors, inferring true biological sequences (ASVs).
4. **Taxonomic Classification**: Assigns taxonomy to inferred sequences using the Silva 138 pre-trained classifier.
5. **Visualization**: Generates interactive taxonomic bar plots to summarize the microbial community composition.
6. **Exporting Results**: Converts visualization files into HTML reports for easy access.

### Usage

```r
sixteenS(input_dir_path = "/path/to/fastq_files")
```

#### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to the directory containing paired-end FASTQ files. The file names must follow the format: `sample_R1.fastq.gz` and `sample_R2.fastq.gz`. |

## Pipeline Steps in Detail

### Step 1: Data Import into QIIME2 Format

Converts paired-end FASTQ files into QIIME2 format using a manifest file. The manifest file contains absolute paths to both forward (R1) and reverse (R2) reads. The imported data is stored as a `.qza` file.

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path aligned_results/sample-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

### Step 2: Quality Control Summary

Generates a summary report of sequence quality to identify trimming parameters. The output is a QIIME2 visualization file (`.qzv`), which can be viewed interactively.

```bash
qiime demux summarize \
  --i-data aligned_results/sample-paired-end-demux.qza \
  --o-visualization aligned_results/sample-paired-end-demux.qzv
```

### Step 3: Denoising with DADA2

Removes sequencing errors and infers amplicon sequence variants (ASVs). Trims reads to 250 bp forward and reverse.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs aligned_results/sample-paired-end-demux.qza \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 250 --p-trunc-len-r 250 \
  --o-table aligned_results/sample-table.qza \
  --o-representative-sequences aligned_results/sample-rep-seqs.qza \
  --o-denoising-stats aligned_results/sample-denoising-stats.qza
```

### Step 4: Denoising Statistics Summary

Converts the denoising statistics into a human-readable format for evaluation.

```bash
qiime metadata tabulate \
  --m-input-file aligned_results/sample-denoising-stats.qza \
  --o-visualization aligned_results/sample-denoising-stats.qzv
```

### Step 5: Taxonomic Classification

Uses the Silva 138 classifier to assign taxonomy to ASVs.

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /home/silva-138-99-nb-classifier.qza \
  --i-reads aligned_results/sample-rep-seqs.qza \
  --o-classification aligned_results/sample-taxonomy.qza
```

### Step 6: Taxonomic Visualization

Generates interactive taxonomic bar plots, showing microbial composition at different taxonomic levels.

```bash
qiime taxa barplot \
  --i-table aligned_results/sample-table.qza \
  --i-taxonomy aligned_results/sample-taxonomy.qza \
  --o-visualization aligned_results/sample-taxa-bar-plots.qzv
```

### Step 7: Exporting Results to HTML

Converts the taxonomic bar plots into an HTML report for easy visualization.

```bash
qiime tools export \
  --input-path aligned_results/sample-taxa-bar-plots.qzv \
  --output-path aligned_results/html/sample-taxa-bar-plots
```

## Output Files

| File | Description |
|------|-------------|
| `aligned_results/sample-paired-end-demux.qza` | Imported QIIME2 artifact containing demultiplexed sequences. |
| `aligned_results/sample-paired-end-demux.qzv` | Visualization of demultiplexed sequences (interactive). |
| `aligned_results/sample-table.qza` | Feature table (ASV abundance per sample). |
| `aligned_results/sample-rep-seqs.qza` | Representative ASV sequences. |
| `aligned_results/sample-denoising-stats.qza` | Denoising statistics from DADA2. |
| `aligned_results/sample-taxonomy.qza` | Taxonomic classifications from the Silva 138 classifier. |
| `aligned_results/sample-taxa-bar-plots.qzv` | Taxonomic bar plots (interactive visualization). |
| `aligned_results/html/sample-taxa-bar-plots/` | Exported HTML report of taxonomic composition. |

## Prerequisites

- **Docker** must be installed and running on your system.
- The input directory must contain paired-end FASTQ files named as:
  - `sample_R1.fastq.gz` (Forward reads)
  - `sample_R2.fastq.gz` (Reverse reads)
- The **Silva 138 classifier file** (`silva-138-99-nb-classifier.qza`) must be available in the Docker container.

## Example

```r
sixteenS(input_dir_path = "/the/input/path")
```

This command runs the 16S rRNA analysis pipeline inside the `repbioinfo/qiime2023:latest` Docker container, mapping the input directory to its appropriate location inside the container.

## Additional Notes for Materials & Methods

- **QIIME2 (version 2023)** was used for all analyses.
- **DADA2** denoising was performed with a truncation length of **250 bp**.
- The **Silva 138 classifier** was used for taxonomy assignment.
- The taxonomic composition was visualized using **QIIME2 bar plots**.
- Results were exported as **HTML reports** for easy visualization.


## ATAC-seq Analysis (`atacSeq` function)

The `atacSeq` function is an R-based frontend for executing the ATAC-seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) analysis pipeline. This function streamlines the execution of the backend pipeline by managing Docker execution directly from R, ensuring reproducibility and ease of use.

### Pipeline Overview

The ATAC-seq analysis pipeline performs the following steps:

1. **Quality Control**: Runs FastQC on raw FASTQ files to assess sequencing quality.
2. **Genome Indexing**: If an index is not found, it generates a Bowtie2 index from the reference genome.
3. **Read Alignment**: Aligns sequencing reads to the reference genome using Bowtie2, allowing soft clipping.
4. **BAM Processing**:
   - Sorts and indexes BAM files using Samtools.
   - Removes mitochondrial (Mt) and plastid (Pt) reads to filter out non-nuclear DNA.
5. **Peak Calling**: Identifies regions of open chromatin using MACS2 (either paired-end or single-end).
6. **BigWig Conversion (Optional)**: Converts BedGraph peak intensity files to BigWig format for visualization in genome browsers.

### Usage

```r
atacSeq(input_dir_path = "/path/to/fastq_files",
        genome_dir_path = "/path/to/fasta_files",
        nThreads = 8)
```

#### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to the directory containing the FASTQ files to be analyzed. |
| `genome_dir_path` | character | Path to the directory containing the FASTA reference genome. The directory must contain a `.fa` or `.fasta` file. |
| `nThreads` | integer (default = 8) | Number of CPU cores for parallelization. |

## Pipeline Steps in Detail

### Step 1: Quality Control (FastQC)

For each input FASTQ file, FastQC generates quality control reports, which are saved in the `results/Quality_ATAC/` directory.

### Step 2: Genome Indexing (Bowtie2 Index)

If an indexed genome is not found in `/genomes/index/`, the pipeline creates one using Bowtie2 (`bowtie2-build`).
Indexing is skipped if an existing Bowtie2 index (`genome_index.1.bt2`) is found.

### Step 3: Read Alignment (Bowtie2)

#### Paired-end data:

```bash
bowtie2 --local --threads 8 -x /genomes/index/genome_index -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz | samtools view -bS - | samtools sort -o sample.sorted.bam
```

#### Single-end data:

```bash
bowtie2 --local --threads 8 -x /genomes/index/genome_index -U sample.fastq.gz | samtools view -bS - | samtools sort -o sample.sorted.bam
```

The resulting BAM file (`sample.sorted.bam`) is indexed:

```bash
samtools index sample.sorted.bam
```

### Step 4: Mitochondrial and Plastid Read Removal

Reads mapped to mitochondrial (Mt) or plastid (Pt) genomes are removed:

```bash
samtools idxstats sample.sorted.bam | cut -f1 | grep -v Mt | grep -v Pt | xargs samtools view -b sample.sorted.bam > sample.sorted.noorg.bam
```

The filtered nuclear BAM file (`sample.sorted.noorg.bam`) is indexed:

```bash
samtools index sample.sorted.noorg.bam
```

### Step 5: Peak Calling (MACS2)

#### Paired-end data:

```bash
macs2 callpeak -t sample.sorted.noorg.bam -q 0.05 --broad -f BAMPE -n sample -B --trackline --outdir results/peaks/
```

#### Single-end data:

```bash
macs2 callpeak -t sample.sorted.noorg.bam -q 0.05 --broad -f BAM -n sample -B --trackline --outdir results/peaks/
```

**Output includes:**  
- NarrowPeak/BroadPeak files (`sample_peaks.broadPeak`)
- BedGraph signal intensity files (`sample_treat_pileup.bdg`)

### Step 6: BigWig Conversion (Optional, for Visualization)

Converts BedGraph to BigWig format for visualization in genome browsers (UCSC, IGV).

```bash
bedGraphToBigWig sample_treat_pileup.bdg genome_sizes.txt sample_treat_pileup.bw
```

The BigWig file (`sample_treat_pileup.bw`) is saved in `results/peaks/`.

## Output Files

| File | Description |
|------|-------------|
| `results/Quality_ATAC/*.html` | FastQC reports for quality control. |
| `results/*.sorted.bam` | Sorted BAM file (aligned reads). |
| `results/*.sorted.bam.bai` | Index file for BAM files. |
| `results/*.sorted.noorg.bam` | Filtered BAM file (removes Mt & Pt reads). |
| `results/*.sorted.noorg.bam.bai` | Index file for filtered BAM. |
| `results/peaks/*.broadPeak` | Peak calling results from MACS2. |
| `results/peaks/*.treat_pileup.bdg` | BedGraph file for signal intensity. |
| `results/peaks/*.treat_pileup.bw` | BigWig file for genome browser visualization. |

## Prerequisites

- **Docker** must be installed and running on your system.
- The input directory must contain FASTQ files for ATAC-seq.
- The genome directory must contain a FASTA file and, if possible, a pre-built Bowtie2 index.

## Example

```r
atacSeq(input_dir_path = "/the/input/path",
        genome_dir_path = "/the/genome/path",
        nThreads = 12)
```

This command runs the ATAC-seq analysis pipeline inside the `repbioinfo/atacseq:latest` Docker container, mapping the input and genome directories to their respective locations inside the container.

## Additional Notes for Materials & Methods

- The **Bowtie2 aligner** was used in local mode to allow soft clipping.
- **Paired-end and single-end read alignment** are both supported.
- **MACS2 peak calling** was performed using a **q-value cutoff of 0.05** and the `--broad` flag.
- **BigWig conversion** was performed for visualization using `bedGraphToBigWig`.


## Bulk RNA-Seq Analysis Functions

The bulk RNA-Seq analysis pipeline provides a fully automated workflow for processing RNA-Seq data, performing alignment, quantification, statistical analysis, and visualization. These R-based frontend functions manage execution inside a Docker container, ensuring reproducibility.

### Pipeline Overview

The bulk RNA-Seq analysis pipeline consists of the following major steps:

1. **Genome Indexing & Read Alignment**  
   - Uses **STAR** to align raw RNA-Seq reads to the reference genome.  
   - Generates sorted BAM files and a gene count matrix.  

2. **Principal Component Analysis (PCA)**  
   - Conducts unsupervised clustering to visualize sample distribution.  

3. **Differential Expression Analysis (DESeq2)**  
   - Identifies significantly differentially expressed genes (DEGs).  
   - Filters genes based on adjusted p-value and log2 fold-change.  
   - Generates DEG tables and Venn diagrams.  

4. **Heatmap Generation**  
   - Visualizes the expression patterns of significant genes.  

5. **Complete Downstream Analysis (CDSA)**  
   - Runs all post-alignment steps in one command.  

---

## 1. Genome Indexing & Read Alignment (`index_align` function)

The `index_align` function performs alignment using **STAR** and generates a gene count matrix.

### Usage

```r
index_align(input_dir_path = "/path/to/fastq_files",
            genome_dir_path = "/path/to/genome_files")
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to the directory containing FASTQ files. |
| `genome_dir_path` | character | Path to the directory containing the FASTA and GTF files. |

### Pipeline Steps

- **Indexing**: If no **STAR** index is found, it creates one (stored in `/genome/star_index/`).
- **Read Trimming**: Adapters are removed using **Cutadapt**.
- **Read Alignment**: Aligns reads to the genome using **STAR** and produces sorted BAM files.
- **Gene Quantification**: Generates a gene count matrix (`gene_count_matrix.csv`).
- **Metadata Creation**: Stores sample group information in `Covariatesstat.csv`.

### Output Files

| File | Description |
|------|-------------|
| `gene_count_matrix.csv` | Count matrix of gene expression levels. |
| `Covariatesstat.csv` | Metadata with sample group labels. |
| `sorted.bam` | Aligned RNA-Seq reads. |

---

## 2. Principal Component Analysis (`pca` function)

The `pca` function performs PCA to assess sample clustering.

### Usage

```r
pca(input_dir_path = "/path/to/results",
    countmatrix_name = "gene_count_matrix.csv",
    metadata_name = "Covariatesstat.csv")
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to directory containing RNA-Seq results. |
| `countmatrix_name` | character | Name of the count matrix file. |
| `metadata_name` | character | Name of the metadata file. |

### Pipeline Steps

- **PCA Calculation**: Transforms gene expression data into principal components.
- **Visualization**: Generates PCA plots to identify sample clustering.

### Output Files

| File | Description |
|------|-------------|
| `gene_count_matrix_pca_plot.png` | PCA plot of RNA-Seq samples. |

---

## 3. Differential Expression Analysis (`deseq2` function)

The `deseq2` function identifies differentially expressed genes (DEGs).

### Usage

```r
deseq2(input_dir_path = "/path/to/results",
       countmatrix_name = "gene_count_matrix.csv",
       metadata_name = "Covariatesstat.csv",
       reference_group = "wt",
       organism = "Drosophilamelanogaster")
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to directory containing RNA-Seq results. |
| `countmatrix_name` | character | Name of the count matrix file. |
| `metadata_name` | character | Name of the metadata file. |
| `reference_group` | character | Name of the control group (e.g., "wt"). |
| `organism` | character | Organism: "Homo sapiens", "Mus musculus", or "Drosophila melanogaster". |

### Pipeline Steps

- **Differential Expression Analysis**: Uses **DESeq2** to compute DEGs.
- **Filtering**: Applies thresholds of **p-adjusted < 0.01** and **log2 fold-change > 2**.
- **Venn Diagram Creation**: Compares DEGs across conditions.

### Output Files

| File | Description |
|------|-------------|
| `DEG_<group>_vs_<reference_group>.csv` | DEG table for each comparison. |
| `filtered_count_matrix.csv` | Filtered gene expression matrix. |
| `venn_diagram.png` | Venn diagram of significant genes. |

---

## 4. Heatmap Generation (`heatmap` function)

The `heatmap` function visualizes gene expression patterns.

### Usage

```r
heatmap(input_dir_path = "/path/to/results",
        countmatrix_name = "filtered_count_matrix.csv",
        metadata_name = "Covariatesstat.csv")
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `input_dir_path` | character | Path to directory containing RNA-Seq results. |
| `countmatrix_name` | character | Name of the filtered count matrix file. |
| `metadata_name` | character | Name of the metadata file. |

### Pipeline Steps

- **Normalization**: Normalizes expression levels across samples.
- **Heatmap Generation**: Uses **pheatmap** to cluster and visualize significant genes.

### Output Files

| File | Description |
|------|-------------|
| `heatmap_filtered.png` | Heatmap of differentially expressed genes. |

---

## 5. Complete Downstream Analysis (`cdsa` function)

The `cdsa` function automates all post-alignment analyses.

### Usage

```r
cdsa(input_dir_path = "/path/to/results",
     countmatrix_name = "gene_count_matrix.csv",
     metadata_name = "Covariatesstat.csv",
     reference_group = "wt",
     organism = "Drosophilamelanogaster")
```

### Pipeline Steps

- Runs **PCA** to visualize sample clustering.
- Performs **DESeq2** Analysis to detect differentially expressed genes.
- Generates **Heatmap** for significant genes.

### Output Files

All PCA, DEG, and heatmap results from previous steps.

---

## Additional Notes for Materials & Methods

- **STAR (v2.7.10a)** was used for genome indexing and read alignment.
- **DESeq2 (v1.34.0)** performed differential expression analysis.
- **PCA** was based on log-transformed gene expression values.
- **DEGs** were filtered at adjusted **p-value < 0.01** and **log2 fold-change > 2**.
- **Heatmaps** were generated using **pheatmap** with row normalization.


## Detect-seq: Genome-Wide Assessment of CBE Off-Targets

The `detectSeq` function is an R-based frontend for executing Detect-seq, a bioinformatics pipeline designed for the unbiased, genome-wide assessment of off-target effects associated with cytosine base editors (CBEs). The function simplifies execution by managing Docker-based execution directly from R, ensuring reproducibility and ease of use.

### Pipeline Overview

The Detect-seq pipeline consists of the following key steps:

1. **Genome Indexing** (if not already available)  
   - Builds **BWA** index for realignment.  
   - Builds **HISAT3N** index with C-to-T conversion handling.  

2. **Read Processing**  
   - Trims sequencing adapters (**Cutadapt**).  
   - Aligns reads to the reference genome (**HISAT3N**).  
   - Extracts low-quality mapped reads for re-alignment with **BWA MEM**.  

3. **Alignment Merging & Deduplication**  
   - Merges HISAT3N and BWA alignments.  
   - Sorts and removes duplicate reads (**Picard**).  

4. **Mutation Detection**  
   - Converts BAM files to **PMAT** format, extracting **C-to-T** conversions.  

5. **Filtering & Visualization**  
   - Filters detected mutations based on user-defined threshold.  
   - Generates **BED** and **WIG** files for downstream analysis.  

---

## Usage

```r
detectSeq(genome_dir_path = "/path/to/genome_files",
          output_dir_path = "/path/to/output_dir",
          fastq_dir_path = "/path/to/fastq_files",
          threshold = 3)
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `genome_dir_path` | character | Path to directory containing FASTA genome file and index (if available). |
| `output_dir_path` | character | Path to directory where output files will be saved. |
| `fastq_dir_path` | character | Path to directory containing FASTQ files. |
| `threshold` | integer | Minimum read count required to retain detected mutations. |

---

## Pipeline Steps in Detail

### Step 1: Genome Indexing

If **BWA** and **HISAT3N** indexes do not exist, the script automatically creates them.

```bash
samtools faidx genome.fa  
bwa index genome.fa  
hisat-3n-build --base-change C,T genome.fa genome_index
```

### Step 2: Read Processing

#### Adapter Trimming (Cutadapt)

Trims adapter sequences and low-quality reads.

```bash
cutadapt -j 0 --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 \
         -a AGATCGGAAGAGCACACGT -A AGATCGGAAGAGCGTCGTG \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz
```

#### Primary Alignment with HISAT3N

Maps reads to CBE-modified genome using **HISAT3N**.

```bash
hisat-3n -x genome_index -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz \
         -p 20 --sensitive --base-change C,T --unique-only --repeat-limit 1000 \
         --no-spliced-alignment -X 700 | samtools view -hb > hisat3n.bam
```

#### Extracting Low-Quality Reads for Re-alignment

Selects reads with mapping quality ≤ 20 for re-alignment.

```bash
samtools view -h hisat3n.bam | awk '$1~"@" || $5 <= 20' | samtools view -hb > low_quality.bam
```

### Step 3: Re-Alignment with BWA MEM

Re-aligns low-quality reads using **BWA MEM**.

```bash
bwa mem genome.fa unmapped_R1.fastq.gz unmapped_R2.fastq.gz \
         -t 20 -M | samtools view -h -b -q 20 -f 3 -F 256 > bwa_realign.bam
```

### Step 4: Merging Alignments & Deduplication

Merges **HISAT3N** and **BWA** alignments, then removes duplicate reads.

```bash
samtools cat -o merged.bam hisat3n.bam bwa_realign.bam
samtools sort -O BAM -o sorted.bam merged.bam
java -jar picard.jar MarkDuplicates I=sorted.bam O=dedup.bam REMOVE_DUPLICATES=true
samtools index dedup.bam
```

### Step 5: Mutation Detection (PMAT Conversion)

Converts BAM files to **PMAT** format, extracting **C-to-T** conversions.

```bash
python bam2pmat.py -i dedup.bam -r genome.fa -o output.pmat -p 20 --mut_type ALL
```

### Step 6: Mutation Filtering & Visualization

Filters detected mutations based on user-defined threshold and generates BED and WIG files.

```bash
awk -v threshold=3 '$9=="CT" && $13>=threshold' output.pmat > filtered.pmat
awk '{print "variableStep chrom="$1" span=1\n"$2" "$9}' filtered.pmat > filtered.wig
```

---

## Output Files

| File | Description |
|------|-------------|
| `dedup.bam` | Final processed BAM file (duplicates removed). |
| `dedup.bam.bai` | Index file for BAM file. |
| `output.pmat` | Raw mutation matrix with detected C-to-T conversions. |
| `filtered.pmat` | Filtered mutation matrix, retaining mutations above threshold. |
| `filtered.bed` | BED file for detected mutations. |
| `filtered.wig` | WIG file for genome browser visualization. |

---

## Prerequisites

- **Docker** must be installed and running on your system.
- The input directory must contain:
  - FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
  - Reference genome (`.fa` or `.fasta`).

---

## Example

```r
detectSeq(genome_dir_path = "/the/genome/dir",
          output_dir_path = "/the/output/dir",
          fastq_dir_path = "/the/fastq/dir",
          threshold = 3)
```

This command runs Detect-seq inside the **repbioinfo/detectseq:latest** Docker container, mapping the genome, output, and FASTQ directories to their appropriate locations inside the container.

---

## Additional Notes for Materials & Methods

- **HISAT3N (v3.0)** was used for primary alignment with C-to-T correction.
- **BWA MEM (v0.7.17)** was used for low-quality read re-alignment.
- **Picard (v2.27.4)** removed duplicate reads.
- **PMAT** files were generated using Detect-seq's `bam2pmat.py` script.
- Mutations were filtered at a threshold of **≥ 3 reads**.



## HTGTS: High-Throughput Genome-wide Translocation Sequencing Analysis

The `HTGTS` function is an R-based frontend for executing the **HTGTS** bioinformatics pipeline, designed to identify and analyze genome-wide translocation events. The function automates the pre-processing, alignment, and translocation detection steps using Docker-based execution, ensuring reproducibility and ease of use.

### Pipeline Overview

The HTGTS pipeline consists of the following key steps:

1. **Library Information Generation**  
   - Extracts sequencing library details from an XML file.  
   - Generates library configuration files required for downstream analysis.  

2. **Read Processing & Demultiplexing**  
   - Identifies specific sequences in FASTQ files.  
   - Separates reads based on predefined primer sequences.  

3. **Alignment & Translocation Mapping**  
   - Aligns reads to the reference genome using **TLPpipeline**.  
   - Identifies potential genome translocations.  

4. **BED File Generation**  
   - Generates **BED** files for downstream analysis and visualization.  

---

## Usage

```r
HTGTS(xml_file_path = "/path/to/HTGTS.xml",
      configType = "HTGTS_mouse",
      data_folder = "/path/to/data",
      fastq1 = "sample_R1.fastq.gz",
      fastq2 = "sample_R2.fastq.gz",
      expInfo = "libseqInfo.txt",
      expInfo2 = "libseqInfo2.txt",
      outDir = "/path/to/output",
      assembly = "mm9")
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `xml_file_path` | character | Path to the XML file containing sequencing library information. |
| `configType` | character | Configuration type (e.g., "HTGTS_mouse", "HTGTS_human"). |
| `data_folder` | character | Path to the directory containing input FASTQ files. |
| `fastq1` | character | Name of the first FASTQ file (R1). |
| `fastq2` | character (optional) | Name of the second FASTQ file (R2, if paired-end). |
| `expInfo` | character | Name of the library information file (`libseqInfo.txt`). |
| `expInfo2` | character | Name of the additional library information file (`libseqInfo2.txt`). |
| `outDir` | character | Directory where output files will be stored. |
| `assembly` | character | Reference genome assembly ("mm9", "mm10", "hg38", "custom"). |

---

## Pipeline Steps in Detail

### Step 1: Library Information Generation

Extracts sequencing library metadata from an XML (Excel) file and generates `libseqInfo.txt` and `libseqInfo2.txt`, which store:

- Restriction enzyme site  
- Primer sequences  
- Breaksite genomic coordinates  

```bash
docker run -v "/data:/Data" repbioinfo/htgts_pipeline_lts_v16 \
       python3 /Algorithm/sample_sheetTolibInfo.py \
       /Data/HTGTS.xlsx /Data/libseqInfo.txt /Data/libseqInfo2.txt HTGTS_mouse
```

### Step 2: Read Processing & Demultiplexing

Identifies specific sequences in FASTQ files and separates reads into correctly formatted sequencing groups.

```bash
GEAT0.2.jar demultiplex -i sample_R1.fastq.gz -o demux_R1.fastq.gz
```

### Step 3: Read Alignment & Translocation Detection

Aligns reads using **TLPpipeline** and identifies potential translocation events in the genome.

```bash
TLPpipeline.pl -i demux_R1.fastq.gz -o alignment_output/ -genome mm9
```

### Step 4: BED File Generation

Creates BED files for visualization of detected translocations.

```bash
bedtools makewindows -b alignment_output/*.bam -w 1000 > translocations.bed
```

---

## Output Files

| File | Description |
|------|-------------|
| `libseqInfo.txt` | Library metadata (restriction sites, primers, genomic coordinates). |
| `libseqInfo2.txt` | Additional library metadata. |
| `alignment_output/*.bam` | Aligned sequencing reads. |
| `translocations.bed` | BED file containing translocation sites. |

---

## Prerequisites

- **Docker** must be installed and running.
- The input directory must contain:
  - **FASTQ files** (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
  - **An XML metadata file** (`HTGTS.xlsx`).

---

## Example

```r
HTGTS(xml_file_path = "/the/path/HTGTS.xlsx",
      configType = "HTGTS_mouse",
      data_folder = "/the/input/data",
      fastq1 = "sample_R1.fastq.gz",
      fastq2 = "sample_R2.fastq.gz",
      expInfo = "libseqInfo.txt",
      expInfo2 = "libseqInfo2.txt",
      outDir = "/the/output/dir",
      assembly = "mm9")
```

This command runs HTGTS inside the **repbioinfo/htgts_pipeline_lts_v16** Docker container, mapping the input/output directories and executing the full pipeline.

---

## Additional Notes for Materials & Methods

- **GEAT0.2** was used for FASTQ demultiplexing.  
- **TLPpipeline (v1.6)** performed read alignment and translocation detection.  
- **BEDTools (v2.30.0)** was used for translocation visualization.  



## SCI: Single-Cell Indexing RNA-seq Pipeline

The `sci_fromfastq` function is an R-based frontend for executing the **SCI-RNA-seq** bioinformatics pipeline, which converts raw FASTQ files into a gene expression matrix, including **UMI counting, quality control, and visualization**. The function automates data processing using Docker-based execution, ensuring reproducibility and ease of use.

### Pipeline Overview

The SCI-RNA-seq pipeline consists of the following key steps:

1. **FASTQ Generation & Preprocessing**  
   - Extracts FASTQ files from raw sequencing data.  
   - Adjusts read names to incorporate barcode and UMI information.  

2. **Poly-A Trimming**  
   - Removes poly-A tails from reads to improve alignment quality.  

3. **Read Alignment**  
   - Aligns reads to a reference genome using **STAR**.  
   - Outputs BAM files and generates alignment statistics.  

4. **Filtering & Sorting BAM Files**  
   - Removes rRNA reads and ambiguously mapped reads.  
   - Assigns reads to genes using BED files.  

5. **UMI Counting & Deduplication**  
   - Computes UMI counts per cell.  
   - Filters cells below the UMI threshold.  

6. **Quality Control & Visualization**  
   - Generates knee plots for UMI distributions.  
   - Produces a final UMI count matrix for downstream analysis.  

---

## Usage

```r
sci_fromfastq(group = "docker",
              folder = "/path/to/working_directory",
              sample.name = "Experiment_01",
              UMI.cutoff = 500)
```

### Parameters

| Parameter       | Type      | Description |
|----------------|----------|-------------|
| `group` | character | `"docker"` or `"sudo"`, depending on user permissions. |
| `folder` | character | Path to working directory containing input files. |
| `sample.name` | character | Name of the experiment/sample being processed. |
| `UMI.cutoff` | integer | Minimum UMI count per cell for inclusion. |

---

## Pipeline Steps in Detail

### Step 1: FASTQ Generation & Preprocessing

Extracts FASTQ files from raw BCL sequencing data (if applicable) and embeds **Read 1 (barcode and UMI)** information into Read 2 names for proper tracking.

```bash
bcl2fastq --runfolder-dir /path/to/sequencing_data --output-dir fastq_output
```

### Step 2: Poly-A Trimming

Removes **poly-A sequences** using **TrimGalore** to improve alignment efficiency.

```bash
trim_galore --paired --output_dir trimmed-fastq sample_R1.fastq.gz sample_R2.fastq.gz
```

### Step 3: Read Alignment (STAR)

Aligns reads to the reference genome using **STAR**, allowing for spliced-mapping of reads.

```bash
STAR --genomeDir /path/to/STAR_index \
     --readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
     --runThreadN 8 --outSAMtype BAM SortedByCoordinate
```

### Step 4: Filtering & Sorting BAM Files

Removes **rRNA reads** and **ambiguously mapped reads**, then sorts BAM files for downstream gene assignment.

```bash
samtools view -hb -q 30 aligned.bam | samtools sort -o filtered_sorted.bam
```

### Step 5: UMI Counting & Deduplication

Assigns reads to genes using a gene annotation BED file and computes UMI counts per sample, removing duplicate UMIs.

```bash
bedtools intersect -a filtered_sorted.bam -b gene_annotations.bed > gene_assigned.bam
umi_tools count --per-gene --gene-tag=XT --stdin=gene_assigned.bam --stdout=UMI_counts.txt
```

### Step 6: Quality Control & Visualization

Generates **knee plots** to visualize the UMI distribution per cell.

```r
Rscript knee-plot.R UMI_counts.txt knee_plot_output.png
```

Outputs a **final UMI count matrix** for downstream analysis.

```bash
gzip -c UMI_counts.txt > final-output/UMI.count.matrix.gz
```

---

## Output Files

| File | Description |
|------|-------------|
| `alignment_report_complete.txt` | Alignment statistics report. |
| `rRNA_report_complete.txt` | rRNA contamination statistics. |
| `UMI_report_complete.txt` | UMI counts per sample. |
| `final-output/rRNA.and.dup.rate.stats` | Final rRNA and duplication rate statistics. |
| `final-output/knee-plots/*.png` | Knee plots visualizing UMI distributions. |
| `prelim.UMI.count.rollup.gz` | Preliminary UMI count matrix. |
| `final-output/cell.annotations` | Filtered cell annotations. |

---

## Prerequisites

- **Docker** must be installed and running.
- The input directory must contain:
  - **FASTQ files** (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
  - **A reference genome** (`.fa`, `.gtf`).

---

## Example

```r
sci_fromfastq(group = "docker",
              folder = "/the/input/data",
              sample.name = "Experiment_01",
              UMI.cutoff = 500)
```

This command runs SCI-RNA-seq inside the **repbioinfo/sci_tomatrix_genome** Docker container, processing raw sequencing data into a gene expression matrix.

---

## Additional Notes for Materials & Methods

- **TrimGalore (v0.6.10)** was used for poly-A trimming.  
- **STAR (v2.7.10a)** performed read alignment to the reference genome.  
- **SAMtools (v1.7)** was used for BAM file processing.  
- **UMI-tools (v1.0.0)** was used for UMI counting and deduplication.  
- **BEDTools (v2.27.1)** assigned reads to genes based on annotation BED files.  
- **Knee plots** were generated using **R ggplot2**.  



## Single Cell RNA-Seq Analysis Pipeline

The `singlecell_*` functions provide an R-based frontend for executing the **Single Cell RNA-seq** pipeline, which processes raw single-cell sequencing data into a **gene expression matrix** and performs **downstream analysis**. The functions automate data processing using Docker-based execution, ensuring reproducibility and ease of use.

### Pipeline Overview

1. **Alignment & Indexing**  
   - Uses **Cell Ranger** to align single-cell reads to a reference genome.  
   - Generates a **gene-cell count matrix**.  

2. **Filtering Mitochondrial & Ribosomal Reads**  
   - Identifies and removes **low-quality cells** based on **rRNA and mitochondrial gene expression**.  

3. **Clustering & Stability Analysis**  
   - Uses **Seurat** to cluster cells and assess cluster stability via bootstrapping.  

4. **Feature Selection (Differential Expression Analysis)**  
   - Identifies **differentially expressed genes (DEGs)** across clusters using **ANOVA, MAST, and edgeR**.  

5. **Enrichment Analysis**  
   - Performs **KEGG/GO pathway enrichment** on differentially expressed genes.  

---

## 1. Alignment & Indexing (`singlecell_alignIndex` function)

Aligns reads and indexes the genome for single-cell RNA-seq experiments.

### Usage

```r
singlecell_alignIndex(input_dir_path = "/path/to/fastq_files",
                      genome_dir_path = "/path/to/genome_files",
                      bamsave = TRUE)
```

### Output Files

| File | Description |
|------|-------------|
| `filtered_feature_bc_matrix/*.mtx` | Filtered gene expression matrix. |
| `filtered_feature_bc_matrix/*.barcodes.tsv.gz` | Modified barcodes including sample names. |
| `filtered_feature_bc_matrix/*.features.tsv.gz` | Gene feature information. |

---

## 2. Filtering Mitochondrial & Ribosomal Reads (`singlecell_mitoRiboFilter` function)

Filters low-quality cells based on rRNA and mitochondrial gene expression.

### Usage

```r
singlecell_mitoRiboFilter(input_file_path = "/path/to/matrix.mtx",
                          mito_min = 0, mito_max = 10,
                          ribo_min = 0, ribo_max = 10,
                          genes_file = "genes.tsv",
                          barcodes_file = "barcodes.tsv")
```

### Output Files

| File | Description |
|------|-------------|
| `filtered_matrix.mtx` | Filtered gene expression matrix. |
| `filtered_features.tsv` | Filtered gene list. |
| `filtered_barcodes.tsv` | Filtered barcode list. |
| `ribosomal_vs_mitochondrial_plot.png` | Quality control plot. |

---

## 3. Clustering & Stability Analysis (`singlecell_clustering` function)

Clusters cells and assesses stability using **Seurat** and **bootstrapping**.

### Usage

```r
singlecell_clustering(input_file_path = "/path/to/matrix.mtx",
                      bootstrap_percentage = 0.1,
                      stability_threshold = 0.8,
                      permutations = 10,
                      genes_file = "genes.tsv",
                      barcodes_file = "barcodes.tsv",
                      resolution = 0.8)
```

### Output Files

| File | Description |
|------|-------------|
| `clustering_stability.output.csv` | Cluster assignments and stability scores. |
| `umap_clusters.png` | UMAP plot of clusters. |
| `umap_stability.png` | UMAP plot of cluster stability. |

---

## 4. Feature Selection (`singlecell_featureSelection` function)

Identifies differentially expressed genes (DEGs) across clusters.

### Usage

```r
singlecell_featureSelection(input_file_path = "/path/to/matrix.mtx",
                            clustering_file = "clustering_stability.output.csv",
                            threshold = 0.8, log2fc = 1, pvalue = 0.05,
                            genes_file = "genes.tsv",
                            barcodes_file = "barcodes.tsv",
                            heatmap = TRUE)
```

### Output Files

| File | Description |
|------|-------------|
| `anova_DE_results.csv` | DEG results for all clusters. |
| `volcano_plot.png` | Volcano plot of DEGs. |
| `heatmap.png` | Heatmap of DEGs. |

---

## 5. Enrichment Analysis (`enrichmentAnalysis` function)

Identifies overrepresented pathways among DEGs.

### Usage

```r
enrichmentAnalysis(input_file_path = "/path/to/anova_DE_results.csv",
                   species = "dmelanogaster",
                   source = "KEGG",
                   separator = ",",
                   max_terms = 20)
```

### Output Files

| File | Description |
|------|-------------|
| `enrichment_plot.pdf` | Bar plot of enriched pathways. |

---

## Prerequisites

- **Docker** must be installed and running.
- The input directory must contain:
  - **FASTQ files** (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
  - **Reference genome** (`.fa`, `.gtf`).
  - **Gene expression matrix** (`.mtx`, `.csv`).

---

## Example Workflow

```r
singlecell_alignIndex("/data/fastq", "/data/genome", TRUE)
singlecell_mitoRiboFilter("/data/matrix.mtx", mito_min = 0, mito_max = 10, ribo_min = 0, ribo_max = 10)
singlecell_clustering("/data/filtered_matrix.mtx", 0.1, 0.8, 10, genes_file = "genes.tsv", barcodes_file = "barcodes.tsv", resolution = 0.8)
singlecell_featureSelection("/data/matrix.mtx", "clustering_stability.output.csv", 0.8, 1, 0.05)
enrichmentAnalysis("/data/anova_DE_results.csv", "dmelanogaster", "KEGG", ",", 20)
```

