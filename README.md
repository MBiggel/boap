# # BOAP: Bacterial ONT Assembly Pipeline

**BOAP** is a Nextflow-based automated pipeline for the assembly, polishing, and quality control of bacterial genomes using Oxford Nanopore Technologies (ONT) Whole Genome Sequencing data.

BOAP takes raw reads (FASTQ or BAM) and produces high-quality, circularized, and polished assemblies with reports on assembly completeness and accuracy.


## Pipeline overview

1.  **Input**: Accepts raw `.fastq.gz` or `.bam` files generated using SUP/HAC (`v5.0` / `v5.2`) basecalling.  
   Inclusion of the basecalling model in the read header is recommended. BAM files are automatically converted to FASTQ.
2.  **Quality Filtering (Nanoq)**: Filters reads based on minimum quality (Q10) and minimum length (1,000 bp).
3.  **Genome Size Estimation**: Uses `Raven` to estimate genome size from a rapid preliminary assembly. Alternatively accepts  user-provided size (e.g., `5m`) to skip raven assembly.
4.  **Downsampling (Filtlong)**: Downsamples high-quality data to a target coverage (default 100x).
5.  **Assembly (Flye)**: Performs *de novo* assembly using Flye with the `--nano-hq` mode.
6.  **Circularization (Dnaapler)**: Re-orients circular chromosomes and plasmids to standardized start positions.
7.  **Polishing (Medaka2)**: Polishes the assembly using ONT's `medaka` with the `--bacteria` model for high consensus accuracy.
8.  **QC & Reporting**: Calculates accuracy via `Alpaqa` and summarizes contig lengths and coverage information.



## Installation & Configuration

### Prerequisites
* [Nextflow](https://www.nextflow.io/) (>=21.10.0)
* [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda)

### Download boap script and install dependencies
Following commands will install python, nanoq, raven-assembler, filtlong, flye, dnaapler, medaka, seqkit, samtools, and scipy into a dedicated boap conda environment:

```bash
git clone https://github.com/MBiggel/boap.git
cd boap
conda env create -f environment.yml
conda env list
```

### Configuration

Important: Before running, update the path to your conda environment in the nextflow.config file, for example

    process.conda = '/home/user/miniconda3/envs/boap'
		
	

##  Usage

### Example commands

#### Run the pipeline on a directory of FASTQ or BAM files
The pipeline will process all valid files in the folder. Genome sizes will be estimated automatically using Raven.

```
nextflow run /path/to/boap.nf -c /path/to/nextflow.config --input /path/to/reads/ --threads 30 
```

#### Run the pipeline on a single input file
Target coverage 80x. Provide the genome size to skip the estimation step. Use --cleanup to automatically delete the work directory after a successful run to save disk space. Include -bg to run the pipeline in the background.
```
nextflow run /path/to/boap.nf -c /path/to/nextflow.config --input "/path/to/reads/sample_01.fastq.gz" --gsize 5.2m --coverage 80 --threads 30 --outdir boap_results --cleanup -bg
```

#### Force a specific medaka2 model
Information on the basecalling model is usually stored in the fastq or bam files and automatically detected by medaka2. If the model information is missing in the input files, specific models can be provided via
```
nextflow run /path/to/boap.nf -c /path/to/nextflow.config --input "*.bam" --force-model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
```
  
  
## Parameters

| Parameter          | Default       | Description                                                                                                                                                                                                                                                                                     |
| ------------------ | ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`          | `input_fastq` | Directory containing `.fastq.gz` / `.bam` files or path to a single file                                                                                                                                                                                                                        |
| `--outdir`         | `boap_output` | Output directory for results                                                                                                                                                                                                                                                                    |
| `--coverage`       | `100`         | Target coverage for downsampling (Filtlong)                                                                                                                                                                                                                                                     |
| `--min_contig_len` | `1000`        | Minimum contig length (Flye)                                                                                                                                                                                                                                                                    |
| `--gsize`          | `null`        | Manual genome size (e.g., `5m`). If not set, calculated automatically via Raven                                                                                                                                                                                                                 |
| `--threads`        | `30`          | Maximum number of threads used for parallel processes                                                                                                                                                                                                                                           |
| `--force_model`    | `false`       | (Optional) Force a specific basecalling model for polishing (e.g., `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`). The pipeline automatically detects the basecalling model from the input FASTQ/BAM headers. This parameter should only be used if the model information is missing in the input |
 



## Output Structure
```
results/
├── assembly/                   # Final polished assemblies
│   ├── sampleA.fasta
│   └── sampleB.fasta
├── summary/                    # QC reports
│   ├── alpaqa_report.tsv       # Statistics on assembly accuracy
│   └── contig_report.tsv       # Contig info (length, coverage, circularity)
└── flye/                       # Intermediate Flye output (assembly graph/logs)
    └── sampleA/
        ├── assembly_graph.gfa
        ├── flye.log
        └── ...
```

The alpaqa report provides detailed information on assembly accuracy. See the [alpaqa](https://github.com/MBiggel/alpaqa) for more details.

## References
**Seqkit**: https://github.com/shenwei356/seqkit
Shen W, Sipos B, and Zhao L (2024). SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta 3, 1–5. doi:10.1002/imt2.191.

**Samtools**: https://github.com/samtools/samtools
Twelve years of SAMtools and BCFtools
Danecek, P, Bonfield, JK, Liddle, J, Marshall, J, Ohan, V, Pollard, MO, Whitwham, A, Keane, T, McCarthy, SA, Davies, RM and Li, H (2021). Twelve years of SAMtools and BCFtools. Gigascience


**Nanoq**: https://github.com/esteinig/nanoq
Steinig E, and Coin L (2022). Nanoq: Ultra-Fast Quality Control for Nanopore Reads. J. Open Source Softw. 7, 2991. doi:10.21105/joss.02991.

**Raven**: https://github.com/lbcb-sci/raven
Vaser R, and Šikić M (2021). Time- and Memory-Efficient Genome Assembly with Raven. Nat. Comput. Sci. 1, 332–336. doi:10.1038/s43588-021-00073-4.

**Filtlong**: https://github.com/rrwick/Filtlong

**Flye**: https://github.com/mikolmogorov/Flye
Kolmogorov M, Yuan J, Lin Y, and Pevzner PA (2019). Assembly of Long, Error-Prone Reads Using Repeat Graphs. Nat. Biotechnol. 37, 540–546. doi:10.1038/s41587-019-0072-8.

**DNAapler**: https://github.com/gbouras13/dnaapler
Bouras G, Grigson SR, Papudeshi B, Mallawaarachchi V, and Roach MJ (2024). Dnaapler: A Tool to Reorient Circular Microbial Genomes. J. Open Source Softw. 9, 5968. doi:10.21105/joss.05968.

**Medaka2**: https://github.com/nanoporetech/medaka

**Alpaqa**: https://github.com/MBiggel/alpaqa

