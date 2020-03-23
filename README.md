# Bioinformatics Things

A curated list of bioinformatics tools, frameworks, libraries, etc. Feel free to [contribute](CONTRIBUTING.md) anything.

## Command Line Tools

- **[parallel](https://www.gnu.org/software/parallel)**:  GNU **parallel** is a shell tool for executing jobs in parallel using one or more computers . [Here](https://www.biostars.org/p/63816/) are some examples
- **[xsv](https://github.com/BurntSushi/xsv)**: A fast CSV command line toolkit written in Rust 
- **[ripgrep](https://crates.io/crates/ripgrep)**: combines the usability of The Silver Searcher with the raw speed of grep
- **[seqtk](https://github.com/lh3/seqtk)**: Toolkit for processing sequences in FASTA/Q formats
- **[bedtools2](https://github.com/arq5x/bedtools2)**: The swiss army knife for genome arithmetic 

## Sequence Preprocessing

- **[fastp](https://github.com/OpenGene/fastp)**:  An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting/merging...)
- **[fastqc](https://github.com/s-andrews/FastQC)**:  A quality control analysis tool for high throughput sequencing data


## Sequence Alignment

- **[bwa](https://github.com/lh3/bwa)**:  Burrow-Wheeler Aligner for short-read alignment
- **[minimap2](https://github.com/lh3/minimap2)**: A versatile pairwise aligner for genomic and spliced nucleotide sequences 
- **[hisat2](https://daehwankimlab.github.io/hisat2)**: Graph-based alignment (Hierarchical Graph FM index) 
- **[STAR](https://github.com/alexdobin/STAR)**:  RNA-seq aligner 
- **[TopHat](http://ccb.jhu.edu/software/tophat/manual.shtml)**: Aligns RNA-Seq reads to a genome in order to identify exon-exon splice junctions.
- **[Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)**: An ultrafast, memory-efficient  short read aligner

## Variant Caller

- **[freebayes](https://github.com/ekg/freebayes)**: Bayesian haplotype-based polymorphism discovery and genotyping.
- **[GATK](https://software.broadinstitute.org/gatk/)**: Variant Discovery in High-Throughput Sequencing Data
- **[samtools/bcftools/htslib](https://github.com/samtools/samtools)**: A suite of tools for manipulating next-generation sequencing data

## Structural Variant Caller

- **[delly](https://github.com/dellytools/delly)**: Structural variant discovery by integrated paired-end and split-read analysis.
- **[lumpy](https://github.com/arq5x/lumpy-sv)**: lumpy: a general probabilistic framework for structural variant discovery.
- **[manta](https://github.com/Illumina/manta)**: Structural variant and indel caller for mapped sequencing data
- **[crest](http://www.stjuderesearch.org/site/lab/zhang)**: maps somatic structural variation in cancer genomes with base-pair resolution. [Here](https://www.nature.com/articles/nmeth.1628) is the paper.
- **[breakdancer](http://breakdancer.sourceforge.net/)**: Genome-wide detection of structural variants from next generation paired-end sequencing reads. [Here](http://www.nature.com/nmeth/journal/v6/n9/abs/nmeth.1363.html) is the paper.

## Copy Number Variant Caller

- **[cnvkit](https://github.com/etal/cnvkit)**: Copy number variant detection from targeted DNA sequencing.
- **[control-freec](http://boevalab.inf.ethz.ch/FREEC/index.html)**: Prediction of copy numbers and allelic content using deep-sequencing data.

## Workflow Language & Manager

- **[nextflow](https://www.nextflow.io)**: A fluent DSL modelled around the UNIX pipe concept, that simplifies  writing parallel and scalable pipelines in a portable manner.
- **[snakemake](https://github.com/snakemake/snakemake)**:  Create **reproducible and scalable** data analyses. Workflows are described via a human readable, Python based language.
- **[WDL](https://github.com/openwdl/wdl)**:  The **Workflow Description Language (WDL)** is a way to specify data processing workflows with a human-readable and writeable syntax.
- **[CWL](https://www.commonwl.org/)**:  An open standard for describing analysis workflows and tools in a way that makes them portable and scalable across a variety of software and hardware environments, from workstations to cluster, cloud, and high performance computing (HPC) environments.
- **[cromwell](https://github.com/broadinstitute/cromwell)**:A Workflow Management System geared towards scientific workflows.
- **[galaxy](https://usegalaxy.org/)** - a popular open-source, web-based platform for data intensive  biomedical research. Has several features, from data analysis to  workflow management to visualization tools.

## Pipelines

- **[Awesome-Pipeline](https://github.com/pditommaso/awesome-pipeline)**: A curated list of awesome pipeline toolkits.
- **[awesome-nextflow](https://github.com/nextflow-io/awesome-nextflow)**: A curated list of nextflow based pipelines.
- **[nf-core](https://nf-co.re/pipelines)**:  A collection of high quality Nextflow pipelines.

## Programming Laguage Library

- **[pysam](https://github.com/pysam-developers/pysam)**: Python wrapper for [samtools](https://github.com/samtools/samtools).
- **[rust-htslib](https://github.com/rust-bio/rust-htslib)**: HTSlib bindings and a high level Rust API for reading and writing BAM files.
- **[bam](https://gitlab.com/tprodanov/bam)**: Rust crate for reading and writing BAM and BGZIP files.

## Data Visualization

- **[d3](https://d3js.org/)**:  Bring data to life with SVG, Canvas and HTML.
- **[echarts](https://echarts.apache.org/zh/index.html)**:  A powerful, interactive charting and visualization library.
- **[matplotlib](https://matplotlib.org/)**:  Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python
- **[ggplot2](https://ggplot2.tidyverse.org/)**: Create Elegant Data Visualisations Using the Grammar of Graphics. [Here](https://ggplot2-book.org/) is a book
- **[plotly](https://plot.ly/)**: Modern Analytic Apps for the Enterprise
- **[bokeh]( https://bokeh.org)**: Publish Sophisticated Dashboards 
- **[antv](https://antv.vision/)**: Liven Data Lively
- **[circos](http://circos.ca/)**: Perl package for circular plots, which are well suited for genomic rearrangements.

## Database
- **[mongodb](https://github.com/mongodb/mongo)**: A general purpose, document-based, distributed database built for modern application
- **[mysql](https://www.mysql.com/)**: Open-Source Relational Database Management System.
- **[leveldb](https://github.com/google/leveldb)**:  LevelDB is a fast  key-value storage library written at Google that provides an ordered  mapping from string keys to string values. 

## Courses

- **[An Introduction to Applied Bioinformatics](http://readiab.org/)**

  >  An **I**ntroduction to **A**pplied **B**ioinformatics (or IAB) is a free, open source interactive text that introduces  readers to core concepts of bioinformatics in the context of their  implementation and application.

- **[bioinformatics](https://github.com/ossu/bioinformatics)**: Path to a free self-taught education in Bioinformatics! 

  > This is a **solid path** for those of you who want to complete a **Bioinformatics** course on your own time, **for free**, with courses from the **best universities** in the World.
  >
  > In our curriculum, we give preference to MOOC (Massive Open Online  Course) style courses because these courses were created with our style  of learning in mind.
  >
  > To become a bioinformatician, you have to learn quite a lot of  science, so be ready for subjects like; Biology, Chemistry, etc...

## Others

- **[MultiQC](http://multiqc.info/)**: Aggregate results from bioinformatics analyses across many samples into a single report.

