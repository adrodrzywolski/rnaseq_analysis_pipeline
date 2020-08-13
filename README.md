# rnaseq analysis pipeline (RAP)
Created by: Adrian Odrzywolski

## Used tools

Simple pipeline to perform rna-seq data preparation.

Pipeline was made using:
* [nextflow](https://github.com/nextflow-io/nextflow) - wrapping whole pipeline,
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - quality control for fastq,
* [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - bbduk as a reads trimmer,
* [hisat2](http://daehwankimlab.github.io/hisat2/) - aligment,
* [samtools](http://www.htslib.org/) - sam/bam manipulation,
* [subread](http://subread.sourceforge.net/) - counting features per exon,
* [qualimap](http://qualimap.bioinfo.cipf.es/) - bam metricies,
* [multiqc](https://github.com/ewels/MultiQC) - final raport.


## Installation

Get nextflow script located on github.

```
git clone https://github.com/adrodrzywolski/rnaseq_analysis_pipeline.git
```
Get docker image:
```
docker pull adrodrzywolski/basic_rnaseq_prep:0.1
```
Get nextflow script:
```
curl -fsSL https://get.nextflow.io | bash
```
## Usage

Basic usage
```bash
./nextflow run rnaseq_analysis_pipeline.nf --reads "path/to/sample/SRR1039509*_[1,2].fastq*" --genome "Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes/1.fa" -with-docker "adrodrzywolski/basic_rnaseq_prep:0.1" --annotation "Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" --output "output"
```

```
Required arguments:
--reads - path to single- or paired-end reads in form of regular exporession. File can be gziped
--genome - path to reference genome file
--annotation - path to annotation file

Optional arguments:
--output - select directory path where all output needs to be stored [default: .]
```

## Output

After run, you'll get files prefixed with fastq file names:
* raw reads fastqc raport (zip,html)
* trimmed reads fastqc raport (zip,html)
* Hisat2 summary (txt)
* qualimap report (whole directory)
* **feature counts (summary)**
* **multiqc report (html)**

## License
The rnaseq pipeline analysis is released under the MIT license.