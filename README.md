# The Problem of Selection Bias in Studies of pre-mRNA Splicing

#### Zachary W. Dwyer and Jeffrey A. Pleiss*

This document contains the data processing and analysis done for "[The Problem of Selection Bias in Studies of pre-mRNA Splicing]()" Custom scripts can be found in the "scripts" folder of this project. Questions and comments can be sent to Zach Dwyer at zwd2@cornell.edu.

### Genome and Annotation Files
The R64-2-1 release of the *Saccharomyces cerevisiae* genome was downloaded from the [Saccharomyces Genome Database](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/). The genome sequence was manually removed from the feature file and saved as sc_feature_R64-2-1.gff (available in resources). The chromosome names of the genome file were manually renamed to match the feature file and saved as sc_genome_R64-2-1.fa.

#### Build HISAT index:
Exon and intron ranges were extracted from sc_feature_R64-2-1.gff. Hisat indexes were built with intron and exon annotations (indexes available in resources/hisat_index, splice site and exon annotatios available in resources).
```
python sc_extract_exons_for_hisat.py sc_feature_R64-2-1.gff > sc_exons.txt
python sc_extract_introns_for_hisat.py sc_feature_R64-2-1.gff > sc_splice_sites.txt

hisat2-build --ss sc_splice_sites.txt --exon sc_exons.txt sc_genome_R64-2-1.fa sc_index_R64-2-1
```

# Figure 1A

## Read Processing

## Adapter Trimming
Reads were trimmed using [fastp](https://github.com/OpenGene/fastp)
```
fastp -w 1 --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT -i raw/WT_A_R1.fastq.gz -I raw/WT_A_R2.fastq.gz -o full/trim/WT_A_R1.fastq.gz -O full/trim/WT_A_R2.fastq.gz --html trim_report/WT_A.html --json trim_report/WT_A.json 2> logs/trim/WT_A.log
```

## Alignment
Trimmed reads were aligned using [hisat2](http://daehwankimlab.github.io/hisat2/manual/)
```
hisat2 -p 8 --new-summary --max-intronlen 2000 --no-unal -x resources/hisat/sc_index_R64-2-1 -1 full/trim/WT_A_R1.fastq.gz -2 full/trim/WT_A_R1.fastq.gz 2> logs/align/WT_A.log | samtools view -bh -q 5 | samtools sort -o full/align/WT_A.bam
```

## Counting
Premature and mature alignments were counted using custom script based off of [HTSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html)
```
python2.7 scripts/feature_count.py --input full/align/WT_A.bam --output full/count/WT_A.txt
```

