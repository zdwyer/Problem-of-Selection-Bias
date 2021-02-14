# The Problem of Selection Bias in Studies of pre-mRNA Splicing

#### Zachary W. Dwyer and Jeffrey A. Pleiss*

This document contains the data processing and analysis done for "[The Problem of Selection Bias in Studies of pre-mRNA Splicing]()" Custom scripts can be found in the "scripts" folder of this project. Questions and comments can be sent to Zach Dwyer at zwd2@cornell.edu.

## Dependencies

| **Package** | **Version** |
|:------------|:------------|
| cowplot     | 1.0.0       |
| DESeq2      | 1.26.0      |
| dplyr       | 0.8.5       |
| ggplot2     | 3.3.0       |
| ggbeeswarm  | 0.6.0       |
| gridextra   | 2.3         |
| reshape2    | 1.4.4       |
| scales      | 1.1.0       |
| tidyr       | 1.0.2       |

# Data Processing

### Genome and Annotation Files
The R64-2-1 release of the *Saccharomyces cerevisiae* genome was downloaded from the [Saccharomyces Genome Database](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/). The genome sequence was manually removed from the feature file and saved as sc_feature_R64-2-1.gff (available in resources). The chromosome names of the genome file were manually renamed to match the feature file and saved as sc_genome_R64-2-1.fa.

#### Build HISAT index:
Exon and intron ranges were extracted from sc_feature_R64-2-1.gff. Hisat indexes were built with intron and exon annotations (indexes available in resources/hisat_index, splice site and exon annotatios available in resources).
```
python sc_extract_exons_for_hisat.py sc_feature_R64-2-1.gff > sc_exons.txt
python sc_extract_introns_for_hisat.py sc_feature_R64-2-1.gff > sc_splice_sites.txt

hisat2-build --ss sc_splice_sites.txt --exon sc_exons.txt sc_genome_R64-2-1.fa sc_index_R64-2-1
```

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

# Analysis

## Set-up R

### Load Libraries

```
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggbeeswarm)
library(scales)
library(cowplot)
library(gridExtra)
```

### Customize ggPlot Theme

```
theme_set(theme_classic())
theme_update(axis.ticks = element_line(color='black'),
             plot.title = element_text(hjust = 0.5))
```

### Load Custom Functions

```
options(scipen=999)
scaleFUN <- function(x) sprintf("%.3f", x)

DESeq_analysis = function(data, sizes, seeds) {

  deseq_results = data.frame(Gene=factor(), Rank=factor(), Type=factor(), baseMean=numeric(), log2FoldChange=numeric(), padj=numeric(), Seed=integer(), Size=integer())

  for(size in sizes) {
    for(seed in seeds) {
      current = data %>% filter(Size==size, Seed==seed) %>%
          melt(id.vars=c("Strain", "Replicate", "Gene", "Rank"), measure.vars=c("Mature", "Premature"), variable.name="Status", value.name="Count") %>%
          unite(Feature, c("Gene", "Rank", "Status"), sep=';') %>%
          unite(Sample, c("Strain", "Replicate"), sep='_') %>%
          select(Feature, Sample, Count) %>%
          dcast(Feature ~ Sample, value.var='Count')
      DEseq_input = current[,-1]
      rownames(DEseq_input) = current[,1]
      DEseq_design = data.frame(row.names = colnames(DEseq_input), Strain = c("WT", "WT", "WT", "prp2-1", "prp2-1", "prp2-1"))
      dds = DESeqDataSetFromMatrix(countData=DEseq_input, colData=DEseq_design, design=~Strain)
      dds$Strain = relevel(dds$Strain, ref='WT')
      
      current_result = tibble::rownames_to_column(as.data.frame(results(DESeq(dds, quiet=TRUE))), "Feature") %>% 
          separate(Feature, c("Gene", "Rank", "Type"), sep=';') %>%
          select(Gene, Rank, Type, baseMean, log2FoldChange, padj) %>%
          mutate(Size = size, Seed= seed)
          
      deseq_results = rbind(deseq_results, current_result)
    }
  }
  return(deseq_results)
}
```

## Figure 1

### Read Auxilary Files
```
curated = read.delim("Currated_List_zwd.txt", header=FALSE, col.names = c("Intron")) %>% mutate(Intron_2 = Intron) %>% separate(Intron_2, c("Gene", "Rank"), sep=";")

intron_properties = read.delim("sc_intron_properties(barrass).txt", header=TRUE) %>% 
                    unite(Intron, c("Gene", "Rank"), sep=';') %>%
                    filter(Intron %in% curated$Intron)
```

### Read Full Dataset
```
full_raw = read.delim("shifted_full.txt", header=TRUE) %>%
                filter(Intron %in% curated$Intron) %>%
                separate(Intron, into=c("Gene", "Rank"), sep=";") %>%
                separate(Strain, into=c("Strain", "Condition"), sep='_')
full_raw$Strain = factor(full_raw$Strain, levels=c("scJS1", "scJS7"))
full_deseq = DESeq_analysis(data=full_raw, sizes=c(0), seeds=c(0))

mpe_expression = full_raw %>% filter(Size==0, Seed==0, Strain=='scJS1') %>%
                    select(Replicate, Gene, Rank, Mature, Premature) %>%
                    group_by(Gene, Rank) %>%
                    summarise(Mature_Average = mean(Mature), Premature_Average=mean(Premature)) %>%
                    unite(Intron, c(Gene, Rank), sep=';')
intron_properties = intron_properties %>% merge(mpe_expression)
```

### Figure 1A

```
fold_change = full_deseq  %>% dcast(Gene+Rank ~ Type, value.var='padj') %>%
mutate(Significant = case_when(is.na(Premature) & is.na(Mature) ~ "NA",
                               Mature <.05 & Premature <.05 ~ "Both",
                               Mature <.05 & (Premature >.05 | is.na(Premature)) ~ "Mature",
                               (Mature >.05 | is.na(Mature)) & Premature <.05 ~ "Premature",
                               (Mature >.05 | is.na(Mature)) & (Premature >.05 | is.na(Premature)) ~ "Neither")) %>%
select(-Premature, -Mature) %>%
unite(Intron, c(Gene,Rank), sep=';') %>%
merge(full_deseq %>% dcast(Gene+Rank~Type, value.var="log2FoldChange") %>% unite(Intron, c(Gene, Rank), sep=';'), by="Intron") %>%
group_by(Significant) %>%
arrange(Premature-Mature) %>%
filter(Significant != "NA") %>%
melt(id.vars=c("Intron", "Significant"), measure.vars=c("Premature", "Mature"), variable.name="Type", value.name="log2FoldChange")

fold_change$Type = factor(fold_change$Type, levels=c("Mature", "Premature"))  
fold_change$Intron = factor(fold_change$Intron, levels=(fold_change %>% filter(Type=='Premature'))$Intron)
fold_change$Significant = factor(fold_change$Significant, levels=c("Both", "Mature", "Premature", "Neither", "NA"))
  
ggplot(fold_change, aes(x=Type, y=Intron, fill=2^log2FoldChange, color=2^log2FoldChange)) +
  facet_grid(Significant~., scales='free_y', space='free_y') +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        legend.position = 'bottom',
        legend.margin = margin(0,0,0,0),
        panel.background = element_blank(),
        panel.spacing = unit(0, 'cm'),
        plot.background = element_blank(),
        plot.margin = margin(0,0,0,0)) +
  scale_fill_gradient2(trans='log2', low='dodgerblue3', mid='grey90', high='#eadf0c', limits=c(.125,8), oob=squish, name=element_blank(), labels=scaleFUN) +
  scale_color_gradient2(trans='log2', low='dodgerblue3', mid='grey90', high='#eadf0c', limits=c(.125,8), oob=squish, name=element_blank(),labels=scaleFUN) +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0)) +
  geom_tile()

fold_change %>% filter(Type=="Premature") %>% 
  group_by(Significant) %>%
  summarise(Count=n())
```
### Read Downsampling Data

```
downsample_raw = read.delim("shifted_specific2.txt", header=TRUE) %>%
                filter(Intron %in% curated$Intron) %>%
                separate(Intron, into=c("Gene", "Rank"), sep=";") %>%
                separate(Strain, into=c("Strain", "Condition"), sep='_')
downsample_raw$Strain = factor(downsample_raw$Strain, levels=c("scJS1", "scJS7"))
downsample_deseq = DESeq_analysis(data=downsample_raw, sizes=c(200000, 400000, 800000), seeds=c(1))
```

### Figure 1B

```
ma_layers = list(
    theme(legend.position = 'none',
        axis.text=element_blank(),
        axis.title=element_blank(),
        plot.margin = margin(0,0,0,0)),
    scale_x_continuous(trans='log10', name="Normalized Counts", breaks = c(.1, 1, 10, 100, 1000, 10000, 100000), limits=c(.1, 100000)),
    scale_color_manual(values=c('black', '#cb181d')),
    geom_hline(yintercept = 1, color='grey', linetype=2),
    geom_point(size=.1))

high_mature = shifted_specific_deseq %>% filter(Seed==1, Size==800000, Type=='Mature', !is.na(padj)) %>% mutate(Significant = padj < .05)
high_mature_ma = ggplot(high_mature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                    ma_layers +
                    scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.001, .004, .016, .064, .25, 1, 2), limits=c(1/1200, 2), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_high_mature.pdf", plot=high_mature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)

high_premature = shifted_specific_deseq %>% filter(Seed==1, Size==800000, Type=='Premature', !is.na(padj)) %>% mutate(Significant = padj < .05)
high_premature_ma = ggplot(high_premature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                      ma_layers +
                      scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.125, 1, 8, .064, 64), limits=c(.125, 64), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_high_premature.pdf", plot=high_premature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)


medium_mature = shifted_specific_deseq %>% filter(Seed==1, Size==400000, Type=='Mature', !is.na(padj)) %>% mutate(Significant = padj < .05)
medium_mature_ma = ggplot(medium_mature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                    ma_layers +
                    scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.001, .004, .016, .064, .25, 1, 2), limits=c(1/1200, 2), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_medium_mature.pdf", plot=medium_mature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)

medium_premature = shifted_specific_deseq %>% filter(Seed==1, Size==400000, Type=='Premature', !is.na(padj)) %>% mutate(Significant = padj < .05)
medium_premature_ma = ggplot(medium_premature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                      ma_layers +
                      scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.125, 1, 8, .064, 64), limits=c(.125, 64), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_medium_premature.pdf", plot=medium_premature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)


low_mature = shifted_specific_deseq %>% filter(Seed==1, Size==200000, Type=='Mature', !is.na(padj)) %>% mutate(Significant = padj < .05)
low_mature_ma = ggplot(low_mature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                    ma_layers +
                    scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.001, .004, .016, .064, .25, 1, 2), limits=c(1/1200, 2), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_low_mature.pdf", plot=low_mature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)

low_premature = shifted_specific_deseq %>% filter(Seed==1, Size==200000, Type=='Premature', !is.na(padj)) %>% mutate(Significant = padj < .05)
low_premature_ma = ggplot(low_premature, aes(x=baseMean, y=2^log2FoldChange, color=Significant)) + 
                      ma_layers +
                      scale_y_continuous(trans='log2', name="Fold Change", breaks = c(.125, 1, 8, .064, 64), limits=c(.125, 64), labels=scaleFUN)
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1B_low_premature.pdf", plot=low_premature_ma, device='pdf', width=2.25,height=2.25, units=c("cm"), useDingbats=FALSE)


counts = shifted_specific_deseq %>% filter(Seed==1, padj<.05) %>%
            group_by(Type, Size) %>%
            summarise(count=n())
```

### Figure 1C

```
boxplot_layers = list(
        facet_grid(~Size, switch='x'),
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none',
        plot.margin = margin(0,0,0,0)),
  scale_color_manual(values=c("#cb181d", "black")),
  geom_quasirandom(size=.05))

boxplot_data = shifted_specific_deseq %>% filter(Seed==1) %>%
    dcast(Gene+Rank+Size ~ Type, value.var='padj') %>%
    mutate(Significant = case_when(is.na(Premature) & Mature < .05 ~ "Significant",
                                   is.na(Mature) & Premature < .05 ~ "Significant",
                                   Mature < .05 & Premature < .05 ~ "Significant",
                                   TRUE ~ "Non-Significant")) %>%
   unite(Intron, c(Gene, Rank), sep=';') %>%
   merge(intron_properties, by='Intron') %>%
   select(Intron, Size, Significant, three_score, intron_length)
   
boxplot_data$Size = factor(boxplot_data$Size, levels=c("800000", "400000", "200000"))
boxplot_data$Significant = factor(boxplot_data$Significant, levels=c("Significant", "Non-Significant"))

three_score_boxplot = ggplot(boxplot_data, aes(x=Significant, y=three_score, color=Significant)) +
  boxplot_layers

intron_length_boxplot =ggplot(boxplot_data, aes(x=Significant, y=intron_length, color=Significant)) +
  boxplot_layers +
  scale_y_continuous(trans='log10', limits=c(50,1100), breaks=c(100, 500, 1000))

intron_length_boxplot

ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1C_intron_length.pdf", plot=intron_length_boxplot, device='pdf', width=5.15, height=2.5, units=c("cm"), useDingbats=FALSE)

wilcox.test(intron_length ~ Significant, data=boxplot_data %>% filter(Size==800000), alternative='two.sided')
wilcox.test(intron_length ~ Significant, data=boxplot_data %>% filter(Size==400000), alternative='two.sided')
wilcox.test(intron_length ~ Significant, data=boxplot_data %>% filter(Size==200000), alternative='two.sided')

wilcox.test(intron_properties$intron_length, (boxplot_data %>% filter(Size==200000, Significant=="Significant"))$intron_length)
```

### Figure 1D

```
quasirandom_layers = list(
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = 'none'),
  scale_alpha_manual(values=c(0,.5)),
  geom_quasirandom(color='#cb181d', size=.001))

fold_change_high = shifted_specific_deseq %>% filter(Seed==1) %>% dcast(Gene+Rank+Size~Type, value.var="log2FoldChange") %>% unite(Intron, c(Gene, Rank), sep=';') %>% filter(Size==800000) %>% select(-Size) %>% mutate(SI_ratio = Premature-Mature)

significant = shifted_specific_deseq %>% filter(Seed==1) %>%
    dcast(Gene+Rank+Size ~ Type, value.var='padj') %>%
    mutate(Significant = case_when(is.na(Premature) & is.na(Mature) ~ "Non-Significant",
                                   is.na(Premature) & Mature < .05 ~ "Significant",
                                   is.na(Mature) & Premature < .05 ~ "Significant",
                                   Mature < .05 & Premature < .05 ~ "Significant",
                                   TRUE ~ "Non-Significant")) %>%
    select(-Premature, -Mature) %>%
    unite(Intron, c(Gene,Rank), sep=';') %>%
    merge(fold_change_high, by="Intron") %>%
    merge(mpe_expression, by="Intron")

significant$Size = factor(significant$Size, levels=c("800000", "400000", "200000"))   

premature_fc = ggplot(significant, aes(x=Size, y=2^Premature, alpha=Significant)) +
               scale_y_continuous(trans='log2', breaks=c(.016, .125, 1, 8, 64)) +
               geom_hline(yintercept = 1, color='grey', linetype=2) +
               quasirandom_layers
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1D_premature_fold_change.pdf", plot=premature_fc, device='pdf', width=2.25, height=2.25, units=c("cm"), useDingbats=FALSE)

wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Premature, (significant %>% filter(Size==400000, Significant=='Significant'))$Premature, alternative='two.sided')
wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Premature, (significant %>% filter(Size==200000, Significant=='Significant'))$Premature, alternative='two.sided')

nrow((significant %>% filter(Size==800000, Significant=='Significant')))
nrow((significant %>% filter(Size==400000, Significant=='Significant')))
nrow((significant %>% filter(Size==200000, Significant=='Significant')))

mature_fc = ggplot(significant, aes(x=Size, y=2^Mature, alpha=Significant)) +
  scale_y_continuous(trans='log2', limits=c(1/128, 2), breaks=c(.008, .032, .125, .5, 1, 2)) +
  geom_hline(yintercept = 1, color='grey', linetype=2) +
  quasirandom_layers

ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1D_mature_fold_change.pdf", plot=mature_fc, device='pdf', width=2.25, height=2.25, units=c("cm"), useDingbats=FALSE)

wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Mature, (significant %>% filter(Size==400000, Significant=='Significant'))$Mature, alternative='two.sided')
wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Mature, (significant %>% filter(Size==200000, Significant=='Significant'))$Mature, alternative='two.sided')

premature_expression = ggplot(significant, aes(x=Size, y=Premature_Average, alpha=Significant)) +
               scale_y_continuous(trans='log10', breaks=c(.1, 1, 10, 100, 1000)) +
               quasirandom_layers
ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1D_premature_expression.pdf", plot=premature_expression, device='pdf', width=2.25, height=2.25, units=c("cm"), useDingbats=FALSE)

wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Premature_Average, (significant %>% filter(Size==400000, Significant=='Significant'))$Premature_Average, alternative='less')
wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Premature_Average, (significant %>% filter(Size==200000, Significant=='Significant'))$Premature_Average, alternative='less')

mature_expression = ggplot(significant, aes(x=Size, y=Mature_Average, alpha=Significant)) +
               scale_y_continuous(trans='log10', breaks=c(.1, 1, 10, 100, 1000, 10000, 100000)) +
               quasirandom_layers

ggsave(filename="~/Documents/Lab/Pleiss/EditorLetter/Figures/Figure1D_mature_expression.pdf", plot=mature_expression, device='pdf', width=2.25, height=2.25, units=c("cm"), useDingbats=FALSE)

wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Mature_Average, (significant %>% filter(Size==400000, Significant=='Significant'))$Mature_Average, alternative='less')
wilcox.test((significant %>% filter(Size==800000, Significant=='Significant'))$Mature_Average, (significant %>% filter(Size==200000, Significant=='Significant'))$Mature_Average, alternative='less')
```

# Figure 2

Gene expression data was obtained from [GSM2535498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2535498)

## Figure 2A

```
rpg = read.delim("Ribosomal_Protein_Genes.txt", header=FALSE, col.names=c("CommonName", "Type", "Source"))

significant = read.delim("supplemental_table_3.txt", header=TRUE) %>%
    select(gene_id=Official.gene.symbol, Chromosome, sigRank=Position.Intron, sigStart=Start.Intron, sigEnd=End.Intron)

fpkm = read.delim("GSM2535498_ctrl48h_13529_4_D223KACXX_genes.fpkm_tracking.txt", header=TRUE) %>% 
  select(gene_id, locus, FPKM) %>% 
  separate(locus, c("first", "stop"), sep=c("-")) %>%
  separate(first, c("chromosome", "start")) %>%
  mutate(Significant=gene_id %in% significant$gene_id) %>%
  mutate(RPG=gene_id %in% rpg$CommonName) %>%
  mutate(Group = case_when(!Significant & !RPG ~ "Non-RPG",
                            !Significant & RPG ~ "RPG",
                             Significant & !RPG ~ "Significant Non-RPG",
                             Significant & RPG ~ "Significant RPG"))

figure2a_dotplot = ggplot(fpkm %>% filter(FPKM > 0), aes(x=Significant, y=FPKM, color=Group, alpha=Group)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(trans='log10', name='FPKM', limits=c(.001, 8000), breaks=c(.1, 10, 1000)) +
  scale_color_manual(values=c("#808080", "black", "#FF7879", "#cb181d")) +
  scale_alpha_manual(values=c(.006, 1, 1, 1)) +
  scale_size_manual(values=c(.01, .1,.1,.1)) +
  geom_jitter(shape=20)

figure2a_boxplot = ggplot(fpkm %>% filter(FPKM > 0), aes(x=Significant, y=FPKM, color=Group)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(trans='log10', name='FPKM', limits=c(.001, 8000), breaks=c(.1, 10, 1000)) +
  scale_color_manual(values=c("#808080", "black", "#FF7879", "#cb181d")) +
  geom_boxplot(outlier.shape = NA)

figure2a = grid.arrange(figure2a_dotplot, figure2a_boxplot, nrow=1)

figure2a
```

#### Significance Testing
```
wilcox.test((fpkm %>% filter(Group=='Non-RPG'))$FPKM, (fpkm %>% filter(Group=='Significant RPG'))$FPKM, alternative='two.sided')
wilcox.test((fpkm %>% filter(Group=='Non-RPG'))$FPKM, (fpkm %>% filter(Group=='Significant Non-RPG'))$FPKM, alternative='two.sided', exact=TRUE)
```

### Figure 2B

```
gene_end = read.delim("UCSC_genes.bed", header=FALSE, col.names=c("Chromosome", "Start", "Stop", "UCSC", "Score", "Orientation", "CDS_Start", "CDS_End", "Blank", "Exon_Count", "Exon_Starts", "Exon_Ends")) %>%
          select(UCSC, Transcript_Stop=Stop)

alias = read.delim("kgSpAlias.txt", header=TRUE) %>% select(UCSC=X.kgID, gene_id=alias)
introns = read.delim("introns.bed", header=FALSE, col.names=c("Chromosome", "Start", "Stop", "Name", "Extra", "Orientation")) %>%
          separate(Name, c("UCSC"), sep='_') %>%
          mutate(Start=Start+1)

sig_introns = significant %>% left_join(introns, by=c("Chromosome"="Chromosome", "sigStart"="Start", "sigEnd"="Stop")) %>%
              unite(Key, c("UCSC", "sigStart", "sigEnd"), sep=';')

rpg_intron = rpg %>% left_join(alias, by=c("CommonName"="gene_id"))

all_introns = introns %>% unite(Key, c("UCSC", "Start", "Stop"), sep=';') %>%
              mutate(Significant = Key %in% sig_introns$Key) %>%
              separate(Key, into=c("UCSC", "Start", "Stop"), sep=';') %>%
              left_join(gene_end, by=c("UCSC"="UCSC")) %>%
              mutate(Distance=Transcript_Stop-as.numeric(Stop)) %>%
              mutate(RPG=UCSC %in% rpg_intron$UCSC) %>%
              mutate(Group = case_when(!Significant & !RPG ~ "Non-RPG",
                            !Significant & RPG ~ "RPG",
                             Significant & !RPG ~ "Significant Non-RPG",
                             Significant & RPG ~ "Significant RPG")) %>%
              arrange(Distance) %>%
              distinct(Chromosome, Start, Stop, .keep_all = TRUE)
              

figure2b_dotplot = ggplot(all_introns, aes(x=Significant, y=Distance, color=Group, alpha=Group)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        #axis.text = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=c("#808080", "black", "#FF7879", "#cb181d")) +
  scale_alpha_manual(values=c(.006, 1, 1, 1)) +
  scale_size_manual(values=c(.01, .1,.1,.1)) +
  geom_jitter(shape=20)

  figure2b_boxplot = ggplot(all_introns, aes(x=Significant, y=Distance, color=Group)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=c("#808080", "black", "#FF7879", "#cb181d")) +
  geom_boxplot(outlier.shape = NA)
  
  figure2b = grid.arrange(figure2b_dotplot, figure2b_boxplot, nrow=1)
  
  figure2b
```  

### Significance Testing

```
wilcox.test((all_introns %>% filter(Group=='Non-RPG'))$Distance, (all_introns %>% filter(Group=='Significant RPG'))$Distance, alternative='two.sided')
wilcox.test((all_introns %>% filter(Group=='Non-RPG'))$Distance, (all_introns %>% filter(Group=='SignificantNon-RPG'))$Distance, alternative='two.sided', exact=TRUE)
```