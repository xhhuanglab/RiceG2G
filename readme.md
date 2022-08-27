# RiceG2G  

### General Introduction 
GWAS in rice is difficult to reach single-gene resolution due to modest linkage disequilibrium decay.In order to systematically prioritize causal genes and variants at trait-associated loci, we developed an integrated genomics approach named RiceG2G. RiceG2G is an integrated genomics approach to systematically prioritize causal genes and variants at trait-associated loci. Which combined genetic associations with gene annotation, transcriptomics, functional genomics and intolerant variants, and produced a well-calibrated score. Briefly speaking, for each gene at a given associated locus, we collected its annotated type , expression pattern in vast spatio-temporal tissues, functional studies of homologs in other plants, variant interpretation, linkage mapping signal in each family and variant-signal correlations among families. Depending on the given trait, a summarized G2G score was then calculated for each candidate gene, and the score rank of candidate genes at an associated locus could reflect the likelihood of any certain gene to be the causal one. Moreover, once a strong candidate gene was determined, the most possible causal QTN could be provided based on variant interpretation, owing to nearly complete information of genetic variants from genome assemblies of parental lines.
### Availability
You can download the source code and the Demo data here:
[RiceG2G](https://github.com/xhhuanglab/RiceG2G)

###  Installation
```
git clone https://github.com/xhhuanglab/RiceG2G
```
### Data preparation
The input file of "RiceG2G" is GWAS sites file output of FarmCPU. Or you can organize your GWAS file into test file format.



The potential impacts of the variants on gene function were quantificationally evaluated for the use in RiceG2G. The intolerance score of indels and SVs in coding regions was based on whether the variant caused frameshift or was located in protein motifs (Lu et al. 2020 Nucleic Acids Res), e.g., setting the score as 25 for frameshift indels. The intolerance score of nonsynonymous SNPs was based on their conservation level, e.g., setting the score as 20 when SIFT=0 (Vaser et al 2016 Nature Protoc). The intolerance score of variants in noncoding regions was based on gene structure (whether located in exon-intron junction, UTR or core promoter regions) and chromatin accessibility (profiled from ATACseq datasets). For example, indels that disrupted exon-intron splicing were scored as 10, and SVs located in chromatin accessibility regions were scored as 5. 