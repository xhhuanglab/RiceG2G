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
![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The test data takes the HZ_Awn_length as an example


### Data preparation
######<font color=#008000>Basic data preparation</font>
**Variant interpretation:** We have prepared all variants found within the NAM population, which need to be downloaded from figshare and called at the  path "/RiceG2G/basic_data" before running.  
**figshare download path:**
https://figshare.com/articles/dataset/RiceG2G_Basicdata3/21117010  
https://figshare.com/articles/dataset/RiceG2G_Basicdata2/21117007  
https://figshare.com/articles/dataset/RiceG2G_Basicdata1/21115783    
**Genome annotation:** "Rice_MSUv7.gff3", "Rice_IRGSP-1.0.gff3".  
**Gene type:** "all.locus_brief_info.7.0" :This file lists the information about the  transposon or retro-transposon elements.  
**house-keeping_gene.list** ：This file lists the information about the house-keeping gene.  
**repressed_geno.list**：H3K27me3-related regions in panicle, root, mature leaf and seedling were downloaded from the RiceENCODE website.  
**Homolog function:** "Basic_Information_Arabidopsis_to_rice.txt", "Basic_Information_Arabidopsis_to_rice_rap.txt", "RAP-MSU_2021-11-11.txt".  
**Expression pattern:** "combineRAPvsMSU_RNAseq.txt".  
######<font color=#008000>Input data preparation</font>
```
cd /RiceG2G
```

**GWAS results:** The input file of "RiceG2G" is GWAS sites file output of FarmCPU. Or you can organize your GWAS file like FarmCPU format.
```
perl ./expand/SelectFarmCPUPeakByCor.pl HZ_Awn_length
perl ./expand/Select_MLM_Peak.pl HZ_Awn_length

```
**WinQTLCart results: **QTL mapping results for 15 populations. You can run the script in the expand file to get. (Note: that this output is run by winqtlcart. You need to put the .map file and .qrt file of winqtlcart into /RiceG2G/data/WinQTLcart.)
```
perl ./expand/WinQTLcart_lod_geno.pl HZ_Awn_length
```
**Select WinQTLCart Peak: **
```
perl ./expand/Select_Winqtlcart_Peak.pl HZ_Awn_length
```

######<font color=#008000>Start RiceG2G</font>
 Collected annotated type , expression pattern in vast spatio-temporal tissues, functional studies of homologs in other plants, variant interpretation, linkage mapping signal in each family and variant-signal correlations among families. 
```
perl ./code/select_FarmCPUpeak_info_RAP.pl HZ_Awn_length young_panicle
perl ./code/RAP_Stpe2_select_FarmCPUpeak_info.pl Awn_length 
```




