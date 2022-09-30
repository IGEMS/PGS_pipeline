# PGS_pipeline
Pipeline to compute polygenic scores

This repository includes information and example codes on how polygenic scores (PGS) are computed in IGEMS Studies. 

The IGEMS PGS workgroup has developed protocols for constructing PGS in a harmonized manner across IGEMS studies. A variety of SNP arrays have been used, but the majority have been typed on the Illumina OmniExpress and Infinium PsychArray BeadChip. Other chips include the Illumina 670k and HumanCoreExome. Data from the various arrays have all been imputed to the 1000 Genomes Project phase 1 version 3 reference panel and/or the Haplotype Reference Consortium panel. All with genome-wide data are of European ancestry. 

PGS have been calculated for IGEMS using a common set of genetic variants (n=952,885) with good imputation quality (INFO score >0.8) and presence (MAF>1%) in all IGEMS studies. To deal with linkage disequilibrium, SBayesR shrinkage has been applied to the GWAS summary results prior to PGS calculations. 

The PGS is then computed using e.g. Plink. For each individual, the number of effect alleles (0, 1, or 2) at each SNP across the genome, weighted by the SNP effect size from the GWAS of the trait. This results in a weighted linear combinations of SNPs, ggregating information across genetic variants associated with a trait (each of which is likely to only have an extremely small effect by itself), including genetic variants that do not achieve standard significance thresholds. 

It is important to note that *any* overlap in inviduals between the discovery (GWAS summary statistics) and target (the IGEMS study) can result in substantial inflation of the predictive ability of the PGS. While the target data likely represent only a small fraction of the discovery sample, that is in fact less important than the proportion of overlapping individuals within the target data (see Wray et al.J Child Psychol Psychiatry 2014. https://doi.org/10.1111/jcpp.12295). Therefore, when IGEMS studies have been part of the GWAS, new summary statistics excluding those studies (i.e., “leave one out”) has first been obtained to avoid sample overlap. 

Tutorial paper (PGS in general): Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020). https://doi.org/10.1038/s41596-020-0353-1 

SBayesR reference: Lloyd-Jones, L.R., Zeng, J., Sidorenko, J. et al. Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nat Commun 10, 5086 (2019). https://doi.org/10.1038/s41467-019-12653-0

Lecture on PGS: https://play.ki.se/media/GWAS_sumstats_PGS/0_lt7j6708
