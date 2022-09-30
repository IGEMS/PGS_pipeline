# PGS_pipeline
Pipeline to compute polygenic scores

This repository includes information and example codes on how polygenic scores (PGS) are computed in IGEMS Studies. 

The IGEMS PGS workgroup has developed protocols for constructing PGS in a harmonized manner across IGEMS studies. A variety of SNP arrays have been used, but the majority have been typed on the Illumina OmniExpress and Infinium PsychArray BeadChip. Other chips include the Illumina 670k and HumanCoreExome. Data from the various arrays have all been imputed to the 1000 Genomes Project phase 1 version 3 reference panel and/or the Haplotype Reference Consortium panel. All with genome-wide data are of European ancestry. 

PGS have been calculated for IGEMS using a common set of genetic variants (n=952,885) with good imputation quality (INFO score >0.8) and presence (MAF>1%) in all IGEMS studies. To deal with linkage disequilibrium, SBayesR shrinkage has been applied to the GWAS summary results prior to PGS calculations. When IGEMS studies have been part of the GWAS, new summary statistics excluding those studies (i.e., “leave one out”) has first been obtained to avoid sample overlap. 

Tutorial paper (PGS in general): Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020). https://doi.org/10.1038/s41596-020-0353-1 

SBayesR reference: Lloyd-Jones, L.R., Zeng, J., Sidorenko, J. et al. Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nat Commun 10, 5086 (2019). https://doi.org/10.1038/s41467-019-12653-0
