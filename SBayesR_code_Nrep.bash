### Code to call SBayes R within the PRS pipeline
### By Ida Karlsson, 20210628
### This code is for when N is reported, a separate code is used when it's not (adding the impute N option)

for i in {21..22}
do
/nfs/AGE/IGEMS/IGEMS_PRS_pipeline/SBayesR_and_RefData/gctb_2.03beta_Linux/gctb \
--sbayes R \
--ldm ../LDmatrix/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${i}_v3_50k.ldm.sparse \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--gwas-summary sumstatfile5.txt \
--seed 123 \
--chain-length 10000 \
--burn-in 2000 \
--robust \
--out-freq 10 \
--out SBayesR_chr_${i} 
done
