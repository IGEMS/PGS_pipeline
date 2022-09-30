################################################################################
# NAME: Template_PRS_SBayesR_STR_20210603.bash
# PURPOSE: Master code for creation of PRS scores in the Swedish Twin Registry data
#
# CREATED: by Ida Karlsson (IK) 20210603
# UPDATED: by IK 20211025 - Update the minor allele when no MAF is available
# (was an error earlier, where the MAF was not correctly updated to represent A1)	
#
# NOTE!! Code should be updated with sample and sumstat info, and saved in the folder
# /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Programs under the name *PHENOTYPE*_*Author*_PRS_SBayesR_STR_*DATE*.bash
# The data will by default be saved in the /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Programs folder,
# and the log in the /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Logs folder, using the same naming convention
#
# To run: nohup bash /nfs/AGE/IGEMS/IGEMS_pheno_PRS/Programs/*PHENOTYPE*_*Author*_PRS_SBayesR_STR_*DATE*.bash &
#
################################################################################
######################## 1. SPECIFY PHENOTYPE AND DATE #########################
################################################################################

### MODIFY
########## 1.1 Specify phenotype (replace BMI with relevant phenotype)
PHENO=BMI

### MODIFY
########## 1.2 Specify today's month and year (e.g. June2021)
DATE=June2021

### MODIFY
########## 1.3 Specify first author of discovery GWAS
AUTHOR=Yengo

### MODIFY
########## 1.4 Specify your username on Vector
USER=idamat

################################################################################
#################### 2. PREPARE FOLDER STRUCTURE AND FILES #####################
################################################################################

### NO MODIFICATIONS TO THIS SECTION, RUN AS IS

### All work will be carried out on the scratch folder, and only files of relevance kept in the /nfs/AGE/IGEMS/IGEMS_Pheno_PRS folder
mkdir /scratch/tmp/$USER
mkdir /scratch/tmp/$USER/IGEMS_PRS
mkdir /scratch/tmp/$USER/IGEMS_PRS/$PHENO
# Move nohup log to work folder, if code is started elsewhere
mv nohup.out /scratch/tmp/$USER/IGEMS_PRS/$PHENO/
# NOTE! If one of these already exists, you'll get a warning message. Just ignore!
cd /scratch/tmp/$USER/IGEMS_PRS/$PHENO

# Assign main directory
HOME_DIR=/nfs/AGE/IGEMS

# File to report Ns and other info 
echo > report.log
echo "PRS:" $PHENO $AUTHOR >> report.log
echo "By:" $USER", "  $DATE >> report.log

################################################################################
######################## 3. UNZIP LD MATRIX, IF NEEDED #########################
################################################################################

# LD matrix was downloaded from GCTB SBayesR webpage:
# https://cnsgenomics.com/software/gctb/#Download (Shrunk sparse matrix)

### USE IF NEEDED
# Please note! Unzipping the LD matrix files takes time, and is not needed if 
# it's already done. Check if the files are already unzipped:
ls /scratch/tmp/$USER/IGEMS_PRS/LDmatrix/ukbEURu_hm3_shrunk_sparse

# If files are already here, no need to do anything else under step 3 - remove
# the below. The list of files should look like this:
# md5sum.txt
# README
# ukbEURu_hm3_all_v3_50k.ldm.sparse.info
# ukbEURu_hm3_chr10_v3_50k.ldm.sparse.bin
# ukbEURu_hm3_chr10_v3_50k.ldm.sparse.info
# ukbEURu_hm3_chr10_v3_50k_sparse.log
# ukbEURu_hm3_chr11_v3_50k.ldm.sparse.bin
# ukbEURu_hm3_chr11_v3_50k.ldm.sparse.info
# ukbEURu_hm3_chr11_v3_50k_sparse.log
# ...

# If you receive the following:
# ls: cannot access /scratch/tmp/$USER/IGEMS_PRS/LDmatrix: No such file or directory
# Or if no files are listed, you need to keep the below code (or run now and then
# remove from the final code, but it only needs to be done once - it takes ~10 min)

# Create a folder for the LD matrix
mkdir /scratch/tmp/$USER/IGEMS_PRS/LDmatrix
# Unzip
unzip /nfs/AGE/IGEMS/IGEMS_PRS_pipeline/SBayesR_and_RefData/ukbEURu_hm3_sparse.zip -d /scratch/tmp/$USER/IGEMS_PRS/LDmatrix/

################################################################################
###################### 4. PREPARE THE SUMMARY STATISTICS #######################
################################################################################

### Please store the original summary statistics in gzip format in the folder 
### /nfs/AGE/IGEMS/IGEMS_pheno_PRS/Data/Original_data/GWAS_sumstats/

### MODIFY
########## 3.1 Specify name of the sumstatfile (change only the file name, not the folder path) 
SUMSTATS=$HOME_DIR/IGEMS_pheno_PRS/Data/Original_data/GWAS_sumstats/Yengo2018_BMI_GIANT_UKB_metaanalysis_noSTR.tbl.gz

########## 3.2 Specify the correct columns in the sumstat file 
# To view the sumstat file:
less -S $SUMSTATS 
# To leave view:
q

### For SBayesR, we need the following columns (.ma format):
### SNP, effect allele, non-effect allele, MAF, beta, se, p-value, number of individuals
### PLEASE NOTE! All these are not always provided in the sumstats (especially the MAF and N). 
### Please leave blank here and see below for how to deal with various issues.

### MODIFY
##### 3.2.1 Assign the columns
SNP=3
EA=4
NEA=5
MAF=11
BETA=6
SE=7
PVAL=8
N=10

### Select one of the below depending on if N is reported in the sumstats (affects the SBayesR code):
# If yes:
NREP=Nrep
# If no:
#NREP=noN
# Also add the N from the GWAS paper here:
#Ngwas=94437 

### MODIFY (depending on sumstat format)
##### 3.2.2 Select columns from the the sumstats
### NOTE! One of these options are selected, depending on whether MAF and N are provided in the sumstats

########## A) If MAF is available:

### A1) If N is also available, use this code: 
zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "freq" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$MAF"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" $'"$N"'}' > sumstatfile.txt

### A2) If N is not available, use this code
#zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "freq" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
#NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$MAF"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$Ngwas"'}' > sumstatfile.txt

########## B) If MAF is NOT available:
### Get MAF from the 1kg reference
#plink --bfile /nfs/AGE/twins.psychchip.data/GRS/Data/DerivedData/ldref/1kgEUR --freq --out 1kG_freq
#awk 'BEGIN {print "RSID" "\t" "MA" "\t" "freq"} NR!=1 {print $2 "\t" $3 "\t" $5}' 1kG_freq.frq > 1kg_MAF.txt

### Then select one of the following, depending on if N is reported in the sumstats
### B1) If N is available, use this code: 
#zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
#NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$N"'}' > sumstatfile_pre.txt
#awk 'NR==FNR {a[$1] = $0; next} ($1) in a {print $0 "\t" a[$1]}' 1kg_MAF.txt sumstatfile_pre.txt > sumstatfile_pre2.txt
#awk 'NR==1 {print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
#NR!=1 {if ($2==$9) print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
#{if ($3==$9) print $1 "\t" $2 "\t" $3 "\t" 1-$10 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' sumstatfile_pre2.txt > sumstatfile.txt

# Report number of SNPs in full GWAS
#echo 'N SNPs in full sumstats, prior to merging with 1kg for MAF' >> report.log
#wc -l < sumstatfile_pre.txt >> report.log

### B2) If N is not available, use this code
#zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$Ngwas"'}' > sumstatfile_pre.txt
#awk 'NR==FNR {a[$1] = $0; next} ($1) in a {print $0 "\t" a[$1]}' 1kg_MAF.txt sumstatfile_pre.txt > sumstatfile_pre2.txt
#awk 'NR==1 {print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
NR!=1 {if ($2==$9) print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
{if ($3==$9) print $1 "\t" $2 "\t" $3 "\t" 1-$10 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' sumstatfile_pre2.txt > sumstatfile.txt

# Report number of SNPs in full GWAS
#echo 'N SNPs in full sumstats, prior to merging with 1kg for MAF' >> report.log
#wc -l < sumstatfile_pre.txt >> report.log


### Check to see it looks as expected
head sumstatfile.txt
# Should have the following format
# RSID    		A1      A2      freq    b       se      p       N
# rs1000000     A       G       0.2219  -0.0012 0.0024  0.6243  579827
# rs10000010    T       C       0.5086  -0.0002 0.0019  0.9318  668280
# rs1000002     T       C       0.4884  -0.0054 0.0020  0.0065  582349
# rs10000023    T       G       0.5817  -0.0031 0.0020  0.1212  582277
# rs1000003     A       G       0.8404  0.0020  0.0028  0.4792  582269
# ...
### NOTE! If e.g lower case letters are used for A1 and A2, or if OR is reported instead of beta, see below for modifications

### NOTE! Depending on the format of the summary statistics, some of the options below may be needed 
##### NOTE! If OR is reported rather than beta, add the following:
# mv sumstatfile.txt sumstat_pre.txt
# awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, log($5), $6, $7, $8}' sumstat_pre.txt > sumstatfile.txt
# head sumstatfile.txt

##### NOTE! If alleles are not in capital letters, add the following:
# Convert alleles to capital letters
# mv sumstatfile.txt sumstat_pre.txt
# awk 'BEGIN{OFS="\t"}{print $1, toupper($2), toupper($3), $4, $5, $6, $7, $8}' sumstat_pre.txt > sumstatfile.txt
# head sumstatfile.txt

################################################################################
############### MODIFICATIONS DONE! Make no changes to the below
##### Remove code not used above (left commented out, optional) and save as specified in the header 
##### Then run as is with the nohup bash command as specified in the header
################################################################################

# Number of SNPs
echo 'N SNPs sumstats from start' >> report.log
wc -l < sumstatfile.txt >> report.log

# Filter out SNPs with alleles different from those in the target data
awk '{print $0 "\t" $2$3 "\t" $3$2}' sumstatfile.txt > sumstatfile2.txt
awk 'NR==FNR {a[$4] = $3; next} ($1) in a {print $0 "\t" a[$1]}' $HOME_DIR/IGEMS_SNPs/SNPs_IGEMS_selected_20210228.txt sumstatfile2.txt > sumstatfile3.txt
awk 'NR==1 {print $0} NR!=1 {if ($9==$11 || $10==$11) print $0}' sumstatfile3.txt > sumstatfile4.txt

# Number of SNPs left after merging with IGEMS SNPs
echo 'N SNPs after merging with IGEMS SNPs and dropping mismatches' >> report.log
wc -l < sumstatfile4.txt >> report.log

# Drop the added columns
awk '{$9="";$10="";$11="";print $0}' OFS='\t' sumstatfile4.txt > sumstatfile5.txt

################################################################################
##### Call the SBayesR script

bash $HOME_DIR/IGEMS_PRS_pipeline/Programs/SBayesR_code_${NREP}.bash
wait

# Concatenate all files
grep 'Name' SBayesR_chr_22.snpRes > header.txt
cat SBayesR_chr_*.snpRes | grep -v -w 'Name' > tmp.txt 
cat header.txt tmp.txt > $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes 
rm tmp.txt

# Number of SNPs in the final sumstatfile
echo 'N SNPs after SBayesR' >> report.log
wc -l < $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes >> report.log

################################################################################
######################## 5. CREATE THE PRS IN PLINK ############################
################################################################################
# Generate file with RSID, effect allele, and beta for scoring
awk '{print $2 "\t" $5 "\t" $8}' $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes > SNP_effects.txt

##### Run Plink to score the PRS, one per target set
# GOSH (NOTE! Duplicates were removed when merging the files)
plink2 \
--pfile $HOME_DIR/GOSH_data_IGEMS_SNPs/GOSH_data_IGEMS_SNPs_20210630 \
--score SNP_effects.txt cols=+scoresums \
--out prs_GOSH

echo 'N check GOSH sample (should be 2214 including the header)' >> report.log
wc -l < prs_GOSH.sscore >> report.log
echo 'N variants processed in the GOSH sample' >> report.log
cat prs_GOSH.log | grep "variants processed." | awk '{print $2}' >> report.log

# TwinGene
plink2 \
--pfile /nfs/projects/STR_GWAS/twingene/TwinGene_HRC_imputation/03_imputed/TwinGene_HRC_imputation \
--rm-dup exclude-mismatch \
--score SNP_effects.txt cols=+scoresums \
--out prs_TG

echo 'N check TwinGene sample (should be 10912 including the header)' >> report.log
wc -l < prs_TG.sscore >> report.log
echo 'N variants processed in the TwinGene sample' >> report.log
cat prs_TG.log | grep "variants processed." | awk '{print $2}' >> report.log

# SALT-Y
plink2 \
--pfile /nfs/projects/STR_GWAS/saltygwas/Data/DerivedData/HRC_imputation_RicopiliQCed_LuYi/SALTY_HRC_imputation.addMZ_PLINK2/SALTY_HRC_imputation.addMZ \
--rm-dup exclude-mismatch \
--score SNP_effects.txt cols=+scoresums \
--out prs_SALTY

echo 'N check SALT-Y sample (should be 6543 including the header)' >> report.log
wc -l < prs_SALTY.sscore >> report.log
echo 'N variants processed in the SALT-Y sample' >> report.log
cat prs_SALTY.log | grep "variants processed." | awk '{print $2}' >> report.log

################################################################################
###################### 7. SAVE THE PRS FILES AND LOGS ##########################
################################################################################

# Copy the PRS files
cp prs_GOSH.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_GOSH_${DATE}.txt
cp prs_TG.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_TG_${DATE}.txt
cp prs_SALTY.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_SALTY_${DATE}.txt

# Copy the report
cp report.log $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_report_${DATE}.log
# Save the main log
cp nohup.out $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_log_${DATE}.out

# Check for errors
# File for error log
echo > $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
echo $PHENO $AUTHOR $DATE >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
grep rror* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
grep fatal nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
grep arning* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
grep ARNING* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log
grep cannot* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_errorcheck_${DATE}.log

# Gzip the sumstats
gzip $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes

################################################################################
############################### END OF FILE ####################################
################################################################################

