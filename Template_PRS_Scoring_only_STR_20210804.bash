################################################################################
# NAME: Template_PRS_Scoring_only_STR_20210804.bash
# PURPOSE: Master code for creation of PRS scores in the Swedish Twin Registry data, 
# when SBayesR has already been done
#
# CREATED: by Ida Karlsson (IK) 20210804
# UPDATED: 	
#
# NOTE!! Code should be updated with sample and sumstat info, and saved in the folder
# /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Programs under the name *PHENOTYPE*_*Author*_PRS_Scoring_only_STR_*DATE*.bash
# The data will by default be saved in the /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Programs folder,
# and the log in the /nfs/AGE/IGEMS/IGEMS_Pheno_PRS/Logs folder, using the same naming convention
#
# NOTE! SBayesR prepared sumstats should be available in the folder
# /nfs/AGE/IGEMS/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/
# Named as: SBayesR_sumstats_PHENO_AUTHOR_DATE.snpRes
#
# To run: nohup bash /nfs/AGE/IGEMS/IGEMS_pheno_PRS/Programs/*PHENOTYPE*_*Author*_PRS_Scoring_only_STR_*DATE*.bash &
#
################################################################################
######################## 1. SPECIFY PHENOTYPE AND DATE #########################
################################################################################

### MODIFY
########## 1.1 Specify phenotype (NOTE! Must be same as name in the SBayesR sumstat file)
PHENO=BMI

### MODIFY
########## 1.2.1 Specify today's month and year (e.g. Jan2021)
DATE=June2021
########## 1.2.2 Specify month and year for SBayesR sumstat name (e.g. Jan2021)
DATE_SBayesR=June2021

### MODIFY
########## 1.3 Specify first author of discovery GWAS
AUTHOR=Yengo

### MODIFY
########## 1.4 Specify your username on Vector
USER=idamat

### NO MODIFICATIONS AFTER THIS SECTION, RUN AS IS

################################################################################
#################### 2. PREPARE FOLDER STRUCTURE AND FILES #####################
################################################################################

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
echo "NOTE! SBayesR already done, this code only does the scoring" >> report.log
echo "" >> report.log

################################################################################

# Number of SNPs in the final sumstatfile, after SBayesR
echo 'N SNPs after SBayesR' >> report.log
wc -l < $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE_SBayesR}.snpRes >> report.log

################################################################################
######################## 3. CREATE THE PRS IN PLINK ############################
################################################################################
# Generate file with RSID, effect allele, and beta for scoring
awk '{print $2 "\t" $5 "\t" $8}' $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/SBayesR_sumstats/SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE_SBayesR}.snpRes > SNP_effects.txt

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
###################### 4. SAVE THE PRS FILES AND LOGS ##########################
################################################################################

# Copy the PRS files
cp prs_GOSH.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_GOSH_${DATE}.txt
cp prs_TG.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_TG_${DATE}.txt
cp prs_SALTY.sscore $HOME_DIR/IGEMS_pheno_PRS/Data/Derived_data/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_STR_SALTY_${DATE}.txt

# Copy the report
cp report.log $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_report_${DATE}.log
# Save the main log
cp nohup.out $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_log_${DATE}.out

# Check for errors
# File for error log
echo > $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
echo $PHENO $AUTHOR $DATE >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
grep rror* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
grep fatal nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
grep arning* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
grep ARNING* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log
grep cannot* nohup.out >> $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_Scoring_only_STR_errorcheck_${DATE}.log

################################################################################
############################### END OF FILE ####################################
################################################################################

