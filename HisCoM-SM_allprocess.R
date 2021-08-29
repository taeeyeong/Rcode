####################################################################################
########################### HisCoM-SM - GMS ###################################
####################### (1) single SNP association ############################
####################################################################################

for (var in names(met2[,3:66])){
  pheno <- data.frame(met2$DIST_ID,met2$DIST_ID, met2[,var])
  colnames(pheno) <- c("FID","IID",var)
  write.table(pheno, paste0("/home2/tyjung/mQTL/Data/Pheno_",var,"_phase6.txt"), row.names = F, quote = F)
}

##Plink code for linear regression (Metabolite~SNP+Age+Sex)
library(multicore)
plink_imp <- function(i){
  var <- names(met2)[i]
  plink_code <- paste0("plink --bfile /home2/tyjung/DeepHisCoM/snp/KARE_Affy5_totalqc --linear --pheno /home2/tyjung/mQTL/Data/Pheno_",var,"_phase6.txt --pheno-name ",var," --covar /home2/tyjung/mQTL/Data/covdata_6.txt --covar-name AGE,SEX --out /home2/tyjung/mQTL/Result/SNP_",var,"_associ_6 --noweb")
  system(plink_code)
}
multi_plinkimp <-mclapply(c(3:18), plink_imp, mc.cores =16) 
multi_plinkimp <-mclapply(c(19:34), plink_imp, mc.cores = 16) 
multi_plinkimp <-mclapply(c(35:50), plink_imp, mc.cores = 16) 
multi_plinkimp <-mclapply(c(51:66), plink_imp, mc.cores = 16) 

## (R) 대사체~SNP association file안에 covariate 효과제거
annot <- read.table("/home2/tyjung/mQTL/HisCoM/annot_PMS_6.set", header = F, stringsAsFactors = F)
head(annot)
dim(annot) # 560 2
length(unique(annot$V1)) #144
length(unique(annot$V2)) #53
names(annot) <- c("group","var")

for(var in unique(annot$var)){
  assoc <- read.table(paste0("/home2/tyjung/mQTL/Result/SNP_",var,"_associ_6.assoc.linear"), header = T, stringsAsFactors = F)
  assoc2 <- assoc[which(assoc$TEST=="ADD"),]
  write.table(assoc2,paste0("/home2/tyjung/mQTL/Clumping/SNP_",var,"_associ_6_only_ADD.assoc.linear"), row.names = F, quote = F)
}


## (Plink2) Clumping 
# Clumped file : containing index SNPs after performing clumping
library(multicore)
plink_imp <- function(i){
  var <- unique(annot$var)[i]
  plink_code <- paste0("/data/bin/plink2 --bfile /home2/tyjung/mQTL/plink/KARE_Affy5_totalqc_627samp --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump /home2/tyjung/mQTL/Clumping/SNP_",var,"_associ_6_only_ADD.assoc.linear --clump-snp-field SNP --clump-field P --out /home2/tyjung/mQTL/Clumping/SNP_",var," --noweb")
  system(plink_code)
}
multi_plinkimp <- mclapply(c(1:27), plink_imp, mc.cores =27)
multi_plinkimp <- mclapply(c(28:53), plink_imp, mc.cores = 26)



library(multicore)
awk_imp <- function(i){
  var <- unique(annot$var)[i]
  # valid.snp: Extract the index SNP
  awk_code1 <- paste0("awk 'NR!=1{print $3}' /home2/tyjung/mQTL/Clumping/SNP_",var,".clumped > /home2/tyjung/mQTL/Clumping/SNP_",var,".valid.snp")
  system(awk_code1)
  # SNP.pvalue: a file containing SNP id and corresponding p-values
  awk_code2 <- paste0("awk '{print $2,$9}' /home2/tyjung/mQTL/Clumping/SNP_",var,"_associ_6_only_ADD.assoc.linear > /home2/tyjung/mQTL/Clumping/",var, "_SNP.pvalue")
  system(awk_code2)
}
multi_awk <- mclapply(c(1:53), awk_imp, mc.cores = 53)

# echo "0.001 0 0.001" >> /home2/tyjung/mQTL/range_list


# profile file : generated PRS file
library(multicore)
plink_imp <- function(i){
  var <- unique(annot$var)[i]
  plink_code <- paste0("/data/bin/plink2 --bfile /home2/tyjung/mQTL/plink/KARE_Affy5_totalqc_627samp --score /home2/tyjung/mQTL/Clumping/SNP_",var,"_associ_6_only_ADD.assoc.linear 2 4 7 header --q-score-range /home2/tyjung/mQTL/range_list /home2/tyjung/mQTL/Clumping/",var,"_SNP.pvalue --extract /home2/tyjung/mQTL/Clumping/SNP_",var,".valid.snp --out /home2/tyjung/mQTL/Clumping/",var," --noweb")
  system(plink_code)
}
multi_plinkimp <- mclapply(c(1:27), plink_imp, mc.cores = 27)
multi_plinkimp <- mclapply(c(28:53), plink_imp, mc.cores = 26)


## datasets for HisCoM 
train_clump <- data.frame(DIST_ID = met3$DIST_ID)
for(var in unique(annot$var)){
  tmp <- read.table(paste0("/home/tyjung/mQTL/PRS_Clumping/",var,".0.001.profile"), header = T, stringsAsFactors = F)
  train_clump[,var] <- tmp$SCORE[match(train_clump$DIST_ID, tmp$FID)]
}

train_clump$PAT <- 0
train_clump$MAT <- 0
train_clump$SEX <- met3$sex[match(train_clump$DIST_ID, met3$DIST_ID)]
train_clump$PHENOTYPE <- met3$status[match(train_clump$DIST_ID, met3$DIST_ID)]
train_clump$PHENOTYPE <- ifelse(train_clump$PHENOTYPE==0,1,2)
table(train_clump$PHENOTYPE)
train_clump2 <- train_clump[,c("DIST_ID","DIST_ID","PAT","MAT","SEX","PHENOTYPE",unique(annot$var))]
names(train_clump2) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE",unique(annot$var))

write.table(train_clump2, "/home/tyjung/mQTL/HisCoM/train_clump_0.001_6.raw", row.names = F, quote = F)


##Covariate file for HisCoM
one_path <- names(which(table(annot$group)==1))
one_var <- annot[which(annot$group %in% one_path),"var"]
cov_pheno <- train_clump2[,c("FID","IID","PHENOTYPE", one_var)]
names(cov_pheno) <- c("FID","IID","PHENOTYPE",one_path)

write.table(cov_pheno,"/home/tyjung/mQTL/HisCoM/ID_pheno_clump_0.001_6.txt", row.names = F, quote = F)


#(Wizard) Clumping data
wisard --expression /home2/tyjung/mQTL/HisCoM/train_clump_0.001_6.raw --geneset /home2/tyjung/mQTL/HisCoM/annot_PMS_6_homosapiens.set --sampvar /home2/tyjung/mQTL/HisCoM/ID_pheno_clump_0.001_6.txt --pharaoh --cname ADRENERGIC_SIGNALING_IN_CARDIOMYOCYTES,AMOEBIASIS,ANTIFOLATE_RESISTANCE,CHAGAS_DISEASE,CIRCADIAN_ENTRAINMENT,CITRATE_CYCLE_TCA_CYCLE,D_ARGININE_AND_D_ORNITHINE_METABOLISM,DOPAMINERGIC_SYNAPSE,ESTROGEN_SIGNALING_PATHWAY,FOXO_SIGNALING_PATHWAY,GASTRIC_ACID_SECRETION,GLUCAGON_SIGNALING_PATHWAY,HUNTINGTON_DISEASE,INSULIN_RESISTANCE,INSULIN_SECRETION,LONG_TERM_DEPRESSION,LONG_TERM_POTENTIATION,MELANOGENESIS,MORPHINE_ADDICTION,NEOMYCIN_KANAMYCIN_AND_GENTAMICIN_BIOSYNTHESIS,OXIDATIVE_PHOSPHORYLATION,PANCREATIC_SECRETION,PARKINSON_DISEASE,PHOSPHOLIPASE_D_SIGNALING_PATHWAY,PHOSPHONATE_AND_PHOSPHINATE_METABOLISM,PROLACTIN_SIGNALING_PATHWAY,REGULATION_OF_ACTIN_CYTOSKELETON,REGULATION_OF_LIPOLYSIS_IN_ADIPOCYTES,RENIN_SECRETION,RIBOFLAVIN_METABOLISM,SELENOCOMPOUND_METABOLISM,SEROTONERGIC_SYNAPSE,SPHINGOLIPID_METABOLISM,SPHINGOLIPID_SIGNALING_PATHWAY,SYNTHESIS_AND_DEGRADATION_OF_KETONE_BODIES,UBIQUINONE_AND_OTHER_TERPENOID_QUINONE_BIOSYNTHESIS,VASCULAR_SMOOTH_MUSCLE_CONTRACTION --pname PHENOTYPE --nperm 100000 --verbose --cv 5 --thread 20 --out /home2/tyjung/mQTL/HisCoM/Result/210722_Clump_0.001/Res_clump_0.001_6 --prorange [1.9608e-100,1000000000]



####################################################################################
#################### HisCoM-SM - GMS calculation ############################
########################### (2) GBLUP ######################################
####################################################################################

map_id <- read.csv("/home2/tyjung/mQTL/Data/Match_Result.csv", header = T, stringsAsFactors = F)
met <- read.csv("/home2/tyjung/mQTL/Data/serrf6_total.csv", header = T, stringsAsFactors = F)
met$DIST_ID <- map_id$DIST_ID[match(met$rand_num,map_id$id)]
sum(is.na(met$DIST_ID)) # 38
met2 <- met[-which(is.na(met$DIST_ID)),]
dim(met2) #653 72
table(met2$status)
met2$status <- gsub(2,1,met2$status) #PreT2D -> T2D
table(met2$status)

# Make GRM file
/home2/tyjung/gcta_1.93.2beta/gcta64 --bfile /home2/tyjung/mQTL/plink/KARE_Affy5_totalqc_627samp --make-grm --out /home2/tyjung/mQTL/GCTA/KARE_Affy5_totalqc --thread-num 10

# To obtain BLUP solutions for the genetic values of individuals

library(multicore)
gcta_imp <- function(i){
  var <- names(met2)[i]
  gcta_code <- paste0("/home2/tyjung/gcta_1.93.2beta/gcta64 --reml --grm /home2/tyjung/mQTL/GCTA/KARE_Affy5_totalqc --pheno /home2/tyjung/mQTL/Data/Pheno_",var,"_phase6.txt --reml-pred-rand --qcovar /home2/tyjung/mQTL/Data/covdata_6.txt --out /home2/tyjung/mQTL/GCTA/",var)
  system(gcta_code)
}
multi_gctaimp <-mclapply(c(3:18), gcta_imp, mc.cores =16) 
multi_gctaimp <-mclapply(c(19:34), gcta_imp, mc.cores = 16) 
multi_gctaimp <-mclapply(c(35:50), gcta_imp, mc.cores = 16) 
multi_gctaimp <-mclapply(c(51:66), gcta_imp, mc.cores = 16) 


# To obtain BLUP solutions for the SNP effects
library(multicore)
gcta_imp <- function(i){
  var <- names(met2)[i]
  gcta_code <- paste0("/home2/tyjung/gcta_1.93.2beta/gcta64 --bfile /home2/tyjung/mQTL/plink/KARE_Affy5_totalqc_627samp --autosome --blup-snp /home2/tyjung/mQTL/GCTA/",var,".indi.blp --out /home2/tyjung/mQTL/GCTA/",var)
  system(gcta_code)
}
multi_gctaimp <-mclapply(c(3:18), gcta_imp, mc.cores =16) 
multi_gctaimp <-mclapply(c(19:34), gcta_imp, mc.cores = 16) 
multi_gctaimp <-mclapply(c(35:50), gcta_imp, mc.cores = 16) 
multi_gctaimp <-mclapply(c(51:66), gcta_imp, mc.cores = 16) 

##annotation file
annot <- read.table("/home2/tyjung/mQTL/HisCoM/annot_PMS_6.set", header = F, stringsAsFactors = F)
head(annot)
dim(annot) # 560 2
length(unique(annot$V1)) #144
length(unique(annot$V2)) #53
names(annot) <- c("group","var")


for(var in unique(annot$var)){
  assoc <- read.table(paste0("/home2/tyjung/mQTL/GCTA/",var,".snp.blp"), header = F, stringsAsFactors = F)
  assoc2 <- assoc[order(abs(assoc$V3), decreasing = T),]
  assoc2 <- assoc2[c(1:(nrow(assoc2)*0.2)),]
  write.table(assoc2,paste0("/home2/tyjung/mQTL/GBLUP_0.2_abs/SNP_",var,"_0.2_abs.snp.blp"), row.names = F, quote = F)
}


library(multicore)
plink_imp <- function(i){
  var <- unique(annot$var)[i]
  plink_code <- paste0("plink --bfile /home2/tyjung/mQTL/plink/KARE_Affy5_totalqc_627samp --score /home2/tyjung/mQTL/GBLUP_0.2_abs/SNP_",var,"_0.2_abs.snp.blp --out /home2/tyjung/mQTL/GBLUP_0.2_abs/",var," --noweb")
  system(plink_code)
}
multi_plinkimp <-mclapply(c(1:18), plink_imp, mc.cores =18) 
multi_plinkimp <-mclapply(c(19:36), plink_imp, mc.cores = 18) 
multi_plinkimp <-mclapply(c(37:53), plink_imp, mc.cores = 17)

## Make dataset for HisCoM 
library(snpStats)
kare_snp <- read.plink(bed = "/home/tyjung/data_snp/KARE_Affy5_totalqc.bed",
                       bim = "/home/tyjung/data_snp/KARE_Affy5_totalqc.bim",
                       fam = "/home/tyjung/data_snp/KARE_Affy5_totalqc.fam")

map_id <- read.csv("/home/tyjung/kare_id_matching/Match_Result.csv", header = T, stringsAsFactors = F)
met <- read.csv("/home/tyjung/kare_id_matching/serrf6_total.csv", header = T, stringsAsFactors = F)
met$DIST_ID <- map_id$DIST_ID[match(met$rand_num,map_id$id)]
sum(is.na(met$DIST_ID)) # 38
met2 <- met[-which(is.na(met$DIST_ID)),]
dim(met2) #653 72
table(met2$status)
met2$status <- gsub(2,1,met2$status) #PreT2D -> T2D
table(met2$status)
met3 <- met2[which(met2$DIST_ID %in% kare_snp$fam$pedigree),]
dim(met3) #627 72


train_gblup <- data.frame(DIST_ID = met3$DIST_ID)
for(var in unique(annot$var)){
  tmp <- read.table(paste0("/home/tyjung/mQTL/GBLUP_0.2_abs/",var,".profile"), header = T, stringsAsFactors = F)
  cat("mean of ",var,": ",mean(tmp$SCORE),"\n")
  train_gblup[,var] <- tmp$SCORE[match(train_gblup$DIST_ID, tmp$FID)]
  
}

train_gblup$PAT <- 0
train_gblup$MAT <- 0
train_gblup$SEX <- met3$sex[match(train_gblup$DIST_ID, met3$DIST_ID)]
train_gblup$PHENOTYPE <- met3$status[match(train_gblup$DIST_ID, met3$DIST_ID)]
train_gblup$PHENOTYPE <- ifelse(train_gblup$PHENOTYPE==0,1,2)
table(train_gblup$PHENOTYPE)
train_gblup2 <- train_gblup[,c("DIST_ID","DIST_ID","PAT","MAT","SEX","PHENOTYPE",unique(annot$var))]
colnames(train_gblup2) <-  c("FID","IID","PAT","MAT","SEX","PHENOTYPE",unique(annot$var))

write.table(train_gblup2, "/home/tyjung/mQTL/HisCoM/train_gblup_0.2_abs_6.raw", row.names = F, quote = F)



##Covariate file for HisCoM
one_path <- names(which(table(annot$group)==1))
one_var <- annot[which(annot$group %in% one_path),"var"]
cov_pheno <- train_gblup2[,c("FID","IID","PHENOTYPE", one_var)]
names(cov_pheno) <- c("FID","IID","PHENOTYPE",one_path)

write.table(cov_pheno,"/home/tyjung/mQTL/HisCoM/ID_pheno_gblup_0.2_abs_6.txt", row.names = F, quote = F)


#gcta data
wisard --expression /home2/tyjung/mQTL/HisCoM/train_gblup_0.2_abs_6.raw --geneset /home2/tyjung/mQTL/HisCoM/annot_PMS_6_homosapiens.set --sampvar /home2/tyjung/mQTL/HisCoM/ID_pheno_gblup_0.2_abs_6.txt --pharaoh --cname ADRENERGIC_SIGNALING_IN_CARDIOMYOCYTES,AMOEBIASIS,ANTIFOLATE_RESISTANCE,CHAGAS_DISEASE,CIRCADIAN_ENTRAINMENT,CITRATE_CYCLE_TCA_CYCLE,D_ARGININE_AND_D_ORNITHINE_METABOLISM,DOPAMINERGIC_SYNAPSE,ESTROGEN_SIGNALING_PATHWAY,FOXO_SIGNALING_PATHWAY,GASTRIC_ACID_SECRETION,GLUCAGON_SIGNALING_PATHWAY,HUNTINGTON_DISEASE,INSULIN_RESISTANCE,INSULIN_SECRETION,LONG_TERM_DEPRESSION,LONG_TERM_POTENTIATION,MELANOGENESIS,MORPHINE_ADDICTION,NEOMYCIN_KANAMYCIN_AND_GENTAMICIN_BIOSYNTHESIS,OXIDATIVE_PHOSPHORYLATION,PANCREATIC_SECRETION,PARKINSON_DISEASE,PHOSPHOLIPASE_D_SIGNALING_PATHWAY,PHOSPHONATE_AND_PHOSPHINATE_METABOLISM,PROLACTIN_SIGNALING_PATHWAY,REGULATION_OF_ACTIN_CYTOSKELETON,REGULATION_OF_LIPOLYSIS_IN_ADIPOCYTES,RENIN_SECRETION,RIBOFLAVIN_METABOLISM,SELENOCOMPOUND_METABOLISM,SEROTONERGIC_SYNAPSE,SPHINGOLIPID_METABOLISM,SPHINGOLIPID_SIGNALING_PATHWAY,SYNTHESIS_AND_DEGRADATION_OF_KETONE_BODIES,UBIQUINONE_AND_OTHER_TERPENOID_QUINONE_BIOSYNTHESIS,VASCULAR_SMOOTH_MUSCLE_CONTRACTION --pname PHENOTYPE --nperm 100000 --verbose --cv 5 --thread 20 --out /home2/tyjung/mQTL/HisCoM/Result/210722_gblup_0.2_abs/Res_gblup_0.2_abs_6 --prorange [1.9608e-120,1000000]


