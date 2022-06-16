order <-  c("participant_id","phenotype","rare_diseases_family_id","SYMBOL","Consequence","participant_type",
            "assembly","genotype","mode_of_inheritance", "segregation_pattern","penetrance","tier",
            "father_affected","mother_affected","full_brothers_affected","full_sisters_affected",
            "participant_phenotypic_sex","X.Uploaded_variation","Location","Allele","IMPACT","Gene","Feature_type",
            "Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
            "Amino_acids","Codons","Existing_variation","STRAND","SYMBOL_SOURCE","HGNC_ID","CANONICAL",
            "MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SIFT","PolyPhen","BLOSUM62",
            "LoFtool","CADD_PHRED","CADD_RAW","SpliceAI_pred_DP_AG","SpliceAI_pred_DP_AL","SpliceAI_pred_DP_DG",
            "SpliceAI_pred_DP_DL","SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL","SpliceAI_pred_SYMBOL","CLIN_SIG","SOMATIC","PHENO","PUBMED","HGVS_OFFSET","AF",
            "AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF",
            "gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF", "gnomAD_OTH_AF",
            "gnomAD_SAS_AF","DISTANCE","FLAGS","DOMAINS","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS",
            "MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","interpretation_cohort_id","interpretation_request_id",
            "sample_id","so_term")


##################
# BIALLELIC -HSP #
#################

unsolved_comphet_lof_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch37/unsolved_comphet_lof_hsp_grch37.csv")
unsolved_comphet_missense_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch37/unsolved_comphet_missense_hsp_grch37.csv")
unsolved_recessive_lof_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch37/unsolved_recessive_lof_hsp_grch37.csv")
unsolved_recessive_LOFmisense_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch37/unsolved_recessive_LOFmissense_hsp_grch37.csv")
unsolved_recessive_missense_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch37/unsolved_recessive_missense_hsp_grch37.csv")

unsolved_comphet_lof_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/comphet_lof_hsp_grch38.csv")
unsolved_comphet_missense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/comphet_missense_hsp_grch38.csv")
unsolved_comphet_LOFmissense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/comphet_LOFmissense_hsp_grch38.csv")
unsolved_recessive_lof_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/recessive_lof_hsp_grch38.csv")
unsolved_recessive_LOFmisense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/recessive_LOFmissense_hsp_grch38.csv")
unsolved_recessive_missense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/biallelic/Grch38/recessive_missense_hsp_grch38.csv")


#Rbind
comphet_lof_hsp <- rbind(unsolved_comphet_lof_hsp_grch37,unsolved_comphet_lof_hsp_grch38)
comphet_lof_hsp <- comphet_lof_hsp[,order]

comphet_missense_hsp <- rbind(unsolved_comphet_missense_hsp_grch37,unsolved_comphet_missense_hsp_grch38)
comphet_missense_hsp <- comphet_missense_hsp[,order]

comphet_LOFmissense_hsp <- unsolved_comphet_LOFmissense_hsp_grch38[,order]

recessive_lof_hsp <- rbind(unsolved_recessive_lof_hsp_grch37, unsolved_recessive_lof_hsp_grch38)
recessive_lof_hsp <- recessive_lof_hsp[,order]

recessive_LOFmissense_hsp <- rbind(unsolved_recessive_LOFmisense_hsp_grch37, unsolved_recessive_LOFmisense_hsp_grch38)
recessive_LOFmissense_hsp <- recessive_LOFmissense_hsp[,order]

recessive_missense_hsp <- rbind(unsolved_recessive_missense_hsp_grch37, unsolved_recessive_missense_hsp_grch38)
recessive_missense_hsp <- recessive_missense_hsp[,order]


###############################
# MONOALLELIC AND DENOVO -HSP #
###############################

unsolved_denovo_lof_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/unsolved_denovo_lof_hsp_grch37.csv")
unsolved_denovo_LOFmissense_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/unsolved_denovo_LOFmissense_hsp_grch37.csv")
unsolved_denovo_missense_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/unsolved_denovo_missense_hsp_grch37.csv")
LOF_monoallelic_father_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/LOF_monoallelic_father_hsp_grch37.csv")
LOF_monoallelic_mother_hsp_grch37 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/LOF_monoallelic_mother_hsp_grch37.csv")
missense_monoallelic_mother_hsp_grch37<- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch37/missense_monoallelic_mother_hsp_grch37.csv")

unsolved_denovo_lof_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/denovo_lof_hsp_grch38.csv")
unsolved_denovo_LOFmissense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/denovo_LOFmissense_hsp_grch38.csv")
unsolved_denovo_missense_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/denovo_missense_hsp_grch38.csv")
LOF_monoallelic_father_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/LOF_monoallelic_father_hsp_grch38.csv")
LOF_monoallelic_mother_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/LOF_monoallelic_mother_hsp_grch38.csv")
missense_monoallelic_father_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/missense_monoallelic_father_hsp_grch38.csv")
missense_monoallelic_mother_hsp_grch38 <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP_output/hsp/monoallelic_denovo/Grch38/missense_monoallelic_mother_hsp_grch38.csv")


#Rbind

denovo_lof_hsp <- rbind(unsolved_denovo_lof_hsp_grch37, unsolved_denovo_lof_hsp_grch38)
denovo_lof_hsp <- denovo_lof_hsp[,order]

denovo_LOFmissense_hsp <- rbind(unsolved_denovo_LOFmissense_hsp_grch37, unsolved_denovo_LOFmissense_hsp_grch38)
denovo_LOFmissense_hsp <- denovo_LOFmissense_hsp[,order]

denovo_missense_hsp <- rbind(unsolved_denovo_missense_hsp_grch37, unsolved_denovo_missense_hsp_grch38)
denovo_missense_hsp <- denovo_missense_hsp[,order]

LOF_monoallelic_father_hsp <- rbind(LOF_monoallelic_father_hsp_grch37, LOF_monoallelic_father_hsp_grch38)
LOF_monoallelic_father_hsp <- LOF_monoallelic_father_hsp[,order]

LOF_monoallelic_mother_hsp <- rbind(LOF_monoallelic_mother_hsp_grch37, LOF_monoallelic_mother_hsp_grch38)
LOF_monoallelic_mother_hsp <- LOF_monoallelic_mother_hsp[,order]

missense_monoallelic_father_hsp <- missense_monoallelic_father_hsp_grch38[,order]

missense_monoallelic_mother_hsp <- rbind(missense_monoallelic_mother_hsp_grch37, missense_monoallelic_mother_hsp_grch38)
missense_monoallelic_mother_hsp <- missense_monoallelic_mother_hsp[,order]

pli <- read.csv("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/pLI_scores.csv")

names(pli)[names(pli) == 'gene'] <- "SYMBOL"

#######

order_2 <-  c("participant_id","phenotype","rare_diseases_family_id","SYMBOL","Consequence","participant_type",
              "assembly","genotype","mode_of_inheritance", "segregation_pattern","penetrance","tier",
              "father_affected","mother_affected","full_brothers_affected","full_sisters_affected",
              "participant_phenotypic_sex","X.Uploaded_variation","Location","Allele","IMPACT","Gene","Feature_type",
              "Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
              "Amino_acids","Codons","Existing_variation","STRAND","SYMBOL_SOURCE","HGNC_ID","CANONICAL",
              "MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SIFT","PolyPhen","BLOSUM62",
              "LoFtool","CADD_PHRED","CADD_RAW","pLI","SpliceAI_pred_DP_AG","SpliceAI_pred_DP_AL","SpliceAI_pred_DP_DG",
              "SpliceAI_pred_DP_DL","SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG",
              "SpliceAI_pred_DS_DL","SpliceAI_pred_SYMBOL","CLIN_SIG","SOMATIC","PHENO","PUBMED","HGVS_OFFSET","AF",
              "AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF",
              "gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF", "gnomAD_OTH_AF",
              "gnomAD_SAS_AF","DISTANCE","FLAGS","DOMAINS","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS",
              "MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","interpretation_cohort_id","interpretation_request_id",
              "sample_id","so_term")


library(dplyr)
comphet_lof_hsp <- left_join(comphet_lof_hsp, pli)[,order_2]
write.csv(comphet_lof_hsp , "comphet_lof_hsp.csv")
comphet_missense_hsp <- left_join(comphet_missense_hsp, pli)[,order_2]
write.csv(comphet_missense_hsp , "comphet_missense_hsp.csv")
comphet_LOFmissense_hsp <- left_join(comphet_LOFmissense_hsp, pli)[,order_2]
write.csv(comphet_LOFmissense_hsp , "comphet_LOFmissense_hsp.csv")

denovo_lof_hsp <- left_join(denovo_lof_hsp, pli)[,order_2]
write.csv(denovo_lof_hsp, "denovo_lof_hsp.csv")
denovo_LOFmissense_hsp <- left_join(denovo_LOFmissense_hsp, pli)[,order_2]
write.csv(denovo_LOFmissense_hsp, "denovo_LOFmissense_hsp.csv")
denovo_missense_hsp <- left_join(denovo_missense_hsp, pli)[,order_2]
write.csv(denovo_missense_hsp, "denovo_missense_hsp.csv")

recessive_lof_hsp <- left_join(recessive_lof_hsp, pli)[,order_2]
write.csv(recessive_lof_hsp, "recessive_lof_hsp.csv")
recessive_LOFmissense_hsp <- left_join(recessive_LOFmissense_hsp, pli)[,order_2]
write.csv(recessive_LOFmissense_hsp, "recessive_LOFmissense_hsp.csv")
recessive_missense_hsp <- left_join(recessive_missense_hsp, pli)[,order_2]
write.csv(recessive_missense_hsp, "recessive_missense_hsp.csv")

LOF_monoallelic_father_hsp <- left_join(LOF_monoallelic_father_hsp, pli)[,order_2]
write.csv(LOF_monoallelic_father_hsp, "LOF_monoallelic_father_hsp.csv")
LOF_monoallelic_mother_hsp <- left_join(LOF_monoallelic_mother_hsp, pli)[,order_2]
write.csv(LOF_monoallelic_mother_hsp, "LOF_monoallelic_mother_hsp.csv")

missense_monoallelic_father_hsp <- left_join(missense_monoallelic_father_hsp, pli)[,order_2]
write.csv(missense_monoallelic_father_hsp, "missense_monoallelic_father_hsp.csv")
missense_monoallelic_mother_hsp <- left_join(missense_monoallelic_mother_hsp, pli)[,order_2]
write.csv(missense_monoallelic_mother_hsp, "missense_monoallelic_mother_hsp.csv")


