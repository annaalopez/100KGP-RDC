#Anna López Carrasco 
# CamBridge: Linking human disease genetics to cell biology

library(Rlabkey)
#"httr" and "jsonlite" required packages. 

#Labkey/main-programme/main-programme_v14_2022-01-27/Bioinformatics/gmc_exit_questionnaire:
#Load data gmc data:
gmc_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="gmc_exit_questionnaire", 
  viewName="", 
  colSelect="participant_id,family_id,interpretation_cohort_id,interpretation_request_id,clinical_report,event_date,case_solved_family,additional_comments,segregation_question,actionability,clinical_utility,reporting_question,variant_group,acmg_classification,assembly,chromosome,position,reference,alternate,gene_name,confirmation_decision,confirmation_outcome,publications", 
  colFilter=NULL, 
  containerFilter=NULL, 
  colNameOpt="rname"
)
#37995 observations and 23 variables.

#Labkey/main-programme/main-programme_v14_2022-01-27/Common/participant:
#Load participant data:
participant_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v12_2021-05-06", 
  schemaName="lists", 
  queryName="participant", 
  viewName="", 
  colSelect="participant_id,gel_case_reference,gel_case_programme_id,latest_case_activity_datetime,participant_sk,programme,handling_gmc_ods_code,handling_gmc,handling_gmc_trust,year_of_birth,participant_ethnic_category,participant_phenotypic_sex,participant_karyotypic_sex,participant_stated_gender,registration_date,registered_at_gmc_ods_code,registered_at_gmc,registered_at_gmc_trust,registered_at_ldp_ods_code,registered_at_ldp_bioinformatics_ods_code,registered_at_ldp_organisation_name,registered_at_ldp_site_name,programme_consent_status,date_of_consent,consent_form,information_sheet,withdrawal_of_consent_date,withdrawal_option,withdrawal_form,rare_diseases_family_sk,rare_diseases_family_id,participant_type,fathers_ethnic_category,fathers_ethnic_category_other,fathers_other_relevant_ancestry,mothers_ethnic_category,mothers_ethnic_category_other,mothers_other_relevant_ancestry,participant_medical_review_date,participant_medical_review_qc_state_code,participant_medical_review_qc_state_description,consanguinity,consanguinity_sk,father_affected,mother_affected,out_of_area_recruitment,penetrance,full_brothers_affected,full_sisters_affected,total_full_brothers,total_full_sisters,biological_relationship_to_proband,other_biological_relationship_to_proband,participant_pipeline_status_id,participant_pipeline_status,reproductive_additional_findings,health_related_additional_findings,normalised_consent_form,duplicated_participant_id", 
  colFilter=NULL, 
  containerFilter=NULL, 
  colNameOpt="rname"
)
#90259 observations and 59 variables

#Filter for "Rare Diseases" programme
table(participant_data[6])
#72955 RD and 17304 cancer

install.packages("dplyr")
library(dplyr)
RD_participant_data <- (participant_data %>% filter(participant_data[6] == "Rare Diseases"))
#72955

#Filter participant data for participant consenting
table(participant_data[23])
consenting_participants <- (RD_participant_data %>% filter(RD_participant_data[23] == "Consenting"))
#72942 consenting participants and 59 variables

#Semi join gmc_data and consenting_participants to return only those participants in gmc_data that are in consenting_participants.
gmc_data_consenting <- semi_join(gmc_data, consenting_participants, by = "participant_id")
#37987 participants and 23 variables

#Labkey/main-programme/main-programme_v12_2021-05-06/Rare Disease/rare_diseases_family:
#Load rare_diseases_family data:
RD_family_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_diseases_family", 
  viewName="", 
  colSelect="rare_diseases_family_sk,rare_diseases_family_id,family_group_type,family_medical_review_date,family_medical_review_qc_state_code,family_medical_review_qc_state_description,multiple_monogenic_likely_code,multiple_monogenic_likely_description,family_non_penetrance_suspected_code,family_non_penetrance_suspected_description,family_medical_review_gmc_ods_code", 
  colFilter=NULL, 
  containerFilter=NULL, 
  colNameOpt="rname"
)
#35002 observations and 11 variables

#Filter RD_family_data for Family Group Type
table(RD_family_data[3])

#Filter only for these Family Group Type:
RD_family_data_filtered <-  filter(family_data, RD_family_data[3] == 'Duo with Mother or Father' | RD_family_data[3] == 'Families with more than three Participants' | RD_family_data[3] == 'Trio with Mother and Father' | RD_family_data[3] == 'Trio with Mother or Father and other Biological Relationship')
#20039 observations

#Rename "Rare Diseases Family Id" column to "Family Id"
RD_family_data_filtered <- rename(RD_family_data_filtered, 'family_id' = 'rare_diseases_family_id')

#Semi join gmc_data_consenting and RD_family_data_filtered
gmc_data_filtered <- semi_join(gmc_data_consenting, RD_family_data_filtered, by = "family_id")
#22186 observations

#Filter gmc_data_filtered for "no", "partially" and "unknown" case solved family
table(gmc_data_filtered[7])
gmc_unsolved <- filter(gmc_data_filtered, gmc_data_filtered[7] == 'no' | gmc_data_filtered[7] == 'partially' | gmc_data_filtered[7] == 'unknown')
#17158 observations

#Labkey/main-programme/main-programme_v14_2022-01-27/Bioinformatics/tiering_data:
#Load tiering data and gene list
tiering_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="tiering_data", 
  viewName="", 
  colSelect="participant_id,rare_diseases_family_id,interpretation_cohort_id,interpretation_request_id,sample_id,phenotype,participant_type,assembly,chromosome,position,reference,alternate,genotype,mode_of_inheritance,segregation_pattern,penetrance,tier,genomic_feature_hgnc,ensembl_id,consequence_type,so_term,father_affected,mother_affected,full_brothers_affected,full_sisters_affected,participant_phenotypic_sex", 
  colFilter=makeFilter(c("genomic_feature_hgnc", "IN", "VPS51;VPS52;VPS53;VPS54;VPS50;EIPR1;EIPR1-IT1;PLEKHA8;VAMP4;VAMP7;C10orf118;SPAST;CHMP1B;IST1;SPG20/SPART;SH3GL1;ESYT1;ATP2B1;ACSL3;ALDH3A2;MBOAT7;RAB5C;RTN4;CCDC47;CCDC136;DIXDC1;KIF5C;MLL5;NAV1;RABEP1;RUFY3;SFPQ;VCAN;ZCCHC7;ZFYVE27;TSG101;VPS25;SNF8;VPS36;VPS28;VPS37B;VPS37A;VPS37C;VPS37D;MVB12A;MVB12B;CHMP4B;CHMP2A;CHMP3;CHMP6;CHMP4A;UBAP1;VPS4A;CHMP4C;CHMP2B;VPS4B;CHMP1A;CHMP5;PDCD6IP;TUBA3E;LEMD2;PDCD6;TUBB2B;RPS27A;TUBB8B;TUBB2A;VTA1;TUBA3D;TUBB1;STAM;TTC19;TUBB4A;PLAA;UBB;TUBA8;USP8;TBK1;PTPN23;WWP1;PLD3;TUBB4B;MITD1;LITAF;WWP2;STAM2;TUBA4A;TUBB3;TUBB8;TUBA1A;TUBAL3;TMED9;LRSAM1;TUBB6;TUBA4B;STAMBP;ROCK1;TUBA1C;PMEL;SPART;UBA52;UBAP1L;TUBA1B;TUBA3C;UBC;ARRDC3;CHMP4BP1;HGS;CD2AP;CXCR4;ARRDC2;CHMP7;DTX3L;AURKB;AP1G2;CD63;CEP55;ARRDC1;CC2D1B;NEDD4;ITCH;IQGAP1;EXOC1;EXOC2;EXOC3;EXOC4;EXOC5;EXOC6;EXOC6B;EXOC7;EXOC8;AAK1;STXBP1;STXBP2;STXBP3;RAB3IP;RAB3IL1;EXOC1;RAB10;RAB12;RAB13;RAB15;RAB3A;RAB3B;RAB3C;RAB3D;RAB40B;RAB40C;RAB8A;RAB8B;EXOC2;EXOC3;EXOC3L1;EXOC3L4;TNFAIP2;ARFGEF1;ARFGEF2;EXOC4;SNAP23;SNAP25;SNAP29;EXOC5;SEC11A;SEC11B;SEC11C;PREB;SEC13;SEC14L1;SEC14L2;SEC14L3;SEC14L4;SEC14L6;EXOC6;EXOC6B;SEC16A;NAPA;NAPB;NSF;GDI1;GDI2;BNIP1;COPG1;COPG2;SEC22A;SEC22B;SEC23A;SEC23B")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
#161034 observations and 26 variables

#Semi join tiering_data and gmc_unsolved
unsolved_tiering_data <- semi_join (tiering_data, gmc_unsolved, by = "participant_id")
#29790 observations

#Filter for unaffected parents ("No" only)
table(unsolved_tiering_data[22])
table(unsolved_tiering_data[23])
unsolved_tiering_data_unaffected <- unsolved_tiering_data %>% filter(unsolved_tiering_data[22] == "No" & unsolved_tiering_data[23] == "No") 
#23776 observations

sec <- c("STXBP1","TXBP2","STXBP3","RAB3IP","RAB3IL1","EXOC1","RAB10","RAB12","RAB13","RAB15","RAB3A","RAB3B","RAB3C","RAB3D","RAB40B","RAB40C","RAB8A","RAB8B","EXOC2","EXOC3","EXOC3L1","EXOC3L4","TNFAIP2","ARFGEF1","ARFGEF2","EXOC4","SNAP23","SNAP25","SNAP29","EXOC5","SEC11A","SEC11B","SEC11C","PREB","SEC13","SEC14L1","SEC14L2","SEC14L3","SEC14L4","SEC14L6","EXOC6","EXOC6B","SEC16A","NAPA","NAPB","NSF","GDI1","GDI2","BNIP1","COPG1","COPG2","SEC22A","SEC22B","SEC23A","SEC23B","SERPINI1")
hsp <- c("VPS51","VPS52","VPS53","VPS54","VPS50","EIPR1","EIPR1-IT1","PLEKHA8","VAMP4","VAMP7","C10orf118","SPAST","CHMP1B","IST1","SPG20/SPART","SH3GL1","ESYT1","ATP2B1","ACSL3","ALDH3A2","MBOAT7","RAB5C","RTN4","CCDC47","CCDC136","DIXDC1","KIF5C","MLL5","NAV1","RABEP1","RUFY3","SFPQ","VCAN","ZCCHC7","ZFYVE27","TSG101","VPS25","SNF8","VPS36","VPS28","VPS37B","VPS37A","VPS37C","VPS37D","MVB12A","MVB12B","CHMP4B","CHMP2A","CHMP3","CHMP6","CHMP4A","UBAP1","VPS4A","CHMP4C","CHMP2B","VPS4B","CHMP1A","CHMP5","PDCD6IP","TUBA3E","LEMD2","PDCD6","TUBB2B","RPS27A","TUBB8B","TUBB2A","VTA1","TUBA3D","TUBB1","STAM","TTC19","TUBB4A","PLAA","UBB","TUBA8","USP8","TBK1","PTPN23","WWP1","PLD3","TUBB4B","MITD1","LITAF","WWP2","STAM2","TUBA4A","TUBB3","TUBB8","TUBA1A","TUBAL3","TMED9","LRSAM1","TUBB6","TUBA4B","STAMBP","ROCK1","TUBA1C","PMEL","SPART","UBA52","UBAP1L","TUBA1B","TUBA3C","UBC","ARRDC3","CHMP4BP1","HGS","CD2AP","CXCR4","ARRDC2","CHMP7","DTX3L","AURKB","AP1G2","CD63","CEP55","ARRDC1","CC2D1B","NEDD4","ITCH","IQGAP1","EXOC1","EXOC2","EXOC3","EXOC4","EXOC5","EXOC6","EXOC6B","EXOC7","EXOC8","AAK1")


###########################################
# CONSEQUENCE TYPE AND SEGREGATION PATTERN #
###########################################

#######
# LOF #
######
table(unsolved_tiering_data_unaffected[20])
LOF <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("feature_truncation","frameshift_variant", "splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost","stop_gained"))

# Simple recessive LOF
unsolved_recessive_lof <- (LOF %>% filter(LOF[15] == 'SimpleRecessive'))
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
recessive_lof <- unsolved_recessive_lof[!duplicated(unsolved_recessive_lof[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
recessive <- recessive_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
recessive <- recessive %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
recessive_lof <- recessive_lof[c(-6,-20)]
# Join 
unsolved_recessive_lof <- left_join(recessive, recessive_lof)
# Remove duplicates
unsolved_recessive_lof <- unsolved_recessive_lof[!duplicated(unsolved_recessive_lof[c(1,2,3,4,5,6)]),]
# Filter for sec and hsp genes
unsolved_recessive_lof_sec <- filter(unsolved_recessive_lof, genomic_feature_hgnc %in% sec)
unsolved_recessive_lof_hsp <- filter(unsolved_recessive_lof, genomic_feature_hgnc %in% hsp)
# Save
write.csv(unsolved_recessive_lof_sec, "unsolved_recessive_lof_sec.csv")
write.csv(unsolved_recessive_lof_hsp, "unsolved_recessive_lof_hsp.csv")
write.csv(unsolved_recessive_lof, "unsolved_recessive_lof.csv")

# DeNovo LOF
unsolved_denovo_lof <- (LOF %>% filter(LOF[15] == 'deNovo'))
table(unsolved_denovo_lof[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
denovo_lof <- unsolved_denovo_lof[!duplicated(unsolved_denovo_lof[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
denovo <- denovo_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
denovo <- denovo %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
denovo_lof <- denovo_lof[c(-6,-20)]
# Join 
unsolved_denovo_lof <- left_join(denovo, denovo_lof)
# Remove duplicates
unsolved_denovo_lof <- unsolved_denovo_lof[!duplicated(unsolved_denovo_lof[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
unsolved_denovo_lof_sec <- filter(unsolved_denovo_lof, genomic_feature_hgnc %in% sec)
unsolved_denovo_lof_hsp <- filter(unsolved_denovo_lof, genomic_feature_hgnc %in% hsp)
#Save 
write.csv(unsolved_denovo_lof_sec, "unsolved_denovo_lof_sec.csv")
write.csv(unsolved_denovo_lof_hsp, "unsolved_denovo_lof_hsp.csv")
write.csv(unsolved_denovo_lof, "unsolved_denovo_lof.csv")

# CompoundHeterozygous LOF
library(tidyverse)
unsolved_comphet_lof <- (LOF %>% filter(LOF[15] == "CompoundHeterozygous")) 
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
comphet_lof <- unsolved_comphet_lof[!duplicated(unsolved_comphet_lof[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
comphet <- comphet_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
comphet <- comphet %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
comphet_lof <- comphet_lof[c(-6,-20)]
# Join 
unsolved_comphet_lof <- left_join(comphet, comphet_lof)
# Remove duplicates
unsolved_comphet_lof <- unsolved_comphet_lof[!duplicated(unsolved_comphet_lof[c(1,2,3,4,5,6)]),]
# Remove participant_id that only appears once
unsolved_comphet_lof <- unsolved_comphet_lof %>%
  group_by(participant_id) %>%
  filter(n()>1)
#Filter for sec and hsp genes
unsolved_comphet_lof_sec <- filter(unsolved_comphet_lof, genomic_feature_hgnc %in% sec)
unsolved_comphet_lof_hsp <- filter(unsolved_comphet_lof, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_comphet_lof_sec, "unsolved_comphet_lof_sec.csv")
write.csv(unsolved_comphet_lof_hsp, "unsolved_comphet_lof_hsp.csv")
write.csv(unsolved_comphet_lof, "unsolved_comphet_lof.csv")

############
# MISSENSE+ #
############

table(unsolved_tiering_data_unaffected[20])
missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("inframe_deletion","inframe_insertion","missense_variant", "stop_gained","stop_lost"))

# Simple recessive missense+ 
unsolved_recessive_missense <- (missense %>% filter(missense[15] == 'SimpleRecessive'))
table(unsolved_recessive_missense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
recessive_missense <- unsolved_recessive_missense[!duplicated(unsolved_recessive_missense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
recessive_miss <- recessive_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
recessive_miss <- recessive_miss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
recessive_missense <- recessive_missense[c(-6,-20)]
# Join 
unsolved_recessive_missense <- left_join(recessive_miss, recessive_missense)
# Remove duplicates
unsolved_recessive_missense <- unsolved_recessive_missense[!duplicated(unsolved_recessive_missense[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
unsolved_recessive_missense_sec <- filter(unsolved_recessive_missense, genomic_feature_hgnc %in% sec)
unsolved_recessive_missense_hsp <- filter(unsolved_recessive_missense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_recessive_missense_sec, "unsolved_recessive_missense_sec.csv")
write.csv(unsolved_recessive_missense_hsp, "unsolved_recessive_missense_hsp.csv")
write.csv(unsolved_recessive_missense, "unsolved_recessive_missense.csv")

# DeNovo missense+
unsolved_denovo_missense <- (missense %>% filter(missense[15] == 'deNovo'))
table(unsolved_denovo_missense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
denovo_missense <- unsolved_denovo_missense[!duplicated(unsolved_denovo_missense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
denovo_miss <- denovo_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
denovo_miss <- denovo_miss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
denovo_missense <- denovo_missense[c(-6,-20)]
# Join 
unsolved_denovo_missense <- left_join(denovo_miss, denovo_missense)
# Remove duplicates
unsolved_denovo_missense <- unsolved_denovo_missense[!duplicated(unsolved_denovo_missense[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
unsolved_denovo_missense_sec <- filter(unsolved_denovo_missense, genomic_feature_hgnc %in% sec)
unsolved_denovo_missense_hsp <- filter(unsolved_denovo_missense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_denovo_missense_sec, "unsolved_denovo_missense_sec.csv")
write.csv(unsolved_denovo_missense_hsp, "unsolved_denovo_missense_hsp.csv")
write.csv(unsolved_denovo_missense, "unsolved_denovo_missense.csv")

# CompoundHeterozygous missense+
unsolved_comphet_missense <- (missense %>% filter(missense[15] == "CompoundHeterozygous"))
table(unsolved_comphet_missense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
comphet_missense <- unsolved_comphet_missense[!duplicated(unsolved_comphet_missense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
comphet_miss <- comphet_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
comphet_miss <- comphet_miss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
comphet_missense <- comphet_missense[c(-6,-20)]
# Join 
unsolved_comphet_missense <- left_join(comphet_miss, comphet_missense)
# Remove duplicates
unsolved_comphet_missense <- unsolved_comphet_missense[!duplicated(unsolved_comphet_missense[c(1,2,3,4,5,6)]),]
# Remove participant_id that only appears once
unsolved_comphet_missense <- unsolved_comphet_missense %>%
  group_by(participant_id) %>%
  filter(n()>1)
#Filter for sec and hsp genes
unsolved_comphet_missense_sec <- filter(unsolved_comphet_missense, genomic_feature_hgnc %in% sec)
unsolved_comphet_missense_hsp <- filter(unsolved_comphet_missense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_comphet_missense_sec, "unsolved_comphet_missense_sec.csv")
write.csv(unsolved_comphet_missense_hsp, "unsolved_comphet_missense_hsp.csv")
write.csv(unsolved_comphet_missense, "unsolved_comphet_missense.csv")


##################
# LOF + MISSENSE #
#################

LOF_missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "start_lost", "stop_gained", "stop_lost"))

# Simple recessive LOF + missense
unsolved_recessive_LOFmissense <- (LOF_missense %>% filter(LOF_missense[15] == 'SimpleRecessive'))
table(unsolved_recessive_LOFmissense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
recessive_LOFmissense <- unsolved_recessive_LOFmissense[!duplicated(unsolved_recessive_LOFmissense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
recessive_LOFmiss <- recessive_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
recessive_LOFmiss <- recessive_LOFmiss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
recessive_LOFmissense <- recessive_LOFmissense[c(-6,-20)]
# Join 
unsolved_recessive_LOFmissense <- left_join(recessive_LOFmiss, recessive_LOFmissense)
# Remove duplicates
unsolved_recessive_LOFmissense <- unsolved_recessive_LOFmissense[!duplicated(unsolved_recessive_LOFmissense[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
unsolved_recessive_LOFmissense_sec <- filter(unsolved_recessive_LOFmissense, genomic_feature_hgnc %in% sec)
unsolved_recessive_LOFmissense_hsp <- filter(unsolved_recessive_LOFmissense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_recessive_LOFmissense_sec, "unsolved_recessive_LOFmissense_sec.csv")
write.csv(unsolved_recessive_LOFmissense_hsp, "unsolved_recessive_LOFmissense_hsp.csv")
write.csv(unsolved_recessive_LOFmissense, "unsolved_recessive_LOFmissense.csv")

# DeNovo LOF + missense
unsolved_denovo_LOFmissense <- (missense %>% filter(missense[15] == 'deNovo'))
table(unsolved_denovo_LOFmissense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
denovo_LOFmissense <- unsolved_denovo_LOFmissense[!duplicated(unsolved_denovo_LOFmissense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
denovo_LOFmiss <- denovo_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
denovo_LOFmiss <- denovo_LOFmiss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
denovo_LOFmissense <- denovo_LOFmissense[c(-6,-20)]
# Join 
unsolved_denovo_LOFmissense <- left_join(denovo_LOFmiss, denovo_LOFmissense)
# Remove duplicates
unsolved_denovo_LOFmissense <- unsolved_denovo_LOFmissense[!duplicated(unsolved_denovo_LOFmissense[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
unsolved_denovo_LOFmissense_sec <- filter(unsolved_denovo_LOFmissense, genomic_feature_hgnc %in% sec)
unsolved_denovo_LOFmissense_hsp <- filter(unsolved_denovo_LOFmissense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_denovo_LOFmissense_sec, "unsolved_denovo_LOFmissense_sec.csv")
write.csv(unsolved_denovo_LOFmissense_hsp, "unsolved_denovo_LOFmissense_hsp.csv")
write.csv(unsolved_denovo_LOFmissense, "unsolved_denovo_LOFmissense.csv")

# CompoundHeterozygous LOF + missense
unsolved_comphet_LOFmissense <- (missense %>% filter(missense[15] == "CompoundHeterozygous"))
table(unsolved_comphet_LOFmissense[20])
# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
comphet_LOFmissense <- unsolved_comphet_LOFmissense[!duplicated(unsolved_comphet_LOFmissense[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
comphet_LOFmiss <- comphet_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequene_type
comphet_LOFmiss <- comphet_LOFmiss %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
comphet_LOFmissense <- comphet_LOFmissense[c(-6,-20)]
# Join 
unsolved_comphet_LOFmissense <- left_join(comphet_LOFmiss, comphet_LOFmissense)
# Remove duplicates
unsolved_comphet_LOFmissense <- unsolved_comphet_LOFmissense[!duplicated(unsolved_comphet_LOFmissense[c(1,2,3,4,5,6)]),]
# Remove participant_id that only appears once
unsolved_comphet_LOFmissense <- unsolved_comphet_LOFmissense %>%
  group_by(participant_id) %>%
  filter(n()>1)
#Filter for sec and hsp genes
unsolved_comphet_LOFmissense_sec <- filter(unsolved_comphet_LOFmissense, genomic_feature_hgnc %in%  sec)
unsolved_comphet_LOFmissense_hsp <- filter(unsolved_comphet_LOFmissense, genomic_feature_hgnc %in% hsp)
#Save
write.csv(unsolved_comphet_LOFmissense_sec, "unsolved_comphet_LOFmissense_sec.csv")
write.csv(unsolved_comphet_LOFmissense_hsp, "unsolved_comphet_LOFmissense_hsp.csv")
write.csv(unsolved_comphet_LOFmissense, "unsolved_comphet_LOFmissense.csv")


#########
#  VEP  #
#########

library(stringr)

#Simple recessive LOF
# SEC GENES
table(unsolved_recessive_lof_sec[12])

unsolved_recessive_lof_sec_grch37 <- unsolved_recessive_lof_sec %>% filter(assembly == "GRCh37")
unsolved_recessive_lof_sec_grch37$base <- paste(as.character(unsolved_recessive_lof_sec_grch37$reference), "/", as.character(unsolved_recessive_lof_sec_grch37$alternate))
unsolved_recessive_lof_sec_grch37_VEP <- unsolved_recessive_lof_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_sec_grch37_VEP_list <- unsolved_recessive_lof_sec_grch37_VEP[!duplicated(unsolved_recessive_lof_sec_grch37_VEP[,2:4]),]
unsolved_recessive_lof_sec_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_lof_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_sec_grch37_VEP_list, "VEP_recessive_lof_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_lof_sec_grch38 <- unsolved_recessive_lof_sec %>% filter(assembly == "GRCh38")
unsolved_recessive_lof_sec_grch38$base <- paste(as.character(unsolved_recessive_lof_sec_grch38$reference), "/", as.character(unsolved_recessive_lof_sec_grch38$alternate))
unsolved_recessive_lof_sec_grch38_VEP <- unsolved_recessive_lof_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_sec_grch38_VEP_list <- unsolved_recessive_lof_sec_grch38_VEP[!duplicated(unsolved_recessive_lof_sec_grch38_VEP[,2:4]),]
unsolved_recessive_lof_sec_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_lof_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_sec_grch38_VEP_list, "VEP_recessive_lof_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_recessive_lof_hsp[12])

unsolved_recessive_lof_hsp_grch37 <- unsolved_recessive_lof_hsp %>% filter(assembly == "GRCh37")
unsolved_recessive_lof_hsp_grch37$base <- paste(as.character(unsolved_recessive_lof_hsp_grch37$reference), "/", as.character(unsolved_recessive_lof_hsp_grch37$alternate))
unsolved_recessive_lof_hsp_grch37_VEP <- unsolved_recessive_lof_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_hsp_grch37_VEP_list <- unsolved_recessive_lof_hsp_grch37_VEP[!duplicated(unsolved_recessive_lof_hsp_grch37_VEP[,2:4]),]
unsolved_recessive_lof_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_lof_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_hsp_grch37_VEP_list, "VEP_recessive_lof_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_lof_hsp_grch38 <- unsolved_recessive_lof_hsp %>% filter(assembly == "GRCh38")
unsolved_recessive_lof_hsp_grch38$base <- paste(as.character(unsolved_recessive_lof_hsp_grch38$reference), "/", as.character(unsolved_recessive_lof_hsp_grch38$alternate))
unsolved_recessive_lof_hsp_grch38_VEP <- unsolved_recessive_lof_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_hsp_grch38_VEP_list <- unsolved_recessive_lof_hsp_grch38_VEP[!duplicated(unsolved_recessive_lof_hsp_grch38_VEP[,2:4]),]
unsolved_recessive_lof_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_lof_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_hsp_grch38_VEP_list, "VEP_recessive_lof_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#Simple recessive missense+
# SEC GENES
table(unsolved_recessive_missense_sec[12])

unsolved_recessive_missense_sec_grch37 <- unsolved_recessive_missense_sec %>% filter(assembly == "GRCh37")
unsolved_recessive_missense_sec_grch37$base <- paste(as.character(unsolved_recessive_missense_sec_grch37$reference), "/", as.character(unsolved_recessive_missense_sec_grch37$alternate))
unsolved_recessive_missense_sec_grch37_VEP <- unsolved_recessive_missense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_sec_grch37_VEP_list <- unsolved_recessive_missense_sec_grch37_VEP[!duplicated(unsolved_recessive_missense_sec_grch37_VEP[,2:4]),]
unsolved_recessive_missense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_missense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_sec_grch37_VEP_list, "VEP_recessive_missense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_missense_sec_grch38 <- unsolved_recessive_missense_sec %>% filter(assembly == "GRCh38")
unsolved_recessive_missense_sec_grch38$base <- paste(as.character(unsolved_recessive_missense_sec_grch38$reference), "/", as.character(unsolved_recessive_missense_sec_grch38$alternate))
unsolved_recessive_missense_sec_grch38_VEP <- unsolved_recessive_missense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_sec_grch38_VEP_list <- unsolved_recessive_missense_sec_grch38_VEP[!duplicated(unsolved_recessive_missense_sec_grch38_VEP[,2:4]),]
unsolved_recessive_missense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_missense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_sec_grch38_VEP_list, "VEP_recessive_missense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_recessive_missense_hsp[12])

unsolved_recessive_missense_hsp_grch37 <- unsolved_recessive_missense_hsp %>% filter(assembly == "GRCh37")
unsolved_recessive_missense_hsp_grch37$base <- paste(as.character(unsolved_recessive_missense_hsp_grch37$reference), "/", as.character(unsolved_recessive_missense_hsp_grch37$alternate))
unsolved_recessive_missense_hsp_grch37_VEP <- unsolved_recessive_missense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_hsp_grch37_VEP_list <- unsolved_recessive_missense_hsp_grch37_VEP[!duplicated(unsolved_recessive_missense_hsp_grch37_VEP[,2:4]),]
unsolved_recessive_missense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_missense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_hsp_grch37_VEP_list, "VEP_recessive_missense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_missense_hsp_grch38 <- unsolved_recessive_missense_hsp %>% filter(assembly == "GRCh38")
unsolved_recessive_missense_hsp_grch38$base <- paste(as.character(unsolved_recessive_missense_hsp_grch38$reference), "/", as.character(unsolved_recessive_missense_hsp_grch38$alternate))
unsolved_recessive_missense_hsp_grch38_VEP <- unsolved_recessive_missense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_hsp_grch38_VEP_list <- unsolved_recessive_missense_hsp_grch38_VEP[!duplicated(unsolved_recessive_missense_hsp_grch38_VEP[,2:4]),]
unsolved_recessive_missense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_missense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_hsp_grch38_VEP_list, "VEP_recessive_missense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#Simple recessive LOF + missense
# SEC GENES
table(unsolved_recessive_LOFmissense_sec[12])

unsolved_recessive_LOFmissense_sec_grch37 <- unsolved_recessive_LOFmissense_sec %>% filter(assembly == "GRCh37")
unsolved_recessive_LOFmissense_sec_grch37$base <- paste(as.character(unsolved_recessive_LOFmissense_sec_grch37$reference), "/", as.character(unsolved_recessive_LOFmissense_sec_grch37$alternate))
unsolved_recessive_LOFmissense_sec_grch37_VEP <- unsolved_recessive_LOFmissense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_sec_grch37_VEP_list <- unsolved_recessive_LOFmissense_sec_grch37_VEP[!duplicated(unsolved_recessive_LOFmissense_sec_grch37_VEP[,2:4]),]
unsolved_recessive_LOFmissense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_sec_grch37_VEP_list, "VEP_recessive_LOFmissense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_LOFmissense_sec_grch38 <- unsolved_recessive_LOFmissense_sec %>% filter(assembly == "GRCh38")
unsolved_recessive_LOFmissense_sec_grch38$base <- paste(as.character(unsolved_recessive_LOFmissense_sec_grch38$reference), "/", as.character(unsolved_recessive_LOFmissense_sec_grch38$alternate))
unsolved_recessive_LOFmissense_sec_grch38_VEP <- unsolved_recessive_LOFmissense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_sec_grch38_VEP_list <- unsolved_recessive_LOFmissense_sec_grch38_VEP[!duplicated(unsolved_recessive_LOFmissense_sec_grch38_VEP[,2:4]),]
unsolved_recessive_LOFmissense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_sec_grch38_VEP_list, "VEP_recessive_LOFmissense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_recessive_LOFmissense_hsp[12])

unsolved_recessive_LOFmissense_hsp_grch37 <- unsolved_recessive_LOFmissense_hsp %>% filter(assembly == "GRCh37")
unsolved_recessive_LOFmissense_hsp_grch37$base <- paste(as.character(unsolved_recessive_LOFmissense_hsp_grch37$reference), "/", as.character(unsolved_recessive_LOFmissense_hsp_grch37$alternate))
unsolved_recessive_LOFmissense_hsp_grch37_VEP <- unsolved_recessive_LOFmissense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_hsp_grch37_VEP_list <- unsolved_recessive_LOFmissense_hsp_grch37_VEP[!duplicated(unsolved_recessive_LOFmissense_hsp_grch37_VEP[,2:4]),]
unsolved_recessive_LOFmissense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_hsp_grch37_VEP_list, "VEP_recessive_LOFmissense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_LOFmissense_hsp_grch38 <- unsolved_recessive_LOFmissense_hsp %>% filter(assembly == "GRCh38")
unsolved_recessive_LOFmissense_hsp_grch38$base <- paste(as.character(unsolved_recessive_LOFmissense_hsp_grch38$reference), "/", as.character(unsolved_recessive_LOFmissense_hsp_grch38$alternate))
unsolved_recessive_LOFmissense_hsp_grch38_VEP <- unsolved_recessive_LOFmissense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_hsp_grch38_VEP_list <- unsolved_recessive_LOFmissense_hsp_grch38_VEP[!duplicated(unsolved_recessive_LOFmissense_hsp_grch38_VEP[,2:4]),]
unsolved_recessive_LOFmissense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_hsp_grch38_VEP_list, "VEP_recessive_LOFmissense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# deNovo LOF
# SEC GENES
table(unsolved_denovo_lof_sec[12])

unsolved_denovo_lof_sec_grch37 <- unsolved_denovo_lof_sec %>% filter(assembly == "GRCh37")
unsolved_denovo_lof_sec_grch37$base <- paste(as.character(unsolved_denovo_lof_sec_grch37$reference), "/", as.character(unsolved_denovo_lof_sec_grch37$alternate))
unsolved_denovo_lof_sec_grch37_VEP <- unsolved_denovo_lof_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_sec_grch37_VEP_list <- unsolved_denovo_lof_sec_grch37_VEP[!duplicated(unsolved_denovo_lof_sec_grch37_VEP[,2:4]),]
unsolved_denovo_lof_sec_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_lof_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_sec_grch37_VEP_list, "VEP_denovo_lof_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_lof_sec_grch38 <- unsolved_denovo_lof_sec %>% filter(assembly == "GRCh38")
unsolved_denovo_lof_sec_grch38$base <- paste(as.character(unsolved_denovo_lof_sec_grch38$reference), "/", as.character(unsolved_denovo_lof_sec_grch38$alternate))
unsolved_denovo_lof_sec_grch38_VEP <- unsolved_denovo_lof_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_sec_grch38_VEP_list <- unsolved_denovo_lof_sec_grch38_VEP[!duplicated(unsolved_denovo_lof_sec_grch38_VEP[,2:4]),]
unsolved_denovo_lof_sec_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_lof_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_sec_grch38_VEP_list, "VEP_denovo_lof_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_denovo_lof_hsp[12])

unsolved_denovo_lof_hsp_grch37 <- unsolved_denovo_lof_hsp %>% filter(assembly == "GRCh37")
unsolved_denovo_lof_hsp_grch37$base <- paste(as.character(unsolved_denovo_lof_hsp_grch37$reference), "/", as.character(unsolved_denovo_lof_hsp_grch37$alternate))
unsolved_denovo_lof_hsp_grch37_VEP <- unsolved_denovo_lof_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_hsp_grch37_VEP_list <- unsolved_denovo_lof_hsp_grch37_VEP[!duplicated(unsolved_denovo_lof_hsp_grch37_VEP[,2:4]),]
unsolved_denovo_lof_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_lof_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_hsp_grch37_VEP_list, "VEP_denovo_lof_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_lof_hsp_grch38 <- unsolved_denovo_lof_hsp %>% filter(assembly == "GRCh38")
unsolved_denovo_lof_hsp_grch38$base <- paste(as.character(unsolved_denovo_lof_hsp_grch38$reference), "/", as.character(unsolved_denovo_lof_hsp_grch38$alternate))
unsolved_denovo_lof_hsp_grch38_VEP <- unsolved_denovo_lof_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_hsp_grch38_VEP_list <- unsolved_denovo_lof_hsp_grch38_VEP[!duplicated(unsolved_denovo_lof_hsp_grch38_VEP[,2:4]),]
unsolved_denovo_lof_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_lof_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_hsp_grch38_VEP_list, "VEP_denovo_lof_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#deNovo missense+
# SEC GENES
table(unsolved_denovo_missense_sec[12])

unsolved_denovo_missense_sec_grch37 <- unsolved_denovo_missense_sec %>% filter(assembly == "GRCh37")
unsolved_denovo_missense_sec_grch37$base <- paste(as.character(unsolved_denovo_missense_sec_grch37$reference), "/", as.character(unsolved_denovo_missense_sec_grch37$alternate))
unsolved_denovo_missense_sec_grch37_VEP <- unsolved_denovo_missense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_sec_grch37_VEP_list <- unsolved_denovo_missense_sec_grch37_VEP[!duplicated(unsolved_denovo_missense_sec_grch37_VEP[,2:4]),]
unsolved_denovo_missense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_missense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_sec_grch37_VEP_list, "VEP_denovo_missense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_missense_sec_grch38 <-  unsolved_denovo_missense_sec %>% filter(assembly == "GRCh38")
unsolved_denovo_missense_sec_grch38$base <- paste(as.character(unsolved_denovo_missense_sec_grch38$reference), "/", as.character(unsolved_denovo_missense_sec_grch38$alternate))
unsolved_denovo_missense_sec_grch38_VEP <- unsolved_denovo_missense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_sec_grch38_VEP_list <- unsolved_denovo_missense_sec_grch38_VEP[!duplicated(unsolved_denovo_missense_sec_grch38_VEP[,2:4]),]
unsolved_denovo_missense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_missense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_sec_grch38_VEP_list, "VEP_denovo_missense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_denovo_missense_hsp[12])

unsolved_denovo_missense_hsp_grch37 <- unsolved_denovo_missense_hsp %>% filter(assembly == "GRCh37")
unsolved_denovo_missense_hsp_grch37$base <- paste(as.character(unsolved_denovo_missense_hsp_grch37$reference), "/", as.character(unsolved_denovo_missense_hsp_grch37$alternate))
unsolved_denovo_missense_hsp_grch37_VEP <- unsolved_denovo_missense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_hsp_grch37_VEP_list <- unsolved_denovo_missense_hsp_grch37_VEP[!duplicated(unsolved_denovo_missense_hsp_grch37_VEP[,2:4]),]
unsolved_denovo_missense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_missense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_hsp_grch37_VEP_list, "VEP_denovo_missense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_missense_hsp_grch38 <- unsolved_denovo_missense_hsp %>% filter(assembly == "GRCh38")
unsolved_denovo_missense_hsp_grch38$base <- paste(as.character(unsolved_denovo_missense_hsp_grch38$reference), "/", as.character(unsolved_denovo_missense_hsp_grch38$alternate))
unsolved_denovo_missense_hsp_grch38_VEP <- unsolved_denovo_missense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_hsp_grch38_VEP_list <- unsolved_denovo_missense_hsp_grch38_VEP[!duplicated(unsolved_denovo_missense_hsp_grch38_VEP[,2:4]),]
unsolved_denovo_missense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_missense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_hsp_grch38_VEP_list, "VEP_denovo_missense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#deNovo LOF + missense
# SEC GENES
table(unsolved_denovo_LOFmissense_sec[12])

unsolved_denovo_LOFmissense_sec_grch37 <- unsolved_denovo_LOFmissense_sec %>% filter(assembly == "GRCh37")
unsolved_denovo_LOFmissense_sec_grch37$base <- paste(as.character(unsolved_denovo_LOFmissense_sec_grch37$reference), "/", as.character(unsolved_denovo_LOFmissense_sec_grch37$alternate))
unsolved_denovo_LOFmissense_sec_grch37_VEP <- unsolved_denovo_LOFmissense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_sec_grch37_VEP_list <- unsolved_denovo_LOFmissense_sec_grch37_VEP[!duplicated(unsolved_denovo_LOFmissense_sec_grch37_VEP[,2:4]),]
unsolved_denovo_LOFmissense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_sec_grch37_VEP_list, "VEP_denovo_LOFmissense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_LOFmissense_sec_grch38 <- unsolved_denovo_LOFmissense_sec %>% filter(assembly == "GRCh38")
unsolved_denovo_LOFmissense_sec_grch38$base <- paste(as.character(unsolved_denovo_LOFmissense_sec_grch38$reference), "/", as.character(unsolved_denovo_LOFmissense_sec_grch38$alternate))
unsolved_denovo_LOFmissense_sec_grch38_VEP <- unsolved_denovo_LOFmissense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_sec_grch38_VEP_list <- unsolved_denovo_LOFmissense_sec_grch38_VEP[!duplicated(unsolved_denovo_LOFmissense_sec_grch38_VEP[,2:4]),]
unsolved_denovo_LOFmissense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_sec_grch38_VEP_list, "VEP_denovo_LOFmissense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_denovo_LOFmissense_hsp[12])

unsolved_denovo_LOFmissense_hsp_grch37 <- unsolved_denovo_LOFmissense_hsp %>% filter(assembly == "GRCh37")
unsolved_denovo_LOFmissense_hsp_grch37$base <- paste(as.character(unsolved_denovo_LOFmissense_hsp_grch37$reference), "/", as.character(unsolved_denovo_LOFmissense_hsp_grch37$alternate))
unsolved_denovo_LOFmissense_hsp_grch37_VEP <- unsolved_denovo_LOFmissense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_hsp_grch37_VEP_list <- unsolved_denovo_LOFmissense_hsp_grch37_VEP[!duplicated(unsolved_denovo_LOFmissense_hsp_grch37_VEP[,2:4]),]
unsolved_denovo_LOFmissense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_hsp_grch37_VEP_list, "VEP_denovo_LOFmissense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_LOFmissense_hsp_grch38 <- unsolved_denovo_LOFmissense_hsp %>% filter(assembly == "GRCh38")
unsolved_denovo_LOFmissense_hsp_grch38$base <- paste(as.character(unsolved_denovo_LOFmissense_hsp_grch38$reference), "/", as.character(unsolved_denovo_LOFmissense_hsp_grch38$alternate))
unsolved_denovo_LOFmissense_hsp_grch38_VEP <- unsolved_denovo_LOFmissense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_hsp_grch38_VEP_list <- unsolved_denovo_LOFmissense_hsp_grch38_VEP[!duplicated(unsolved_denovo_LOFmissense_hsp_grch38_VEP[,2:4]),]
unsolved_denovo_LOFmissense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_hsp_grch38_VEP_list, "VEP_denovo_LOFmissense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# Comphet LOF
# SEC GENES
table(unsolved_comphet_lof_sec[12])

unsolved_comphet_lof_sec_grch38 <- unsolved_comphet_lof_sec %>% filter(assembly == "GRCh38")
unsolved_comphet_lof_sec_grch38$base <- paste(as.character(unsolved_comphet_lof_sec_grch38$reference), "/", as.character(unsolved_comphet_lof_sec_grch38$alternate))
unsolved_comphet_lof_sec_grch38_VEP <- unsolved_comphet_lof_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_sec_grch38_VEP_list <- unsolved_comphet_lof_sec_grch38_VEP[!duplicated(unsolved_comphet_lof_sec_grch38_VEP[,2:4]),]
unsolved_comphet_lof_sec_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_lof_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_sec_grch38_VEP_list, "VEP_comphet_lof_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_comphet_lof_hsp[12])

unsolved_comphet_lof_hsp_grch37 <- unsolved_comphet_lof_hsp %>% filter(assembly == "GRCh37")
unsolved_comphet_lof_hsp_grch37$base <- paste(as.character(unsolved_comphet_lof_hsp_grch37$reference), "/", as.character(unsolved_comphet_lof_hsp_grch37$alternate))
unsolved_comphet_lof_hsp_grch37_VEP <- unsolved_comphet_lof_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_hsp_grch37_VEP_list <- unsolved_comphet_lof_hsp_grch37_VEP[!duplicated(unsolved_comphet_lof_hsp_grch37_VEP[,2:4]),]
unsolved_comphet_lof_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_lof_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_hsp_grch37_VEP_list, "VEP_comphet_lof_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_lof_hsp_grch38 <- unsolved_comphet_lof_hsp %>% filter(assembly == "GRCh38")
unsolved_comphet_lof_hsp_grch38$base <- paste(as.character(unsolved_comphet_lof_hsp_grch38$reference), "/", as.character(unsolved_comphet_lof_hsp_grch38$alternate))
unsolved_comphet_lof_hsp_grch38_VEP <- unsolved_comphet_lof_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_hsp_grch38_VEP_list <- unsolved_comphet_lof_hsp_grch38_VEP[!duplicated(unsolved_comphet_lof_hsp_grch38_VEP[,2:4]),]
unsolved_comphet_lof_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_lof_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_hsp_grch38_VEP_list, "VEP_comphet_lof_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Comphet missense+
# SEC GENES
table(unsolved_comphet_missense_sec[12])

unsolved_comphet_missense_sec_grch37 <- unsolved_comphet_missense_sec %>% filter(assembly == "GRCh37")
unsolved_comphet_missense_sec_grch37$base <- paste(as.character(unsolved_comphet_missense_sec_grch37$reference), "/", as.character(unsolved_comphet_missense_sec_grch37$alternate))
unsolved_comphet_missense_sec_grch37_VEP <- unsolved_comphet_missense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_sec_grch37_VEP_list <- unsolved_comphet_missense_sec_grch37_VEP[!duplicated(unsolved_comphet_missense_sec_grch37_VEP[,2:4]),]
unsolved_comphet_missense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_missense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_sec_grch37_VEP_list, "VEP_comphet_missense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_missense_sec_grch38 <- unsolved_comphet_missense_sec %>% filter(assembly == "GRCh38")
unsolved_comphet_missense_sec_grch38$base <- paste(as.character(unsolved_comphet_missense_sec_grch38$reference), "/", as.character(unsolved_comphet_missense_sec_grch38$alternate))
unsolved_comphet_missense_sec_grch38_VEP <- unsolved_comphet_missense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_sec_grch38_VEP_list <- unsolved_comphet_missense_sec_grch38_VEP[!duplicated(unsolved_comphet_missense_sec_grch38_VEP[,2:4]),]
unsolved_comphet_missense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_missense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_sec_grch38_VEP_list, "VEP_comphet_missense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_comphet_missense_hsp[12])

unsolved_comphet_missense_hsp_grch37 <- unsolved_comphet_missense_hsp %>% filter(assembly == "GRCh37")
unsolved_comphet_missense_hsp_grch37$base <- paste(as.character(unsolved_comphet_missense_hsp_grch37$reference), "/", as.character(unsolved_comphet_missense_hsp_grch37$alternate))
unsolved_comphet_missense_hsp_grch37_VEP <- unsolved_comphet_missense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_hsp_grch37_VEP_list <- unsolved_comphet_missense_hsp_grch37_VEP[!duplicated(unsolved_comphet_missense_hsp_grch37_VEP[,2:4]),]
unsolved_comphet_missense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_missense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_hsp_grch37_VEP_list, "VEP_comphet_missense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_missense_hsp_grch38 <- unsolved_comphet_missense_hsp %>% filter(assembly == "GRCh38")
unsolved_comphet_missense_hsp_grch38$base <- paste(as.character(unsolved_comphet_missense_hsp_grch38$reference), "/", as.character(unsolved_comphet_missense_hsp_grch38$alternate))
unsolved_comphet_missense_hsp_grch38_VEP <- unsolved_comphet_missense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_hsp_grch38_VEP_list <- unsolved_comphet_missense_hsp_grch38_VEP[!duplicated(unsolved_comphet_missense_hsp_grch38_VEP[,2:4]),]
unsolved_comphet_missense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_missense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_hsp_grch38_VEP_list, "VEP_comphet_missense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Comphet LOF + missense
# SEC GENES
table(unsolved_comphet_LOFmissense_sec[12])

unsolved_comphet_LOFmissense_sec_grch37 <- unsolved_comphet_LOFmissense_sec %>% filter(assembly == "GRCh37")
unsolved_comphet_LOFmissense_sec_grch37$base <- paste(as.character(unsolved_comphet_LOFmissense_sec_grch37$reference), "/", as.character(unsolved_comphet_LOFmissense_sec_grch37$alternate))
unsolved_comphet_LOFmissense_sec_grch37_VEP <- unsolved_comphet_LOFmissense_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_sec_grch37_VEP_list <- unsolved_comphet_LOFmissense_sec_grch37_VEP[!duplicated(unsolved_comphet_LOFmissense_sec_grch37_VEP[,2:4]),]
unsolved_comphet_LOFmissense_sec_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_sec_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_sec_grch37_VEP_list, "VEP_comphet_LOFmissense_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_LOFmissense_sec_grch38 <- unsolved_comphet_LOFmissense_sec %>% filter(assembly == "GRCh38")
unsolved_comphet_LOFmissense_sec_grch38$base <- paste(as.character(unsolved_comphet_LOFmissense_sec_grch38$reference), "/", as.character(unsolved_comphet_LOFmissense_sec_grch38$alternate))
unsolved_comphet_LOFmissense_sec_grch38_VEP <- unsolved_comphet_LOFmissense_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_sec_grch38_VEP_list <- unsolved_comphet_LOFmissense_sec_grch38_VEP[!duplicated(unsolved_comphet_LOFmissense_sec_grch38_VEP[,2:4]),]
unsolved_comphet_LOFmissense_sec_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_sec_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_sec_grch38_VEP_list, "VEP_comphet_LOFmissense_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# HSP GENES
table(unsolved_comphet_LOFmissense_hsp[12])

unsolved_comphet_LOFmissense_hsp_grch37 <- unsolved_comphet_LOFmissense_hsp %>% filter(assembly == "GRCh37")
unsolved_comphet_LOFmissense_hsp_grch37$base <- paste(as.character(unsolved_comphet_LOFmissense_hsp_grch37$reference), "/", as.character(unsolved_comphet_LOFmissense_hsp_grch37$alternate))
unsolved_comphet_LOFmissense_hsp_grch37_VEP <- unsolved_comphet_LOFmissense_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_hsp_grch37_VEP_list <- unsolved_comphet_LOFmissense_hsp_grch37_VEP[!duplicated(unsolved_comphet_LOFmissense_hsp_grch37_VEP[,2:4]),]
unsolved_comphet_LOFmissense_hsp_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_hsp_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_hsp_grch37_VEP_list, "VEP_comphet_LOFmissense_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_LOFmissense_hsp_grch38 <- unsolved_comphet_LOFmissense_hsp %>% filter(assembly == "GRCh38")
unsolved_comphet_LOFmissense_hsp_grch38$base <- paste(as.character(unsolved_comphet_LOFmissense_hsp_grch38$reference), "/", as.character(unsolved_comphet_LOFmissense_hsp_grch38$alternate))
unsolved_comphet_LOFmissense_hsp_grch38_VEP <- unsolved_comphet_LOFmissense_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_hsp_grch38_VEP_list <- unsolved_comphet_LOFmissense_hsp_grch38_VEP[!duplicated(unsolved_comphet_LOFmissense_hsp_grch38_VEP[,2:4]),]
unsolved_comphet_LOFmissense_hsp_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_hsp_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_hsp_grch38_VEP_list, "VEP_comphet_LOFmissense_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

###############
# MONOALLELIC #
###############

###################
# Affected father #
###################

# Filter affected parents (father)
affected_parents_father <- c((unsolved_tiering_data[22] == "Yes" & unsolved_tiering_data[23] == "No"))
unsolved_tiering_data_affected_father <- unsolved_tiering_data %>% filter(affected_parents_father)
table(unsolved_tiering_data_affected_father[22])
table(unsolved_tiering_data_affected_father[23])

# Filter for monoallelic (father)
table(unsolved_tiering_data_affected_father[14])
monoallelic_father <- filter(unsolved_tiering_data_affected_father, mode_of_inheritance %in% c("monoallelic_not_imprinted", "xlinked_monoallelic"))

# Father LOF
LOF_father <- filter(monoallelic_father, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant", "splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost"))

# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
LOF_father <- LOF_father[!duplicated(LOF_father[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
l_father<- LOF_father %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequence type
lo_father<- l_father %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
loss_father <- LOF_father[c(-6,-20)]
# Join 
LOF_monoallelic_father <- left_join(lo_father, loss_father)
# Remove duplicates
LOF_monoallelic_father <- LOF_monoallelic_father[!duplicated(LOF_monoallelic_father[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
LOF_monoallelic_father_sec <- filter(LOF_monoallelic_father, genomic_feature_hgnc %in% sec)
LOF_monoallelic_father_hsp <- filter(LOF_monoallelic_father, genomic_feature_hgnc %in% hsp)


# Father missense
missense_father <- filter(monoallelic_father, consequence_type %in% c("inframe_deletion","inframe_insertion","missense_variant", "stop_gained","stop_lost"))

# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
missense_father <- missense_father[!duplicated(missense_father[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
m_father<- missense_father %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequence type
mi_father<- m_father %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
miss_father <- missense_father[c(-6,-20)]
# Join 
missense_monoallelic_father <- left_join(mi_father, miss_father)
# Remove duplicates
missense_monoallelic_father <- missense_monoallelic_father[!duplicated(missense_monoallelic_father[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
missense_monoallelic_father_sec <- filter(missense_monoallelic_father, genomic_feature_hgnc %in% sec)
missense_monoallelic_father_hsp <- filter(missense_monoallelic_father, genomic_feature_hgnc %in% hsp)


###################
# Affected mother #
###################
# Filter affected parents (mother)
affected_parents_mother <- c((unsolved_tiering_data[22] == "No" & unsolved_tiering_data[23] == "Yes"))
unsolved_tiering_data_affected_mother <- unsolved_tiering_data %>% filter(affected_parents_mother)

# Filter for monoallelic (mother)
table(unsolved_tiering_data_affected_mother[14])
monoallelic_mother <- filter(unsolved_tiering_data_affected_mother, mode_of_inheritance %in% c("monoallelic_not_imprinted", "xlinked_monoalle"))
# Remove "deNovo" segregation pattern
monoallelic_mother <- filter(monoallelic_mother, monoallelic_mother[15] != "deNovo")

# Mother LOF
LOF_mother <- filter(monoallelic_mother, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant", "splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost"))

# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
LOF_mother <- LOF_mother[!duplicated(LOF_mother[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
l_mother <- LOF_mother %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequence type
lo_mother <- l_mother %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
loss_mother <- LOF_mother[c(-6,-20)]
# Join 
LOF_monoallelic_mother <- left_join(lo_mother, loss_mother)
# Remove duplicates
LOF_monoallelic_mother <- LOF_monoallelic_mother[!duplicated(LOF_monoallelic_mother[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
LOF_monoallelic_mother_sec <- filter(LOF_monoallelic_mother, genomic_feature_hgnc %in% sec)
LOF_monoallelic_mother_hsp <- filter(LOF_monoallelic_mother, genomic_feature_hgnc %in% hsp)


# Mother missense
missense_mother <- filter(monoallelic_mother, consequence_type %in% c("inframe_deletion","inframe_insertion","missense_variant", "stop_gained","stop_lost"))

# Remove duplicates based on participant_id, position, reference, alternate,phenotype and consequence type
missense_mother <- missense_mother[!duplicated(missense_mother[c(1,6,10,11,12,20)]),]
# Join consequence type based on participant_id, position, reference, alternate and phenotype
m_mother <- missense_mother %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))
# Join phenotype based on participant_id, position, reference, alternate and consequence type
mi_mother <- m_mother %>%
  group_by(participant_id, position, reference, alternate, consequence_type) %>%
  summarise(phenotype = paste(phenotype, collapse = " | "))
# Remove consequence type and phenotype columns
miss_mother <- missense_mother[c(-6,-20)]
# Join 
missense_monoallelic_mother <- left_join(mi_mother, miss_mother)
# Remove duplicates
missense_monoallelic_mother <- missense_monoallelic_mother[!duplicated(missense_monoallelic_mother[c(1,2,3,4,5,6)]),]
#Filter for sec and hsp genes
missense_monoallelic_mother_sec <- filter(missense_monoallelic_mother, genomic_feature_hgnc %in% sec)
missense_monoallelic_mother_hsp <- filter(missense_monoallelic_mother, genomic_feature_hgnc %in% hsp)


###################
# VEP MONOALLELIC #
###################

###################
# Affected father #
###################

# SEC GENES
#Loss of function
table(LOF_monoallelic_father_sec[12])
LOF_monoallelic_father_sec_grch37 <- LOF_monoallelic_father_sec %>% filter(assembly == "GRCh37")
LOF_monoallelic_father_sec_grch37$base <- paste(as.character(LOF_monoallelic_father_sec_grch37$reference), "/", as.character(LOF_monoallelic_father_sec_grch37$alternate))
LOF_monoallelic_father_sec_grch37_VEP <- LOF_monoallelic_father_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_father_sec_grch37_VEP_list <- LOF_monoallelic_father_sec_grch37_VEP[!duplicated(LOF_monoallelic_father_sec_grch37_VEP[,2:4]),]
LOF_monoallelic_father_sec_grch37_VEP_list$base <- str_replace_all(LOF_monoallelic_father_sec_grch37_VEP_list$base, " ", "")
write.table(LOF_monoallelic_father_sec_grch37_VEP_list, "VEP_LOF_monoallelic_father_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

LOF_monoallelic_father_sec_grch38 <- LOF_monoallelic_father_sec %>% filter(assembly == "GRCh38")
LOF_monoallelic_father_sec_grch38$base <- paste(as.character(LOF_monoallelic_father_sec_grch38$reference), "/", as.character(LOF_monoallelic_father_sec_grch38$alternate))
LOF_monoallelic_father_sec_grch38_VEP <- LOF_monoallelic_father_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_father_sec_grch38_VEP_list <- LOF_monoallelic_father_sec_grch38_VEP[!duplicated(LOF_monoallelic_father_sec_grch38_VEP[,2:4]),]
LOF_monoallelic_father_sec_grch38_VEP_list$base <- str_replace_all(LOF_monoallelic_father_sec_grch38_VEP_list$base, " ", "")
write.table(LOF_monoallelic_father_sec_grch38_VEP_list, "VEP_LOF_monoallelic_father_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Missense
table(missense_monoallelic_father_sec[12])
missense_monoallelic_father_sec_grch37 <- missense_monoallelic_father_sec %>% filter(assembly == "GRCh37")
missense_monoallelic_father_sec_grch37$base <- paste(as.character(missense_monoallelic_father_sec_grch37$reference), "/", as.character(missense_monoallelic_father_sec_grch37$alternate))
missense_monoallelic_father_sec_grch37_VEP <- missense_monoallelic_father_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_father_sec_grch37_VEP_list <- missense_monoallelic_father_sec_grch37_VEP[!duplicated(missense_monoallelic_father_sec_grch37_VEP[,2:4]),]
missense_monoallelic_father_sec_grch37_VEP_list$base <- str_replace_all(missense_monoallelic_father_sec_grch37_VEP_list$base, " ", "")
write.table(missense_monoallelic_father_sec_grch37_VEP_list, "VEP_missense_monoallelic_father_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

missense_monoallelic_father_sec_grch38 <- missense_monoallelic_father_sec_sec %>% filter(assembly == "GRCh38")
missense_monoallelic_father_sec_grch38$base <- paste(as.character(missense_monoallelic_father_sec_grch38$reference), "/", as.character(missense_monoallelic_father_sec_grch38$alternate))
missense_monoallelic_father_sec_grch38_VEP <- missense_monoallelic_father_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_father_sec_grch38_VEP_list <- missense_monoallelic_father_sec_grch38_VEP[!duplicated(missense_monoallelic_father_sec_grch38_VEP[,2:4]),]
missense_monoallelic_father_sec_grch38_VEP_list$base <- str_replace_all(missense_monoallelic_father_sec_grch38_VEP_list$base, " ", "")
write.table(missense_monoallelic_father_sec_grch38_VEP_list, "VEP_missense_monoallelic_father_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#HSP GENES
table(LOF_monoallelic_father_hsp[12])
LOF_monoallelic_father_hsp_grch37 <- LOF_monoallelic_father_hsp %>% filter(assembly == "GRCh37")
LOF_monoallelic_father_hsp_grch37$base <- paste(as.character(LOF_monoallelic_father_hsp_grch37$reference), "/", as.character(LOF_monoallelic_father_hsp_grch37$alternate))
LOF_monoallelic_father_hsp_grch37_VEP <- LOF_monoallelic_father_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_father_hsp_grch37_VEP_list <- LOF_monoallelic_father_hsp_grch37_VEP[!duplicated(LOF_monoallelic_father_hsp_grch37_VEP[,2:4]),]
LOF_monoallelic_father_hsp_grch37_VEP_list$base <- str_replace_all(LOF_monoallelic_father_hsp_grch37_VEP_list$base, " ", "")
write.table(LOF_monoallelic_father_hsp_grch37_VEP_list, "VEP_LOF_monoallelic_father_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

LOF_monoallelic_father_hsp_grch38 <- LOF_monoallelic_father_hsp %>% filter(assembly == "GRCh38")
LOF_monoallelic_father_hsp_grch38$base <- paste(as.character(LOF_monoallelic_father_hsp_grch38$reference), "/", as.character(LOF_monoallelic_father_hsp_grch38$alternate))
LOF_monoallelic_father_hsp_grch38_VEP <- LOF_monoallelic_father_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_father_hsp_grch38_VEP_list <- LOF_monoallelic_father_hsp_grch38_VEP[!duplicated(LOF_monoallelic_father_hsp_grch38_VEP[,2:4]),]
LOF_monoallelic_father_hsp_grch38_VEP_list$base <- str_replace_all(LOF_monoallelic_father_hsp_grch38_VEP_list$base, " ", "")
write.table(LOF_monoallelic_father_hsp_grch38_VEP_list, "VEP_LOF_monoallelic_father_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Missense
table(missense_monoallelic_father_hsp[12])
missense_monoallelic_father_hsp_grch37 <- missense_monoallelic_father_hsp %>% filter(assembly == "GRCh37")
missense_monoallelic_father_hsp_grch37$base <- paste(as.character(missense_monoallelic_father_hsp_grch37$reference), "/", as.character(missense_monoallelic_father_hsp_grch37$alternate))
missense_monoallelic_father_hsp_grch37_VEP <- missense_monoallelic_father_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_father_hsp_grch37_VEP_list <- missense_monoallelic_father_hsp_grch37_VEP[!duplicated(missense_monoallelic_father_hsp_grch37_VEP[,2:4]),]
missense_monoallelic_father_hsp_grch37_VEP_list$base <- str_replace_all(missense_monoallelic_father_hsp_grch37_VEP_list$base, " ", "")
write.table(missense_monoallelic_father_hsp_grch37_VEP_list, "VEP_missense_monoallelic_father_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

missense_monoallelic_father_hsp_grch38 <- missense_monoallelic_father_hsp %>% filter(assembly == "GRCh38")
missense_monoallelic_father_hsp_grch38$base <- paste(as.character(missense_monoallelic_father_hsp_grch38$reference), "/", as.character(missense_monoallelic_father_hsp_grch38$alternate))
missense_monoallelic_father_hsp_grch38_VEP <- missense_monoallelic_father_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_father_hsp_grch38_VEP_list <- missense_monoallelic_father_hsp_grch38_VEP[!duplicated(missense_monoallelic_father_hsp_grch38_VEP[,2:4]),]
missense_monoallelic_father_hsp_grch38_VEP_list$base <- str_replace_all(missense_monoallelic_father_hsp_grch38_VEP_list$base, " ", "")
write.table(missense_monoallelic_father_hsp_grch38_VEP_list, "VEP_missense_monoallelic_father_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


###################
# Affected mother #
###################

# SEC GENES
#Loss of function
table(LOF_monoallelic_mother_sec[12])
LOF_monoallelic_mother_sec_grch37 <- LOF_monoallelic_mother_sec %>% filter(assembly == "GRCh37")
LOF_monoallelic_mother_sec_grch37$base <- paste(as.character(LOF_monoallelic_mother_sec_grch37$reference), "/", as.character(LOF_monoallelic_mother_sec_grch37$alternate))
LOF_monoallelic_mother_sec_grch37_VEP <- LOF_monoallelic_mother_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_mother_sec_grch37_VEP_list <- LOF_monoallelic_mother_sec_grch37_VEP[!duplicated(LOF_monoallelic_mother_sec_grch37_VEP[,2:4]),]
LOF_monoallelic_mother_sec_grch37_VEP_list$base <- str_replace_all(LOF_monoallelic_mother_sec_grch37_VEP_list$base, " ", "")
write.table(LOF_monoallelic_mother_sec_grch37_VEP_list, "VEP_LOF_monoallelic_mother_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

LOF_monoallelic_mother_sec_grch38 <- LOF_monoallelic_mother_sec %>% filter(assembly == "GRCh38")
LOF_monoallelic_mother_sec_grch38$base <- paste(as.character(LOF_monoallelic_mother_sec_grch38$reference), "/", as.character(LOF_monoallelic_mother_sec_grch38$alternate))
LOF_monoallelic_mother_sec_grch38_VEP <- LOF_monoallelic_mother_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_mother_sec_grch38_VEP_list <- LOF_monoallelic_mother_sec_grch38_VEP[!duplicated(LOF_monoallelic_mother_sec_grch38_VEP[,2:4]),]
LOF_monoallelic_mother_sec_grch38_VEP_list$base <- str_replace_all(LOF_monoallelic_mother_sec_grch38_VEP_list$base, " ", "")
write.table(LOF_monoallelic_mother_sec_grch38_VEP_list, "VEP_LOF_monoallelic_mother_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Missense
table(missense_monoallelic_mother_sec[12])
missense_monoallelic_mother_sec_grch37 <- missense_monoallelic_mother_sec %>% filter(assembly == "GRCh37")
missense_monoallelic_mother_sec_grch37$base <- paste(as.character(missense_monoallelic_mother_sec_grch37$reference), "/", as.character(missense_monoallelic_mother_sec_grch37$alternate))
missense_monoallelic_mother_sec_grch37_VEP <- missense_monoallelic_mother_sec_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_mother_sec_grch37_VEP_list <- missense_monoallelic_mother_sec_grch37_VEP[!duplicated(missense_monoallelic_mother_sec_grch37_VEP[,2:4]),]
missense_monoallelic_mother_sec_grch37_VEP_list$base <- str_replace_all(missense_monoallelic_mother_sec_grch37_VEP_list$base, " ", "")
write.table(missense_monoallelic_mother_sec_grch37_VEP_list, "VEP_missense_monoallelic_mother_sec_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

missense_monoallelic_mother_sec_grch38 <- missense_monoallelic_mother_sec_sec %>% filter(assembly == "GRCh38")
missense_monoallelic_mother_sec_grch38$base <- paste(as.character(missense_monoallelic_mother_sec_grch38$reference), "/", as.character(missense_monoallelic_mother_sec_grch38$alternate))
missense_monoallelic_mother_sec_grch38_VEP <- missense_monoallelic_mother_sec_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_mother_sec_grch38_VEP_list <- missense_monoallelic_mother_sec_grch38_VEP[!duplicated(missense_monoallelic_mother_sec_grch38_VEP[,2:4]),]
missense_monoallelic_mother_sec_grch38_VEP_list$base <- str_replace_all(missense_monoallelic_mother_sec_grch38_VEP_list$base, " ", "")
write.table(missense_monoallelic_father_sec_grch38_VEP_list, "VEP_missense_monoallelic_mother_sec_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#HSP GENES
table(LOF_monoallelic_mother_hsp[12])
LOF_monoallelic_mother_hsp_grch37 <- LOF_monoallelic_mother_hsp %>% filter(assembly == "GRCh37")
LOF_monoallelic_mother_hsp_grch37$base <- paste(as.character(LOF_monoallelic_mother_hsp_grch37$reference), "/", as.character(LOF_monoallelic_mother_hsp_grch37$alternate))
LOF_monoallelic_mother_hsp_grch37_VEP <- LOF_monoallelic_mother_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_mother_hsp_grch37_VEP_list <- LOF_monoallelic_mother_hsp_grch37_VEP[!duplicated(LOF_monoallelic_mother_hsp_grch37_VEP[,2:4]),]
LOF_monoallelic_mother_hsp_grch37_VEP_list$base <- str_replace_all(LOF_monoallelic_mother_hsp_grch37_VEP_list$base, " ", "")
write.table(LOF_monoallelic_mother_hsp_grch37_VEP_list, "VEP_LOF_monoallelic_mother_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

LOF_monoallelic_mother_hsp_grch38 <- LOF_monoallelic_mother_hsp %>% filter(assembly == "GRCh38")
LOF_monoallelic_mother_hsp_grch38$base <- paste(as.character(LOF_monoallelic_mother_hsp_grch38$reference), "/", as.character(LOF_monoallelic_mother_hsp_grch38$alternate))
LOF_monoallelic_mother_hsp_grch38_VEP <- LOF_monoallelic_mother_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
LOF_monoallelic_mother_hsp_grch38_VEP_list <- LOF_monoallelic_mother_hsp_grch38_VEP[!duplicated(LOF_monoallelic_mother_hsp_grch38_VEP[,2:4]),]
LOF_monoallelic_mother_hsp_grch38_VEP_list$base <- str_replace_all(LOF_monoallelic_mother_hsp_grch38_VEP_list$base, " ", "")
write.table(LOF_monoallelic_mother_hsp_grch38_VEP_list, "VEP_LOF_monoallelic_mother_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Missense
table(missense_monoallelic_mother_hsp[12])
missense_monoallelic_mother_hsp_grch37 <- missense_monoallelic_mother_hsp %>% filter(assembly == "GRCh37")
missense_monoallelic_mother_hsp_grch37$base <- paste(as.character(missense_monoallelic_mother_hsp_grch37$reference), "/", as.character(missense_monoallelic_mother_hsp_grch37$alternate))
missense_monoallelic_mother_hsp_grch37_VEP <- missense_monoallelic_mother_hsp_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_mother_hsp_grch37_VEP_list <- missense_monoallelic_mother_hsp_grch37_VEP[!duplicated(missense_monoallelic_mother_hsp_grch37_VEP[,2:4]),]
missense_monoallelic_mother_hsp_grch37_VEP_list$base <- str_replace_all(missense_monoallelic_mother_hsp_grch37_VEP_list$base, " ", "")
write.table(missense_monoallelic_mother_hsp_grch37_VEP_list, "VEP_missense_monoallelic_mother_hsp_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

missense_monoallelic_mother_hsp_grch38 <- missense_monoallelic_mother_hsp %>% filter(assembly == "GRCh38")
missense_monoallelic_mother_hsp_grch38$base <- paste(as.character(missense_monoallelic_mother_hsp_grch38$reference), "/", as.character(missense_monoallelic_mother_hsp_grch38$alternate))
missense_monoallelic_mother_hsp_grch38_VEP <- missense_monoallelic_mother_hsp_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
missense_monoallelic_mother_hsp_grch38_VEP_list <- missense_monoallelic_mother_hsp_grch38_VEP[!duplicated(missense_monoallelic_mother_hsp_grch38_VEP[,2:4]),]
missense_monoallelic_mother_hsp_grch38_VEP_list$base <- str_replace_all(missense_monoallelic_mother_hsp_grch38_VEP_list$base, " ", "")
write.table(missense_monoallelic_mother_hsp_grch38_VEP_list, "VEP_missense_monoallelic_mother_hsp_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
