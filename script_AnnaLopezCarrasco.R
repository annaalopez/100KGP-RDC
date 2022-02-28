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
#160536 observations and 26 variables

#Semi join tiering_data and gmc_unsolved
unsolved_tiering_data <- semi_join (tiering_data, gmc_unsolved, by = "participant_id")
#29696 observations

#Filter for unaffected parents ("No" only)
table(unsolved_tiering_data[22])
table(unsolved_tiering_data[23])
unsolved_tiering_data_unaffected <- unsolved_tiering_data %>% filter(unsolved_tiering_data[22] == "No" & unsolved_tiering_data[23] == "No") 
#23716 observations

#Filter for LOF
LOF <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant","splice_acceptor_variant","splice_donor_variant", "splive_region_variant","start_lost"))
#3981

#Simple recessive LOF
unsolved_recessive_lof <- (LOF %>% filter(LOF[15] == "SimpleRecessive"))
table(unsolved_recessive_lof[20])
#120 0bservations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
recessive <- unsolved_recessive_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
recessive_lof <- unsolved_recessive_lof[!duplicated(unsolved_recessive_lof[c(1,6,10,11,12)]),]

#Remove consequence type column
recessive_lof <- recessive_lof[-20]

#Join
unsolved_comphet_lof <- merge(recessive, recessive_lof)
write.csv(unsolved_recessive_lof, "unsolved_recessive_lof.csv")
#93 observations

#DeNovo LOF
unsolved_denovo_lof <- (LOF %>% filter(LOF[15] == "deNovo"))
table(unsolved_denovo_lof[20])
#72 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
denovo <- unsolved_denovo_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
denovo_lof <- unsolved_denovo_lof[!duplicated(unsolved_denovo_lof[c(1,6,10,11,12)]),]

#Remove phenotype column
denovo_lof <- comphet_lof[-20]

#Join
unsolved_denovo_lof <- merge(denovo, denovo_lof)
write.csv(unsolved_denovo_lof, "unsolved_denovo_lof.csv")
#57 observations

#CompoundHeterozygous LOF
unsolved_comphet_lof <- (LOF %>% filter(LOF[15] == "CompoundHeterozygous"))
#

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
comphet <- unsolved_comphet_lof %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
comphet_lof <- unsolved_comphet_lof[!duplicated(unsolved_comphet_lof[c(1,6,10,11,12)]),]
#Remove phenotype column
comphet_lof <- comphet_lof[-20]

#Join
unsolved_comphet_lof <- merge(comphet, comphet_lof)

#Remove participant_id that only appears once
unsolved_comphet_lof <- unsolved_comphet_lof %>%
  group_by(participant_id) %>%
  filter(n()>1)
write.csv(unsolved_comphet_lof, "unsolved_comphet_lof.csv")
#205 observations

#Filter for missense+
table(unsolved_tiering_data_unaffected[20])
missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("inframe_deletion", "inframe_insertion", "missense_Variant", "start_gained", "stop_lost"))

#Simple recessive missense+
unsolved_recessive_missense <- (missense %>% filter(missense[15] == "SimpleRecessive"))
#172 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
recessive_miss <- unsolved_recessive_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
recessive_missense <- unsolved_recessive_missense[!duplicated(unsolved_recessive_missense[c(1,6,10,11,12)]),]

#Remove consequence type column
recessive_missense <- recessive_missense[-20]

#Join
unsolved_comphet_missense <- merge(recessive_miss, recessive_missense)
write.csv(unsolved_recessive_missense, "unsolved_recessive_missense.csv")
#166 observations

#DeNovo missense+ 
unsolved_denovo_missense <- (missense %>% filter(missense[15] == "deNovo"))
#109 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
denovo_miss <- unsolved_denovo_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
denovo_missense<- unsolved_denovo_missense[!duplicated(unsolved_denovo_missense[c(1,6,10,11,12)]),]

#Remove consequence type column
denovo_missense <- denovo_missense[-20]

#Join
unsolved_denovo_missense <- merge(denovo_miss, denovo_missense)
write.csv(unsolved_denovo_missense, "unsolved_denovo_missense.csv")
#

#CompoundHeterozygous missense+
unsolved_comphet_missense <- (missense %>% filter(missense[15] == "CompoundHeterozygous"))
#568 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
comphet_miss <- unsolved_comphet_missense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
comphet_missense <- unsolved_comphet_missense[!duplicated(unsolved_comphet_missense[c(1,6,10,11,12)]),]
#Remove phenotype column
comphet_missense <- comphet_missense[-20]

#Join
unsolved_comphet_missense <- merge(comphet_miss, comphet_missense)

#Remove participant_id that only appears once
unsolved_comphet_missense <- unsolved_comphet_missense %>%
  group_by(participant_id) %>%
  filter(n()>1)
write.csv(unsolved_comphet_missense, "unsolved_comphet_missense.csv")
#474 observations

#Filter for LOF + missense
LOF_missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant","splice_acceptor_variant","splice_donor_variant", "splive_region_variant","start_lost", "missense_variant"))

#Simple recessive LOF + missense
unsolved_recessive_LOFmissense <- (LOF_missense %>% filter(LOF_missense[15] == 'SimpleRecessive'))
#292 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
recessive_LOFmiss <- unsolved_recessive_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
recessive_LOFmissense <- unsolved_recessive_LOFmissense[!duplicated(unsolved_recessive_LOFmissense[c(1,6,10,11,12)]),]
#Remove phenotype column
recessive_LOFmissense <- recessive_LOFmissense[-20]

#Join
unsolved_recessive_LOFmissense <- merge(recessive_LOFmiss, recessive_LOFmissense)
write.csv(unsolved_recessive_LOFmissense, "unsolved_recessive_LOFmissense.csv")
#221 observations

#deNovo LOF + missense
unsolved_denovo_LOFmissense <- (LOF_missense %>% filter(LOF_missense[15] == 'deNovo'))
#

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
denovo_LOFmiss <- unsolved_denovo_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
denovo_LOFmissense <- unsolved_denovo_LOFmissense[!duplicated(unsolved_denovo_LOFmissense[c(1,6,10,11,12)]),]

#Remove phenotype column
denovo_LOFmissense <- denovo_LOFmissense[-20]

#Join
unsolved_denovo_LOFmissense <- merge(denovo_LOFmiss, denovo_LOFmissense)
write.csv(unsolved_denovo_LOFmissense, "unsolved_denovo_LOFmissense.csv")
#1101 observations

#CompoundHeterozygous LOF + missense
unsolved_comphet_LOFmissense <- (LOF_missense %>% filter(LOF_missense[15] == "CompoundHeterozygous"))
#568 observations

#Join consequence_type names based on participant_id, position, reference, alternate and phenotype
comphet_LOFmiss <- unsolved_comphet_LOFmissense %>%
  group_by(participant_id, position, reference, alternate, phenotype) %>%
  summarise(consequence_type = paste(consequence_type, collapse = " | "))

#Remove duplicates based on participant_id, position, reference and alternate
comphet_LOFmissense <- unsolved_comphet_LOFmissense[!duplicated(unsolved_comphet_LOFmissense[c(1,6,10,11,12)]),]
#Remove phenotype column
comphet_LOFmissense <- comphet_LOFmissense[-20]

#Join
unsolved_comphet_LOFmissense <- merge(comphet_LOFmiss, comphet_LOFmissense)

#Remove participant_id that only appears once
unsolved_comphet_LOFmissense <- unsolved_comphet_LOFmissense %>%
  group_by(participant_id) %>%
  filter(n()>1)
write.csv(unsolved_comphet_LOFmissense, "unsolved_comphet_LOFmissense.csv")
#474 observations

#########
#  VEP  #
#########

library(stringr)

#Simple recessive LOF
table(unsolved_recessive_lof[12])

unsolved_recessive_lof_grch37 <- filter(unsolved_recessive_lof, unsolved_recessive_lof[8] == "GRCh37")
unsolved_recessive_lof_grch37$base <- paste(as.character(unsolved_recessive_lof_grch37$reference), "/", as.character(unsolved_recessive_lof_grch37$alternate))
unsolved_recessive_lof_grch37_VEP <- unsolved_recessive_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_grch37_VEP_list <- unsolved_recessive_lof_grch37_VEP[!duplicated(unsolved_recessive_lof_grch37_VEP[,2:4]),]
unsolved_recessive_lof_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_grch37_VEP_list, "VEP_recessive_lof_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_recessive_lof_grch38 <- filter(unsolved_recessive_lof, unsolved_recessive_lof[8] == "GRCh38")
unsolved_recessive_lof_grch38$base <- paste(as.character(unsolved_recessive_lof_grch38$reference), "/", as.character(unsolved_recessive_lof_grch38$alternate))
unsolved_recessive_lof_grch38_VEP <- unsolved_recessive_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_grch38_VEP_list <- unsolved_recessive_lof_grch38_VEP[!duplicated(unsolved_recessive_lof_grch38_VEP[,2:4]),]
unsolved_recessive_lof_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_grch38_VEP_list, "VEP_recessive_lof_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#Simple recessive missense+
table(unsolved_recessive_missense[12])

unsolved_recessive_missense_grch37 <- filter(unsolved_recessive_missense, unsolved_recessive_missense[8] == "GRCh37")
unsolved_recessive_missense_grch37$base <- paste(as.character(unsolved_recessive_missense_grch37$reference), "/", as.character(unsolved_recessive_missense_grch37$alternate))
unsolved_recessive_missense_grch37_VEP <- unsolved_recessive_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_grch37_VEP_list <- unsolved_recessive_missense_grch37_VEP[!duplicated(unsolved_recessive_missense_grch37_VEP[,2:4]),]
unsolved_recessive_missense_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_grch37_VEP_list, "VEP_recessive_missense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_recessive_missense_grch38 <- filter(unsolved_recessive_missense, unsolved_recessive_missense[8] == "GRCh38")
unsolved_recessive_missense_grch38$base <- paste(as.character(unsolved_recessive_missense_grch38$reference), "/", as.character(unsolved_recessive_missense_grch38$alternate))
unsolved_recessive_missense_grch38_VEP <- unsolved_recessive_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_grch38_VEP_list <- unsolved_recessive_missense_grch38_VEP[!duplicated(unsolved_recessive_missense_grch38_VEP[,2:4]),]
unsolved_recessive_missense_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_grch38_VEP_list, "VEP_recessive_missense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#Simple recessive LOF + missense
table(unsolved_recessive_LOFmissense[12])

unsolved_recessive_LOFmissense_grch37 <- filter(unsolved_recessive_LOFmissense, unsolved_recessive_LOFmissense[8] == "GRCh37")
unsolved_recessive_LOFmissense_grch37$base <- paste(as.character(unsolved_recessive_LOFmissense_grch37$reference), "/", as.character(unsolved_recessive_LOFmissense_grch37$alternate))
unsolved_recessive_LOFmissense_grch37_VEP <- unsolved_recessive_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_grch37_VEP_list <- unsolved_recessive_LOFmissense_grch37_VEP[!duplicated(unsolved_recessive_LOFmissense_grch37_VEP[,2:4]),]
unsolved_recessive_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_grch37_VEP_list, "VEP_recessive_LOFmissense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_recessive_LOFmissense_grch38 <- filter(unsolved_recessive_LOFmissense, unsolved_recessive_LOFmissense[8] == "GRCh38")
unsolved_recessive_LOFmissense_grch38$base <- paste(as.character(unsolved_recessive_LOFmissense_grch38$reference), "/", as.character(unsolved_recessive_LOFmissense_grch38$alternate))
unsolved_recessive_LOFmissense_grch38_VEP <- unsolved_recessive_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_grch38_VEP_list <- unsolved_recessive_LOFmissense_grch38_VEP[!duplicated(unsolved_recessive_LOFmissense_grch38_VEP[,2:4]),]
unsolved_recessive_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_grch38_VEP_list, "VEP_recessive_LOFmissense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#DeNovo LOF
table(unsolved_denovo_lof[12])

unsolved_denovo_lof_grch37 <- filter(unsolved_denovo_lof, unsolved_denovo_lof[8] == "GRCh37")
unsolved_denovo_lof_grch37$base <- paste(as.character(unsolved_denovo_lof_grch37$reference), "/", as.character(unsolved_denovo_lof_grch37$alternate))
unsolved_denovo_lof_grch37_VEP <- unsolved_denovo_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_grch37_VEP_list <- unsolved_denovo_lof_grch37_VEP[!duplicated(unsolved_denovo_lof_grch37_VEP[,2:4]),]
unsolved_denovo_lof_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_grch37_VEP_list, "VEP_denovo_lof_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_denovo_lof_grch38 <- filter(unsolved_denovo_lof, unsolved_denovo_lof[8] == "GRCh38")
unsolved_denovo_lof_grch38$base <- paste(as.character(unsolved_denovo_lof_grch38$reference), "/", as.character(unsolved_denovo_lof_grch38$alternate))
unsolved_denovo_lof_grch38_VEP <- unsolved_denovo_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_grch38_VEP_list <- unsolved_denovo_lof_grch38_VEP[!duplicated(unsolved_denovo_lof_grch38_VEP[,2:4]),]
unsolved_denovo_lof_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_grch38_VEP_list, "VEP_denovo_lof_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#DeNovo missense+
table(unsolved_denovo_missense[12])

unsolved_denovo_missense_grch37 <- filter(unsolved_denovo_missense, unsolved_denovo_missense[8] == "GRCh37")
unsolved_denovo_missense_grch37$base <- paste(as.character(unsolved_denovo_missense_grch37$reference), "/", as.character(unsolved_denovo_missense_grch37$alternate))
unsolved_denovo_missense_grch37_VEP <- unsolved_denovo_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_grch37_VEP_list <- unsolved_denovo_missense_grch37_VEP[!duplicated(unsolved_denovo_missense_grch37_VEP[,2:4]),]
unsolved_denovo_missense_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_grch37_VEP_list, "VEP_denovo_missense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_denovo_missense_grch38 <- filter(unsolved_denovo_missense, unsolved_denovo_missense[8] == "GRCh38")
unsolved_denovo_missense_grch38$base <- paste(as.character(unsolved_denovo_missense_grch38$reference), "/", as.character(unsolved_denovo_missense_grch38$alternate))
unsolved_denovo_missense_grch38_VEP <- unsolved_denovo_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_grch38_VEP_list <- unsolved_denovo_missense_grch38_VEP[!duplicated(unsolved_denovo_missense_grch38_VEP[,2:4]),]
unsolved_denovo_missense_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_grch38_VEP_list, "VEP_denovo_missense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#DeNovo LOF + missense
table(unsolved_denovo_LOFmissense[12])

unsolved_denovo_LOFmissense_grch37 <- filter(unsolved_denovo_LOFmissense, unsolved_denovo_LOFmissense[8] == "GRCh37")
unsolved_denovo_LOFmissense_grch37$base <- paste(as.character(unsolved_denovo_LOFmissense_grch37$reference), "/", as.character(unsolved_denovo_LOFmissense_grch37$alternate))
unsolved_denovo_LOFmissense_grch37_VEP <- unsolved_denovo_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_grch37_VEP_list <- unsolved_denovo_LOFmissense_grch37_VEP[!duplicated(unsolved_denovo_LOFmissense_grch37_VEP[,2:4]),]
unsolved_denovo_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_grch37_VEP_list, "VEP_denovo_LOFmissense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_denovo_LOFmissense_grch38 <- filter(unsolved_denovo_LOFmissense, unsolved_denovo_LOFmissense[8] == "GRCh38")
unsolved_denovo_LOFmissense_grch38$base <- paste(as.character(unsolved_denovo_LOFmissense_grch38$reference), "/", as.character(unsolved_denovo_LOFmissense_grch38$alternate))
unsolved_denovo_LOFmissense_grch38_VEP <- unsolved_denovo_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_grch38_VEP_list <- unsolved_denovo_LOFmissense_grch38_VEP[!duplicated(unsolved_denovo_LOFmissense_grch38_VEP[,2:4]),]
unsolved_denovo_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_grch38_VEP_list, "VEP_denovo_LOFmissense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#CompoundHeterozygous LOF
table(unsolved_comphet_lof[12])

unsolved_comphet_lof_grch37 <- filter(unsolved_comphet_lof, unsolved_comphet_lof[8] == "GRCh37")
unsolved_comphet_lof_grch37$base <- paste(as.character(unsolved_comphet_lof_grch37$reference), "/", as.character(unsolved_comphet_lof_grch37$alternate))
unsolved_comphet_lof_grch37_VEP <- unsolved_comphet_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_grch37_VEP_list <- unsolved_comphet_lof_grch37_VEP[!duplicated(unsolved_comphet_lof_grch37_VEP[,2:4]),]
unsolved_comphet_lof_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_grch37_VEP_list, "VEP_comphet_lof_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_comphet_lof_grch38 <- filter(unsolved_comphet_lof, unsolved_comphet_lof[8] == "GRCh38")
unsolved_comphet_lof_grch38$base <- paste(as.character(unsolved_comphet_lof_grch38$reference), "/", as.character(unsolved_comphet_lof_grch38$alternate))
unsolved_comphet_lof_grch38_VEP <- unsolved_comphet_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_grch38_VEP_list <- unsolved_comphet_lof_grch38_VEP[!duplicated(unsolved_comphet_lof_grch38_VEP[,2:4]),]
unsolved_comphet_lof_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_grch38_VEP_list, "VEP_comphet_lof_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#CompoundHeterozygous missense+
table(unsolved_comphet_missense[12])

unsolved_comphet_missense_grch37 <- filter(unsolved_comphet_missense, unsolved_comphet_missense[8] == "GRCh37")
unsolved_comphet_missense_grch37$base <- paste(as.character(unsolved_comphet_missense_grch37$reference), "/", as.character(unsolved_comphet_missense_grch37$alternate))
unsolved_comphet_missense_grch37_VEP <- unsolved_comphet_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_grch37_VEP_list <- unsolved_comphet_missense_grch37_VEP[!duplicated(unsolved_comphet_missense_grch37_VEP[,2:4]),]
unsolved_comphet_missense_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_grch37_VEP_list, "VEP_comphet_missense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_comphet_missense_grch38 <- filter(unsolved_comphet_missense, unsolved_comphet_missense[8] == "GRCh38")
unsolved_comphet_missense_grch38$base <- paste(as.character(unsolved_comphet_missense_grch38$reference), "/", as.character(unsolved_comphet_missense_grch38$alternate))
unsolved_comphet_missense_grch38_VEP <- unsolved_comphet_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_grch38_VEP_list <- unsolved_comphet_missense_grch38_VEP[!duplicated(unsolved_comphet_missense_grch38_VEP[,2:4]),]
unsolved_comphet_missense_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_grch38_VEP_list, "VEP_comphet_missense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#CompoundHeterozygous LOF + missense
table(unsolved_comphet_LOFmissense[12])

unsolved_comphet_LOFmissense_grch37 <- filter(unsolved_comphet_LOFmissense, unsolved_comphet_LOFmissense[8] == "GRCh37")
unsolved_comphet_LOFmissense_grch37$base <- paste(as.character(unsolved_comphet_LOFmissense_grch37$reference), "/", as.character(unsolved_comphet_LOFmissense_grch37$alternate))
unsolved_comphet_LOFmissense_grch37_VEP <- unsolved_comphet_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_grch37_VEP_list <- unsolved_comphet_LOFmissense_grch37_VEP[!duplicated(unsolved_comphet_LOFmissense_grch37_VEP[,2:4]),]
unsolved_comphet_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_grch37_VEP_list, "VEP_comphet_LOFmissense_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

unsolved_comphet_LOFmissense_grch38 <- filter(unsolved_comphet_LOFmissense, unsolved_comphet_LOFmissense[8] == "GRCh38")
unsolved_comphet_LOFmissense_grch38$base <- paste(as.character(unsolved_comphet_LOFmissense_grch38$reference), "/", as.character(unsolved_comphet_LOFmissense_grch38$alternate))
unsolved_comphet_LOFmissense_grch38_VEP <- unsolved_comphet_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_grch38_VEP_list <- unsolved_comphet_LOFmissense_grch38_VEP[!duplicated(unsolved_comphet_LOFmissense_grch38_VEP[,2:4]),]
unsolved_comphet_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_grch38_VEP_list, "VEP_comphet_LOFmissense_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = c("Chromosome", "Start", "End", "Base", "Strand"))

#Filter affected parents (only father)
affected_parents_father <- c((unsolved_tiering_data[22] == "Yes" & unsolved_tiering_data[23] == "No"))
unsolved_tiering_data_affected_father <- unsolved_tiering_data %>% filter(affected_parents_father)
table(unsolved_tiering_data_affected_father[22])
table(unsolved_tiering_data_affected_father[23])

#Filter for monoallelic
table(unsolved_tiering_data_affected_father[14])
monoallelic_father <- filer(unsolved_tiering_data_affected_father, mode_of_inheritance %in% c("monoallelic_not_imprented", "xlinked_monoallelic"))
table(monoallelic_father[14])
table(monoallelic_father[15])
table(monoallelic_father[20])

#Filter affected parents (only mother)
affected_parents_mother <- c((unsolved_tiering_data[23] == "Yes" & unsolved_tiering_data[22] == "No"))
unsolved_tiering_data_affected_mother <- unsolved_tiering_data %>% filter(affected_parents_mother)
table(unsolved_tiering_data_affected_mother[22])
table(unsolved_tiering_data_affected_mother[23])

#Filter for monoallelic
table(unsolved_tiering_data_affected_mother[14])
monoallelic_mother <-  filer(unsolved_tiering_data_affected_mother, mode_of_inheritance %in% c("monoallelic_not_imprented", "xlinked_monoallelic"))

#Remove "deNovo" segregation pattern
monoallelic_mother <- filter(monoallelic_mother, monoallelic_mother[15] != "deNovo")
table(monoallelic_mother[14])
table(monoallelic_mother[15])
table(monoallelic_mother[20])

# Father LOF
LOF_father <- filter(monoallelic_father, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant","splice_acceptor_variant","splice_donor_variant", "splive_region_variant","start_lost"))
#218 observations

# Father missense
missense_father <- filter(monoallelic_father, consequence_type %in% c("inframe_deletion", "inframe_insertion", "missense_Variant", "start_gained", "stop_lost"))
#438 observations

# Mother LOF
LOF_mother <- filter(monoallelic_mother, consequence_type %in% c("NMD_transcript_variant","feature_truncation","frameshift_variant","splice_acceptor_variant","splice_donor_variant", "splive_region_variant","start_lost"))
#279 observations

# Mother missense 
missense_mother <- filter(monoallelic_mother, consequence_type %in% c("inframe_deletion", "inframe_insertion", "missense_Variant", "start_gained", "stop_lost"))
#563 observations



