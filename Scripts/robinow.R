library(Rlabkey)

# Load data from gmc_exit_questionnaire:
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

# Load data from Participants:
participant_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="participant", 
  viewName="", 
  colSelect="participant_id,gel_case_reference,gel_case_programme_id,latest_case_activity_datetime,participant_sk,programme,handling_gmc_ods_code,handling_gmc,handling_gmc_trust,year_of_birth,participant_ethnic_category,participant_phenotypic_sex,participant_karyotypic_sex,participant_stated_gender,registration_date,registered_at_gmc_ods_code,registered_at_gmc,registered_at_gmc_trust,registered_at_ldp_ods_code,registered_at_ldp_bioinformatics_ods_code,registered_at_ldp_organisation_name,registered_at_ldp_site_name,programme_consent_status,date_of_consent,consent_form,information_sheet,withdrawal_of_consent_date,withdrawal_option,withdrawal_form,rare_diseases_family_sk,rare_diseases_family_id,participant_type,fathers_ethnic_category,fathers_ethnic_category_other,fathers_other_relevant_ancestry,mothers_ethnic_category,mothers_ethnic_category_other,mothers_other_relevant_ancestry,participant_medical_review_date,participant_medical_review_qc_state_code,participant_medical_review_qc_state_description,consanguinity,consanguinity_sk,father_affected,mother_affected,out_of_area_recruitment,penetrance,full_brothers_affected,full_sisters_affected,total_full_brothers,total_full_sisters,biological_relationship_to_proband,other_biological_relationship_to_proband,participant_pipeline_status_id,participant_pipeline_status,reproductive_additional_findings,health_related_additional_findings,normalised_consent_form,duplicated_participant_id", 
  colFilter=NULL, 
  containerFilter=NULL, 
  colNameOpt="rname"
)

table(participant_data[6])
library(dplyr)

# Filter for Rare Diseases
RD_participant_data <- (participant_data %>% filter(participant_data[6] == "Rare Diseases"))

# Filter participant data for participant consenting
table(RD_participant_data[23])
consenting_participants <- (RD_participant_data %>% filter(RD_participant_data[23] == "Consenting"))

# Semi join gmc_data and consenting_participants to return only those participants in gmc_data that are in consenting_participants
gmc_data_consenting <- semi_join(gmc_data,consenting_participants, by = "participant_id")

# Load data from rare_diseases_family
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

# Filter RD_family_data for Family Group Type
# Know the names of the Family Group Type:
table(RD_family_data[3])
RD_family_data_filtered <- filter(RD_family_data, RD_family_data[3] == 'Duo with Mother or Father' | RD_family_data[3] == 'Families with more than three Participants' | RD_family_data[3] == 'Trio with Mother and Father' | RD_family_data[3] == 'Trio with Mother or Father and other Biological Relationship')
table(RD_family_data_filtered[3])

# Rename "Rare Diseases Family Id" column to "Family Id"
names(RD_family_data_filtered)[names(RD_family_data_filtered) == 'rare_diseases_family_id'] <- "family_id"

# Semi join gmc_data_consenting and RD_family_data_filtered
gmc_data_filtered <- semi_join(gmc_data_consenting, RD_family_data_filtered, by = "family_id")

# Filter gmc_data_filtered for "unsolved", "partially solved" or "unknown" participants
table(gmc_data_filtered[7])
gmc_unsolved <- filter(gmc_data_filtered, gmc_data_filtered[7] == 'no' | gmc_data_filtered[7] == 'partially' | gmc_data_filtered[7] == 'unknown')

# Load data from tiereing_data and gene list
tiering_data <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="tiering_data", 
  viewName="", 
  colSelect="participant_id,rare_diseases_family_id,interpretation_cohort_id,interpretation_request_id,sample_id,phenotype,participant_type,assembly,chromosome,position,reference,alternate,genotype,mode_of_inheritance,segregation_pattern,penetrance,tier,genomic_feature_hgnc,ensembl_id,consequence_type,so_term,father_affected,mother_affected,full_brothers_affected,full_sisters_affected,participant_phenotypic_sex", 
  colFilter=makeFilter(c("genomic_feature_hgnc", "IN", "ABL1;ABL2;ANKRD6;CCDC88A;CCDC88C;CELSR1;CELSR3;CSNK1D;CSNK1E;CSNK2A1;CTHRC1;CTNNB1;DAAM1;DAAM2;DAB2;DACT1;DVL1;DVL2;DVL3;FGD1;FGFR1;FGFR2;FGFR3;FGFR4;FRZB;FZD1;FZD2;FZD3;FZD4;FZD5;FZD6;FZD7;FZD8;FZD9;FZD10;GPC3;GPC4;GSK3B;IFT80;INPPL;ITCH;KIF26B;KLHL12;KLHL20;LRP5;LRP6;MKS1;MLLT3;MYOC;NXN;PLEKHA4;PLK1;PORCN;PRICKLE1;PRICKLE2;PTK7;PTPN11;RAC1;RAC3;RHOA;RNF213;ROCK1;ROCK2;ROR1;ROR2;RNF43;RPGRIP1L;RSPO1;RSPO3;RYK;SCRIBBLE;SEC63;SH3PXD2B;SFRP1;SFRP2;SFRP4;SFRP5;SMURF1;SMURF2;TMEM67;TXNDC2;VANGL1;VANGL2;WNT1;WNT2;WNT2B;WNT3;WNT3A;WNT5A;WNT5B;WNT6;WNT7A;WNT7B;WNT8A;WNT8B;WNT10A;WNT10B;WNT11;WNT14;WNT15;WNT16;WWP1;WWP2;ZNRF3")),
  containerFilter=NULL, 
  colNameOpt="rname"
)

# Semi join tiering_data and gmc_unsolved
unsolved_tiering_data <- semi_join(tiering_data, gmc_unsolved, by = "participant_id")

#Filter for unaffected parents
table(unsolved_tiering_data[22])
table(unsolved_tiering_data[23])

# Taking into account "No" 
unsolved_tiering_data_unaffected <- (unsolved_tiering_data %>% filter(unsolved_tiering_data[22] == "No" & unsolved_tiering_data[23] == "No"))

# robinow <- c("ABL1","ABL2","ANKRD6","CCDC88A", "CCDC88C", "CELSR1", "CELSR3", "CSNK1D", "CSNK1E", "CSNK2A1", "CTHRC1", "CTNNB1", "DAAM1", "DAAM2", "DAB2", "DACT1", "DVL1", "DVL2", "DVL3", "FGD1", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FRZB", "FZD1", "FZD2", "FZD3", "FZD4", "FZD2", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10", "GPC3", "GPC4", "GSK3B", "IFT80", "INPPL", "ITCH", "KIF26B", "KLHL12", "KLHL20", "LRP5", "LRP6", "MKS1", "MLLT3", "MYOC", "NXN", "PLEKHA4", "PLK1", "PORCN", "PRICKLE1", "PRICKLE2", "PTK7", "PTPN11", "RAC1", "RAC3", "RHOA", "RNF213", "ROCK1", "ROCK2", "ROR1", "ROR2", "RNF43", "RPGRIP1L", "RSPO1", "RSPO3", "RYK", "SCRIBBLE", "SEC63", "SH3PXD2B", "SFRP1", "SFRP2", "SFRP4", "SFRP5", "SMURF1", "SMURF2", "TMEM67", "TXNDC2", "VANGL1", "VANGL2", "WNT1","WNT2", "WNT2B", "WNT3", "WNT3A", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT10A", "WNT10B", "WNT11", "WNT14", "WNT15", "WNT16", "WWP1", "WWP2", "ZNRF3")

#Filter for LOF
table(unsolved_tiering_data_unaffected[20])
LOF <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("feature_truncation","frameshift_variant", "splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost","stop_gained"))

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
unsolved_recessive_lof <- merge(recessive, recessive_lof)

# Remove duplicates
unsolved_recessive_lof <- unsolved_recessive_lof[!duplicated(unsolved_recessive_lof[c(1,2,3,4,5,6)]),]

# Save
write.csv(unsolved_recessive_lof, "unsolved_recessive_lof_robinow.csv")

# DeNovo LOF
unsolved_denovo_lof <- (LOF %>% filter(LOF[15] == 'deNovo'))
denovo_lof <- unsolved_denovo_lof[!duplicated(unsolved_denovo_lof[c(1,6,10,11,12,20)]),]

# Save
write.csv(denovo_lof, "unsolved_denovo_lof_robinow.csv")

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
unsolved_comphet_lof <- merge(comphet, comphet_lof)
# Remove duplicates
unsolved_comphet_lof <- unsolved_comphet_lof[!duplicated(unsolved_comphet_lof[c(1,2,3,4,5,6)]),]

# Remove participant_id that only appears once
unsolved_comphet_lof <- unsolved_comphet_lof %>%
  group_by(participant_id) %>%
  filter(n()>1)

#Save
write.csv(unsolved_comphet_lof, "unsolved_comphet_lof_robinow.csv")

# Filter for missense+ 
table(unsolved_tiering_data_unaffected[20])
missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("inframe_deletion","inframe_insertion","missense_variant"))

# Simple recessive missense+ 
unsolved_recessive_missense <- (missense %>% filter(missense[15] == 'SimpleRecessive'))

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
unsolved_recessive_missense <- merge(recessive_miss, recessive_missense)
# Remove duplicates
unsolved_recessive_missense <- unsolved_recessive_missense[!duplicated(unsolved_recessive_missense[c(1,2,3,4,5,6)]),]

#Save
write.csv(unsolved_recessive_missense, "unsolved_recessive_missense_robinow.csv")

# DeNovo missense+
unsolved_denovo_missense <- (missense %>% filter(missense[15] == 'deNovo'))

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
unsolved_denovo_missense <- merge(denovo_miss, denovo_missense)
# Remove duplicates
unsolved_denovo_missense <- unsolved_denovo_missense[!duplicated(unsolved_denovo_missense[c(1,2,3,4,5,6)]),]

#Save
write.csv(unsolved_denovo_missense, "unsolved_denovo_missense_robinow.csv")

# CompoundHeterozygous missense+
unsolved_comphet_missense <- (missense %>% filter(missense[15] == "CompoundHeterozygous"))

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
unsolved_comphet_missense <- merge(comphet_miss, comphet_missense)
# Remove duplicates
unsolved_comphet_missense <- unsolved_comphet_missense[!duplicated(unsolved_comphet_missense[c(1,2,3,4,5,6)]),]

# Remove participant_id that only appears once
unsolved_comphet_missense <- unsolved_comphet_missense %>%
  group_by(participant_id) %>%
  filter(n()>1)

#Save
write.csv(unsolved_comphet_missense, "unsolved_comphet_missense_robinow.csv")

# LOF + missense
LOF_missense <- filter(unsolved_tiering_data_unaffected, consequence_type %in% c("feature_truncation","frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "start_lost", "stop_gained", "stop_lost"))

# Simple recessive LOF + missense
unsolved_recessive_LOFmissense <- (LOF_missense %>% filter(LOF_missense[15] == 'SimpleRecessive'))

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
unsolved_recessive_LOFmissense <- merge(recessive_LOFmiss, recessive_LOFmissense)
# Remove duplicates
unsolved_recessive_LOFmissense <- unsolved_recessive_LOFmissense[!duplicated(unsolved_recessive_LOFmissense[c(1,2,3,4,5,6)]),]
#Save
write.csv(unsolved_recessive_LOFmissense, "unsolved_recessive_LOFmissense_robinow.csv")

# DeNovo LOF + missense
unsolved_denovo_LOFmissense <- (missense %>% filter(missense[15] == 'deNovo'))

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
unsolved_denovo_LOFmissense <- merge(denovo_LOFmiss, denovo_LOFmissense)
# Remove duplicates
unsolved_denovo_LOFmissense <- unsolved_denovo_LOFmissense[!duplicated(unsolved_denovo_LOFmissense[c(1,2,3,4,5,6)]),]

#Save
write.csv(unsolved_denovo_LOFmissense, "unsolved_denovo_LOFmissense_robinow.csv")

# CompoundHeterozygous LOF + missense
unsolved_comphet_LOFmissense <- (missense %>% filter(missense[15] == "CompoundHeterozygous"))

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
unsolved_comphet_LOFmissense <- merge(comphet_LOFmiss, comphet_LOFmissense)
# Remove duplicates
unsolved_comphet_LOFmissense <- unsolved_comphet_LOFmissense[!duplicated(unsolved_comphet_LOFmissense[c(1,2,3,4,5,6)]),]

# Remove participant_id that only appears once
unsolved_comphet_LOFmissense <- unsolved_comphet_LOFmissense %>%
  group_by(participant_id) %>%
  filter(n()>1)

#Save
write.csv(unsolved_comphet_LOFmissense, "unsolved_comphet_LOFmissense_robinow.csv")


#########
#  VEP  #
#########

library(stringr)

#Simple recessive LOF
table(unsolved_recessive_lof[12])
unsolved_recessive_lof_grch37 <- unsolved_recessive_lof %>% filter(assembly == "GRCh37")
write.csv(unsolved_recessive_lof_grch37, "unsolved_recessive_lof_robinow_grch37.csv")

unsolved_recessive_lof_grch37$base <- paste(as.character(unsolved_recessive_lof_grch37$reference), "/", as.character(unsolved_recessive_lof_grch37$alternate))
unsolved_recessive_lof_grch37_VEP <- unsolved_recessive_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_grch37_VEP_list <- unsolved_recessive_lof_grch37_VEP[!duplicated(unsolved_recessive_lof_grch37_VEP[,2:4]),]
unsolved_recessive_lof_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_grch37_VEP_list, "VEP_recessive_lof_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_lof_grch38 <- unsolved_recessive_lof %>% filter(assembly == "GRCh38")
write.csv(unsolved_recessive_lof_grch38, "unsolved_recessive_lof_robinow_grch38.csv")

unsolved_recessive_lof_grch38$base <- paste(as.character(unsolved_recessive_lof_grch38$reference), "/", as.character(unsolved_recessive_lof_grch38$alternate))
unsolved_recessive_lof_grch38_VEP <- unsolved_recessive_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_lof_grch38_VEP_list <- unsolved_recessive_lof_grch38_VEP[!duplicated(unsolved_recessive_lof_grch38_VEP[,2:4]),]
unsolved_recessive_lof_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_lof_grch38_VEP_list, "VEP_recessive_lof_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#DeNovo LOF
table(unsolved_denovo_lof[8])
unsolved_denovo_lof_grch37 <- unsolved_denovo_lof %>% filter(assembly == "GRCh37")
write.csv(unsolved_denovo_lof_grch37, "unsolved_denovo_lof_robinow_grch37.csv")

unsolved_denovo_lof_grch37$base <- paste(as.character(unsolved_denovo_lof_grch37$reference), "/", as.character(unsolved_denovo_lof_grch37$alternate))
unsolved_denovo_lof_grch37_VEP <- unsolved_denovo_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_grch37_VEP_list <- unsolved_denovo_lof_grch37_VEP[!duplicated(unsolved_denovo_lof_grch37_VEP[,2:4]),]
unsolved_denovo_lof_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_grch37_VEP_list, "VEP_denovo_lof_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_lof_grch38 <- unsolved_denovo_lof %>% filter(assembly == "GRCh38")
write.csv(unsolved_denovo_lof_grch38, "unsolved_denovo_lof_robinow_grch38.csv")

unsolved_denovo_lof_grch38$base <- paste(as.character(unsolved_denovo_lof_grch38$reference), "/", as.character(unsolved_denovo_lof_grch38$alternate))
unsolved_denovo_lof_grch38_VEP <- unsolved_denovo_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_lof_grch38_VEP_list <- unsolved_denovo_lof_grch38_VEP[!duplicated(unsolved_denovo_lof_grch38_VEP[,2:4]),]
unsolved_denovo_lof_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_lof_grch38_VEP_list, "VEP_denovo_lof_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Comphet LOF
table(unsolved_comphet_lof[12])
unsolved_comphet_lof_grch37 <- unsolved_comphet_lof %>% filter(assembly == "GRCh37")
write.csv(unsolved_comphet_lof_grch37, "unsolved_comphet_lof_robinow_grch37.csv")

unsolved_comphet_lof_grch37$base <- paste(as.character(unsolved_comphet_lof_grch37$reference), "/", as.character(unsolved_comphet_lof_grch37$alternate))
unsolved_comphet_lof_grch37_VEP <- unsolved_comphet_lof_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_grch37_VEP_list <- unsolved_comphet_lof_grch37_VEP[!duplicated(unsolved_comphet_lof_grch37_VEP[,2:4]),]
unsolved_comphet_lof_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_lof_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_grch37_VEP_list, "VEP_comphet_lof_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_lof_grch38 <- unsolved_comphet_lof %>% filter(assembly == "GRCh38")
write.csv(unsolved_comphet_lof_grch38, "unsolved_comphet_lof_robinow_grch38.csv")

unsolved_comphet_lof_grch38$base <- paste(as.character(unsolved_comphet_lof_grch38$reference), "/", as.character(unsolved_comphet_lof_grch38$alternate))
unsolved_comphet_lof_grch38_VEP <- unsolved_comphet_lof_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_lof_grch38_VEP_list <- unsolved_comphet_lof_grch38_VEP[!duplicated(unsolved_comphet_lof_grch38_VEP[,2:4]),]
unsolved_comphet_lof_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_lof_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_lof_grch38_VEP_list, "VEP_comphet_lof_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Simple recessive missense
table(unsolved_recessive_missense[12])
unsolved_recessive_missense_grch37 <- unsolved_recessive_missense %>% filter(assembly == "GRCh37")
write.csv(unsolved_recessive_missense_grch37, "unsolved_recessive_missense_robinow_grch37.csv")

unsolved_recessive_missense_grch37$base <- paste(as.character(unsolved_recessive_missense_grch37$reference), "/", as.character(unsolved_recessive_missense_grch37$alternate))
unsolved_recessive_missense_grch37_VEP <- unsolved_recessive_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_grch37_VEP_list <- unsolved_recessive_missense_grch37_VEP[!duplicated(unsolved_recessive_missense_grch37_VEP[,2:4]),]
unsolved_recessive_missense_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_grch37_VEP_list, "VEP_recessive_missense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_missense_grch38 <- unsolved_recessive_missense %>% filter(assembly == "GRCh38")
write.csv(unsolved_recessive_missense_grch38, "unsolved_recessive_missense_robinow_grch38.csv")

unsolved_recessive_missense_grch38$base <- paste(as.character(unsolved_recessive_missense_grch38$reference), "/", as.character(unsolved_recessive_missense_grch38$alternate))
unsolved_recessive_missense_grch38_VEP <- unsolved_recessive_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_missense_grch38_VEP_list <- unsolved_recessive_missense_grch38_VEP[!duplicated(unsolved_recessive_missense_grch38_VEP[,2:4]),]
unsolved_recessive_missense_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_missense_grch38_VEP_list, "VEP_recessive_missense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#DeNovo missense
table(unsolved_denovo_missense[12])
unsolved_denovo_missense_grch37 <- unsolved_denovo_missense %>% filter(assembly == "GRCh37")
write.csv(unsolved_denovo_missense_grch37, "unsolved_denovo_missense_robinow_grch37.csv")

unsolved_denovo_missense_grch37$base <- paste(as.character(unsolved_denovo_missense_grch37$reference), "/", as.character(unsolved_denovo_missense_grch37$alternate))
unsolved_denovo_missense_grch37_VEP <- unsolved_denovo_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_grch37_VEP_list <- unsolved_denovo_missense_grch37_VEP[!duplicated(unsolved_denovo_missense_grch37_VEP[,2:4]),]
unsolved_denovo_missense_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_grch37_VEP_list, "VEP_denovo_missense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_missense_grch38 <- unsolved_denovo_missense %>% filter(assembly == "GRCh38")
write.csv(unsolved_denovo_missense_grch38, "unsolved_denovo_missense_robinow_grch38.csv")

unsolved_denovo_missense_grch38$base <- paste(as.character(unsolved_denovo_missense_grch38$reference), "/", as.character(unsolved_denovo_missense_grch38$alternate))
unsolved_denovo_missense_grch38_VEP <- unsolved_denovo_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_missense_grch38_VEP_list <- unsolved_denovo_missense_grch38_VEP[!duplicated(unsolved_denovo_missense_grch38_VEP[,2:4]),]
unsolved_denovo_missense_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_missense_grch38_VEP_list, "VEP_denovo_missense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#Comphet missense
table(unsolved_comphet_missense[12])
unsolved_comphet_missense_grch37 <- unsolved_comphet_missense %>% filter(assembly == "GRCh37")
write.csv(unsolved_comphet_missense_grch37, "unsolved_comphet_missense_robinow_grch37.csv")

unsolved_comphet_missense_grch37$base <- paste(as.character(unsolved_comphet_missense_grch37$reference), "/", as.character(unsolved_comphet_missense_grch37$alternate))
unsolved_comphet_missense_grch37_VEP <- unsolved_comphet_missense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_grch37_VEP_list <- unsolved_comphet_missense_grch37_VEP[!duplicated(unsolved_comphet_missense_grch37_VEP[,2:4]),]
unsolved_comphet_missense_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_missense_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_grch37_VEP_list, "VEP_comphet_missense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_missense_grch38 <- unsolved_comphet_missense %>% filter(assembly == "GRCh38")
write.csv(unsolved_comphet_missense_grch38, "unsolved_comphet_missense_robinow_grch38.csv")

unsolved_comphet_missense_grch38$base <- paste(as.character(unsolved_comphet_missense_grch38$reference), "/", as.character(unsolved_comphet_missense_grch38$alternate))
unsolved_comphet_missense_grch38_VEP <- unsolved_comphet_missense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_missense_grch38_VEP_list <- unsolved_comphet_missense_grch38_VEP[!duplicated(unsolved_comphet_missense_grch38_VEP[,2:4]),]
unsolved_comphet_missense_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_missense_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_missense_grch38_VEP_list, "VEP_comphet_missense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#Simple recessive LOF + missense
table(unsolved_recessive_LOFmissense[12])
unsolved_recessive_LOFmissense_grch37 <- unsolved_recessive_LOFmissense %>% filter(assembly == "GRCh37")
write.csv(unsolved_recessive_LOFmissense_grch37, "unsolved_recessive_LOFmissense_robinow_grch37.csv")

unsolved_recessive_LOFmissense_grch37$base <- paste(as.character(unsolved_recessive_LOFmissense_grch37$reference), "/", as.character(unsolved_recessive_LOFmissense_grch37$alternate))
unsolved_recessive_LOFmissense_grch37_VEP <- unsolved_recessive_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_grch37_VEP_list <- unsolved_recessive_LOFmissense_grch37_VEP[!duplicated(unsolved_recessive_LOFmissense_grch37_VEP[,2:4]),]
unsolved_recessive_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_grch37_VEP_list, "VEP_recessive_LOFmissense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_recessive_LOFmissense_grch38 <- unsolved_recessive_LOFmissense %>% filter(assembly == "GRCh38")
write.csv(unsolved_recessive_LOFmissense_grch38, "unsolved_recessive_LOFmissense_robinow_grch38.csv")

unsolved_recessive_LOFmissense_grch38$base <- paste(as.character(unsolved_recessive_LOFmissense_grch38$reference), "/", as.character(unsolved_recessive_LOFmissense_grch38$alternate))
unsolved_recessive_LOFmissense_grch38_VEP <- unsolved_recessive_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_recessive_LOFmissense_grch38_VEP_list <- unsolved_recessive_LOFmissense_grch38_VEP[!duplicated(unsolved_recessive_LOFmissense_grch38_VEP[,2:4]),]
unsolved_recessive_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_recessive_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_recessive_LOFmissense_grch38_VEP_list, "VEP_recessive_LOFmissense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#DeNovo LOF + missense
table(unsolved_denovo_LOFmissense[12])
unsolved_denovo_LOFmissense_grch37 <- unsolved_denovo_LOFmissense %>% filter(assembly == "GRCh37")
write.csv(unsolved_denovo_LOFmissense_grch37, "unsolved_denovo_LOFmissense_robinow_grch37.csv")

unsolved_denovo_LOFmissense_grch37$base <- paste(as.character(unsolved_denovo_LOFmissense_grch37$reference), "/", as.character(unsolved_denovo_LOFmissense_grch37$alternate))
unsolved_denovo_LOFmissense_grch37_VEP <- unsolved_denovo_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_grch37_VEP_list <- unsolved_denovo_LOFmissense_grch37_VEP[!duplicated(unsolved_denovo_LOFmissense_grch37_VEP[,2:4]),]
unsolved_denovo_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_grch37_VEP_list, "VEP_denovo_LOFmissense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_denovo_LOFmissense_grch38 <- unsolved_denovo_LOFmissense %>% filter(assembly == "GRCh38")
write.csv(unsolved_denovo_LOFmissense_grch38, "unsolved_denovo_LOFmissense_robinow_grch38.csv")

unsolved_denovo_LOFmissense_grch38$base <- paste(as.character(unsolved_denovo_LOFmissense_grch38$reference), "/", as.character(unsolved_denovo_LOFmissense_grch38$alternate))
unsolved_denovo_LOFmissense_grch38_VEP <- unsolved_denovo_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_denovo_LOFmissense_grch38_VEP_list <- unsolved_denovo_LOFmissense_grch38_VEP[!duplicated(unsolved_denovo_LOFmissense_grch38_VEP[,2:4]),]
unsolved_denovo_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_denovo_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_denovo_LOFmissense_grch38_VEP_list, "VEP_denovo_LOFmissense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#Comphet LOF + missense
table(unsolved_comphet_LOFmissense[12])
unsolved_comphet_LOFmissense_grch37 <- unsolved_comphet_LOFmissense %>% filter(assembly == "GRCh37")
write.csv(unsolved_comphet_LOFmissense_grch37, "unsolved_comphet_LOFmissense_robinow_grch37.csv")

unsolved_comphet_LOFmissense_grch37$base <- paste(as.character(unsolved_comphet_LOFmissense_grch37$reference), "/", as.character(unsolved_comphet_LOFmissense_grch37$alternate))
unsolved_comphet_LOFmissense_grch37_VEP <- unsolved_comphet_LOFmissense_grch37 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_grch37_VEP_list <- unsolved_comphet_LOFmissense_grch37_VEP[!duplicated(unsolved_comphet_LOFmissense_grch37_VEP[,2:4]),]
unsolved_comphet_LOFmissense_grch37_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_grch37_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_grch37_VEP_list, "VEP_comphet_LOFmissense_robinow_grch37", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

unsolved_comphet_LOFmissense_grch38 <- unsolved_comphet_LOFmissense %>% filter(assembly == "GRCh38")
write.csv(unsolved_comphet_LOFmissense_grch38, "unsolved_comphet_LOFmissense_robinow_grch38.csv")

unsolved_comphet_LOFmissense_grch38$base <- paste(as.character(unsolved_comphet_LOFmissense_grch38$reference), "/", as.character(unsolved_comphet_LOFmissense_grch38$alternate))
unsolved_comphet_LOFmissense_grch38_VEP <- unsolved_comphet_LOFmissense_grch38 %>% mutate(start = position, end = position + (nchar(alternate) + nchar(reference) -2), strand = 1) %>% select(chromosome, start, end, base, strand)
unsolved_comphet_LOFmissense_grch38_VEP_list <- unsolved_comphet_LOFmissense_grch38_VEP[!duplicated(unsolved_comphet_LOFmissense_grch38_VEP[,2:4]),]
unsolved_comphet_LOFmissense_grch38_VEP_list$base <- str_replace_all(unsolved_comphet_LOFmissense_grch38_VEP_list$base, " ", "")
write.table(unsolved_comphet_LOFmissense_grch38_VEP_list, "VEP_comphet_LOFmissense_robinow_grch38", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


