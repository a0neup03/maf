

clinical <- GDCquery_clinic(project = "TCGA-CHOL", type = "clinical")

tcgaProjList<-c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-UCEC","TCGA-UCS","TCGA-READ","TCGA-COAD","TCGA-LIHC","TCGA-HNSC","TCGA-ESCA","TCGA-PRAD","TCGA-STAD","TCGA-THCA","TCGA-LUAD","TCGA-LUSC","TCGA-KIRC","TCGA-KIRP","TCGA-KICH")
tcgaProjList<-c("TCGA-BLCA","TCGA-CHOL")

#mdList will have all data frames for rcgaProjList
mdListDF<-data.frame()
for(s in tcgaProjList){
  mdList<-c(mdListDF,getjoinedBiospcClinc(s))
  if(dim(mdListDF)[1]<1){
    mdListDF<-getjoinedBiospcClinc(s)
  }else{
    print("joining")
    temp<-getjoinedBiospcClinc(s)
    mdListDF<-dplyr::bind_rows(mdListDF,temp)  
  }
  
  
  
}

require(tidyr)

#######function
getjoinedBiospcClinc<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  
  #rename all cols from clinical table with suffix clinical
  colnames(clinicalBRCA)<- paste0("clinical.",colnames(clinicalBRCA))
  
  #expand biospecimen data in the order portions, portions.analytes, portions.analytes.aliquots
  
  # toUnpack<-c("portions", "portions.analytes", "portions.analytes.aliquots")
  
  ##only "portions" column is present.. need a way to unlist the "portions" column 
  
  # for(s in toUnpack){
  #  biospecimenBRCA<-tidyr::expand(biospecimenBRCA,s)
  # }
  
  
  #add patient barcode to biospecimen data
  biospecimenBRCA<- biospecimenBRCA %>% mutate(clinical.bcr_patient_barcode=substr(submitter_id,1,nchar(as.character(submitter_id))-4))
  #join clinical and biospecimen
  finalJoined<-plyr::join(clinicalBRCA,biospecimenBRCA,by="clinical.bcr_patient_barcode")
  return(finalJoined)
}

#download clinical and biospecimen data by Projname
printDatadim<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  print(paste("clinical dim:",dim(clinicalBRCA)))
  print(paste("biospc dim:",dim(biospecimenBRCA)))
}


#remove cols with all NA values
naCols<-colnames(mdListDF)[sapply(mdListDF, function(x)all(is.na(x)))]
mdListDFNONA<-mdListDF[,!(colnames(mdListDF) %in% naCols)]
#keep rows with RNA samples only
#mdListDFRNA<-mdListDF%>%filter(portions.analytes.analyte_type_id == "R")
#mdList<-mdListDFNONA

##not requirded
ulMD<-unlist(mdList)
mdJoined<-rbindlist(unlist(mdListDFNONA))
n1<-colnames(t)
n2<-colnames(BRCAMetadata)
n3<-colnames(mdListDF)
setdiff(n2,n1)

##filter using the columns
##age=== clinical.age_at_index
##gender==clinical.gender
##age==clinical.age_at_diagnosis

CGA_Clinical.tidy <- mdListDF %>%
  dplyr::rename(
    Age = clinical.age_at_diagnosis, Gender = clinical.gender,
    Tumor_Sample_Barcode = sample_id #Smoking_history = tobacco_smoking_history,
    #Smoking_indicator = tobacco_smoking_history_indicator,
  ) %>%
  filter(sample_type %in% c(
    "Solid Tissue Normal", "Primary Tumor", "Metastatic", "Recurrent Tumor",
    "Primary Blood Derived Cancer - Peripheral Blood"
  )) %>% # Additional - New Primary, Additional Metastatic, FFPE Scrolls total 14 sample removed
  mutate(Gender = case_when(
    Gender == "FEMALE" ~ "Female",
    Gender == "MALE" ~ "Male",
    TRUE ~ NA_character_
  ), Tumor_stage = case_when(
    clinical.ajcc_pathologic_stage == "Stage 0" ~ "0",
    clinical.ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "I",
    clinical.ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "II",
    clinical.ajcc_pathologic_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
    clinical.ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "IV",
    clinical.ajcc_pathologic_stage == "Stage X" ~ "X",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    Gender = factor(Gender, levels = c("Male", "Female")),
    Tumor_stage = factor(Tumor_stage, levels = c("0", "I", "II", "III", "IV", "X"))
  )

if (!file.exists("results/TCGA_tidy_Clinical.RData")) {
  dir.create("results", showWarnings = FALSE)
  save(TCGA_Clinical.tidy, file = "results/TCGA_tidy_Clinical.RData")
}


###things to do..
select sample_id with G quadruplex mutations


##########
#getting clinical data using RTCGA

BiocInstaller::biocLite(c("RTCGA.mRNA","RTCGA.clinical"))
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)

##create clinical data and extract for cancers
clin<- survivalTCGA(BRCA.clinical,OV.clinical,GBM.clinical,
                    extract.cols = c("admin.disease_code","patient.drugs.drug.therapy_types.therapy_type"))

brca_clin_orig<-survivalTCGA(BRCA.clinical,
                             extract.cols = c("patient.gender","patient.race","patient.ethnicity","patient.days_to_birth",
                                              "patient.vital_status","patient.drugs.drug.therapy_types.therapy_type",
                                              "patient.stage_event.pathologic_stage",
                                              "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
                                              "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n",
                                              "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m"))
df<-as.data.frame(BRCA.clinical)
##search columns with specific sub strings
df1 <- colnames(df[ , grepl("disease_code", names(df), perl = TRUE ) ])

brca_clin_orig %>% sapply(function(x) sum(is.na(x)))

#check number of rows with missing data

na_rows<-brca_clin_orig %>% apply(MARGIN = 1, FUN=function(x) sum(is.na(x)))
sum(na_rows>0)


brca_clin_orig %>% sapply(summary) 
brca_clin_orig %>% filter(times<0) %>% nrow()


##clean up helper functions
clean_pathologic_stage<-function(x) {
  x %>% str_replace_all(c(
    "stage iv[a-d]*"="stage4",
    "stage [i]{3}[a-d]*"="stage3",
    "stage i{2}[a-d]*"="stage2",
    "stage i{1}[a-d]*"="stage1",
    "stage x"="stageX"))
}


clean_pathologyTstage<-function(x) {
  x %>% str_replace_all(c(
    "\\s*tx"="tx",
    "\\s*t1[a-z]*"="t1",
    "\\s*t2[a-z]*"="t2",
    "\\s*t3[a-z]*"="t3",
    "\\s*t4[a-z]*"="t4"))
  
}


clean_pathologyMstage<-function(x) {
  x %>% str_replace_all(c("cm0\\s+\\(i[\\+,-]\\)"="cm0"))
  
}

clean_pathologyNstage<- function(x) {
  x %>% str_trim(.) %>%
    str_replace_all(c(
      "nx"="nx",
      "n1[a-z]*"="n1",
      "n2[a-z]*"="n2",
      "n3[a-z]*"="n3",
      "n4[a-z]*"="n4",
      "n0"="n0",
      "n0\\s+\\([a-z]*[\\+|-]\\)"="n0"))
}

##transform the data
library(stringr)
library(forcats)
brca_clin <- brca_clin_orig %>% 
  rename(gender=patient.gender,
         race=patient.race,
         ethnicity=patient.ethnicity,
         vital_status=patient.vital_status,
         therapy_type=patient.drugs.drug.therapy_types.therapy_type,
         pathologic_stage=patient.stage_event.pathologic_stage,
         pathologyTstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_t,
         pathologyNstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_n,
         pathologyMstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_m) %>%
  filter(gender=="female") %>%
  filter(times >0) %>% 
  mutate(age=abs(as.numeric(patient.days_to_birth))/365,
                 therapy_type=ifelse(is.na(therapy_type),"No Info",therapy_type),
                 therapy_type=fct_lump(therapy_type,3),
                 race=fct_lump(race,2),
                 pathologic_stage=str_trim(pathologic_stage),
                 pathologyTstage=str_trim(pathologyTstage),
                 pathologyNstage=str_trim(pathologyNstage),
                 pathologyMstage=str_trim(pathologyMstage),
         pathologic_stage=clean_pathologic_stage(pathologic_stage),
         pathologyTstage=clean_pathologyTstage(pathologyTstage),
         pathologyNstage=clean_pathologyNstage(pathologyNstage),
         pathologyMstage=clean_pathologyMstage(pathologyMstage),
                 years_to_event=times/365,
                 agecat=cut(age,breaks=c(0,40,60,Inf),labels=c("young","middle","old"))
  )

##convert specified columns from character to factor type
convert_to_factor <-c("ethnicity","pathologic_stage","pathologyTstage","pathologyMstage","pathologyNstage")
brca_clin <-brca_clin %>% mutate_at(convert_to_factor,factor)

sapply(brca_clin,class)
dim(brca_clin)

##remove unnecessary columns
brca_clin <-brca_clin %>% dplyr::select(-starts_with("patient"))





https://github.com/ramaanathan/SurvivalAnalysis/blob/master/Survival_Analysis_TCGA.R




grep(c("disease_code","drug.therapy"), names(BRCA.clinical))
clin<-
brca_clin_orig<-survivalTCGA(BRCA.clinical,extract.cols=)

