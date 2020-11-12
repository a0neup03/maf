

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