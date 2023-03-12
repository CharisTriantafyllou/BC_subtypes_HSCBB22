
##############################################################################################################
################################################# START ######################################################
##############################################################################################################

### Define data directory ###

data.dir <- "./Datasets/"


##############################################################################################################
##################################### COMBINE CLINICAL TABLES ################################################
##############################################################################################################

### Load data ###

clinical.patient <- readr::read_delim(paste0(data.dir,"data_clinical_patient.txt"))
clinical.sample <- readr::read_delim(paste0(data.dir,"data_clinical_sample.txt"))


### Filter out descriptions ###

clinical.patient <- as.data.frame(clinical.patient)
idx <- which(clinical.patient[,1]=="PATIENT_ID")
colnames(clinical.patient) <- clinical.patient[idx,]
clinical.patient <- clinical.patient[-c(1:idx),]

clinical.sample <- as.data.frame(clinical.sample)
idx <- which(clinical.sample[,1]=="PATIENT_ID")
colnames(clinical.sample) <- clinical.sample[idx,]
clinical.sample <- clinical.sample[-c(1:idx),]


### find common columns and patient IDs

common.cols <- intersect(colnames(clinical.patient),colnames(clinical.sample)) 
common.patients <- intersect(clinical.sample$PATIENT_ID,clinical.patient$PATIENT_ID)


### Rename rows of dataframes

rownames(clinical.patient) <- clinical.patient$PATIENT_ID
rownames(clinical.sample) <- clinical.sample$PATIENT_ID


### first merge unique columns for the same patients

clinical.extended <- cbind(clinical.sample[common.patients,!colnames(clinical.sample) %in% common.cols],
                           clinical.patient[common.patients,!colnames(clinical.patient) %in% common.cols])


### add common columns
# It is expected to be only PATIENT ID

clinical.extended <- cbind(clinical.sample[common.patients,common.cols],clinical.extended)


### Restore column names in the new data frame 
# If there is only one common column then the name is lost so we have to restore it

if(length(common.cols)==1){
  colnames(clinical.extended)[1] <- common.cols
  message(paste0("common column is ",common.cols))
}

### Clean memory space

gc()


##############################################################################################################
######################### MODIFY AND SELECT CLINICAL VALUES FOR SURVIVAL ANALYSIS ############################
##############################################################################################################

### Use a simple name

clinical <- clinical.extended


### Create and add a couple of extra variables

clinical$T_M_B_NON_SYNONYMOUS_LOG10 <- log10(as.numeric(clinical$TMB_NONSYNONYMOUS))
clinical$M_S_I_SCORE_MANTIS_LOG10 <- log10(as.numeric(clinical$MSI_SCORE_MANTIS))


### Simplify this variable ...
#About Stages
#browseURL("https://www.cancer.org/cancer/breast-cancer/understanding-a-breast-cancer-diagnosis/stages-of-breast-cancer.html")
#browseURL("https://www.biostars.org/p/181470/")

clinical$PATH_STAGE <- clinical$AJCC_PATHOLOGIC_TUMOR_STAGE
clinical$PATH_STAGE[grep("I{3}",clinical$PATH_STAGE)] <- 3
clinical$PATH_STAGE[grep("I{2}",clinical$PATH_STAGE)] <- 2
clinical$PATH_STAGE[grep("I{1}",clinical$PATH_STAGE)] <- 1
clinical$PATH_STAGE[!clinical$PATH_STAGE %in% c(1,2,3)] <- NA
clinical$PATH_STAGE <- as.numeric(clinical$PATH_STAGE)

### Correct the values of this variable

clinical$SUBTYPE <- gsub("^BRCA_","",clinical$SUBTYPE)


### Transform from months to years

clinical$DFS_MONTHS <- round(as.numeric(clinical$DFS_MONTHS)/12,2) #AFTER THAT WE HAVE TO RENAME THE COLUMN
clinical$DFS_STATUS <- as.numeric(unlist(lapply(clinical$DFS_STATUS,
                                                function(x){stringr::str_split(x,":")[[1]][1]})))


### Transform from months to years
clinical$OS_MONTHS <- round(as.numeric(clinical$OS_MONTHS)/12,2) #AFTER THAT WE HAVE TO RENAME THE COLUMN
clinical$OS_STATUS <- as.numeric(unlist(lapply(clinical$OS_STATUS,
                                                function(x){stringr::str_split(x,":")[[1]][1]})))


### Transform from months to years

clinical$DSS_MONTHS <- round(as.numeric(clinical$DSS_MONTHS)/12,2) #AFTER THAT WE HAVE TO RENAME THE COLUMN
clinical$DSS_STATUS <- as.numeric(unlist(lapply(clinical$DSS_STATUS,
                                               function(x){stringr::str_split(x,":")[[1]][1]})))


### Transform from months to years

clinical$PFS_MONTHS <- round(as.numeric(clinical$PFS_MONTHS)/12,2) #AFTER THAT WE HAVE TO RENAME THE COLUMN
clinical$PFS_STATUS <- as.numeric(unlist(lapply(clinical$PFS_STATUS,
                                               function(x){stringr::str_split(x,":")[[1]][1]})))


### Select variable to proceed with ...

clinical <- clinical[,c("PATIENT_ID","DFS_STATUS","DFS_MONTHS","OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS",
                        "PFS_STATUS","PFS_MONTHS","SUBTYPE","AGE","SEX","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_STAGING_EDITION","ETHNICITY",
                        "RACE","PATH_M_STAGE","PATH_N_STAGE","PATH_T_STAGE","PATH_STAGE","TUMOR_TYPE",
                        "ANEUPLOIDY_SCORE","MSI_SCORE_MANTIS","M_S_I_SCORE_MANTIS_LOG10","MSI_SENSOR_SCORE",
                        "TMB_NONSYNONYMOUS","T_M_B_NON_SYNONYMOUS_LOG10","RADIATION_THERAPY",
                        "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT")]


### Rename columns that their values underwent changes

colnames(clinical)[1:9] <- c("ID","DFS_EVENT","DFS_YEARS","OS_EVENT","OS_YEARS",
                             "DSS_EVENT","DSS_YEARS","PFS_EVENT","PFS_YEARS")


##############################################################################################################
################################################## SAVE ######################################################
##############################################################################################################

new.dir <- paste0("./ProcessedDatasets/")
dir.create(new.dir,recursive = T)

#write
readr::write_delim(clinical,paste0(new.dir,"custom_clinical.txt"))
readr::write_excel_csv(clinical,paste0(new.dir,"custom_clinical.csv"))


##############################################################################################################
################################################## END #######################################################
##############################################################################################################

