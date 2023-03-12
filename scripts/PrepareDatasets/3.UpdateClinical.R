
##############################################################################################################
############################################# UPDATE CLINICAL ################################################
##############################################################################################################

### LOAD CLINICAL

filepath <- "../PrepareDatasets/ProcessedDatasets/custom_clinical.txt"
clinical <- as.data.frame(readr::read_delim(filepath))


### LOAD sample info regarding the germline mutation analysis

filepath <- "../PrepareDatasets/Datasets/samples.txt"
germ_analysis_samples <- as.data.frame(readr::read_delim(filepath))
germ_analysis_samples <- colnames(germ_analysis_samples)
germ_analysis_samples <- germ_analysis_samples[-c(1:9)]
germ_analysis_patients <- unlist(lapply(germ_analysis_samples,function(x) paste(strsplit(x,"-",fixed = T)[[1]][1:3],collapse = "-")))
filtered.clinical <- clinical[clinical$ID %in% germ_analysis_patients,]


### LOAD Hereditary cases
filepath <- "../PrepareDatasets/ProcessedDatasets/hereditary_BC_cases.txt"
hereditary_BC <- as.data.frame(readr::read_delim(filepath))
sporadic_clinical <- filtered.clinical[!filtered.clinical$ID %in% hereditary_BC$bcr_patient_barcode,]


##############################################################################################################
########################################## SAVE SPORADIC CLINICAL ############################################
##############################################################################################################

new.dir <- paste0("./ProcessedDatasets/")
dir.create(new.dir,recursive = T)

#write
readr::write_delim(sporadic_clinical ,paste0(new.dir,"sporadic_clinical.txt"))
readr::write_excel_csv(sporadic_clinical ,paste0(new.dir,"sporadic_clinical.csv"))


##############################################################################################################
#################################################### END #####################################################
##############################################################################################################

