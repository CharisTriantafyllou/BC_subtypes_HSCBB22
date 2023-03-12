
##############################################################################################################
################################################# START ######################################################
##############################################################################################################

### Load packages
library(AnnotationDbi)
library(org.Hs.eg.db)


##############################################################################################################
#################################### DISCOVER HEREDITARY CANCER CASES ########################################
##############################################################################################################


### Find genes involved in hereditary Breast Cancer
# List downloaded from: "https://testguide.labmed.uw.edu/public/view/BROCA"

filepath <- "../PrepareDatasets/Datasets/broca_list.xlsx"
broca_genes <- openxlsx::read.xlsx(filepath)
broca_breast <- broca_genes[grep("Breast",broca_genes[,3],fixed=T),]
predisposing_genes <- broca_breast$Gene


### Find the appropriate annonation key - it is ALIAS most probably

sum(predisposing_genes %in% keys(org.Hs.eg.db,"ALIAS"))
sum(predisposing_genes %in% keys(org.Hs.eg.db,"SYMBOL"))


### Convert ALIAS IDs to SYMBOL IDs

predisposing_gene_symbol <- mapIds(org.Hs.eg.db,predisposing_genes,"SYMBOL","ALIAS")


### Investigate Hereditary Cases using information about sporadic mutation in TCGA data
#Downloaded from
#browseURL("https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Germline-AWG")

filepath <- "../PrepareDatasets/Datasets/PCA_pathVar_integrated_filtered_adjusted.tsv"
germline_mut <- as.data.frame(readr::read_delim(filepath,delim = "\t"))
tcga_brca_germ <- germline_mut[germline_mut$cancer %in% "BRCA",] #select BRCA cases only


### Find the appropriate annonation key - it is ALIAS most probably

sum(tcga_brca_germ$HUGO_Symbol %in% keys(org.Hs.eg.db,"ALIAS"))
sum(tcga_brca_germ$HUGO_Symbol %in% keys(org.Hs.eg.db,"SYMBOL"))


### Convert to SYMBOL

tcga_brca_germ$gene_Symbol <- mapIds(org.Hs.eg.db,tcga_brca_germ$HUGO_Symbol,"SYMBOL","ALIAS")


################################## SELECT SAMPLES ################################

#just a test
#predisposing_gene_symbol %in% unique(tcga_brca_germ$gene_Symbol)

hereditary_BC <- tcga_brca_germ[tcga_brca_germ$HUGO_Symbol %in% predisposing_gene_symbol,]

#just a test - how many genes:12
length(unique(hereditary_BC$HUGO_Symbol))

hereditary <- unique(hereditary_BC$bcr_patient_barcode)

##############################################################################################################
################################################## SAVE ######################################################
##############################################################################################################

new.dir <- paste0("./ProcessedDatasets/")
dir.create(new.dir,recursive = T)

#write
readr::write_delim(hereditary_BC ,paste0(new.dir,"hereditary_BC_cases.txt"))
readr::write_excel_csv(hereditary_BC ,paste0(new.dir,"hereditary_BC_cases.csv"))
rm(hereditary_BC)

##############################################################################################################
################################################### END ######################################################
##############################################################################################################
