
##############################################################################################################
################################################### START ####################################################
##############################################################################################################

### Load packages

library(genefilter)
library(biomaRt)


##############################################################################################################
################################################# LOAD DATA ##################################################
##############################################################################################################

### load RNA expression from cbioportal PanCancer 2018 dataset

file.path <- "./Datasets/data_RNA_Seq_v2_expression_median.txt"
rna.df <- as.data.frame(readr::read_delim(file.path))


### load processed clinical data 

file.path <- "./ProcessedDatasets/sporadic_clinical.txt"
clinical.df <- as.data.frame(readr::read_delim(file.path))


##############################################################################################################
#################################### PROCESS AND FILTER DATA #################################################
##############################################################################################################

### Keep gene annotation and gene expression data on different objects

rna.anno <- rna.df[,grep("^TCGA",colnames(rna.df),invert = T)]
rna.exp <- rna.df[,grep("^TCGA",colnames(rna.df))]


### Add patient information to label data (there are no duplicates expected in this data set)

rna.samples <- colnames(rna.exp)
rna.patients <- gsub("-01$","",rna.samples)
colnames(rna.exp) <- rna.patients


### Keep the same patients in clinical and expression dataset to filter out hereditary cases

patients.filt <- intersect(clinical.df$ID,colnames(rna.exp))
clinical.df <- clinical.df[match(patients.filt,clinical.df$ID),]
rna.exp.new <- rna.exp[,match(patients.filt,colnames(rna.exp))]
rna.df.new <- cbind(rna.anno,rna.exp.new)


### Keep the highest gene expression value per gene (ENTREZ ID) in case of moultiple values

filt.tmp <- as.data.frame(cbind(rna.df.new$Entrez_Gene_Id,rowSds(rna.exp.new)))
idx <- order(filt.tmp[,2],decreasing = T) #make ordered index
rna.df.new <- rna.df.new[idx,] #order
rna.df.new <- rna.df.new[!duplicated(rna.df.new$Entrez_Gene_Id),] #filter duplicates
rownames(rna.df.new) <- rna.df.new$Entrez_Gene_Id #Rename rows


### Filter based on gene type

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version=107) 
genes <- getBM(attributes=c('entrezgene_id','entrezgene_accession','gene_biotype',
                            'transcript_biotype','external_gene_name'),
               filters='entrezgene_id',
               values = rna.df.new$Entrez_Gene_Id,
               mart=ensembl) #use biomart to collect annotation information

lncRNAs <- genes[genes$gene_biotype %in% "lncRNA",] #filter biomart result
lncRNAs_entrez <- unique(lncRNAs$entrezgene_id) #store ENTREZ IDs of validated lncRNA genes
genes_to_keep <- lncRNAs_entrez #create a vector of genes to keep
readr::write_delim(as.data.frame(genes_to_keep),"ProcessedDatasets/genes_to_keep.txt") #save gene list
rna.df.new <- rna.df.new[match(genes_to_keep,rownames(rna.df.new)),] #filter expression data set


##############################################################################################
############################### COMBINE CLINICAL WITH EXPRESSION #############################
##############################################################################################

### Provide path and label information
dir="./ProcessedDatasets/"
my_id="tcga_brca_rna_non_transformed_custom"


### Use simple names

clinical <- clinical.df #clinical
ex <- rna.df.new[,-c(1:2)] #expression data without gene annotation

### Merge and save

if(mean(clinical$ID==rownames(t(ex)))==1){
  message("Rownames match, merging clinical and expression data ..")
  combined <- as.data.frame(cbind(clinical,t(ex)))
  message("Saving clinical and expression data ..")
  readr::write_delim(combined,paste0(dir,my_id,"_pheno-expression.txt"),delim ="\t",quote = "none")
  readr::write_excel_csv(combined,paste0(dir,my_id,"_pheno-expression.csv"))
  message("Merged clinical and expression data tables are ready")
}

##############################################################################################
############################################# END ############################################
##############################################################################################

