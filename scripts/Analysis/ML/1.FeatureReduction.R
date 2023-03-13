
##############################################################################################################
################################################### START ####################################################
##############################################################################################################

### Load packages

library(genefilter)
library(caret)
library(glmnet)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RColorBrewer)
library(ROCR)
library(precrec)


### Load data

dataset <- as.data.frame(readr::read_delim("../../PrepareDatasets/ProcessedDatasets/tcga_brca_rna_non_transformed_custom_pheno-expression.txt"))


### Define paths and labels

id <- "TCGAsporadicBC"
plot.dir <- paste0("plots/exploratory/")
dir.create(plot.dir,recursive = T)


##############################################################################################################
########################################### PREPROCESSING ####################################################
##############################################################################################################

### Define new subtype labels 

data <- dataset #use shorter name
clin <- data[,grep("[[:alpha:]]",colnames(data))] # filter out genes - gene labels are numerical (ENTREZ IDs)

data <- data[clin$SUBTYPE %in% c("LumA","LumB","Basal","Her2"),]
clin <- clin[clin$SUBTYPE %in% c("LumA","LumB","Basal","Her2"),]

clin$SUBTYPE[clin$SUBTYPE %in% c("LumA","LumB")] <- "Auspicious" #ER +
clin$SUBTYPE[clin$SUBTYPE %in% c("Basal","Her2")] <- "Ominous" #ER - 


### Prepare expression dataset

exp <- data[,grep("[[:alpha:]]",colnames(data),invert = T)]
rownames(exp) <- data$ID

exp0 <- exp
exp <- log2(exp+1)
#exp <- exp[,!colSums(exp == 0) > 0.1*dim(exp)[1]] #remove those with 0 in more than 10% of samples
t_exp <- t(exp)
t_exp0 <- t(exp0)


##############################################################################################################
######################################## EXPLORATORY PLOTS ###################################################
##############################################################################################################

png(paste0(plot.dir,id,"_","30_var_samples.png"),width = 4000, height = 3000, res=450)
BiocGenerics::boxplot(t_exp0[,order(apply(t_exp0,2,sd),decreasing = T)[1:30]],las=2)
title(paste0(id," 30 samples\n with higher SD"))
dev.off()

png(paste0(plot.dir,id,"_","30_var_samples_log2.png"),width = 4000, height = 3000, res=450)
BiocGenerics::boxplot(t_exp[,order(apply(t_exp,2,sd),decreasing = T)[1:30]],las=2)
title(paste0(id," 30 samples\n with higher SD"))
dev.off()

png(paste0(plot.dir,id,"_","30_var_genes.png"),width = 4000, height = 3000, res=450)
BiocGenerics::boxplot(exp[,order(apply(exp,2,sd),decreasing = T)[1:30]],las=2)
title(paste0(id," 30 genes\n with higher SD"))
dev.off()


### Define group labels

group <- as.factor(clin$SUBTYPE)
group = relevel(group, ref="Auspicious")


### PCA

pca = prcomp(exp)
smr <- summary(pca)
var.prop <- round(smr$importance[2,]*100,2)


png(paste0(plot.dir,id,"_","first_pca.png"),width = 4000, height = 3000, res=600)
# Take PC1 and PC2 for the plot
plot(pca$x[,1:2],col=group, pch=19,
     xlab=paste0("PC1: ",var.prop[1],"%"),
     ylab=paste0("PC2:",var.prop[2],"%"),
     main=paste0("PCA using all features (n=",ncol(exp),")"))
# include a legend for points
legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
dev.off()


##############################################################################################################
############################################ STANDARDIZE #####################################################
##############################################################################################################

# Make it into a matrix object and standardize
d_mat0 = as.matrix(exp)
d_mat=scale(d_mat0)

# As before, we want this to be a factor
d_resp = group

png(paste0(plot.dir,id,"_","30_var_genes_unscaled.png"),width = 4000, height = 3000, res=450)
BiocGenerics::boxplot(d_mat0[,order(apply(d_mat0,2,sd),decreasing = T)[1:30]],las=2)
title(paste0(id," 30 genes\n with higher SD unscaled"))
dev.off()

png(paste0(plot.dir,id,"_","30_var_genes_scaled.png"),width = 4000, height = 3000, res=450)
BiocGenerics::boxplot(d_mat[,order(apply(d_mat,2,sd),decreasing = T)[1:30]],las=2)
title(paste0(id," 30 genes\n with higher SD scaled"))
dev.off()


d_mat=d_mat0


##############################################################################################################
##################################### SPLIT DATASET TO TRAIN #################################################
##############################################################################################################

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(1)
train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]


##############################################################################################################
###################################### HYPERPARAMETER TUNING #################################################
##############################################################################################################

id <- "TCGAsporadicBC"
alphas <- seq(0,1,0.1)

set.seed(1)
foldid <- sample(1:10, size = length(y_train), replace = TRUE)

plot.dir <- paste0("plots/cross_validation/")
dir.create(plot.dir,recursive = T)

res.dir <- paste0("results/cross_validation/")
dir.create(res.dir,recursive = T)

for(a in alphas){
  cv <- cv.glmnet(x_train, y_train,foldid = foldid, alpha = a,type="class",family="binomial")
  saveRDS(cv,paste0(res.dir,"alpha_",a,"_",id,".RDS"))
  png(paste0(plot.dir,id,"_","cross-validation","_alpha_",a,".png"),width = 4000, height = 3000, res=450)
  plot(cv)
  title(paste0(id," alpha=",a),line=-2)
  dev.off()
}


rds_list <- list.files(res.dir)[grep(".RDS",list.files(res.dir))]
cv_alpha_stats <- matrix(0,length(rds_list),6)
cv_alpha_stats <- as.data.frame(cv_alpha_stats)
colnames(cv_alpha_stats) <- c("Alpha","Lamba","Index","MinClassError","SE","Nonzero")

for(rds in rds_list){
  idx<- which(rds_list == rds)
  cv <- readRDS(paste0(res.dir,rds))
  res_coef = coef(cv, s="lambda.min")
  stats <- as.data.frame(print(cv))
  A <- strsplit(rds,"_",fixed=T)[[1]][2]
  stat_sum <- c(A,unlist(unname(stats[1,])))
  cv_alpha_stats[idx,] <- stat_sum
}

#rds <- "alpha_1_TCGAsporadicBC.RDS"
#cv <- readRDS(paste0(res.dir,rds))
#res_coef = coef(cv, s="lambda.min")

cv_alpha_stats <- as.data.frame(cv_alpha_stats)
cv_alpha_stats <- as.data.frame(apply(cv_alpha_stats, 2, as.numeric))
cv_alpha_stats <- cv_alpha_stats[order(cv_alpha_stats$Alpha),]

readr::write_delim(cv_alpha_stats,paste0(res.dir,id,"_cv_alpha_tuning_stats.txt"))
readr::write_excel_csv(cv_alpha_stats,paste0(res.dir,id,"_cv_alpha_tuning_stats.csv"))

cv_alpha_stats_less_predictors_top3 <- cv_alpha_stats[order(cv_alpha_stats$Nonzero,decreasing = F),][1:3,]
cv_alpha_stats_best <- cv_alpha_stats_less_predictors_top3[order(cv_alpha_stats_less_predictors_top3$MinClassError,decreasing = F),][1,]
best_alpha <- as.numeric(cv_alpha_stats_best[,"Alpha"])

res.dir <- "results/final_model/"
dir.create(res.dir)


##############################################################################################################
########################################### MODEL FITTING ####################################################
##############################################################################################################

### Proceed with best alpha
tfit = glmnet(
  x = x_train,
  y = y_train,
  alpha = best_alpha,
  family = "binomial"
)

png(paste0(res.dir,id,"_","best_model_coefs.png"),width = 4000, height = 3000, res=450)
plot(tfit)
dev.off()

res = cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = best_alpha,
  foldid = foldid,
  family = "binomial",
  type="class")

png(paste0(res.dir,id,"_","best_model_lambda_cross-validation.png"),width = 4000, height = 3000, res=450)
plot(res)
title(paste0(id," alpha=",best_alpha),line=-2)
dev.off()


##############################################################################################################
################################# FIRST EVALUATION OF MODEL PERFORMANCE ######################################
##############################################################################################################

# Test/Make prediction on test dataset
res$lambda.min
y_pred = predict(res, newx=x_test, type="class", s="lambda.min")

confusion_matrix = table(y_pred, y_test)
confusionMatrix(confusion_matrix)

# Evaluation statistics
print(confusion_matrix)
print(paste0("Sensitivity: ",sensitivity(confusion_matrix)))
print(paste0("Specificity: ",specificity(confusion_matrix)))
print(paste0("Precision: ",precision(confusion_matrix)))

res_coef = coef(res, s="lambda.min") # the "coef" function returns a sparse matrix
dim(res_coef)

head(res_coef)

# get coefficients with non-zero values
res_coef = res_coef[res_coef[,1] != 0,]
# note how performing this operation changed the type of the variable
head(res_coef)

# remove first coefficient as this is the intercept, a variable of the model itself
res_coef = res_coef[-1]

relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes


relevant_gene_names <- mapIds(org.Hs.eg.db,keys = relevant_genes,column = "SYMBOL",keytype = "ENTREZID" )

gene_list <- cbind(relevant_gene_names,relevant_genes)

readr::write_excel_csv(as.data.frame(gene_list),paste0(res.dir,id,"_","selected_features.csv"))


###### PCA2 


pca.slc = prcomp(exp[,relevant_genes])
smr.slc <- summary(pca.slc)
var.prop.slc <- round(smr.slc$importance[2,]*100,2)

png(paste0(res.dir,id,"_","selected_feature_pca.png"),width = 4000, height = 3000, res=600)
# Take PC1 and PC2 for the plot
plot(pca.slc$x[,1:2],col=group, pch=19,
     xlab=paste0("PC1: ",var.prop.slc[1],"%"),
     ylab=paste0("PC2:",var.prop.slc[2],"%"),
     main=paste0("PCA using resulting features (n=",length(relevant_genes),")"))
# include a legend for points
legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
dev.off()


##############################################################################################################
############################# SPOT EXPRESSION PATTERNS OF INFORMATIVE GENES ##################################
##############################################################################################################


# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)



library(gplots)
# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  #RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)
dev.off()
png(paste0(res.dir,id,"_","heatmap.png"),width = 3000, height = 4000, res=600)
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  #RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)
dev.off()
# Extract the hierarchical cluster from heatmap to class "hclust"
hc = as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters = cutree(hc, k=2)
table(clusters)

cluster1_entrez <- names(which(clusters == 1))
cluster1_symbol <- mapIds(org.Hs.eg.db,keys = cluster1_entrez,column = "SYMBOL",keytype = "ENTREZID" )
cluster2_entrez <- names(which(clusters == 2))
cluster2_symbol <- mapIds(org.Hs.eg.db,keys = cluster2_entrez,column = "SYMBOL",keytype = "ENTREZID" )

cluster1 <- cbind(cluster1_symbol,cluster1_entrez)
cluster2 <- cbind(cluster2_symbol,cluster2_entrez)
colnames(cluster1) <- c("SYMBOL","ENTREZ")
colnames(cluster2) <- c("SYMBOL","ENTREZ")

readr::write_delim(as.data.frame(cluster1),paste0(res.dir,id,"_cluster1.txt"))
readr::write_excel_csv(as.data.frame(cluster1),paste0(res.dir,id,"_cluster1.csv"))
readr::write_delim(as.data.frame(cluster2),paste0(res.dir,id,"_cluster2.txt"))
readr::write_excel_csv(as.data.frame(cluster2),paste0(res.dir,id,"_cluster2.csv"))


##############################################################################################################
################################ SECOND EVALUATION OF MODEL PERFORMANCE ######################################
##############################################################################################################

cfit <-res
cnf <- confusion.glmnet(cfit, newx = x_test, newy = y_test)
print(cnf)

saveRDS(cnf,paste0(res.dir,id,"_Confussion_matrix.RDS"))


# Run model over training dataset
lasso.model <- cv.glmnet(x = x_train, y = y_train,family = 'binomial',foldid = foldid,alpha=best_alpha, type.measure = 'auc')

# Apply model to testing dataset
lasso.prob <- predict(lasso.model,type="response",newx = x_test, s = 'lambda.min')
pred <- prediction(lasso.prob, y_test)

#ROC
precrec_obj <- evalmod(scores = pred@predictions, labels = pred@labels)

png(paste0(res.dir,id,"_ROC.png"),width = 4000, height = 3000, res=600)
autoplot(precrec_obj)
dev.off()

roc.res <- as.data.frame(attr(precrec_obj,"auc"))
readr::write_delim(roc.res,paste0(res.dir,id,"_roc.results.txt"))
readr::write_excel_csv(roc.res,paste0(res.dir,id,"_roc.results.csv"))





