options(warn=-1)

library("PRROC")
library("glmnet")
library("iDRW")
library("igraph")
set.seed(42)

args <- commandArgs(TRUE)

label <- as.character(args[1])
methylome_file <- as.character(args[2])
proteome_file <- as.character(args[3])
transcriptome_file <- as.character(args[4])
clinical_file <- as.character(args[5])
train_samples_file <- as.character(args[6])
test_samples_file <- as.character(args[7])

################################Load KEGG ###############################
load(file="/tools/iDRW/data/directGraph.KEGGgraph.rda")
load(file="/tools/iDRW//data/pathSet.KEGGgraph.rda")

m <- directGraph 
g <- directGraph
p <- directGraph

gene_delim <- c('m.', 'g.', 'p.') # genes from Methyl-seq methylation(m), RNA-Seq gene expression(g), Protein(p) profile

V(m)$name <- paste(gene_delim[1],V(m)$name,sep="")
V(g)$name <-paste(gene_delim[2],V(g)$name,sep="")
V(p)$name <-paste(gene_delim[3],V(p)$name,sep="")

gcm <- (m %du% g) %du% p

########################Preprocessing######################################
methylome_raw <- t(read.table(methylome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
transcriptome_raw <- t(read.table(transcriptome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
proteome_raw <- t(read.table(proteome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))

colnames(methylome_raw) = paste(gene_delim[1], colnames(methylome_raw), sep="")
colnames(transcriptome_raw) = paste(gene_delim[2], colnames(transcriptome_raw), sep="")
colnames(proteome_raw) = paste(gene_delim[3], colnames(proteome_raw), sep="")

clinical <- read.table(clinical_file,sep='\t',row.names=1,header=TRUE)
clinical[,label] <- list(sapply(clinical[,label], function(x) if (x==1) 1 else 0))

train_samples_df <- read.table(train_samples_file)
test_samples_df <- read.table(test_samples_file)

train_samples <- train_samples_df$"V1"
test_samples <- test_samples_df$"V1"

samples_common_omics <- Reduce(intersect,
                               list(rownames(methylome_raw), rownames(proteome_raw), rownames(transcriptome_raw), rownames(clinical)))

# Train dataset
samples_common_omics_train <- intersect(train_samples, samples_common_omics)

methylome_train=methylome_raw[samples_common_omics_train,]
transcriptome_train=transcriptome_raw[samples_common_omics_train,]
proteome_train=proteome_raw[samples_common_omics_train,]

## Test dataset
samples_common_omics_test <- intersect(test_samples, samples_common_omics)

methylome_test=methylome_raw[samples_common_omics_test,]
transcriptome_test=transcriptome_raw[samples_common_omics_test,]
proteome_test=proteome_raw[samples_common_omics_test,]

################iDRW for sample representation as pathway activation#####################################
train_data <- list(methylome_train,transcriptome_train,proteome_train)
test_data <- list(methylome_test,transcriptome_test,proteome_test)
                                
res_train <- get.iDRWP(x=train_data, y=clinical[samples_common_omics_train,], globalGraph=gcm, pathSet=pathSet, class.outcome=label,covs=NULL, family="binomial", Gamma=0.3, Corr=FALSE)

pa_train <- res_train[[1]]
vertexWeight_train <- res_train[[2]]
vertexZP_train <- res_train[[3]]
                                
pa_test <- get.iDRWP_test(x=test_data, vertexWeight_train, vertexZP_train, pathSet)

######################LASSO prediction model######################
# Train
x_train = pa_train$pathActivity[samples_common_omics_train,]
y_train = as.vector(clinical[samples_common_omics_train,label])
x_test = pa_test$pathActivity[samples_common_omics_test,]
y_test = as.vector(clinical[samples_common_omics_test,label])

cv_model <- cv.glmnet(x_train, y_train, alpha=0, nfolds=5, type.measure='class', family='binomial')
best_lambda <- cv_model$lambda.1se

best_model <- glmnet(x_train, y_train, alpha = 0, lambda = best_lambda, family='binomial')

# Test
y_predicted <- predict(best_model, s=best_lambda, newx=x_test, type='response')
                                
#model <- glmnet(x_train, y_train, type.measure='class', family='binomial')
#best_lambda <- model$lambda.min
#
#fitted_model <- glmnet(x_train, y_train, family='binomial')
#
### Test
#y_predicted <- predict(fitted_model, s=best_lambda, newx=x_test, type='response')

# Evaluation
fg = y_predicted[y_test==1,]
#bg = y_predicted[y_test==0,]

# ROC Curve    
#roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = y_predicted, weights.class0 = y_test, curve = T)
# PR Curve
#pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr <- pr.curve(scores.class0 = y_predicted, weights.class0 = y_test, curve = T)

file <- file('iDRW_performance.txt')
writeLines(c(paste0("AUROC: ",roc$auc), paste0("AUPRC: ",pr$auc.integral)), file)
close(file)
