options(warn=-1)

library("mixOmics",verbose=FALSE)
library("PRROC")

set.seed(42)

args <- commandArgs(TRUE)

label <- as.character(args[1])
methylome_file <- as.character(args[2])
proteome_file <- as.character(args[3])
transcriptome_file <- as.character(args[4])
clinical_file <- as.character(args[5])
train_samples_file <- as.character(args[6])
test_samples_file <- as.character(args[7])

methylome_raw <- t(read.table(methylome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
proteome_raw <- t(read.table(proteome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
transcriptome_raw <- t(read.table(transcriptome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))

subgroup <- read.table(clinical_file,sep='\t',row.names=1,header=TRUE)

train_samples_df <- read.table(train_samples_file)
test_samples_df <- read.table(test_samples_file)

train_samples <- train_samples_df$"V1"
test_samples <- test_samples_df$"V1"

samples_common_omics <- Reduce(intersect, 
                               list(rownames(methylome_raw), rownames(proteome_raw), rownames(transcriptome_raw)))
# Train model
samples_common_omics_train <- intersect(train_samples, samples_common_omics)

X_train <- list(
                methylome=methylome_raw[samples_common_omics_train,], 
                proteome=proteome_raw[samples_common_omics_train,],
                transcriptome=transcriptome_raw[samples_common_omics_train,])
Y_train <- as.factor(subgroup[as.character(samples_common_omics_train),label])

list.keepX <- list(methylome=c(50,2), proteome=c(50,2), transcriptome=c(50,2))

MyResult.diablo <- block.splsda(X_train, Y_train, keepX=list.keepX)


## Test
samples_common_omics_test <- intersect(test_samples, samples_common_omics)

X_test <- list(
               methylome=methylome_raw[samples_common_omics_test,], 
               proteome=proteome_raw[samples_common_omics_test,],
               transcriptome=transcriptome_raw[samples_common_omics_test,])
Y_test <- as.factor(subgroup[as.character(samples_common_omics_test),label])

Mypredict.diablo <- predict(MyResult.diablo, newdata=X_test)

a = Mypredict.diablo$class$centroids.dist$methylome[,1] == 2
b = Mypredict.diablo$class$centroids.dist$transcriptome[,1] == 2
c = Mypredict.diablo$class$centroids.dist$proteome[,1] == 2

probs = (a+b+c)/3

fg <- probs[Y_test == 2]
bg <- probs[Y_test == 1]

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

file <- file('DIABLO_performance.txt')
writeLines(c(paste0("AUROC: ",roc$auc), paste0("AUPRC: ",pr$auc.integral)), file)
close(file)
