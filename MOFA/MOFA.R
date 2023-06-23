options(warn=-1)

library("reticulate")
library("PRROC")
library("glmnet")
predict <- predict.glmnet
library("MOFA2")
library("igraph")
set.seed(42)

library(reticulate)
use_python("/data/project/dabin/conda_env/asthma_R/bin/python", required=TRUE)
use_condaenv("/data/project/dabin/conda_env/asthma_R", required=TRUE)
args <- commandArgs(TRUE)

label <- as.character(args[1])
methylome_file <- as.character(args[2])
proteome_file <- as.character(args[3])
transcriptome_file <- as.character(args[4])
clinical_file <- as.character(args[5])
train_samples_file <- as.character(args[6])
test_samples_file <- as.character(args[7])


#########Preprocessing###########
methylome_raw <- t(read.table(methylome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
transcriptome_raw <- t(read.table(transcriptome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))
proteome_raw <- t(read.table(proteome_file,row.names=1,header=TRUE,check.names=FALSE,sep='\t'))

clinical <- read.table(clinical_file,sep='\t',row.names=1,header=TRUE)
clinical[,label] <- list(sapply(clinical[,label], function(x) if (x==1) 1 else 0))

train_samples_df <- read.table(train_samples_file)
test_samples_df <- read.table(test_samples_file)

train_samples <- train_samples_df$"V1"
test_samples <- test_samples_df$"V1"

samples_common_omics <- Reduce(intersect,
                               list(rownames(methylome_raw), rownames(proteome_raw), rownames(transcriptome_raw), rownames(clinical)))

# Dataset
methylome=t(methylome_raw[samples_common_omics,])
transcriptome=t(transcriptome_raw[samples_common_omics,])
proteome=t(proteome_raw[samples_common_omics,])

################ MOFA for sample representation as factor #####################################
data <- list(view_1=methylome, view_2=transcriptome, view_3=proteome)
MOFAobject <-create_mofa(data)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("trained_model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE) 

model <- load_model(outfile)
factors <- get_factors(model, factors="all")

######################LASSO prediction model######################
# Train
samples_common_omics_train <- intersect(train_samples, samples_common_omics)
samples_common_omics_test <- intersect(test_samples, samples_common_omics)

x_train = factors$group1[samples_common_omics_train,]
y_train = as.vector(clinical[samples_common_omics_train,label])
x_test = factors$group1[samples_common_omics_test,]
y_test = as.vector(clinical[samples_common_omics_test,label])

cv_model <- cv.glmnet(x_train, y_train, alpha = 0, nfolds=5, type.measure='class', family='binomial')
best_lambda <- cv_model$lambda.1se

best_model <- glmnet(x_train, y_train, alpha = 0, lambda = best_lambda, family='binomial')

# Test
y_predicted <- predict(best_model, s=best_lambda, newx=x_test, type='response')

# Evaluation
fg = y_predicted[y_test==1,]
bg = y_predicted[y_test==0,]

# ROC Curve    
#roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = y_predicted, weights.class0 = y_test, curve = T)
# PR Curve
#pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr <- pr.curve(scores.class0 = y_predicted, weights.class0 = y_test, curve = T)
file <- file('MOFA_performance.txt')
writeLines(c(paste0("AUROC: ",roc$auc), paste0("AUPRC: ",pr$auc.integral)), file)
close(file)

