get.genes.sig <- function(X, Y, class.outcome, covs=NULL, family="cox") {

  library(survival)
  library(stats)

  cat('Calculating gene signature...')

  Survdata <- data.frame(X, Y)

  stats <- matrix(NA,nrow=ncol(X),ncol=2)
  rownames(stats) <- colnames(X)
  colnames(stats) <- c("Score","p-value")

  if(!is.null(covs)) add_strata <- paste("+strata(", paste(covs, collapse = ","), ")")
  else add_strata <- NULL

  for(i in 1 : ncol(X)){
    if(family == "cox"){
      # calculate the signature of each gene using Cox PH model
      res.coxph <- coxph(as.formula(paste("Surv(time, status)~ ", colnames(Survdata)[i], add_strata)), Survdata)
      stats[i,] <- summary(res.coxph)$coefficients[c(4,5)]

    } else if(family == "binomial") {
      # calculate the signature of each gene using t-test
      stats_tmp <- t.test(X[Y[,class.outcome]==0,i], X[Y[,class.outcome]==1,i], var.equal=TRUE)
      stats[i,1] <- stats_tmp$statistic
      stats[i,2] <- stats_tmp$p.value

    } else if(family == "multinomial") {
      res.aov <- aov(as.formula(paste(colnames(Survdata)[i], "~ as.factor(",class.outcome,")")), data=Survdata)
      stats[i,] <- as.numeric(summary(res.aov)[[1]][1,4:5])
    }
  }

  cat('Done\n')
  return(stats)
}
