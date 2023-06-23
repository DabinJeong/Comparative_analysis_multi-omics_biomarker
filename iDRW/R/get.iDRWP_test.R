#' @title get iDRW-generated pathway profile
#'
#' @description This function infers pathway activities across samples using integrative directed random walk approach (iDRW) on multi-layered gene-gene graph.
#'
#' @param x list of feature matrices (genomic profiles). The length of \code{x} >= 1. Each \code{x} should be samples as rows and genes (features) as columns.
#' @param y clinical matrix that contains response variable. The type of response variable can be either survival time or binary/multinomial response.
#' @param globalGraph directed multi-layered gene-gene graph.
#' @param pathSet The set of pathways.
#' @param class.outcome The variable name of response variable. For \code{class.outcome} is survival time \code{time}, two columns of survival time (\code{time}) and event status (\code{status}) should be provided.
#' @param covs A vector of variable names of covariates (optional).
#' @param family The type of response variable - \code{cox} (default) / \code{binomial} / \code{multinomial}.
#' @param Gamma The restart probability of random walks (default: \code{0.3}).
#' @param Corr A logical variable to indicate between-layer interactions are defined by correlation between genes (default: \code{FALSE}).
#'
#' @return A list of following elements.
#' \itemize{
#' \item pathActivity : A pathway profile inferred by iDRW (samples x pathways).
#' \item sigGenes : A list of the set of significant genes within pathways.
#' \item w : A final weight vector of genes by the product of directed random walk with restart on graphs.
#' \item pathSet : A list of the set of genes within pathways.
#' \item vertexZP : A data frame containing statistic score (\code{Score}) and p-value of the statistical test (\code{p-value}) for each gene.
#' }
#'
#' @export
#'
get.iDRWP_test <-
  function(x, vertexWeight, x_sig, pathSet) {
    x_norm <- list(0)                                                                                          
                                                                                                                     
    for(i in 1:length(x)) {                                                                                          
      # z-normalize gene profile                                                                                     
      x_norm[[i]] <- apply(x[[i]], 2, function(k) (k - mean(k)) / sd(k) ^ as.logical(sd(k)))                         
    }
    x <- Reduce(cbind, x_norm)   

    # pathway activity inference method
    pA <- getPathActivity(x = x, pathSet = pathSet, w = vertexWeight, vertexZP = x_sig)

    return(pA)
  }
