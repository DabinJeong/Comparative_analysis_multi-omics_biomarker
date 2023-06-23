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
#' @examples
#' library(iDRW)
#' library(igraph)
#'
#' g <- directGraph
#' c <- directGraph
#' m <- directGraph
#'
#' gene_delim <- c('g.', 'c.', 'm.') # genes from RNA-Seq gene expression(g), CNV(c), Methylation(m) profile
#'
#' V(g)$name <- paste(gene_delim[1],V(g)$name,sep="")
#' V(c)$name <-paste(gene_delim[2],V(c)$name,sep="")
#' V(m)$name <-paste(gene_delim[3],V(m)$name,sep="")
#'
#' gcm <- (g %du% c) %du% m
#'
#' class.outcome <- "time"
#' covs <- c("age", "gender", "stageT", "stageN", "stageM")
#' family <- "cox"
#'
#' pa <- get.iDRWP(x=list(exp, cna, methyl), y=clinical, globalGraph=gcm, pathSet=pathSet, class.outcome=class.outcome,
#'                 covs=covs, family=family, Gamma=0.3, Corr=FALSE)
#'
#'
#'
#' @export
#'
get.iDRWP <-
  function(x, y, globalGraph, pathSet, class.outcome, covs=NULL, family="cox", Gamma=0.3, Corr=FALSE) {

    x_norm <- list(0)
    x_sig <- list(0)
    gene_weight <- list(0)

    for(i in 1:length(x)) {

      # z-normalize gene profile
      x_norm[[i]] <- apply(x[[i]], 2, function(k) (k - mean(k)) / sd(k) ^ as.logical(sd(k)))

      # calculate the p-value of each gene
      x_sig[[i]] <- get.genes.sig(X=x_norm[[i]], Y=y, class.outcome=class.outcome, covs=covs, family=family)

      cat('dimension of x ', dim(x_norm[[i]]), '\n')

      # initialize gene weights
      geneWeight <- -log(x_sig[[i]][ ,2]+2.2e-16)
      geneWeight[which(is.na(geneWeight))] <- 0
      gene_weight[[i]] <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))

      cat('gene weight initialization Done\n')
    }

    # assign initial weights to the pathway graph
    W0 <- getW0(gene_weight, globalGraph)

    cat('Performing directed random walk...\n')

    # get adjacency matrix of the (integrated) gene-gene graph
    W = getW(G = globalGraph, x = x_norm, Corr=Corr)

    # perform DRW on gene-gene graph
    vertexWeight <- DRW(W = W, p0 = W0, gamma = Gamma)
    names(vertexWeight) <- names(W0)

    x <- Reduce(cbind, x_norm)
    x_sig <- Reduce(rbind, x_sig)

    # pathway activity inference method
    pA <- getPathActivity(x = x, pathSet = pathSet, w = vertexWeight, vertexZP = x_sig)

    return(list(pA, vertexWeight, x_sig))
  }
