getW0 <- function(gene_weight, globalGraph){
    
    # initial gene weight vector W0
    Vertexs <- V(globalGraph)
    W0 <- rep(0, length(Vertexs))
    names(W0) <- Vertexs$name
    
    for(p in 1:length(gene_weight)) {
      for(i in 1:length(W0)){
        idx <- which(names(gene_weight[[p]]) == names(W0[i]))
        if(length(idx) > 0){
          W0[i] <- gene_weight[[p]][idx]
        }
      }
    }
    
    W0 <- W0/sum(W0)
    
    return(W0)
  }
