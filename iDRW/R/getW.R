getW <- function(G, x, Corr=FALSE) {

  # number of layers
  L = length(x)

  if(L > 1) {
    # within-layer interactions
    if(!Corr) G <- as.matrix(get.adjacency(G))

    # between-layer interactions
    for(i in 1:(L-1)) {
      for (j in (i+1):L) {
        G <- add_edges_btw(G=G, from=x[[i]], to=x[[j]],reverse = T, Corr=Corr)
      }
    }

    if(Corr) W <- as.matrix(get.adjacency(G))
    else W <- G

  } else {
    W <- as.matrix(get.adjacency(G))
  }

  cat('number of nodes', dim(W), '\n')
  cat('number of edges', sum(W), '\n')

  cat('Adjacency matrix W complete')

  return(W)
}

add_edges_btw <- function(G, from, to, reverse=FALSE, Corr=FALSE, theta=0.5) {

  if(!Corr) {
    nE <- sum(G)

    from <- intersect(colnames(from), colnames(G))
    to <- intersect(colnames(to), colnames(G))

    G[from, to] <- 1
    if(reverse) G[to, from] <- 1

    cat('added',sum(G)-nE,'edges\n')

  } else {
    corr_r <- get_corr(from,to,val='r')
    edge_list <- data.frame(from=rownames(corr_r$ft)[which(abs(corr_r$ft)>=theta, arr.ind = T)[,1]],
                            to=colnames(corr_r$ft)[which(abs(corr_r$ft)>=theta, arr.ind = T)[,2]])
    if(reverse) {
      edge_list <- rbind(edge_list, data.frame(from=rownames(corr_r$tf)[which(abs(corr_r$tf)>=theta, arr.ind = T)[,1]],
                                               to=colnames(corr_r$tf)[which(abs(corr_r$tf)>=theta, arr.ind = T)[,2]]))
    }

    nE <- gsize(G)
    G <- add.edges(G, get.edgelist(graph.data.frame(edge_list, directed = T)))
    cat('added',gsize(G)-nE,'edges\n')
  }

  return(G)
}

get_corr <- function(from, to, val='p'){
  corr <- rcorr(x=from,y=to,type = "pearson")
  if(val=='p'){
    corr_p <- corr$P
    return(list(ft=corr_p[colnames(from),colnames(to)],tf=corr_p[colnames(to),colnames(from)]))

  } else if(val=='r') {
    corr_r <- corr$r
    return(list(ft=corr_r[colnames(from),colnames(to)],tf=corr_r[colnames(to),colnames(from)]))
  }
}


