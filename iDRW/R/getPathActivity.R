getPathActivity <-
function(x, pathSet, w, vertexZP){

  cat('Infer pathway activities...')

	# Infer pathway activity
	pathActivity <- c()
	sigGenes <- vector("list", length = length(pathSet))  # The significant genes to compute pathway activity
	names(sigGenes) <- names(pathSet)

    # vertex is weighted with DRW method
    for (i in 1 : length(pathSet)){
        Vpathwayi <- pathSet[[i]]
        if (length(Vpathwayi) > 0){
          n <- 0    # the number of differential genes in ith pathway
          pathActivity_tmp <- matrix(nrow=dim(x)[1],ncol=1,data=0)
          sigGenesi <- c()
          Idx_pathwayi <- c()
          for (j in 1 : length(Vpathwayi)){
            Idxs <- which(substring(colnames(x),3)==Vpathwayi[j])
            for (Idx in Idxs){
                if (length(Idx) > 0){
                  if ( colnames(x)[Idx] %in% names(w)){
                    idx <- which(vertexZP[colnames(x)[Idx],"p-value"] < 0.05)
                    if(length(idx) > 0){
                      for (k in 1:length(idx)) {
                        pathActivity_tmp <- pathActivity_tmp + sign(vertexZP[colnames(x)[Idx[idx[k]]],"Score"]) * w[colnames(x)[Idx[idx[k]]]] * x[,Idx[idx[k]]]
                        n <- n + 1
                        Idx_pathwayi <- rbind(Idx_pathwayi,Idx[idx[k]])
                        sigGenesi <- c(sigGenesi, colnames(x)[Idx[idx[k]]])
                      }
                    }
                  }
                }
            }
          }
          if(n > 0){
            pathActivity_tmp <- pathActivity_tmp / sqrt(sum(w[colnames(x)[Idx_pathwayi]]^2))
            colnames(pathActivity_tmp) <- names(pathSet)[i]
            pathActivity <- cbind(pathActivity, pathActivity_tmp)
            sigGenes[[i]] <- sigGenesi
          }
        }
    }

	rownames(pathActivity) <- rownames(x)

	return(list(pathActivity=pathActivity, sigGenes=sigGenes, w=w, pathSet=pathSet, vertexZP=vertexZP))
}
