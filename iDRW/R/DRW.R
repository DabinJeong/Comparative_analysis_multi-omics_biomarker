DRW <- function(W, p0, gamma=0.3){

    # Add Ground Node, construct new adjacent matrix
    newrow <- matrix(1,1,dim(W)[2])
    rownames(newrow) <- c("GN")
    W1 <- rbind(W,newrow)
    newcol <- matrix(1,dim(W1)[1],1)
    colnames(newcol) <- c("GN")
    WGN <- cbind(W1,newcol)   # adjacency matrix after adding ground node
    p0 <- t(as.matrix(p0/sum(p0)))

    # The initial probability of the ground node is 0.
    p0 <- cbind(p0,0)
    colnames(p0)[dim(p0)[2]] = "GN"

    PT <- p0

    k <- 0
    delta <- 1

    # reverse the direction of the edges
    WGN <- t(WGN)

    # normalization
    for (i in 1:dim(WGN)[1]){
      sumr <- sum(WGN[i,])
      if(sumr == 0){
        WGN[i,] <-numeric(length=length(WGN[i,]))
      }
      if(sumr > 0){
        WGN[i,] <- WGN[i,]/sumr
      }
    }
    WGN <- t(WGN)

    # iteration
    while(delta > 1e-10){
      PT1 <- (1-gamma)*WGN
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*p0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      cat(delta," ")
      PT <- PT4
      k <- k + 1
    }
    cat('converged\n')

    PT <- t(PT)
    rownames(PT) <- NULL
    res <- drop(PT[1:(dim(PT)[1]-1)])

    return(res)
  }
