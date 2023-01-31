



distf <- function(A, B = NULL, method, ncores = NULL, p = NULL) {
  library(parallel)

  A <- as.matrix(A)
  if (is.null(B)) {
    B <- A
  } else {
    B <- as.matrix(B)
  }
  
  nrowmax <- max(nrow(A), nrow(B))
  
  if (is.null(ncores)) {
    if (nrowmax > 100){
      numCores <- 1
      
    } else {
      numCores <- 1
    }
    
  }

  cl <- makeCluster(numCores, type = "PSOCK")
  clusterEvalQ(cl, library(fastDist))
  clusterExport(cl, varlist = c("B"))
  
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregestry"))
  }
  
  q <- nrow(A) %/% numCores
  r <- nrow(A) %% numCores
  
  my_list1 <- lapply(split(A,                    # Split matrix into list
                    c(rep(1:numCores, each = q ), 
                      rep(numCores,   each = r ))
                    ), matrix, ncol = ncol(A))
  
  
  
  as.matrix(do.call(rbind.data.frame, parLapply(cl, my_list1, fdist, method = method, B = B)))
}

fdist <- function(A, B = NULL, method, p = NULL) {
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregestry"))
  }
  A <- as.matrix(A)
  if (is.null(B)) {
    B <- A
  } else {
    B <- as.matrix(B)
  }
  if (is.na(fdistregistry$get_entry(method)$p)) {
    result <- fdistregistry$get_entry(method)$fun(A, B)
  } else {
    result <- fdistregistry$get_entry(method)$fun(A, B, p)
  }
  
  return(result)
}
