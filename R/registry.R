library(registry)
fdistregistry <- registry()

fdistregistry$set_field("method", 
                        type = "character", 
                        is_key = TRUE,
                        index_FUN = match_partial_ignorecase)
fdistregistry$set_field("fun", 
                        type = "function", 
                        is_key = FALSE)
fdistregistry$set_field("p", 
                        type = "numeric", 
                        is_key = FALSE)

Rcpp::sourceCpp("src/fastDist.cpp")
fdistregistry$set_entry(method = "euclidean",  
                        fun    = euclidean,
                        p      = NA)
fdistregistry$set_entry(method = "manhattan",  
                        fun    = manhattan,
                        p      = NA)
fdistregistry$set_entry(method = "minkowski",  
                        fun    = minkowski,
                        p      = 2)


