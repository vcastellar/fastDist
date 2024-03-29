library(registry)
fdistregistry <- registry()

fdistregistry$set_field("method", 
                        type = "character", is_key = TRUE,
                        index_FUN = match_partial_ignorecase)
fdistregistry$set_field("fun", 
                        type = "function", is_key = FALSE)
fdistregistry$set_field("p", 
                        type = "numeric", is_key = FALSE)
fdistregistry$set_field("description", 
                         type = "character", is_key = FALSE)

Rcpp::sourceCpp("src/fastDist.cpp")
fdistregistry$set_entry(method = "euclidean",  
                        fun    = .euclidean,
                        description = "distancia euclidea")

fdistregistry$set_entry(method = "manhattan",  
                        fun    = .manhattan,
                        description = "distancia de manhattan")

fdistregistry$set_entry(method = "minkowski",  
                        fun    = .minkowski,
                        p      = 2,
                        description = "distancia de Minkowski")

fdistregistry$set_entry(method = "canberra",  
                        fun    = .canberra,
                        description = "distancia de canberra")


fdistregistry$set_entry(method = "supremum",  
                        fun    = .supremum,
                        description = "distancia del supremo")

fdistregistry$set_entry(method = "mahalanobis",  
                        fun    = .mahalanobis,
                        description = "distancia de Mahalanobis")





