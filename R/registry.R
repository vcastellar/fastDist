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

fdistregistry$set_entry(method = "euclidean",  
                        fun    = .euclidean,
                        description = "Euclidean distance")

fdistregistry$set_entry(method = "manhattan",  
                        fun    = .manhattan,
                        description = "Manhattan distance")

fdistregistry$set_entry(method = "minkowski",  
                        fun    = .minkowski,
                        p      = 2,
                        description = "Minkowski distance")

fdistregistry$set_entry(method = "correlation",  
                        fun    = .correlation,
                        description = "correlation distance")

fdistregistry$set_entry(method = "cosine",  
                        fun    = .cosine,
                        description = "cosine distance")

fdistregistry$set_entry(method = "canberra",  
                        fun    = .canberra,
                        description = "Canberra distance")


fdistregistry$set_entry(method = "supremum",  
                        fun    = .supremum,
                        description = "supremum distance")

fdistregistry$set_entry(method = "mahalanobis",  
                        fun    = .mahalanobis,
                        description = "Mahalanobis distance")





