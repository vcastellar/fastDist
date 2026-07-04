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

fdistregistry$set_entry(method = "squared_euclidean",
                        fun    = .squared_euclidean,
                        description = "squared Euclidean distance")

fdistregistry$set_entry(method = "bray_curtis",
                        fun    = .braycurtis,
                        description = "Bray-Curtis distance")

fdistregistry$set_entry(method = "hellinger",
                        fun    = .hellinger,
                        description = "Hellinger distance")

fdistregistry$set_entry(method = "chi_squared",
                        fun    = .chisquared,
                        description = "chi-squared distance")

fdistregistry$set_entry(method = "jensen_shannon",
                        fun    = .jensenshannon,
                        description = "Jensen-Shannon distance")

fdistregistry$set_entry(method = "haversine",
                        fun    = .haversine,
                        description = "Haversine (great-circle) distance")

fdistregistry$set_entry(method = "standardized_euclidean",
                        fun    = .standardized_euclidean,
                        description = "standardized Euclidean distance")

fdistregistry$set_entry(method = "spearman",
                        fun    = .spearman,
                        description = "Spearman correlation distance")

fdistregistry$set_entry(method = "mahalanobis",
                        fun    = .mahalanobis,
                        description = "Mahalanobis distance")

fdistregistry$set_entry(method = "hamming",
                        fun    = .hamming,
                        description = "Hamming distance")

fdistregistry$set_entry(method = "jaccard",
                        fun    = .jaccard,
                        description = "Jaccard distance")

fdistregistry$set_entry(method = "gower",
                        fun    = .gower,
                        description = "Gower distance")




