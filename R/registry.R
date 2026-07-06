#' Registry of Distance Backends
#'
#' A registry object (created with the \pkg{registry} package) holding the
#' distance backends
#' available to [fdist()]. Each entry stores the `method` name (the lookup
#' key), the C++ backend function `fun`, a short human-readable
#' `description`, and `params`: a named list with the default value of every
#' extra parameter that method accepts (used unless overridden through
#' `...` in [fdist()]). Methods that take no extra parameter simply store
#' an empty list.
#'
#' Currently registered extra parameters are:
#'
#' - `p`: exponent of the Minkowski distance (defaults to `2`).
#' - `radius`: sphere radius used by the Haversine distance (defaults to
#'   `6371`, the mean Earth radius in km).
#' - `base`: logarithm base used by the Jensen-Shannon distance (defaults to
#'   `exp(1)`, i.e. natural log).
#' - `threshold`: value above which a feature is considered "present" when
#'   binarizing data for the Jaccard (defaults to `0`) and Hamming (defaults
#'   to `NA`, meaning raw values are compared instead of binarized) distances.
#' - `regularize`: ridge added to the diagonal of the covariance matrix used
#'   by the Mahalanobis distance before inversion (defaults to `0`).
#' - `weights`: per-feature scale used by the standardized Euclidean
#'   distance instead of the sample variance of `A` (defaults to `NULL`).
#' - `cov`: precomputed covariance matrix used by the Mahalanobis distance
#'   instead of the sample covariance of `A` (defaults to `NULL`).
#'
#' The registry is mainly useful to discover which methods are available and
#' which extra parameters they accept:
#'
#' - `fdistregistry$get_entry_names()` returns the method names accepted by
#'   the `method` argument of [fdist()].
#' - `fdistregistry$get_entry("euclidean")` returns the full entry of a
#'   given method.
#' - `fdistregistry$get_entry("haversine")$params` returns the extra
#'   parameters (with their defaults) accepted by a given method.
#'
#' @format A `registry` object with fields `method`, `fun`, `description`
#'   and `params`.
#'
#' @examples
#' # methods accepted by fdist()
#' fdistregistry$get_entry_names()
#'
#' # description of every registered method
#' vapply(fdistregistry$get_entry_names(),
#'        function(m) fdistregistry$get_entry(m)$description,
#'        character(1))
#'
#' # extra parameters (and their defaults) accepted by a specific method
#' fdistregistry$get_entry("haversine")$params
#' fdistregistry$get_entry("mahalanobis")$params
#'
#' @seealso [fdist()]
#' @export
fdistregistry <- registry::registry()

fdistregistry$set_field("method",
                        type = "character", is_key = TRUE,
                        index_FUN = registry::match_partial_ignorecase)
fdistregistry$set_field("fun",
                        type = "function", is_key = FALSE)
fdistregistry$set_field("description",
                        type = "character", is_key = FALSE)
fdistregistry$set_field("params",
                        is_key = FALSE, default = list())

fdistregistry$set_entry(method = "euclidean",
                        fun    = .euclidean,
                        description = "Euclidean distance")

fdistregistry$set_entry(method = "manhattan",
                        fun    = .manhattan,
                        description = "Manhattan distance")

fdistregistry$set_entry(method = "minkowski",
                        fun    = .minkowski,
                        params = list(p = 2),
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
                        params = list(base = exp(1)),
                        description = "Jensen-Shannon distance")

fdistregistry$set_entry(method = "haversine",
                        fun    = .haversine,
                        params = list(radius = 6371),
                        description = "Haversine (great-circle) distance")

fdistregistry$set_entry(method = "standardized_euclidean",
                        fun    = .standardized_euclidean,
                        params = list(weights = NULL),
                        description = "standardized Euclidean distance")

fdistregistry$set_entry(method = "spearman",
                        fun    = .spearman,
                        description = "Spearman correlation distance")

fdistregistry$set_entry(method = "mahalanobis",
                        fun    = .mahalanobis,
                        params = list(cov = NULL, regularize = 0),
                        description = "Mahalanobis distance")

fdistregistry$set_entry(method = "hamming",
                        fun    = .hamming,
                        params = list(threshold = NA_real_),
                        description = "Hamming distance")

fdistregistry$set_entry(method = "jaccard",
                        fun    = .jaccard,
                        params = list(threshold = 0),
                        description = "Jaccard distance")

fdistregistry$set_entry(method = "gower",
                        fun    = .gower,
                        description = "Gower distance")
