#' Probabilistic Pathway Scoring
#'
#' Computes probabilistic activity scores for predefined biological pathways
#' based on count data, using either a negative binomial (NB) or zero-inflated
#' negative binomial (ZINB) model. Supports multithreading and batch correction.
#' @param counts A numeric matrix or data frame of gene expression counts
#' @param method Character string specifying the statistical model to use.
#'   Options are:
#'   \describe{
#'     \item{"nb"}{Negative Binomial model (default).}
#'     \item{"zinb"}{Zero-Inflated Negative Binomial model.}
#'   }
#' @param ncores Integer specifying the number of CPU cores to use for parallel
#'   computation. Default is 4.
#' @param batch Optional vector or factor indicating batch assignment of samples,
#'   used for batch-effect correction. Default is \code{NULL}.
#' @param offset Optional numeric matrix of the same dimension as \code{counts},
#'   representing model offsets (e.g., library size normalization factors).
#'   Default is \code{NULL}.
#' @param maxit Integer specifying the maximum number of iterations allowed
#'   for model fitting. Default is 30.
#' @param PathwaySet Character vector specifying which KEGG (or other database)
#'   pathway categories to include in scoring. Default includes:
#'   \itemize{
#'     \item "Metabolism"
#'     \item "Genetic Information Processing"
#'     \item "Environmental Information Processing"
#'     \item "Cellular Processes"
#'     \item "Organismal Systems"
#'     \item "Human Diseases"
#'     \item "Drug Development"
#'   }
#'
#' @returns A list or data frame containing:
#'   \describe{
#'     \item{\code{ScoreMatrix}}{A matrix of pathway activity scores per sample.}
#'     \item{\code{ModelParams}}{Estimated model parameters for each gene.}
#'     \item{\code{Convergence}}{Logical vector indicating model convergence.}
#'   }
#' @importFrom future plan
#' @importFrom doFuture registerDoFuture
#' @importFrom foreach foreach %dopar%
#' @importFrom glmmTMB glmmTMB glmmTMBControl sigma nbinom2
#' @importFrom graphite pathways pathwayGraph
#' @importFrom KEGGREST keggGet
#' @importFrom clusterProfiler bitr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr inner_join group_by mutate ungroup summarise %>%
#' @importFrom igraph page.rank
#' @export
#' @examples
#' \dontrun{
#' data(example_counts)
#' result <- ProbScorer(counts = example_counts,
#'                      method = "nb",
#'                      ncores = 8,
#'                      batch = NULL,
#'                      maxit = 50)
#' head(result$ScoreMatrix)
#' }


ProbScorer <- function(counts,
                       method = c("nb", "zinb"),
                       ncores = 4,
                       batch = NULL,
                       offset = NULL,
                       maxit = 30,
                       PathwaySet = c("Metabolism",
                                      "Genetic Information Processing",
                                      "Environmental Information Processing",
                                      "Cellular Processes",
                                      "Organismal Systems",
                                      "Human Diseases",
                                      "Drug Development"),
                       memory_size = 4 * 1024^3
                       ) {

  # Input validation
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("counts must be a matrix or Matrix object")
  }

  if (is.null(rownames(counts))) {
    stop("counts must have row names (gene names)")
  }

  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("Cell_", seq_len(ncol(counts)))
  }

  # Match method argument
  method <- match.arg(method)

  # Set future globals max size to handle large objects in parallel
  old_option <- options(future.globals.maxSize = memory_size)
  on.exit(options(old_option), add = TRUE)

  # Get pathway weights (assuming this function exists)
  if (!exists("get_pathway_weights")) {
    stop("get_pathway_weights function not found. Please make sure it's available in your environment.")
  }

  weighs <- get_pathway_weights(PathwaySet = PathwaySet)

  if (length(weighs) == 0) {
    stop("No pathway weights found for the specified PathwaySet")
  }

  # Get all genes across pathways
  all_genes <- unique(unlist(lapply(weighs, names)))

  # Filter counts to include only genes present in pathways
  counts_filtered <- counts[rownames(counts) %in% all_genes, , drop = FALSE]

  if (nrow(counts_filtered) == 0) {
    stop("No genes from the counts matrix are present in the specified pathways")
  }

  message(sprintf("Using %d genes out of %d total genes for pathway scoring",
                  nrow(counts_filtered), nrow(counts)))

  # Calculate offset if not provided
  if (is.null(offset)) {
    offset <- colSums(as.matrix(counts_filtered))
  }

  # Fit the probabilistic model
  message(sprintf("Fitting %s model using %d cores...", method, ncores))
  cmf_mat <- fit_and_cmf(
    counts = counts_filtered,
    method = method,
    ncores = ncores,
    maxit = maxit,
    batch = batch,
    offset = offset
  )

  # Initialize result matrix
  res_mat <- matrix(
    NA,
    nrow = length(weighs),
    ncol = ncol(cmf_mat),
    dimnames = list(names(weighs), colnames(cmf_mat))
  )

  # Calculate pathway scores
  message("Calculating pathway activity scores...")
  for (i in seq_along(weighs)) {
    pathway_name <- names(weighs)[i]
    score_weights <- weighs[[i]]
    genes_in_pathway <- names(score_weights)

    # Find intersection of pathway genes and expressed genes
    genes_use <- intersect(genes_in_pathway, rownames(cmf_mat))

    if (length(genes_use) == 0) {
      warning(sprintf("No genes from pathway '%s' found in expression matrix", pathway_name))
      next
    }

    # Normalize weights for the available genes
    weights_use <- score_weights[genes_use]
    weights_use <- weights_use / sum(weights_use)

    # Extract expression submatrix
    submat <- cmf_mat[genes_use, , drop = FALSE]

    # Calculate pathway score: weighted sum of expressions
    pathway_score <- as.numeric(crossprod(weights_use, submat))

    # Store result
    res_mat[i, ] <- pathway_score
  }

  # Remove pathways with no valid scores
  valid_pathways <- apply(res_mat, 1, function(x) !all(is.na(x)))
  res_mat <- res_mat[valid_pathways, , drop = FALSE]

  if (nrow(res_mat) == 0) {
    stop("No valid pathway scores could be calculated")
  }

  message(sprintf("Successfully calculated scores for %d pathways across %d cells",
                  nrow(res_mat), ncol(res_mat)))

  return(res_mat)
}

