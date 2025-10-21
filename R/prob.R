

fit_gene_counts_parallel <- function(counts,
                                     method = c("nb", "zinb"),
                                     ncores = 4,
                                     batch = NULL,
                                     offset = NULL,
                                     maxit = 30) {
  method <- match.arg(method)
  stopifnot(all(counts >= 0))
  print("Transformation Begins")
  n_genes <- nrow(counts)
  genes   <- rownames(counts)
  ncell   <- ncol(counts)
  
  # --- External Functions ---
  safe_extract <- function(fit, want_zi, gene, ncell) {
    # --- mu_hat ---
    mu_hat <- tryCatch(as.numeric(predict(fit, type = "response")),
                       error = function(e) NULL)
    
    # 检测 mu_hat 是否异常
    if (is.null(mu_hat) || all(is.na(mu_hat)) || all(is.infinite(mu_hat))) {
      warning(sprintf("%s: invalid fitted values (NULL, NA, or Inf) — removed from results.", gene))
      return(NULL)
    }
    
    # --- theta_hat ---
    theta_hat <- tryCatch(as.numeric(glmmTMB::sigma(fit)),
                          error = function(e) NA_real_)
    
    # --- pi_hat (only for ZINB) ---
    pi_hat <- if (want_zi) {
      tryCatch(as.numeric(predict(fit, type = "zprob")[1]),
               error = function(e) NA_real_)
    } else NA
    
    # --- Final result ---
    list(gene = gene, mu = mu_hat, theta = theta_hat, pi = pi_hat)
  }
  
  
  # --- Parallel setting（future + doFuture）---
  doFuture::registerDoFuture()
  old_plan <- future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(old_plan), add = TRUE)
  
  # --- preprocessing offset/batch ---
  log_offset <- if (!is.null(offset)) log(as.numeric(offset)) else rep(0, ncell)
  batch_fac  <- if (!is.null(batch)) factor(batch) else factor(rep(1, ncell))
  
  # --- Define the formula in advance ---
  form_nb  <- if (is.null(batch)) y ~ 1 else y ~ batch
  form_zi  <- ~ 1
  
  # --- Main loop ---
  results <- foreach::foreach(i = seq_len(n_genes),
                              .packages = "glmmTMB") %dopar% {
                                y <- as.numeric(counts[i, ])
                                gene_name <- genes[i]
                                
                                # Skip all-zero genes
                                if (all(y == 0)) {
                                  warning(sprintf("%s gene expression was removed: all zero counts.", gene_name))
                                  return(NULL)
                                }
                                
                                # Check sparsity: too few nonzero entries
                                nonzero_n <- sum(y > 0)
                                if (nonzero_n < 3) {
                                  warning(sprintf("%s gene expression was removed: too sparse for reliable regression convergence.", gene_name))
                                  return(NULL)
                                }
                                
                                dat <- data.frame(y = y, batch = batch_fac, log_offset = log_offset)
                                
                                fit <- tryCatch({
                                  fit_obj <- suppressWarnings(   # 捕获但不打断警告
                                    if (method == "nb") {
                                      glmmTMB::glmmTMB(
                                        form_nb,
                                        family = glmmTMB::nbinom2(link = "log"),
                                        ziformula = ~0,
                                        offset = dat$log_offset,
                                        data = dat,
                                        control = glmmTMB::glmmTMBControl(
                                          optimizer = stats::optim,
                                          optArgs = list(method = "BFGS", maxit = maxit)
                                        )
                                      
                                      )
                                    } else {
                                      glmmTMB::glmmTMB(
                                        form_nb,
                                        ziformula = form_zi,
                                        family = glmmTMB::nbinom2(link = "log"),
                                        offset = dat$log_offset,
                                        data = dat,
                                        control = glmmTMB::glmmTMBControl(
                                          optimizer = stats::optim,
                                          optArgs = list(method = "BFGS", maxit = maxit)
                                        )
                                      )
                                    }
                                  )
                                  
                                  # 成功返回
                                  fit_obj
                                  
                                }, error = function(e) {
                                  warning(sprintf("%s: model fitting failed (%s)", gene_name, conditionMessage(e)))
                                  return(FALSE)
                                })
                                
                                
                                #print(paste0(round(i / n_genes * 100, 2), "%"))
                                
                                if (isFALSE(fit) || is.null(fit)) {
                                  warning(sprintf("%s gene expression was removed: regression failed to converge.", gene_name))
                                  return(NULL)
                                } else {
                                  safe_extract(fit, want_zi = (method == "zinb"), gene = gene_name, ncell = ncell)
                                }
                                
                              }
  results <- results[!vapply(results, function(x) {
    is.null(x) || (is.list(x) && all(vapply(x, is.null, logical(1))))
  }, logical(1))]
  print("All genes transformed")
  results
}


##  Sub-function: Extract mu and theta from the fitting results
extract_mu_theta <- function(fit_list) {
  mu_mat <- do.call(rbind, lapply(fit_list, function(x) x$mu))
  rownames(mu_mat) <- sapply(fit_list, function(x) x$gene)

  theta_vec <- sapply(fit_list, function(x) x$theta)
  names(theta_vec) <- sapply(fit_list, function(x) x$gene)

  return(list(mu_mat = mu_mat, theta_vec = theta_vec))
}

## Sub-function: Calculate the CMF matrix
counts_to_cmf <- function(counts_mat, mu_mat, theta_vec) {
  if (!all(dim(counts_mat) == dim(mu_mat))) {
    stop("counts_mat and mu_mat must have the same dimensions")
  }
  if (length(theta_vec) != nrow(counts_mat)) {
    stop("theta_vec length must equal number of rows in counts_mat")
  }

  theta_mat <- matrix(theta_vec, nrow = nrow(counts_mat), ncol = ncol(counts_mat))

  # q_adj: Non-zero "count - 1", zero remains "0"
  q_adj <- pmax(counts_mat - 1, 0)
  
  
  
  prob_mat <- pnbinom(q = q_adj, mu = mu_mat, size = theta_mat)
  
  fix_nan_prob <- function(prob_mat) {
    # --- Get row names ---
    gene_names <- rownames(prob_mat)

    # --- Detect NaN counts per row ---
    nan_counts <- apply(prob_mat, 1, function(x) sum(is.nan(x)))
    bad_rows <- which(nan_counts > 0)
    
    # --- Report and fix ---
    if (length(bad_rows) > 0) {
      message("!  The following genes contain NaN values (likely due to very low expression in some batches):")
      for (i in bad_rows) {
        message(sprintf("  %s: %d NaN values detected — partial batch expression too low, probability set to 0.",
                        gene_names[i], nan_counts[i]))
      }
      # Replace NaN with 0
      prob_mat[is.nan(prob_mat)] <- 0
    }
    
    # --- Return cleaned matrix ---
    return(prob_mat)
  }
  
  
  prob_mat <- fix_nan_prob(prob_mat)
  
  # The probability of "count = 0" is forcibly set to 0.
  prob_mat[counts_mat == 0] <- 0

  return(prob_mat)
}

# Main function
fit_and_cmf <- function(counts, method = c("nb", "zinb"), ncores, batch, offset, maxit) {
  method <- match.arg(method)
  
  # --- Ensure that counts is a matrix ---
  if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix")) {
    counts <- as.matrix(counts)
  }
  if (!is.matrix(counts)) stop("counts must be a matrix")
  
  ncell <- ncol(counts)
  
  # --- check batch ---
  if (!is.null(batch)) {
    if (length(batch) != ncell) {
      warning("Length of batch not equal to number of cells, ignoring batch")
      batch <- NULL
    }
  }
  
  # --- check offset ---
  if (!is.null(offset)) {
    if (length(offset) != ncell) {
      warning("Length of offset not equal to number of cells, ignoring offset")
      offset <- NULL
    }
  }
  
  # --- Parallel fitting ---
  fit_list <- fit_gene_counts_parallel(counts, method = method, ncores = ncores,
                                       batch = batch, offset = offset, maxit = maxit)
  
  # --- Check if any fits succeeded ---
  if (length(fit_list) == 0) {
    stop("No genes were successfully fitted. Please check input data or parameters.")
  }
  
  # --- Get successfully fitted gene names ---
  fitted_genes <- vapply(fit_list, function(x) x$gene, character(1))
  fitted_genes <- fitted_genes[!is.na(fitted_genes) & fitted_genes != ""]
  
  # --- Filter counts to keep only those genes ---
  missing_genes <- setdiff(rownames(counts), fitted_genes)
  if (length(missing_genes) > 0) {
    message("The following genes were removed because fitting failed: ",
            paste(missing_genes, collapse = ", "))
  }
  
  counts <- counts[fitted_genes, , drop = FALSE]

  # --- Extract parameters ---
  res <- extract_mu_theta(fit_list)
  
  # --- Calculate CMF probability matrix ---
  cmf_mat <- counts_to_cmf(counts, res$mu_mat, res$theta_vec)
  message(sprintf("A total of %d genes were used for scoring.", nrow(cmf_mat)))
  
  return(cmf_mat)
}
