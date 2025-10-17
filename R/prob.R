

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
    mu_hat <- tryCatch(as.numeric(predict(fit, type = "response")),
                       error = function(e) rep(NA_real_, ncell))
    theta_hat <- tryCatch(as.numeric(glmmTMB::sigma(fit)),
                          error = function(e) NA_real_)
    pi_hat <- if (want_zi) {
      tryCatch(as.numeric(predict(fit, type = "zprob")[1]),
               error = function(e) NA_real_)
    } else NA
    list(gene = gene, mu = mu_hat, theta = theta_hat, pi = pi_hat)
  }

  # --- Parallel setting（future + doFuture）---
  doFuture::registerDoFuture()
  # Windows/Cross-platform recommendation multisession；Linux/macOS or multicore
  old_plan <- future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(old_plan), add = TRUE)

  # --- preprocessing offset/batch ---
  log_offset <- if (!is.null(offset)) log(as.numeric(offset)) else rep(0, ncell)
  batch_fac  <- if (!is.null(batch)) factor(batch) else factor(rep(1, ncell))

  # --- Define the formula in advance ---
  form_nb  <- if (is.null(batch)) y ~ 1 else y ~ batch
  form_zi  <- ~ 1

  # --- Main loop (keeping the foreach unchanged)---
  results <- foreach::foreach(i = seq_len(n_genes),
                              .packages = "glmmTMB") %dopar% {
                                y <- as.numeric(counts[i, ])

                                if (all(y == 0)) {
                                  return(list(gene = genes[i],
                                              mu = rep(NA_real_, ncell),
                                              theta = NA_real_,
                                              pi = if (method == "zinb") rep(NA_real_, ncell) else NA))
                                }

                                dat <- data.frame(y = y, batch = batch_fac, log_offset = log_offset)

                                fit <- tryCatch({
                                  if (method == "nb") {
                                    glmmTMB::glmmTMB(
                                      form_nb,
                                      family = glmmTMB::nbinom2(link = "log"),
                                      ziformula = ~0,
                                      offset = dat$log_offset,
                                      data = dat,
                                      control = glmmTMB::glmmTMBControl(
                                        optimizer = stats::optim,
                                        optArgs = list(method = "BFGS")
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
                                        optArgs = list(method = "BFGS")
                                      )
                                    )
                                  }
                                }, error = function(e) NULL)

                                print(paste(i/n_genes*100, "%",sep = ""))

                                if (is.null(fit)) {
                                  list(gene = genes[i],
                                       mu = rep(NA_real_, ncell),
                                       theta = NA_real_,
                                       pi = if (method == "zinb") rep(NA_real_, ncell) else NA)
                                } else {
                                  safe_extract(fit, want_zi = (method == "zinb"), gene = genes[i], ncell = ncell)
                                }
                              }
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

  # The probability of "count = 0" is forcibly set to 0.
  prob_mat[counts_mat == 0] <- 0

  return(prob_mat)
}

## Main function: Integrating all steps
fit_and_cmf <- function(counts, method = c("nb", "zinb"), ncores = 4, batch = NULL, offset = NULL, maxit = 30) {
  method <- match.arg(method)

  # Ensure that counts is a matrix
  if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix")) {
    counts <- as.matrix(counts)
  }
  if (!is.matrix(counts)) stop("counts must be a matrix")

  ncell <- ncol(counts)

  # check batch
  if (!is.null(batch)) {
    if (length(batch) != ncell) {
      warning("Length of batch not equal to number of cells, ignoring batch")
      batch <- NULL
    }
  }

  # check offset
  if (!is.null(offset)) {
    if (length(offset) != ncell) {
      warning("Length of offset not equal to number of cells, ignoring offset")
      offset <- NULL
    }
  }

  # Parallel fitting
  fit_list <- fit_gene_counts_parallel(counts, method = method, ncores = ncores,
                                       batch = batch, offset = offset, maxit = 30)

  # extract parameter
  res <- extract_mu_theta(fit_list)

  # Calculate the CMF probability matrix
  cmf_mat <- counts_to_cmf(counts, res$mu_mat, res$theta_vec)

  return(cmf_mat)
}


