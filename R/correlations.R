#' Plot correlation statistics in a dataset
#'
#' @param m A matrix of binary observations.
#'
#' @return A matrix of correlation statistics.
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' feature_correlations(data)
#' @export
feature_correlations <- function(m) {
  X = m
  feature_names = colnames(X)
  n <- nrow(X)
  p <- ncol(X)
  
  res_mat <- matrix(NA, p, p)
  
  for(i in 1:p){
    for(j in 1:p){
      xi <- X[,i]
      xj <- X[,j]
      
      n11 <- sum(xi == 1 & xj == 1)
      n1. <- sum(xi == 1)
      n.1 <- sum(xj == 1)
      
      E11 <- n1. * n.1 / n
      
      res_mat[i,j] <- (n11 - E11) / sqrt(E11 + 1e-12)
    }
  }
  
  diag(res_mat) = NA
  rownames(res_mat) = feature_names
  colnames(res_mat) = feature_names
  
  return(res_mat)
}

#' Plot correlation statistics in a dataset
#'
#' @param m A matrix of binary observations.
#'
#' @return A ggplot2 object.
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' plot_correlations(data)
#' @export
plot_correlations <- function(m) {
  feature_names = colnames(m)
  res_mat = feature_correlations(m)
  df <- reshape2::melt(res_mat)
  colnames(df) <- c("row", "col", "value")
  
  df$row <- factor(df$row, levels = feature_names)
  df$col <- factor(df$col, levels = feature_names)
  
  return(ggplot2::ggplot(df, ggplot2::aes(col, row, fill = value)) +
           ggplot2::geom_tile() +
           ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red") +
           ggplot2::theme_minimal() +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
           panel.grid = ggplot2::element_blank()) + 
           ggplot2::labs(x = NULL, y=NULL, fill = "Residual")
  )
}

#' Phylogenetic correlation statistics in a dataset
#'
#' @param m A matrix of binary observations. Row names should correspond to tree tips.
#' @param tree A phylogenetic tree.
#' @param pre.prune Boolean (default TRUE) whether to preprocess by removing non-significant features from uncorrected correlation first
#'
#' @return A matrix of likelihood-ratio-test statistics for dependence vs independence.
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' colnames(data) = c("a", "b", "c")
#' rownames(data) = tree$tip.label
#' phylo_assoc_matrix(data, tree)
#' @export
phylo_assoc_matrix <- function(m, 
                               tree,
                               pre.prune = TRUE) {
  
  if(pre.prune == TRUE) {
    cm = feature_correlations(m)
    cdf <- reshape2::melt(cm)
    sigdf = cdf[!is.na(cdf$value) & cdf$value > 2,]
    if(nrow(sigdf) == 0) {
      return(NULL)
    }
    sigfeats = unique(sigdf$Var1)
    sigindices = which(colnames(m) %in% sigfeats)
    sigm = m[,sigindices]
    X = sigm
  } else {
    X = m
  }
  
  p <- ncol(X)
  assoc <- matrix(NA, p, p)
  rownames(assoc) <- colnames(X)
  colnames(assoc) <- colnames(X)
  
  # convert to corHMM format once
  tip_names <- rownames(X)
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      
      dat <- data.frame(
        species = tip_names,
        trait1 = X[, i],
        trait2 = X[, j]
      )
      
      # ensure alignment
      dat <- dat[match(tree$tip.label, dat$species), ]
      
      # build joint states
      dat$state <- paste0(dat$trait1, dat$trait2)
      
      dat_use <- data.frame(species = dat$species, state = dat$state)
      
      # drop missing if any
      dat_use <- na.omit(dat_use)
      
      # independent model
      fit_ind <- try(
        corHMM::corHMM(tree, dat_use, rate.cat = 1, model = "ER"),
        silent = TRUE
      )
      
      # dependent model
      fit_dep <- try(
        corHMM::corHMM(tree, dat_use, rate.cat = 1, model = "ARD"),
        silent = TRUE
      )
      
      if (inherits(fit_ind, "try-error") || inherits(fit_dep, "try-error")) {
        next
      }
      
      lrt <- 2 * (fit_dep$loglik - fit_ind$loglik)
      
      assoc[i, j] <- lrt
      assoc[j, i] <- lrt
    }
  }
  
  diag(assoc) <- 0
  return(assoc)
}
