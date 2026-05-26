#' Get correlation statistics of pairs vs independent features
#'
#' @param m A matrix of binary observations. Each row should correspond to the ith tree tip observation.
#' @param tree A phylogenetic tree linking observations.
#'
#' @return A named list containing dAICs and AICs for independent and coupled model fits
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' phylo_correlations(data, tree)
#' @export
phylo_correlations = function(m, tree) {
  L = ncol(m)
  pair_mat = ind_mat = matrix(NA, ncol=L, nrow=L)
  single_vec = rep(NA, L)
  for(i in 1:L) {
    m.0 = as.matrix(m[,i], col=1)
    rownames(m.0) = NULL
    colnames(m.0) = NULL
    fit = hypermk::mk_infer_phylogenetic(m.0, tree, flux.samples = 1)
    single_vec[i] = fit$fitted_mk$AIC
  }
  for(i in 1:(L-1)) {
    for(j in (i+1):L) {
      m.0 = m[,c(i,j)]
      rownames(m.0) = NULL
      colnames(m.0) = NULL
      fit = hypermk::mk_infer_phylogenetic(m.0, tree, flux.samples = 1)
      pair_mat[i,j] = pair_mat[j,i] = fit$fitted_mk$AIC
      ind_mat[i,j] = ind_mat[j,i] = single_vec[i]+single_vec[j]
    }
  }
  dAICs = pair_mat-ind_mat
  colnames(dAICs) = rownames(dAICs) = colnames(m)
  return(list(dAICs=dAICs,
              pair_mat=pair_mat,
              ind_mat=ind_mat))
}