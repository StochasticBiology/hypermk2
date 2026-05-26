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
