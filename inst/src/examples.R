expt = 2
if(expt == 1 | expt == 2) {
  # 50 observations of L=10 random samples -- kind of OK? 81 states
  # 20 observations of L=22 -- 39 states, fitted in a minute
  L = 5
  n = 20
  tree = ape::rtree(n)
  if(expt == 1) {
    m = matrix(rbinom(n*L, 1, 0.5), ncol=L, nrow=n)
  } else {
    m = matrix(rep(c(1,0,0,0,0,
                     1,1,0,0,0,
                     1,1,1,0,0,
                     1,1,1,1,0,
                     1,1,1,1,1), n/L), ncol=L, nrow=n, byrow = TRUE)
  }
}
# pset will store lists of possible ancestral states
if(expt == 3) {
  load(system.file("data", "tb-test-tree.RData", package = "hypermk2"))
  load(system.file("data", "tb-test-df.RData", package = "hypermk2"))
  m = matrix(ncol = ncol(this.df)-1)
  n = length(tree$tip.label)
  for(i in 1:n) {
    r = which(this.df[,1]==tree$tip.label[i])
    m = rbind(m, as.matrix(this.df[r,2:ncol(this.df)]))
  }
  m = m[2:nrow(m),]
}
if(expt == 4) {
  load(system.file("data", "kp-test-tree.RData", package = "hypermk2"))
  load(system.file("data", "kp-test-data.RData", package = "hypermk2"))
}

hypermk2(m, tree)
