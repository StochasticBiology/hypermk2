#' Get statistics of a null model fit of independent features
#'
#' @param m A matrix of binary observations. Each row should correspond to the ith tree tip observation.
#' @param tree A phylogenetic tree linking observations.
#' @param ... other parameters to pass to hypermk2
#'
#' @return A named list containing total and feature-specific likelihoods and AICs
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' hypermk2_independent(data, tree)
#' @export
hypermk2_independent = function(m,
                                tree,
                                ...) {
  res.df = data.frame()
  for(i in 1:ncol(m)) {
    this.f = matrix(m[,i], ncol=1)
    if(sum(this.f) != 0 & sum(this.f) != nrow(m)) {
      this.fit = hypermk2(this.f, tree, nwalker = 1, ...)
      res.df = rbind(res.df, data.frame(feature=i, loglik = this.fit$fitted_mk$loglikelihood, AIC = this.fit$fitted_mk$AIC))
    }
  }
  #  res.df = rbind(res.df, data.frame(feature=0, loglik=sum(res.df$loglik), AIC=sum(res.df$AIC)))
  return(list(loglik =sum(res.df$loglik),
              AIC = sum(res.df$AIC),
              by.feature = res.df))
}

#' Use HyperMk2 to fit a reversible evolutionary accumulation model
#'
#' @param m A matrix of binary observations. Each row should correspond to the ith tree tip observation.
#' @param tree A phylogenetic tree linking observations.
#' @param reversible Boolean (default TRUE) whether to allow reversible transitions
#' @param nwalker Integer (default 10000), the number of random walkers to simulate on the inferred transition network to sample fluxes
#' @param force.origin Boolean (default FALSE), whether to force the root of the tree to have state 0^L
#' @param compare.null Boolean (default FALSE), whether to compare a null model of independent characters
#'
#' @return A named list containing the fitted Mk model object, inferred fluxes between states, the number of features, set of transitions in the reduced space, and feature names
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' hypermk2(data, tree)
#' @export
hypermk2 = function(m,
                    tree,
                    reversible = TRUE,
                    nwalker=10000,
                    force.origin = FALSE,
                    compare.null = FALSE) {
  verbose = FALSE
  n = length(tree$tip.label)
  L = ncol(m)
  
  bin_to_dec <- function(x) {
    sum(x * 2^((length(x)-1):0))
  }
  binS_to_dec <- function(x) {
    x = as.numeric(strsplit(x, "")[[1]])
    sum(x * 2^((length(x)-1):0))
  }
  message("Building reduced state space...")
  if(reversible == FALSE) {
    mstr = apply(m, 1, paste0, collapse = "")
    arb = hyperdags::simplest_arborescence(mstr)
    g = arb$rewired.graph
    edges = igraph::as_data_frame(g, what = "edges")
    colnames(edges) = c("From", "To")
    trans = as.data.frame(apply(edges, c(1,2), binS_to_dec))
  }
  if(reversible == TRUE) {
    state.set = build_states(tree, m)
    sample.states = sample_states(tree, state.set)
    trans = sample.states$edges
    if(force.origin == TRUE) {
      trans = rbind(trans, data.frame(From=0,
                                      To=unique(sample.states$edges$From)))
    }
    
    trans = trans[trans$From != trans$To,]
  }
  
  trans.df = trans
  # try to pull this together into inference for the Mk model
  # relabel states in the transition set
  stateset = unique(c(trans$From, trans$To))
  statetrans = data.frame()
  for(i in 1:nrow(trans)) {
    this.from = which(stateset == trans$From[i])
    this.to = which(stateset == trans$To[i])
    statetrans = rbind(statetrans, data.frame(From=this.from, To=this.to))
  }
  statetrans = unique(statetrans)
  zero.state = which(stateset == 0)
  
  # relabel observations on the original dataset
  stateobs = c()
  mdec = as.vector(apply(m, 1, bin_to_dec))
  for(i in 1:nrow(m)) {
    this.obs = which(stateset == bin_to_dec(m[i,]))
    stateobs = c(stateobs, this.obs)
  }
  
  # use this relabelling to motivate the indexed Q-matrix
  Q = matrix(0, nrow=length(stateset), ncol=length(stateset))
  for(i in 1:nrow(statetrans)) {
    Q[statetrans$From[i],statetrans$To[i]] = i
  }
  
  
  trees = tree
  Nstates = length(stateset)
  tip_states = stateobs
  rate_model = Q
  message(Nstates, " states, ", length(which(as.vector(Q) != 0)), " transitions ", length(stateobs), " tips")
  message("Fitting Mk2 model...")
  if(force.origin == TRUE) {
    root_prior = rep(0, Nstates)
    root_prior[zero.state] = 1
    fit.model = castor::fit_mk(trees = trees, Nstates = Nstates,
                               root_prior = root_prior,
                               tip_states = stateobs, rate_model = Q)
  } else {
    fit.model = castor::fit_mk(trees = trees, Nstates = Nstates,
                               tip_states = stateobs, rate_model = Q)
  }
  
  # get the set of supported transitions
  model.df = hypermk::mk_pull_transitions(fit.model)
  model.df = model.df[model.df$From != model.df$To,]
  model.df$From = model.df$From+1
  model.df$To = model.df$To+1
  
  message("simulating walkers")
  # simulate walkers
  mk.rev.df = hypermk::mk_pull_transitions(fit.model, reversible = TRUE)
  mk.rev.df = mk.rev.df[mk.rev.df$From != mk.rev.df$To, ]
  mk.rev.df$Rate[mk.rev.df$Rate == Inf] = 10 * max(mk.rev.df$Rate[mk.rev.df$Rate != Inf])
  mk.rev.df$From = mk.rev.df$From +1
  mk.rev.df$To = mk.rev.df$To + 1
  mk.rev.df$Flux = 0
  
  DecToBinS <- function(x, len) {
    s = c()
    for(j in (len-1):0)
    {
      if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
    }
    return(paste(s, collapse=""))
  }
  
  statesbin = sapply(stateset, DecToBinS, len=L)
  onecounts = as.vector(sapply(statesbin, stringr::str_count, "1"))
  for (walk in 1:nwalker) {
    if(force.origin == TRUE) {
      state = zero.state
    } else {
      state = mk.rev.df$From[sample(1:nrow(mk.rev.df), 1)]
    }
    for (t in 1:(2 * L)) {
      trans = fit.model$transition_matrix[state, ]
      trans[state] = 0
      if (sum(trans) == 0) {
        break
      }
      this.out = sample(1:Nstates, size = 1, prob = trans)
      ref = which(mk.rev.df$From == state & mk.rev.df$To == this.out)
      mk.rev.df$Flux[ref] = mk.rev.df$Flux[ref] + 1
      state = this.out
    }
  }
  mk.rev.df
  
  r.df = mk.rev.df
  r.df$From = stateset[r.df$From]
  r.df$To = stateset[r.df$To]
  
  r.df$FromS = sapply(r.df$From, DecToBinS, len=L)
  r.df$ToS = sapply(r.df$To, DecToBinS, len=L)
  
  if(compare.null == FALSE) {
    hyperfit = list(mk2_fluxes = r.df,
                    fitted_mk = fit.model,
                    trans = trans.df,
                    L = L,
                    force.origin = force.origin,
                    feature.names = colnames(m))
  } else {
    message("Fitting null model...")
    null.fit = hypermk2_independent(m, tree, reversible=reversible, force.origin=force.origin)
    hyperfit = list(mk2_fluxes = r.df,
                    fitted_mk = fit.model,
                    trans = trans.df,
                    L = L,
                    force.origin = force.origin,
                    feature.names = colnames(m),
                    null.fit = null.fit)
  }
  
  return(hyperfit)
}

