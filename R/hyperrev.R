# get the set of binary vectors where every disagreed bit between x and y is expanded into both 0 and 1
expand_disagreements <- function(x, y) {
  stopifnot(length(x) == length(y))

  # For each position, keep either the single value (if equal)
  # or both values (if different)
  choices <- Map(function(a, b) {
    if (a == b) a else c(a, b)
  }, x, y)

  # Generate all combinations
  combos <- expand.grid(choices)

  # Convert rows to vectors
  result <- split(as.matrix(combos), seq(nrow(combos)))
  lapply(result, as.numeric)
}

closest_by_hamming <- function(lst, target) {
  stopifnot(all(sapply(lst, length) == length(target)))

  # Compute Hamming distances
  dists <- sapply(lst, function(v) sum(v != target))

  # Return sublist with minimum distance
  lst[dists == min(dists)]
}

# given two lists of binary vectors, return a list of binary vectors
# that minimise the sum of hamming distances to an element of the two
# lists. in other words the set of i for which argmax_j(i,list1[j])
# + argmax_k(i, list2[k]) takes a minimum value. also return the
# subsets of the two lists containing the js and ks
min_sum_hamming_three_lists <- function(lst1, lst2) {
  n <- length(lst1[[1]])

  # matrices
  m1 <- do.call(rbind, lst1)
  m2 <- do.call(rbind, lst2)

  # --- find optimal i's (bitwise) ---
  choices <- lapply(seq_len(n), function(i) {
    vals1 <- unique(m1[, i])
    vals2 <- unique(m2[, i])

    cost <- function(bit) {
      min(abs(vals1 - bit)) + min(abs(vals2 - bit))
    }

    costs <- sapply(c(0,1), cost)
    c(0,1)[costs == min(costs)]
  })

  combos <- expand.grid(choices)
  candidates <- split(as.matrix(combos), seq(nrow(combos))) |> lapply(as.numeric)

  # helper
  hdist <- function(a, b) sum(a != b)

  # --- collect witnesses ---
  witness1 <- list()
  witness2 <- list()

  for (i in candidates) {
    d1 <- apply(m1, 1, hdist, b = i)
    d2 <- apply(m2, 1, hdist, b = i)

    witness1 <- c(witness1, lst1[d1 == min(d1)])
    witness2 <- c(witness2, lst2[d2 == min(d2)])
  }

  # remove duplicates (by string encoding)
  dedup <- function(lst) {
    keys <- sapply(lst, paste0, collapse = "")
    lst[!duplicated(keys)]
  }

  list(
    i_list = dedup(candidates),
    lst1_witness = dedup(witness1),
    lst2_witness = dedup(witness2)
  )
}

# given two lists of binary vectors, return a list of binary vectors
# that minimise the sum of irreversible 0->1 changes to an element of the two
# lists. in other words the set of i for which argmax_j(i,list1[j])
# + argmax_k(i, list2[k]) takes a minimum value. also return the
# subsets of the two lists containing the js and ks
min_irreversible_ancestor <- function(lst1, lst2) {
  # helper to deduplicate list of vectors
  dedup_vecs <- function(lst) {
    keys <- sapply(lst, paste0, collapse = "")
    lst[!duplicated(keys)]
  }

  results <- list()
  scores <- c()
  idx <- 1

  for (i in seq_along(lst1)) {
    for (j in seq_along(lst2)) {
      v1 <- lst1[[i]]
      v2 <- lst2[[j]]

      anc <- pmin(v1, v2)   # bitwise AND
      score <- sum(anc)     # maximize number of 1's

      results[[idx]] <- list(
        i = anc,
        j = v1,
        k = v2,
        score = score
      )
      scores[idx] <- score
      idx <- idx + 1
    }
  }

  # keep optimal
  max_score <- max(scores)
  best <- results[scores == max_score]

  # collect outputs
  i_list <- lapply(best, `[[`, "i")
  w1 <- lapply(best, `[[`, "j")
  w2 <- lapply(best, `[[`, "k")

  list(
    i_list = dedup_vecs(i_list),
    lst1_witness = dedup_vecs(w1),
    lst2_witness = dedup_vecs(w2),
    score = max_score
  )
}

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
                    force.origin = FALSE) {
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

  if(reversible == FALSE) {
    mstr = apply(m, 1, paste0, collapse = "")
    arb = hyperdags::simplest_arborescence(mstr)
    g = arb$rewired.graph
    edges = igraph::as_data_frame(g, what = "edges")
    colnames(edges) = c("From", "To")
    trans = as.data.frame(apply(edges, c(1,2), binS_to_dec))
  }
  if(reversible == TRUE) {
  pset = list()
  for(i in 1:(n+tree$Nnode)) {
    if(i <= n) {
      pset[[i]] = list(m[i,])
    } else {
      pset[[i]] = list()
    }
  }

  # iterate over tree
  change = TRUE
  while(change == TRUE) {
    change = FALSE
    # go over internal nodes
    for(i in 1:tree$Nnode) {
      noderef = i+n
      # if we don't have guesses for this node yet
      if(length(pset[[noderef]]) == 0) {
        change = TRUE
        d = phangorn::Children(tree, noderef)
        if(d[1] > d[2]) { tmp = d[2]; d[2] = d[1]; d[1] = tmp; }
        # if both daughters are tips, they both have single states
        if(d[1] <= n & d[2] <= n) {
          x1 = pset[[d[1]]][[1]]
          x2 = pset[[d[2]]][[1]]
          # this node gets the set of disagreements
          pset[[noderef]] = expand_disagreements(x1, x2)
          if(verbose == TRUE) {
            cat("Node ", noderef, " has tip daughters ", d[1], " ", d[2], " so gets\n ")
            print(pset[[d[1]]])
            print("(+)")
            print(pset[[d[2]]])
            print("(=)")
            print(pset[[noderef]])
          }
        } else if(d[1] <= n & d[2] > n) {
          # one daughter a node, the other a tip
          # this has ceased to matter for the algorithm
          # it's just like two node daughters
          if(length(pset[[d[2]]]) > 0) {
            tmp = min_sum_hamming_three_lists(pset[[d[1]]], pset[[d[2]]])
              pset[[d[1]]] = tmp$lst1_witness
              pset[[d[2]]] = tmp$lst2_witness
              pset[[noderef]] = tmp$i_list

            if(verbose == TRUE) {
              cat("Node ", noderef, " has tip daughter ", d[1], " and node daughter ", d[2], " so gets\n")
              print(pset[[d[1]]])
              print("(+)")
              print(pset[[d[2]]])
              print("(=)")
              print(pset[[noderef]])
            }
          }
        } else if(d[1] > n & d[2] > n) {
          # two node daughters
          if(length(pset[[d[1]]]) > 0 & length(pset[[d[2]]]) > 0) {
            # this node gets the lowest hamming distance set
            # the daughter sets get reduced to the witnesses for this
            # minimum set
                tmp = min_sum_hamming_three_lists(pset[[d[1]]], pset[[d[2]]])
              pset[[d[1]]] = tmp$lst1_witness
              pset[[d[2]]] = tmp$lst2_witness
              pset[[noderef]] = tmp$i_list
            if(verbose == TRUE) {
              cat("Node ", noderef, " has node daughters ", d[1], " ", d[2], " so gets\n")
              print(pset[[d[1]]])
              print("(+)")
              print(pset[[d[2]]])
              print("(=)")
              print(pset[[noderef]])
            }
          }
        }
      }
    }
  }

  # if there are any sets remaining, pull one element
  for(i in 1:length(pset)) {
    pset[[i]] = pset[[i]][[1]]
  }


  # now produce a set of transitions
  trans = data.frame()
  for(i in 1:tree$Nnode) {
    noderef = i+n
    d = phangorn::Children(tree, noderef)
    # transition set
    trans = rbind(trans, data.frame(From=bin_to_dec(pset[[noderef]]),
                                    To=bin_to_dec(pset[[d[[1]]]])))
    trans = rbind(trans, data.frame(From=bin_to_dec(pset[[noderef]]),
                                    To=bin_to_dec(pset[[d[[2]]]])))
    if(force.origin == TRUE) {
      trans = rbind(trans, data.frame(From=0,
                                      To=bin_to_dec(pset[[noderef]])))
    }
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
  if(FALSE) {
    fit.model = castor::fit_mk(trees = trees, Nstates = Nstates,
                               tip_states = stateobs, rate_model = "ARD")
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

  hyperfit = list(mk2_fluxes = r.df,
                  fitted_mk = fit.model,
                  trans = trans.df,
                  L = L,
                  force.origin = force.origin,
                  feature.names = colnames(m))

  return(hyperfit)
}

