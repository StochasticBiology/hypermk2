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


#' Build a cheap approximation to reduced state space, avoiding Fitch-like computation
#'
#' @param m A matrix of binary observations. Each row should correspond to the ith tree tip observation.
#' @param tree A phylogenetic tree linking observations.
#' @param force.origin Boolean (default FALSE), whether to force the root of the tree to have state 0^L
#'
#' @return A dataframe with transitions between (decimal) states in the reduced space
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' cheap_transition_set(data, tree)
#' @export
cheap_transition_set = function(m,
                                tree,
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
  
  return(trans)
}
