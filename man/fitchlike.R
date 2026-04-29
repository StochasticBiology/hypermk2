#' Fitch-like algorithm to build possible minimum-evolution state sets over a tree
#'
#' @param tree A phylogenetic tree linking observations.
#' @param tip_states A numeric matrix with 0s and 1s for tip observations. Row i should correspond to observation at tip i.
#'
#' @return A list of lists. Each element corresponds to a tip or node in the tree. Each subelement corresponds to a possible state at that tip, along with a matrix describing the indices of "witness" states at the descendant vertices that correspond to the given state at this vertex.
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' build_states(tree, data)
#' @export
build_states <- function(tree, tip_states) {

  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  total_nodes <- Ntip + Nnode
  
  states <- vector("list", total_nodes)
  
  # --- helper: merge identical parent states ---
  merge_states <- function(state_list) {
    keys <- sapply(state_list, function(x) paste0(x$state, collapse = ""))
    ukeys <- unique(keys)
    
    out <- vector("list", length(ukeys))
    
    for (k in seq_along(ukeys)) {
      idx <- which(keys == ukeys[k])
      
      # merge pairs
      pairs <- do.call(rbind, lapply(state_list[idx], function(x) x$pairs))
      pairs <- unique(pairs)
      
      out[[k]] <- list(
        state = state_list[[idx[1]]]$state,
        pairs = pairs
      )
    }
    out
  }
  
  # --- tips ---
  for (i in seq_len(Ntip)) {
    states[[i]] <- list(
      list(state = tip_states[i, ], pairs = NULL)
    )
  }
  
  # --- postorder traversal ---
  tree_post <- ape::reorder.phylo(tree, "postorder")
  internal_nodes <- unique(tree_post$edge[,1])
  internal_nodes <- internal_nodes[internal_nodes > Ntip]
  
  hdist <- function(a, b) sum(a != b)
  
  for (node in internal_nodes) {
    children <- phangorn::Children(tree_post, node)
    left <- children[1]
    right <- children[2]
    
    left_states <- states[[left]]
    right_states <- states[[right]]
    
    candidates <- list()
    idx <- 1
    
    for (i in seq_along(left_states)) {
      for (j in seq_along(right_states)) {
        
        v1 <- left_states[[i]]$state
        v2 <- right_states[[j]]$state
        
        # --- Fitch-style bitwise possibilities ---
        choices <- lapply(seq_along(v1), function(k) {
          if (v1[k] == v2[k]) v1[k] else c(0,1)
        })
        
        combos <- expand.grid(choices)
        
        for (r in seq_len(nrow(combos))) {
          parent_state <- as.numeric(combos[r, ])
          
          candidates[[idx]] <- list(
            state = parent_state,
            pair  = c(i, j)
          )
          idx <- idx + 1
        }
      }
    }
    
    # --- compute costs over pairs ---
    costs <- sapply(candidates, function(x) {
      i_vec <- x$state
      pair <- x$pair
      
      y <- left_states[[pair[1]]]$state
      z <- right_states[[pair[2]]]$state
      
      hdist(i_vec, y) + hdist(i_vec, z)
    })
    
    min_cost <- min(costs)
    best <- candidates[costs == min_cost]
    
    # --- convert to (state, pairs) structure ---
    structured <- lapply(best, function(x) {
      list(
        state = x$state,
        pairs = matrix(x$pair, nrow = 1)
      )
    })
    
    # --- merge duplicate parent states ---
    states[[node]] <- merge_states(structured)
  }
  
  states
}

#' Sample a single minimum-evolution instance from the Fitch-like algorithm output
#'
#' @param tree A phylogenetic tree linking observations.
#' @param states A list of lists describing possible states and witnesses, output from build_states
#'
#' @return A named list containing (a) a list of binary vectors, with each vector corresponding to the sampled state at that vertex in the tree; (b) a dataframe of From-To transitions on the tree (in decimal representation)
#' @examples
#' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rtree(3)
#' state.set = build_states(tree, data)
#' sample_states(tree, state.set)
#' @export
sample_states <- function(tree, states) {
  
  Ntip <- length(tree$tip.label)
  
  # robust root detection
  root <- setdiff(tree$edge[,1], tree$edge[,2])[1]
  
  sampled <- vector("list", length(states))
  
  recurse <- function(node, state_idx) {
    s <- states[[node]][[state_idx]]
    sampled[[node]] <<- s$state
    
    if (node <= Ntip) return()
    
    children <- phangorn::Children(tree, node)
    
    pairs <- s$pairs
    if(nrow(pairs) == 1) {
      pick <- pairs[1,]
    } else {
      pick <- pairs[sample(nrow(pairs), 1), ]
    }  
    recurse(children[1], pick[1])
    recurse(children[2], pick[2])
  }
  
  root_idx <- sample(seq_along(states[[root]]), 1)
  recurse(root, root_idx)
  
  # --- helper: binary vector -> decimal ---
  bin_to_dec <- function(v) {
    sum(v * 2^((length(v)-1):0))
  }
  
  # --- build edge dataframe ---
  edge_df <- data.frame(
    From = integer(nrow(tree$edge)),
    To   = integer(nrow(tree$edge))
  )
  
  for (e in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[e, 1]
    child  <- tree$edge[e, 2]
    
    edge_df$From[e] <- bin_to_dec(sampled[[parent]])
    edge_df$To[e]   <- bin_to_dec(sampled[[child]])
  }
  
  list(
    states = sampled,
    edges  = edge_df
  )
}