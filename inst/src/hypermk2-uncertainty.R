tree = ape::stree(4, type="balanced")
data = matrix(c(1,0,0,
                NA,NA,NA,
                0,1,0,
                0,0,1), byrow = TRUE, ncol=3, nrow=4)
colnames(data) = c("a", "b", "c")
rownames(data) = tree$tip.label
hyperinf::plot_hyperinf_data(data, tree)
states = build_states(data, tree)
states = cheap_transition_set(data, tree)
states
set.seed(1)
sample_states(tree, states, expand.uncertainty = TRUE)

fit = hypermk2(data, tree)

##### linear example

set.seed(1)
sim.dyn = hyperdags::simulate_accumulation(32, 4, dynamics="linear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
colnames(m) = 1:4
fit = hypermk2(m, tree)

hyperinf::plot_hyperinf(fit)

fits = fits2 = list()
for(i in 1:5) {
  mnew = m
  if(i > 1) {
  mnew[sample(length(mnew), (i-1)*0.1*length(mnew))] = NA
  }
  fits[[i]] = hypermk2(mnew, tree)
  fits2[[i]] = hypermk2(mnew, tree, expand.uncertainty=FALSE)
}
hyperinf::plot_hyperinf_comparative(fits, style="full")
hyperinf::plot_hyperinf_comparative(fits2, style="full")

hyperinf::plot_hyperinf_ordering_matrices(fits, type="relative")

####### bilinear example

set.seed(5)
sim.dyn = hyperdags::simulate_accumulation(50, 6, dynamics="bilinear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
colnames(m) = 1:6
fit = hypermk2(m, tree)

hyperinf::plot_hyperinf(fit)
hyperinf::plot_hyperinf_data(m, tree)

fits = fits2 = list()
for(i in 1:5) {
  mnew = m
  if(i > 1) {
    mnew[sample(length(mnew), (i-1)*0.1*length(mnew))] = NA
  }
  fits[[i]] = hypermk2(mnew, tree, expand.uncertainty=FALSE)
  fits2[[i]] = hypermk2(mnew, tree, expand.uncertainty=FALSE)
}
hyperinf::plot_hyperinf_data(mnew, tree)
hyperinf::plot_hyperinf_comparative(fits, style="full")
hyperinf::plot_hyperinf_comparative(fits2, style="full")

hyperinf::plot_hyperinf_ordering_matrices(fits, type="relative")
