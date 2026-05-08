library(hyperdags)
library(hyperinf)
library(ggraph)
library(ggpubr)

# simply returns a binary (character string) of length len from a decimal
DecToBinS <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

set.seed(5)
sim.dyn = simulate_accumulation(50, 6, dynamics="bilinear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
plot_hyperinf_data(m, tree)
fit = hyperinf(m, tree, method="hypermk2")
fit.r = hyperinf(m, tree, method="hypermk2", reversible = TRUE)
fit.0 = hyperinf(m, tree, method="hypermk2", reversible = TRUE, compare.null = TRUE)

my.labels = function(x) {
  return(sapply(as.numeric(x), DecToBinS, len=6))
}

ggarrange(plot_hyperinf_data(m, tree),
                  ggarrange(
                    plot_hyperinf(fit) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
                  plot_hyperinf(fit.r) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
                  nrow = 2, labels=c("B", "C")),
          nrow = 1, labels=c("A", ""))

ggpubr::ggarrange(plot_hyperinf(fit) + ggraph::geom_node_text(ggplot2::aes(label=name)),
                  plot_hyperinf(fit.r) + ggraph::geom_node_text(ggplot2::aes(label=name)) )

hyperinf_AIC(fit)   
hyperinf_AIC(fit.r) 
fit.0$null.fit$AIC

sim.dyn.lin = simulate_accumulation(50, 6, dynamics="linear")
tree.lin = sim.dyn.lin$my.tree
m.lin <- do.call(rbind, sim.dyn.lin$x)[1:length(sim.dyn.lin$my.tree$tip.label),]
rownames(m.lin) = tree.lin$tip.label
plot_hyperinf_data(m.lin, tree.lin)
fit.lin = hyperinf(m.lin, tree.lin, method="hypermk2")
fit.lin.r = hyperinf(m.lin, tree.lin, method="hypermk2", reversible = TRUE)
fit.lin.0 = hyperinf(m.lin, tree.lin, method="hypermk2", reversible = TRUE, compare.null = TRUE)


ggpubr::ggarrange(plot_hyperinf(fit.lin) + ggraph::geom_node_text(ggplot2::aes(label=name)),
                  plot_hyperinf(fit.lin.r) + ggraph::geom_node_text(ggplot2::aes(label=name)) )

hyperinf_AIC(fit.lin)   
hyperinf_AIC(fit.lin.r) 



set.seed(5)
sim.dyn = simulate_accumulation(20, 4, dynamics="linear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
plot_hyperinf_data(m, tree)
fit.mk = hyperinf(m, tree, method="hypermk", reversible = TRUE)
fit.mk2 = hyperinf(m, tree, method="hypermk2", reversible = TRUE)
fit.mk2.0 = hyperinf(m, tree, method="hypermk2", reversible = TRUE, compare.null= TRUE)
fit.mk.ir = hyperinf(m, tree, method="hypermk", reversible = FALSE)
fit.mk2.ir = hyperinf(m, tree, method="hypermk2", reversible = FALSE)

my.labels = function(x) {
  return(sapply(as.numeric(x), DecToBinS, len=4))
}

ggarrange(plot_hyperinf_data(m, tree),
ggarrange(
  plot_hyperinf(fit.mk) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
  plot_hyperinf(fit.mk2) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
  plot_hyperinf(fit.mk.ir) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
  plot_hyperinf(fit.mk2.ir) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
  nrow=2, ncol=2, labels=c("Bi", "Ci", "ii", "ii")
), nrow=1, widths=c(1,1.5), labels=c("A", "") )

lapply(list(fit.mk, fit.mk2, fit.mk.ir, fit.mk2.ir), hyperinf_AIC)
fit.mk2.0$null.fit$AIC


ggarrange(plot_hyperinf_data(m, tree),
          ggarrange(
            plot_hyperinf(fit.mk2) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
            plot_hyperinf(fit.mk) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
            nrow = 2, labels=c("B", "C")),
          nrow = 1, labels=c("A", ""))

hyperinf_AIC(fit.mk)
hyperinf_AIC(fit.mk2)
plot_hyperinf(fit.mk)
plot_hyperinf(fit.mk2)
