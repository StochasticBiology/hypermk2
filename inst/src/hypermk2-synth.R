library(hyperdags)
library(hyperinf)
library(ggraph)
library(ggpubr)

sf = 2

# simply returns a binary (character string) of length len from a decimal
DecToBinS <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

#### bilinear dynamics, compare reversible and irreversible Mk2

set.seed(5)
sim.dyn = simulate_accumulation(50, 6, dynamics="bilinear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
plot_hyperinf_data(m, tree)
fit = fit.r = fit.0 = list()
for(i in 1:10) {
  fit[[i]] = hyperinf(m, tree, method="hypermk2", reversible = FALSE)
  fit.r[[i]] = hyperinf(m, tree, method="hypermk2")
  fit.0[[i]] = hyperinf(m, tree, method="hypermk2", compare.null = TRUE)
}

# visualise transition graphs from different instances
plot.comp.1 = ggarrange(
  plot_hyperinf_data(m, tree),
  plot_hyperinf_comparative(fit, style="full"),
  plot_hyperinf_comparative(fit.r, style="full", bend=1),
  labels=c("A", "B", "C"), nrow=1
)
png("plot-comp-1.png", width=900*sf, height=300*sf, res=72*sf)
print(plot.comp.1)
dev.off()

my.labels = function(x) {
  return(sapply(as.numeric(x), DecToBinS, len=6))
}

plot.out.1 = ggarrange(plot_hyperinf_data(m, tree),
                       ggarrange(
                         plot_hyperinf(fit[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
                         plot_hyperinf(fit.r[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
                         nrow = 2, labels=c("B", "C")),
                       nrow = 1, labels=c("A", ""))

png("plot-out-1.png", width=500*sf, height=500*sf, res=72*sf)
print(plot.out.1)
dev.off()

ggpubr::ggarrange(plot_hyperinf(fit[[1]]) + ggraph::geom_node_text(ggplot2::aes(label=name)),
                  plot_hyperinf(fit.r[[1]]) + ggraph::geom_node_text(ggplot2::aes(label=name)) )

hyperinf_AIC(fit[[1]])   
hyperinf_AIC(fit.r[[1]]) 
fit.0[[1]]$null.fit$AIC

###### linear dynamics, compare Mk and Mk2

set.seed(5)
sim.dyn = simulate_accumulation(20, 4, dynamics="linear")
tree = sim.dyn$my.tree
m <- do.call(rbind, sim.dyn$x)[1:length(sim.dyn$my.tree$tip.label),]
rownames(m) = tree$tip.label
plot_hyperinf_data(m, tree)
fit.mk = hyperinf(m, tree, method="hypermk")
fit.mk.ir = hyperinf(m, tree, method="hypermk", reversible = FALSE)

fit.mk2 = fit.mk2.0 = fit.mk2.ir = list()
for(i in 1:10) {
  fit.mk2[[i]] = hyperinf(m, tree, method="hypermk2")
  fit.mk2.0[[i]] = hyperinf(m, tree, method="hypermk2", compare.null= TRUE)
  fit.mk2.ir[[i]] = hyperinf(m, tree, method="hypermk2", reversible = FALSE)
}

my.labels = function(x) {
  return(sapply(as.numeric(x), DecToBinS, len=4))
}

ggarrange(plot_hyperinf_data(m, tree),
          ggarrange(
            plot_hyperinf(fit.mk) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
            plot_hyperinf(fit.mk2[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
            plot_hyperinf(fit.mk.ir) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
            plot_hyperinf(fit.mk2.ir[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
            nrow=2, ncol=2, labels=c("Bi", "Ci", "ii", "ii")
          ), nrow=1, widths=c(1,1.5), labels=c("A", "") )

lapply(list(fit.mk, fit.mk2[[1]], fit.mk.ir, fit.mk2.ir[[1]]), hyperinf_AIC)
fit.mk2.0$null.fit$AIC

# visualise transition graphs from different instances
plot.comp.2 = ggarrange(plot_hyperinf_data(m, tree),
                        plot_hyperinf_comparative(fit.mk2, style="full"),
                        plot_hyperinf_comparative(fit.mk2.ir, style="full"),
                        labels=c("A", "B", "C"), nrow=1
)

png("plot-comp-2.png", width=900*sf, height=300*sf, res=72*sf)
print(plot.comp.2)
dev.off()

plot.out.2 = ggarrange(plot_hyperinf_data(m, tree),
                       ggarrange(
                         plot_hyperinf(fit.mk) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
                         plot_hyperinf(fit.mk2[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
                         plot_hyperinf(fit.mk.ir) + geom_node_text(size=3, angle=45, aes(label=my.labels(name)))+  coord_cartesian(clip = "off"),
                         plot_hyperinf(fit.mk2.ir[[1]]) + geom_node_text(size=3, angle=45, aes(label=my.labels(name))) +  coord_cartesian(clip = "off"),
                         nrow = 2, ncol = 2, labels=c("Bi", "Ci", "ii", "ii")),
                       nrow = 1, labels=c("A", ""))

png("plot-out-2.png", width=800*sf, height=400*sf, res=72*sf)
print(plot.out.2)
dev.off()

hyperinf_AIC(fit.mk)
hyperinf_AIC(fit.mk2[[1]])
plot_hyperinf(fit.mk)
plot_hyperinf(fit.mk2[[1]])
