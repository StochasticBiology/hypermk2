library(hyperinf)

run.inference = TRUE
n.samples = 1

load("kp-test-data.Rdata")
load("kp-test-tree.Rdata")

rownames(m) = tree$tip.label

fit.hmk2 = fit.hmk21 = list()
if(run.inference == TRUE) {
  
  for(i in 1:n.samples) {
  fit.hmk2[[i]] = hyperinf(m, tree, reversible=TRUE)
  fit.hmk21[[i]] = hyperinf(m, tree, reversible=TRUE, force.origin=TRUE)
  }
  
  load("kp-hypertraps-fit.Rdata")
  res.tmp$feature.names = res.tmp$featurenames
  fit.ht = res.tmp
  
  fit.set = list(fit.hmk2, fit.hmk21, fit.ht)
  save(fit.set, file="fits-kp-test-hmk2-hmk21-ht-many.Rdata")
}

load("fits-kp-test-hmk2-hmk21-ht.Rdata")
fit.hmk2 = fit.set[[1]]
fit.hmk21 = fit.set[[2]]
fit.ht = fit.set[[3]]

ordering_matrix(fit.ht)
plot.kp.1 = plot_hyperinf_ordering_matrices(list(fit.hmk2, fit.ht), 
                                            feature.names = substr(fit.hmk2$feature.names, 1,5),
                                            expt.names=c("Mk2", "TraPS"), 
                                            type="relative")

custom_label <- function(x) {
  ifelse(x == "16", "{Bla_chr}",
         ifelse(x == "24", "{Bla_chr,SHV_m}", ""))
}

plot.kp.2 = plot_hyperinf(fit.hmk2, threshold = 0.02,
                          feature.names = substr(fit.hmk2$feature.names, 1,5)) + 
  ggraph::geom_node_text(ggplot2::aes(label=custom_label(name)))

mred = m
colnames(mred) = substr(fit.hmk2$feature.names, 1,5)
kp.plot.set = ggpubr::ggarrange(plot_hyperinf_data(mred, tree, font.size=3) ,
                                ggpubr::ggarrange(plot.kp.2, plot.kp.1, nrow=2, labels=c("B", "C")), 
                                labels=c("A", ""), widths=c(1,2), nrow=1)

sf = 2
png("kp-plot-set.png", width=800*sf, height=600*sf, res=72*sf)
print(kp.plot.set)
dev.off()



