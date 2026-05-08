library(hyperinf)
library(ggpubr)

library(readxl)
library(phytools)
library(hypertrapsct)

system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls")
o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")

# get phylogeny linking isolates
# Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt")
tree = read.tree("41588_2014_BFng2878_MOESM34_ESM.txt")

##########
##### Curate data

# extract isolate ID and resistance profiles to our ten drugs (discarding mutation info)
# remove any incomplete profiles and recast as binary strings
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL", "PRO")
df = o.df[,which(colnames(o.df) %in% col.interest)]
missing.rows = unique(which(df == ".", arr.ind=TRUE)[,1])
final.df = df[-missing.rows,]
final.df = final.df %>%
  mutate(across(-Isolate, ~ ifelse(. == "R", 1, 0)))
final.df = as.data.frame(final.df)

##########
##### Evolutionary accumulation modelling
library(hyperinf)

set.seed(123)
small.df = final.df[sample(1:nrow(final.df), 50),]
prune.tree = keep.tip(tree, small.df$Isolate)
plot_hyperinf_data(small.df, prune.tree)

fit.tb.small = hyperinf(small.df, prune.tree, method="hypermk2")
fit.hmm = hyperinf(small.df, prune.tree)
fit.ht = hyperinf(small.df, prune.tree, method = "hypertraps")

# simply returns a binary (character string) of length len from a decimal
DecToBinS <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}
plot_hyperinf(fit.tb.small) + 
  ggraph::geom_node_text(ggplot2::aes(label=name))

DecToBinS(800, fit.tb.small$L)

plot_hyperinf_data(small.df, prune.tree)

###
#fit.tb.small.1 = hyperinf(small.df, prune.tree, method="hypermk2")
fit.tb.small.1 = fit.tb.small

#save(fit.tb.small.1, file="fit-tb-small-1.Rdata")
plot.tb.small.1 = plot_hyperinf(fit.tb.small.1, threshold = 0.05) + 
  ggraph::geom_node_text(ggplot2::aes(label=name))
plot.tb.small.2 = plot_hyperinf_ordering_matrices(list(fit.tb.small.1, fit.hmm, fit.ht),
                                                  expt.names = c("Mk2", "HMM", "TraPS"))
plot.tb.small.3 = plot_hyperinf_comparative(list(fit.tb.small.1, fit.hmm, fit.ht), 
                          style = "full",
                          expt.names = c("Mk2", "HMM", "TraPS"))
tb.trellis = ggarrange(plot_hyperinf_data(small.df, prune.tree),
          plot.tb.small.1, plot.tb.small.2,
          plot.tb.small.3, labels=c("A", "B", "C", "D"), 
          heights=c(1.5,2))

sf = 2
png("tb-trellis-mk2.png", width=600*sf, height=500*sf, res=72*sf)
print(tb.trellis)
dev.off()

#save(fit.ht, file="~/Dropbox/Documents/2026_Projects/Jakob/fitht.Rdata")
ordering_matrix(fit.ht)
###########

plot_hyperinf_ordering_matrices(fit.tb.small.1)
plot_hyperinf_ordering_matrices(fit.hmm)
plot_hyperinf_ordering_matrices(fit.ht)

######
fit.tb = hyperinf(final.df, tree, method="hypermk2")
save(fit.tb, file="fit-tb-full.Rdata")
plot_hyperinf(fit.tb, threshold = 0.25)                
