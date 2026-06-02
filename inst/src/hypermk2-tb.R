library(hyperinf)
library(ggpubr)
library(readxl)
library(phytools)
library(hypertrapsct)

run.inference = TRUE
pull.data = TRUE

if(pull.data == TRUE) {
  system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls")
  # get phylogeny linking isolates
  # Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
  system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt")
}

o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")
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

set.seed(123)
small.df = final.df[sample(1:nrow(final.df), 50),]
prune.tree = keep.tip(tree, small.df$Isolate)
plot_hyperinf_data(small.df, prune.tree)

fit.tb.small = list()
if(run.inference == TRUE) {
  for(i in 1:10) {
    fit.tb.small[[i]] = hyperinf(small.df, prune.tree, method="hypermk2")
  }
  fit.hmm = hyperinf(small.df, prune.tree)
  fit.ht = hyperinf(small.df, prune.tree, method = "hypertraps")
  
  fits = list(fit.tb.small = fit.tb.small,
              fit.hmm = fit.hmm,
              fit.ht = fit.ht)
  save(fits, file="hypermk2-tb-fits-many.Rdata")
}

load("hypermk2-tb-fits-many.Rdata")
fit.tb.small = fits$fit.tb.small
fit.hmm = fits$fit.hmm
fit.ht = fits$fit.ht

# simply returns a binary (character string) of length len from a decimal
DecToBinS <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}
plot_hyperinf(fit.tb.small[[1]]) + 
  ggraph::geom_node_text(ggplot2::aes(label=name))

plot_hyperinf_data(small.df, prune.tree)

plot_hyperinf_comparative(fit.tb.small)

fit.tb.small.1 = fit.tb.small[[1]]

plot.tb.small.1 = plot_hyperinf(fit.tb.small.1, threshold = 0.05) + 
  ggraph::geom_node_text(ggplot2::aes(label=name))

plot.tb.small.1a = plot_hyperinf_comparative(fit.tb.small[1:5], 
                                            style = "full", threshold = 0)

plot.tb.small.2 = plot_hyperinf_ordering_matrices(list(fit.tb.small.1, fit.hmm, fit.ht),
                                                  expt.names = c("Mk2", "HMM", "TraPS"),
                                                  type="relative")
plot.tb.small.2a = plot_hyperinf_ordering_matrices(list(fit.tb.small.1, fit.hmm, fit.ht),
                                                  expt.names = c("Mk2", "HMM", "TraPS"),
                                                  type="absolute")
plot.tb.small.3 = plot_hyperinf_comparative(list(fit.tb.small.1, fit.hmm, fit.ht), 
                                            style = "full",
                                            expt.names = c("Mk2", "HMM", "TraPS"))
tb.trellis = ggarrange(plot_hyperinf_data(small.df, prune.tree),
                       plot.tb.small.1a, plot.tb.small.2,
                       plot.tb.small.3, labels=c("A", "B", "C", "D"), 
                       heights=c(1.5,2))

sf = 2
png("tb-trellis-mk2.png", width=600*sf, height=500*sf, res=72*sf)
print(tb.trellis)
dev.off()
             
plot.tb.small.3a = plot_hyperinf_comparative(c(fit.tb.small[1:5], list(fit.hmm, fit.ht)), 
                                            style = "full",
                                            expt.names = c(rep("Mk2", length(fit.tb.small[1:5])), 
                                                           "HMM", "TraPS"))
