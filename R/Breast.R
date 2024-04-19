library(INLA)
library(scico)
library(foreach)
library(rstan)
source("../scripts/SPOTS/helpers.R")
source("../scripts/utils.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/indepMCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")
source("../INLA/spotCCAR.R")
source("../INLA/spotMCCAR.R")


##### Predictions for breast cancer
breast = readSpotscancerlist("~/Documents/postdoc/MCAR/data/spots/cancer/")
aar = c("fbh", "fbl", "mac2", "unknown", "lymph", "mac1")
#W = cancerNbhdMat(breast)
genepairs = list(
  "CD4" = "Cd4",
  "CD8a" = "Cd8a",
  "CD366" = "Havcr2",
  "CD279" = "Pdcd1",
  "Ly-6g" = "Ly6g",
  "Ly-6C" = "Ly6c1",
  "CD19" = "Cd19",
  "Cd45" = "Ptprc",
  "CD25" = "Il2ra",
  "CD11c" = "Itgax",
  "F4/80" = "Adgre1",
  "I-A/I-E" = "H2-Ab1", 
  "NK-1.1" = "Klrb1c",
  "Ly-6A/E" = "Ly6a",
  "CD274" = "Cd274",
  "CD86" = "Cd86",
  "CD192 (CCR2)" = "Ccr2",
  "CD326" = "Epcam",
  "CD38" = "Cd38",
  "IgD" = "Ighd",
  "CD140a" = "Pdgfra",
  "CD11a" = "Itgal",
  "P2X7R" = "P2rx7",
  "CD1d" = "Cd1d1",
  "Notch 4" = "Notch4",
  "CD31" = "Pecam1",
  "Podoplanin" = "Pdpn",
  "CD45R/B220" = "Ptprc",
  "CD11b" = "Itgam",
  "CD202b" = "Tek"
)

# based on Fig. 6e
genes = list(lymph = c("B2m", "Retnla", "Postn", "Gabpa", "Spp1", "Cd74", "Ccr2", "Apoe", "Cd14", "Igha",
                       "Ighg2c", "Jchain", "Ctsb", "Chil1", "Gas6"),
             mac2 =  c("Cd74", "H2-Ab1", "H2-Aa", "Ccl8", "Apoe", "C1qa", "Spp1", "Irf7"),
             mac1 = c("Trf", "Ikbkb", "Epcam"),
             fbl = c("Krt7", "Krt8", "Krt18"),
             fbh = c("Col4a1", "Col4a2", "Col5a1", "Acta2", "Vegfa", "Vim", "Mmp2", "Bsg"))


as.data.frame(t(breast$RNA)) %>%
  select(all_of(unname(sapply(genepairs, c)))) %>%
  colSums()

as.data.frame(t(breast$RNA)) %>%
  select(all_of(genes$fbh)) %>%
  colSums()

# TODO: choose based on correlations instead
prot = "CD117"
#(prot = rownames(breast$Protein)[15])
rnaids = which(rownames(breast$RNA) %in% unname(unlist(sapply(genes, c))))
rnaids = which(rownames(breast$RNA) %in% genes$fbl)
test = apply(breast$RNA[rnaids,], 1, function(x){cor(x, breast$Protein[which(rownames(breast$Protein) == prot),])})
sort(abs(test), T)

df = data.frame("spot" = dimnames(breast$Protein)[[2]], 
                "prot" = unname(breast$Protein[which(rownames(breast$Protein) == prot),]), 
                "size_prot" = unname(colSums(breast$Protein)) / median(colSums(breast$Protein)), 
                "size_rna" = unname(colSums(breast$RNA)) / median(colSums(breast$RNA)),
                "idx" = 1:ncol(breast$Protein)) %>%
  full_join(., breast$coords, by = "spot") %>%
  full_join(breast$AAR, by = "spot") %>%
  select(!AARs)

rownames(breast$RNA[rnaids,])
rna = data.frame(t(breast$RNA[rnaids, ]), row.names = c())
names(rna) = paste("rna", 1:length(rnaids), sep = "")

df = cbind(df, rna)

df %>% ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = prot)) +
  scale_color_scico(palette = "vik")

rnaidx = 1:length(genes$lymph)
rna = t(breast$RNA[which(rownames(breast$RNA) %in% names(sort(abs(test), T)[rnaidx])),])

rna = as.data.frame(t(breast$RNA)) %>% select(all_of(genes))
rownames(rna) = NULL
rna = as.data.frame(t(breast$RNA)) %>%
  select(all_of(names(sort(test, T)[rnaidx])))
names(rna) = paste("rna", 1:length(genes$lymph), sep = "")
df = cbind(df, rna)

# For the paired approach
rna = data.frame(t(breast$RNA), row.names = c()) %>% select(all_of(c("Ly6a", "Bsg")))
names(rna) = paste("rna", 1:ncol(rna), sep = "")
df = cbind(df, rna)
# transform protein to the same scale as rna pair
#df$prot = as.integer(df$prot * sqrt(var(df$rna1)/var(df$prot)))

df %>% ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = rna1)) +
  scale_color_scico(palette = "vik")

#idx = sample(1:2000, 100, F)
dat = predData("CD11c", NULL, breast, aar, geneindex = sample(1:2000, 200, F))
df = dat$df
set.seed(1); pred_idx = sample(1:nrow(df), 200)
set.seed(12); pred_idx = sample(which(df$imagecol <= 200), 200)
pred_idx = which(df$imagecol >= 100 & df$imagecol<= 200 & df$imagerow <= 550 & df$imagerow >= 435)

df %>%
  mutate(istest = ifelse(1:nrow(df) %in% pred_idx, 1, 0)) %>%
  ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = istest))

df_pred = df[pred_idx,]
df$prot[pred_idx] = NA

coordsmat = cbind(df$imagerow, df$imagecol)
W = matrix(0, nrow(coordsmat), nrow(coordsmat))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    cp = crossprod(coordsmat[i,] - coordsmat[j,])
    if(cp > 0 && cp < 55){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

# go through all genes
my.cluster <- parallel::makeCluster(5)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:(sum(str_detect(names(df), "^rna[0-9]*$"))+1)) %dopar% {
  if(i == 1){
    try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = F, family = c("poisson", "poisson")))
  } else{
    #try(spotsInla(df, W, "prot", paste("rna",1:(i-1),sep=""), aar = aar, neighbors = T))
    #try(spotsInla(df, W, "prot", rnas[[i-1]], aar = aar, neighbors = T))
    try(spotsInla(df, W, "prot", paste("rna", i-1, sep = ""), aar = aar, neighbors = T, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred, T) %>% arrange((dic))
displayResults(models, df, df_pred, T)

models[[3]] = try(spotsInla(df, W, "prot", c("rna1", "rna2"), aar = aar, neighbors = F, ind = F, family = c("poisson", "poisson")))

models[[1]]$summary.hyperpar
models[[3]]$summary.hyperpar
prot

# ordered way
cores = 3
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:cores) %dopar% {
  if(i == 1){
    try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = T, family = c("poisson", "poisson")))
  } else{
    try(spotsInla(df, W, "prot", paste("rna",1:(i-1),sep=""), aar = aar, neighbors = T, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred, F) %>% arrange((dic))

models[[11]]$summary.hyperpar
displayResults(models2, df, df_pred, F)


df_pred = df_pred %>%
  mutate(lambda0 = models[[1]]$summary.fitted.values$mean[pred_idx],
         lambda1 = models[[2]]$summary.fitted.values$mean[pred_idx])

df_pred %>%
  pivot_longer(cols = c(prot,lambda0, lambda1)) %>%
  ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = value)) +
  facet_grid(name ~ .) +
  scale_color_scico(palette = "vik")

df_pred %>%
  mutate(res1 = scale(prot - lambda0),
         res2 = scale(prot - lambda1)) %>%
  pivot_longer(cols = c(res1,res2)) %>%
  ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = value)) +
  facet_grid(name ~ .) +
  scale_color_scico(palette = "vik")
  


rbind(displayResults(models2, df, df_pred, F) %>% mutate(neighbors = F),
      displayResults(models, df, df_pred, T) %>% mutate(neighbors = T)) %>%
  arrange((dic))

# Find the first one, move on to the next
idx = c(1:5, 7:length(genes))
my.cluster <- parallel::makeCluster(5)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models2 = foreach(i = idx) %dopar% {
  cat("iter = ", i, "\n", sep = "")
  try(spotsInla(df, W, "prot", c("rna6", paste("rna", i-1, sep = "")), aar = aar, neighbors = T))
}
parallel::stopCluster(cl = my.cluster)



# just 1 rna but trying all
cores = 2
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:cores) %dopar% {
  cat("iter = ", i, "\n", sep = "")
  if(i == 1){
    try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = T, family = c("poisson", "poisson")))
  } else{
    try(spotsInla(df, W, "prot", paste("rna",i-1,sep=""), aar = aar, neighbors = T,family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

loc = "/Users/lukar818/Documents/postdoc/research/RnaProt/plots/bcancer/"

#top_preds = names(sort(abs(test), T)[1:3])
models = models[c(1, 3)]
test = models
models = models[c(1,9)]

models$preds = c("Vim")
models$preds = dat$top_preds[1]
models$df = df_pred
prot
#save(models, file = paste(loc, "F480_cancer.RData"))
save(models, file = paste(loc, "CD117_cancer_unpaired.RData", sep = ""))

### STAN BREAST
ccar = try(spotsInla(df, W, "prot", "rna1", aar = aar, neighbors = T))

D = diag(colSums(W))
n = nrow(W)
W_n = sum(W)/2
W_sparse = matrix(0, nrow = W_n, ncol =2)
counter = 1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(W[i,j] == 1){
      W_sparse[counter,1] = i
      W_sparse[counter,2] = j
      counter = counter + 1
    }
  }
}
D_sparse = colSums(W)
lambda = eigen(D-W, T)$values

x = cbind(1, df %>% select( fbl, mac1, mac2, unknown, lymph) %>% as.matrix())

## STAN
sp_d <- list(n = nrow(W),
             p = ncol(x),
             N = 2*nrow(df),
             y = c(df$rna1, df$prot),
             X = x,
             W = W,
             W_n = sum(W)/2,
             W_sparse = W_sparse,
             D_sparse = D_sparse,
             lambda = lambda)

sp_fit <- stan('CCAR.stan', data = sp_d, iter = 2e4, warmup = 2000, thin = 5, chains = 1)
posterior <- extract(fit)