library(INLA)
library(tidyverse)
library(scico)
library(foreach)
source("../scripts/highplex/helpers.R")
source("../scripts/utils.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")
source("../INLA/spotCCAR.R")
source("../INLA/spotMCCAR.R")

loc = "~/Documents/postdoc/MCAR/data/highplexdata/"
saveloc = "/Users/lukar818/Documents/postdoc/research/RnaProt/plots/highplex/"

## human skin
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*skin_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## human thymus
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*thymus_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## human tonsil
protein = read.table(paste(loc, "GSM6578071_humantonsil_protein.tsv.gz", sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## mouse kdieny
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*kidney_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

names(protein) = str_replace(names(protein), "^(C{1}D{1}[0-9]+[A-z]*)\\.+.*$", "\\1") %>% 
  str_replace(., "[..]*[A-Z]{7,}", "")

protein$size = rowSums(protein[, 1:(ncol(protein)-2)]) / median(rowSums(protein[, 1:(ncol(protein)-2)]))


#max(rowSums(protein[, 1:(ncol(protein)-2)]) / mean(rowSums(protein[, 1:(ncol(protein)-2)])))
## human skin
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*skin_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## human thymus
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*thymus_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## human tonsil
rna = read.table(paste(loc,"GSM6578062_humantonsil_RNA.tsv.gz", sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## mouse kidney
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*kidney_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)


names(rna) = str_replace(names(rna), "^CD(.*)$", "Cd\\1")
rna_top = names(sort(colSums(rna[,!colnames(rna) %in% c("x", "y", "size")]), T)[1:250])
rna$size = rowSums(rna[, 1:(ncol(rna)-2)]) / median(rowSums(rna[, 1:(ncol(rna)-2)]))

###
names(sort(colSums(protein[,!colnames(protein) %in% c("x", "y")]), T))[1:60]

protein %>% select(x,y,CD223) %>%
  ggplot() + geom_point(aes(x = x, y = y, color = CD223))+
  scale_color_scico(palette = "vik")

top_preds[which(!(top_preds %in% c("(Intercept)")))][1:20]

rna %>% select(x,y,LAGE3) %>%
  ggplot() + geom_point(aes(x = x, y = y, color = LAGE3))+
  scale_color_scico(palette = "vik")

mdat = full_join(rna %>% select(all_of(rna_top), x, y), 
                 protein %>% select(x,y,CD107a,size), 
                 by = c("x","y"), suffix = c("rna", "prot")) %>%
  mutate(id = 1:nrow(.))

coords = cbind(mdat$x, mdat$y)
W = matrix(0, nrow(coords), nrow(coords))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    x = coords[i,]
    if(crossprod(coords[i,] - coords[j,]) == 1){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

protform = as.formula(paste("CD107a~ ", paste(names(mdat)[1:length(rna_top)], collapse= "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(id, model = m)), 
                 data = mdat, family = "poisson", offset = log(size))

top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
top_preds = top_preds[which(!(top_preds %in% c("(Intercept)")))][1:5]

df = full_join(rna %>% select(all_of(c("x", "y", c("LAMP1",top_preds[1:4]), "size"))), 
               protein %>% select(x,y,CD107a,size), 
               by = c("x","y"), suffix = c("rna", "prot")) %>%
  mutate(idx = 1:nrow(.))
names(df) = c("x", "y", paste("rna", 1:5, sep = ""), "size_rna", "prot", "size_prot", "idx")

df %>% ggplot() + geom_point(aes(x = x, y = y, color = prot)) + scale_color_scico(palette = "vik")

coords = cbind(df$x, df$y)
W = matrix(0, nrow(coords), nrow(coords))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    x = coords[i,]
    if(crossprod(coords[i,] - coords[j,]) == 1){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

#set.seed(12); pred_idx =sample(which(df$x >= 15 & df$x <= 40 & df$y <= 42 & df$y >= 27), 169)
set.seed(12); pred_idx =sample(1:nrow(df), as.integer(nrow(df)*0.1))
df %>% mutate(istest = ifelse((1:nrow(.)) %in% pred_idx, 1, 0)) %>% ggplot() + geom_point(aes(x = x, y = y, color = istest))
  
pred_df = df[pred_idx,]
df$prot[pred_idx] = NA

cores = 4
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:cores) %dopar% {
  cat("iter = ", i, "\n", sep = "")
  if(i == 1){
    try(hpInla(df, W, "prot", character(), neighbors = F, family = c("poisson", "poisson")))
  } else{
    try(hpInla(df, W, "prot", paste("rna",1:(i-1),sep=""), neighbors = F, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

data.frame(
  k = sapply(models, function(x){(nrow(x$summary.hyperpar)-2)}),
  dic = sapply(models, function(x){x$dic$dic}),
  rmse = sapply(models, function(x){sqrt(mean((pred_df$prot - x$summary.fitted.values$mean[pred_idx])^2))}),
  rmsetrain = sapply(models, function(x){sqrt(mean((df$prot[-pred_idx] - x$summary.fitted.values$mean[-pred_idx])^2))}),
  cpo = sapply(models, function(x){sum(log(x$cpo$cpo[-pred_idx]))}),
  waic = sapply(models, function(x){x$waic$waic})
)

models[[1]]$summary.hyperpar
models[[2]]$summary.hyperpar
models[[3]]$summary.hyperpar
models[[4]]$summary.hyperpar

test = try(hpInla(df, W, "prot", "rna3", neighbors = F, family = c("poisson", "poisson")))

sqrt(mean((pred_df$prot - test$summary.fitted.values$mean[pred_idx])^2))
which.max((pred_df$prot - models[[2]]$summary.fitted.values$mean[pred_idx])^2)
ggplot() + geom_density(aes(x = (pred_df$prot - models[[2]]$summary.fitted.values$mean[pred_idx])^2))

(pred_df$prot[48]-models[[2]]$summary.fitted.values$mean[pred_idx[48]])^2

pred_df$prot
models[[2]]$summary.fitted.values$mean[pred_idx]
(38-18)^2

pred_df[48,]
which.max(models[[1]]$summary.fitted.values$mean)
models[[1]]$summary.fitted.values$mean[1601]

models$df = df
models$pred_df = pred_df
models$preds = c(top_preds[1:4])

save(models, file = paste(saveloc, "CD223_skin.RData", sep = ""))
