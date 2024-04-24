## Script prediction on the spots spleen
library(INLA)
library(scico)
library(tidyverse)
library(foreach)
library(rstan)

source("../scripts/SPOTS/helpers.R")
source("../scripts/highplex/helpers.R")
source("../scripts/utils.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")

# CD4 prediction example
spots = readSpotsSpleen("~/Documents/postdoc/MCAR/data/spots/spleen/")
aar = c("pulp", "bf", "mz", "pals")
names(spots)[2] = "Protein"
rownames(spots$Protein)
dat = predData("F4-80", NULL, spots, aar, geneindex = 1:200, genepair = "Adgre1")

df = dat$df
pred_idx = which(df$imagecol > 290 & df$imagecol < 390 & df$imagerow < 310 & df$imagerow > 210)
df_pred = df[pred_idx,]
df$prot[pred_idx] = NA

coordsmat = cbind(df$imagerow, df$imagecol)
W = matrix(0, nrow(coordsmat), nrow(coordsmat))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    cp = crossprod(coordsmat[i,] - coordsmat[j,])
    if(cp > 0 && cp < 45){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

my.cluster <- parallel::makeCluster(sum(str_detect(names(df), "^rna[0-9]*$"))+1)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:(sum(str_detect(names(df), "^rna[0-9]*$"))+1)) %dopar% {
  if(i == 1){
    #try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = T, family = c("poisson", "poisson")))
    try(hpInla(df, W, "prot", character(),  neighbors = T, family = c("poisson", "poisson")))
  } else{
    #try(spotsInla(df, W, "prot", paste("rna",1:(i-1),sep=""), aar = aar, neighbors = T, family = c("poisson", "poisson")))
    try(hpInla(df, W, "prot", paste("rna",1:(i-1),sep=""), neighbors = T, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred, T)

m = spotsInla(df, W, "prot", "rna1", aar = aar, neighbors = T, family = c("poisson", "poisson"))

load("/Users/lukar818/Documents/postdoc/research/RnaProt/plots/cortext/CD4_cortex.RData")

temp = model

model$ccar = m

model = model[c(1,5,2,3,4)]
names(model)
save(model, file = "/Users/lukar818/Documents/postdoc/research/RnaProt/plots/cortext/CD4_cortex.RData")


## GLMS
# cd4glm
protform <- as.formula(paste(paste("prot", "~"), paste(aar[-1],collapse = "+")))

glmf480 <- inla(protform, data = df, family = "poisson",
               offset = log(size_prot),
               control.predictor = list(compute = TRUE, link = 1))

summary(glmf480)

glms = list("cd4" = glmcd4, "f480" = glmf480)

save(glms, file = "/Users/lukar818/Documents/postdoc/research/RnaProt/plots/other/glmsSpleen.RData")
