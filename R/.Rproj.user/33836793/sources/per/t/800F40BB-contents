## Script prediction on the spots spleen
library(INLA)
library(scico)
library(tidyverse)
library(foreach)
library(rstan)

source("../scripts/SPOTS/helpers.R")
source("../scripts/utils.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")
source("../INLA/spotCCAR.R")

# CD4 prediction example
spots = readSpotsSpleen("~/Documents/postdoc/MCAR/data/spots/spleen/")
aar = c("pulp", "bf", "mz", "pals")
names(spots)[2] = "Protein"
dat = predData("CD4", NULL, spots, aar, geneindex = 1:200, genepair = "Cd4")
df = dat$df

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

x = cbind(1, df %>% select(bf, mz, pals) %>% as.matrix())

#m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
#rnaform <- as.formula(paste(paste("rna1", "~") , paste(aar[-1],collapse = "+")))
#l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, 
#              family = "poisson", offset = log(size_rna))

ccar = try(spotsInla(df, W, "prot", "rna1", aar = aar, neighbors = F))

ccar$summary.fixed
ccar$summary.hyperpar
ccar$summary.random$idx$mean

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
             lambda = lambda,
             log_size_rna = df$size_rna,
             log_size_prot = df$size_prot)

init_fun <- function(lcar, ccar){
  list(list(b2 = ccar$summary.fixed$mean, 
            theta2 = ccar$summary.hyperpar$mean[2],
            tau2 = exp(ccar$summary.hyperpar$mean[2]),
            nu2 = ccar$summary.hyperpar$mean[1],
            eta0 = ccar$summary.hyperpar$mean[3],
            alpha2 = plogis(ccar$summary.hyperpar$mean[1]),
            phi2 = ccar$summary.random$idx$mean,
            b1 = lcar$summary.fixed$mean, 
            theta1 = lcar$summary.hyperpar$mean[2], 
            tau1 = exp(lcar$summary.hyperpar$mean[2]),
            nu1 = lcar$summary.hyperpar$mean[1],
            alpha1 = plogis(lcar$summary.hyperpar$mean[1]),
            phi1 = lcar$summary.random$idx$mean))
}

sp_fit <- stan('CCAR.stan', data = sp_d, iter = 1000, warmup = 500, thin = 1, chains = 1, seed = 123)
posterior <- extract(fit)

?plogis

