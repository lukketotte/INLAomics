library(INLA)
source("../scripts/SPOTS/helpers.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")

fig_pairs = list("CD3" = "Cd3e",
                 "F480" = "Adgre1",
                 "CD163" = "Cd163",
                 "CD29" = "Itgb1",
                 "CD68" = "Cd68",
                 "IgM" = "Ighm",
                 "CD38" = "Cd38",
                 "MadCAM1" = "Madcam1",
                 "EpCAM" = "Epcam",
                 "CD11b" = "Itgam",
                 "CD105" = "Eng",
                 "CD31" = "Pecam1",
                 "CD20" = "Ms4a1",
                 "CD169" = "Siglec1",
                 "IgD" = "Ighd",
                 "CD4" = "Cd4",
                 "CD8" = "Cd8a",
                 "CD19" = "Cd19",
                 "B220" = "Ptprc")

## SPLEEN
df = SpotsProteinData("~/Documents/postdoc/MCAR/data/spots/spleen/", fig_pairs)

#Visualize the AARs
df %>%
  mutate(AARs = factor(AARs, levels = c("PALS","Marginal zone","Red pulp", "B follicle"))) %>%
  ggplot() + geom_point(aes(x = imagecol, y = imagerow, color = AARs)) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(text = element_text(size = 14), 
        legend.text.align = 0,
        legend.title.align=0.2,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=10),
        #panel.spacing.y = unit(0.8, "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")
#ggsave("/Users/lukar818/Documents/postdoc/research/RnaProt/plots/aarcolor.png", width = 9, height = 5, dpi=700)
ggsave(paste("/Users/lukar818/Documents/postdoc/research/RnaProt/plots/pdfplots/spleen", "/annotation.pdf", sep = ""), width = 9, height = 5.2, dpi=600)


protid = 16
preds = c("Cd3e", "Adgre1", "Cd163")
k = length(preds)
aar = c("pulp", "bf", "mz", "pals")
X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))

k = length(preds)
X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))
mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                  "idx" = 1:(k*nrow(df)), 
                  "size" = rep(df$size_rna, k))
mdat = cbind(mdat, X)
names(mdat)[4:ncol(mdat)] = paste(aar, rep(paste("_", 1:k, sep = ""), each = 4), sep = "")
rnaform = as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse= "+"), "-1"))

spotsInla = function(df, W, protein, preds, aar = c("pulp", "bf", "mz", "pals")){
  k = length(preds)
  if(k == 1){
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
    rnaform <- as.formula(paste(paste(preds, "~") , paste(aar[2:4],collapse = "+")))
    l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, 
                  family = "poisson", offset = log(size_rna))
    
    m <- inla.CCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
  } else {
    X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))
    mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                      "idx" = 1:(k*nrow(df)), 
                      "size" = rep(df$size_rna, k))
    mdat = cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] = paste(aar, rep(paste("_", 1:k, sep = ""), each = 4), sep = "")
    rnaform = as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse= "+"), "-1"))
    
    m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
    m.car <- inla(update(rnaform, .~. + f(idx, model = m)), 
                  data = mdat, family = "poisson", offset = log(size))
    
    m <- inla.MCCAR.model(W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k), 
                          k = k, alpha.min = 0, alpha.max = 1)
  }
  
  protform <- as.formula(paste(paste(protein, "~"), paste(aar[2:4],collapse = "+")))
  mc.car <- inla(update(protform,.~. + f(idx, model = m)), 
                 data = df, family = "poisson",
                 offset = log(size_prot),
                 control.compute = list(dic = TRUE))
  
  return(mc.car)
}

## Create neighborhood matrix
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

# Recreate CD3 row
PROT = 4
protform = as.formula(paste(paste(names(fig_pairs)[PROT], "~"), paste(c("bf", "mz", "pals",unname(sapply(fig_pairs, c))), collapse = "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(idx, model = m)), data = df, family = "poisson", offset = log(size_prot))
top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
top_preds = top_preds[which(!(top_preds %in% c("(Intercept)","bf","pals","mz", fig_pairs[PROT][[1]])))]

# Protein-gene pair
preds = fig_pairs[[PROT]]
ccar = spotsInla(df, W, names(fig_pairs)[PROT], preds, aar)

## add to cd4
loc = "~/Documents/postdoc/MCAR/r/spotscorr/"

# Protein | 2 genes
preds = c(fig_pairs[[PROT]], top_preds[1])
c2car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c2car$dic$dic < ccar$dic$dic

# Protein | 3 genes
preds = c(fig_pairs[[PROT]], top_preds[1:2])
c3car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c3car$dic$dic < c2car$dic$dic

# Protein | 4 genes
preds = c(fig_pairs[[PROT]], top_preds[1:3])
c4car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c4car$dic$dic < c3car$dic$dic

# Protein | 5 genes
preds = c(fig_pairs[[PROT]], top_preds[1:4])
c4car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# DIC goes up after 5 so we stop.
c2car$summary.hyperpar

### As a loop
models = list()
for(i in 1:6){
  if(i == 1){
    preds = fig_pairs[[PROT]]
    models[[i]] = try(spotsInla(df, W, names(fig_pairs)[PROT], preds))
  } else {
    preds = c(fig_pairs[[PROT]], top_preds[1:iter])
    models[[i]] = spotsInla(df, W, names(fig_pairs)[PROT], preds)
  }
}
# choose best filtering out those which errored
models = models[which(sapply(models, class) == "inla")]
best = models[which.min(sapply(models, function(x){x$dic$dic}))][[1]]
# number of genes
k = (nrow(best$summary.hyperpar) - 2)/2

models[[1]] = car
models[[2]] = ccar
models[[4]] = ccar3

# check which converged
models = models[which(sapply(models, class) == "inla")]
best = models[which.min(sapply(models, function(x){x$dic$dic}))][[1]]
k = (nrow(best$summary.hyperpar) - 2)/2

# timeout: class(r) == "try-error" && length(grep("timeout", geterrmessage())) 
# error: class(r) == "try-error" && !length(grep("timeout", geterrmessage()))

### As a while-loop
maxiters = 6
flag = T
iter = 1

c(fig_pairs[[PROT]], top_preds[1:2])

while(flag && iter <= maxiters){
  if(iter  == 1){
    preds = fig_pairs[[PROT]]
    car = spotsInla(df, W, names(fig_pairs)[PROT], preds)
  } else{
    preds = c(fig_pairs[[PROT]], top_preds[1:iter])
    temp = try(spotsInla(df, W, names(fig_pairs)[PROT], preds))
  }
  
  if(class(temp) == "inla"){
    # check if DIC lowered
    if(car$dic$dic > temp$dic$dic){
      car = temp
    } else {
      dicUp = iter
    }
  }
}
