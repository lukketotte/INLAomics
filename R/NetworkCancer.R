library(INLA)
library(scico)
source("../scripts/SPOTS/helpers.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/spotCCAR.R")
source("../INLA/MCCAR.R")
source("../INLA/spotMCCAR.R")


## Protein to rna DB
db = hgnc::import_hgnc_dataset("https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt")
names(db)
db$name[str_detect(db$name, "^CD.*$")]
db$name[str_detect(db$name, "^Ig.*$")]
which(db$name == "CD3 delta subunit of T-cell receptor complex")
as.data.frame(db[2499,])
names(db)
db$locus_group

inla.setOption(inla.timeout = 10*60)

# See figure 2e of Spots
breast_genes = list(fibhigh = c("Col4a1", "Col4a2", "Col5a1", "Acta2", "Vegfa", "Vim", "Mmp2", "Bsg"),
                    fiblow = c("Krt7", "Krt8", "Krt18"),
                    mac2 = c("Cd74", "H2-Ab1", "H2-Aa", "Ccl8", "Apoe", "C1qa", "Spp1", "Irf7"),
                    mac1 = c("Trf", "Ikbkb", "Epcam"),
                    lymph = c("B2m", "Retnla", "Postn", "Gapdh", "Ccr2", "Apoe", "Cd14", "Igha",
                              "Ighg2c", "Jchain", "Ctsb", "Chil1", "Gas6"))
genes = unname(unlist(breast_genes))

df = SpotsCancerData("~/Documents/postdoc/MCAR/data/spots/cancer/", genes) %>%
  select(!c(AARs, spot)) %>%
  mutate(idx = 1:nrow(.))

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

names(df)

PROT = "CD86"
aar = c("fbh", "fbl", "mac2", "unknown", "lymph", "mac1")
protform = as.formula(paste(paste(PROT, "~"), 
                            paste(c(aar[-1], names(df)[names(df) %in% str_replace(genes, "-", "_")]), 
                                  collapse = "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(idx, model = m)), data = df[-pred_idx,], family = "poisson", offset = log(size_prot))
top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
(top_preds = top_preds[which(!(top_preds %in% c("(Intercept)", aar[-1])))])

df %>%
  ggplot() + geom_point(aes(x = imagerow, y = imagecol, color = Krt7)) +
  scale_color_scico(palette = "vik") + theme_bw() + scale_x_reverse() +
  scale_y_reverse() + xlab("") + ylab("")

set.seed(1234); pred_idx = sample(1:nrow(df), 198) 
df_pred = df[pred_idx,]
df[pred_idx, names(df)== PROT] = NA

car = spotsInla(df, W, PROT, character(), aar = aar, family = c("poisson", "poisson"))
ccar = spotsInla(df, W, PROT, "Krt7", aar = aar, family = c("poisson", "poisson"))
ccar2 = spotsInla(df,W,PROT, c("Acta2", "Chil1"), aar = aar, neighbors = F, family = c("poisson", "nbinomial"))
ccar3 = spotsInla(df,W,PROT, top_preds[1:3], aar = aar, family = c("poisson", "nbinomial"))

car$dic$dic
ccar$dic$dic
summary(car)

car$summary.hyperpar
ccar$summary.hyperpar
ccar2$summary.hyperpar
-mean(log(car$cpo$cpo[-pred_idx]))
-mean(log(ccar$cpo$cpo[-pred_idx]))
-mean(log(ccar2$cpo$cpo[-pred_idx]))

sqrt(mean((df_pred[,names(df_pred) == PROT] - car$summary.fitted.values$mean[pred_idx])^2))
sqrt(mean((df_pred[,names(df_pred) == PROT] - ccar$summary.fitted.values$mean[pred_idx])^2))
sqrt(mean((df_pred[,names(df_pred) == PROT] - ccar2$summary.fitted.values$mean[pred_idx])^2))

sqrt(mean((df[-pred_idx,names(df) == PROT] - car$summary.fitted.values$mean[-pred_idx])^2))
sqrt(mean((df[-pred_idx,names(df) == PROT] - ccar$summary.fitted.values$mean[-pred_idx])^2))
sqrt(mean((df[-pred_idx,names(df) == PROT] - ccar2$summary.fitted.values$mean[-pred_idx])^2))

# logic should be ccar_0,ccar_1, ..., ccar_k

m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
rnaform <- as.formula(paste(paste(top_preds[1], "~") , paste(aars[-1],collapse = "+")))
l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, 
              family = "zeroinflatedpoisson1", offset = log(size_rna))

l.car$summary.hyperpar


m <- inla.CCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
mc.car <- inla(update(protform,.~. + f(idx, model = m)), 
               data = df, family = "poisson",
               offset = log(size_prot),
               control.compute = list(dic = TRUE))


k = 2
naars = length(aar)
preds = top_preds[1:2]
X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))

mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                  "idx" = 1:(k*nrow(df)), 
                  "size" = rep(df$size_rna, k))
mdat = cbind(mdat, X)

ncol(mdat)
paste(aar, rep(paste("_", 1:k, sep = ""), each = naars), sep = "")
mdat
names(mdat)[4:ncol(mdat)] = paste(aar, rep(paste("_", 1:k, sep = ""), each = naars), sep = "")
rnaform = as.formula(paste("rna ~", paste(names(mdat)[naars:(ncol(mdat))], collapse= "+"), "-1"))

m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
m.car <- inla(update(rnaform, .~. + f(idx, model = m)), 
              data = mdat, family = "poisson", offset = log(size))



m <- inla.MCCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = matrix(m.car$summary.random$idx$mean, ncol = k), k = 2)
mc2.car <- inla(update(protform,.~. + f(idx, model = m)), 
                data = df, family = "poisson",
                offset = log(size_prot),
                control.compute = list(dic = TRUE))

mc2.car$summary.hyperpar
mc2.car$dic$dic
mc.car$dic$dic
c.car$dic$dic

car = spotsInla(df,W,"CD4",character(), aar = aar)

car$dic$dic
c.car$dic$dic

mc.car2 = spotsInla(df,W,"CD4", top_preds[1], aar = aar, family = c("zeroinflatedpoisson1", "poisson"))

mc.car2$dic$dic

mc.car2$summary.hyperpar
mc.car$summary.hyperpar
