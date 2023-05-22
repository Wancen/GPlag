library(plgp)
library(mvtnorm)
library(tidyverse)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(pbapply)
library(GenomicRanges)
library(InteractionSet)
source("TMGGP_functions.R")
# epCounts %>% mcols() %>% colnames()

## Load libraries for accessing counts
library(DelayedArray)
library(SummarizedExperiment)

## Load expectedCounts
# expectedCounts <- readRDS("data/h5/expectedCounts.rds")

# ## Update path to count data
# path(assay(expectedCounts)) <- "data/h5/expectedCounts.h5"

# # ## Access count data
# obs_hic = mcols(epCounts)[,183:189]  |> as.matrix()
# M <- as(obs_hic, "HDF5Matrix")
# exp_hic = counts(expectedCounts)[,1:7]  |> as.matrix()
# expectedCounts_all = 1 + (obs_hic - exp_hic) / (exp_hic + 1)
# epCounts$hic_max = rowMax(expectedCounts_all)
# epCounts$expectedCounts <- counts(expectedCounts)[,"LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic"]


# ## 22813 pairs left
# epCounts$distance = pairdist(epCounts)
# epCounts_nonstatic <- epCounts[epCounts$anchor1.K27_cluster != "static" & epCounts$anchor2.RNA_cluster!="static" & epCounts$distance < 250000]
# # save(epCounts_nonstatic, shiftprior, theta, file = "LIMA_binnedEPCounts_nonstatic.rda")
# # load("LIMA_binnedEPCounts_nonstatic.rda")

# # epCounts_nonstatic
# peak<-epCounts_nonstatic %>% mcols() %>%
#   as.data.frame() %>%
#   dplyr::select(anchor1.K27_m0000_VST:anchor1.K27_m0360_VST) %>% t() %>% `colnames<-`(epCounts_nonstatic$anchor1.K27_name)

# gene<-epCounts_nonstatic %>% mcols() %>%
#   as.data.frame() %>%
#   dplyr::select(anchor2.RNA_m0000_VST:anchor2.RNA_m0360_VST) %>% t() %>% `colnames<-`(epCounts_nonstatic$anchor2.RNA_hgnc)

# dat <- rbind(gene,peak) %>% `colnames<-`(paste(colnames(gene),colnames(peak),sep = "-"))

# ## second filter
# peak_max = apply(peak,2,max)
# peak_index = peak_max > median(peak_max)
# gene_max = apply(gene,2,max)
# gene_index = gene_max > median(gene_max)

# filter2 = peak_index & gene_index
# dat_filter = dat[,filter2]
# epCounts_filter <- epCounts_nonstatic[filter2]

load(file = "LIMA_binnedEPCounts_nonstatic.rda")
peak<-epCounts_filter %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(anchor1.K27_m0000_VST:anchor1.K27_m0360_VST) %>% t() %>% `colnames<-`(paste0(seqnames(anchors(epCounts_filter, type="first"))," ",start(anchors(epCounts_filter, type="first")),"-", end(anchors(epCounts_filter, type="first"))))

gene<-epCounts_filter %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(anchor2.RNA_m0000_VST:anchor2.RNA_m0360_VST) %>% t() %>% `colnames<-`(epCounts_filter$anchor2.RNA_hgnc)

dat_filter <- rbind(gene,peak) %>% `colnames<-`(paste(colnames(peak),colnames(gene),sep = "-"))


t=c(0,0.5,1,1.5,2,4,6)
n=length(t)

tinput <- c(t,t) %>% as.matrix()
delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
group.index <- rep(c(1,2),each=n)
group.size <- rle(group.index)
meanGene = colMeans(dat_filter[1:n,])
meanGeneMat <- matrix(rep(meanGene, each = 2 * n), nrow = 2 * n)
Y.shift <- dat_filter - meanGeneMat

## use ccf to derive shift initialization
negative_index = floor(10*log10(n/2))+1
shiftprior <- pbsapply(1:ncol(Y.shift), function(i){
  cc <- ccf(Y.shift[1:n,i],Y.shift[(n+1):(2*n),i],plot=FALSE)
  t[cc$lag[negative_index+which.max(cc$acf[negative_index:(2*negative_index-1)])-1]+1]
},cl=18)
table(shiftprior)

# Run GP first to derive mean b as initialization
# library(GauPro)
eps <- sqrt(.Machine$double.eps)
# theta <- pbsapply(1:ncol(peak), function(i){
#   y = peak[,i] %>% unlist() |> unname()
#   y <- y-mean(y)
#   gp <- GauPro(X=t, Z=y)
#   theta <- gp$theta
#   # gpi <- newGPsep(t %>% as.matrix(), y, d = 1, g = 0.1 * var(y), dK = TRUE)
#   # mle <- mleGPsep(gpi, param="d", tmin=0.1, tmax=10)
#   # theta <- 1/mle$d
#   return(min(theta,10))
# })
# summary(theta)
# hist(theta)

library(mlegp)
theta <- pbsapply(1:ncol(dat_filter), function(i){
  y = dat_filter[(n+1):(2*n),i] %>% unlist() |> unname()
  y <- y-mean(y)
  invisible(capture.output(gp <- mlegp(matrix(t,ncol = 1), y, verbose = 0)))
  theta <- gp$beta
  # gpi <- newGPsep(t %>% as.matrix(), y, d = 1, g = 0.1 * var(y), dK = TRUE)
  # mle <- mleGPsep(gpi, param="d", tmin=0.1, tmax=10)
  # theta <- 1/mle$d
  return(min(theta,10))
}, cl=18)

# which(colnames(dat_filter)=="NAV3-RR_peak_44733")
save(epCounts_filter, shiftprior, theta, file = "LIMA_binnedEPCounts_nonstatic.rda")


res2 <- pbsapply(1:ncol(dat_filter),function(i){
  s.init <- shiftprior[i]
  b.init <- theta[i]
  out <- optim(c(softplusinv(b.init), softplusinv(1), s.init, 0, 0.01 * var(Y.shift[,i])),
          nl4, method="L-BFGS-B",
          lower = c(-5,-10,0,-7,eps),
          upper = c(20,Inf,2,7,0.1*var(Y.shift[,i])),
          t = tinput, Y = Y.shift[,i], delta = delta, group.size = group.size)
  # res[i,]<- out$par
  bhat = softplus(out$par[1])
  ahat = softplus(out$par[2])
  shat = out$par[3]
  muhat = out$par[4]
  ghat = out$par[5]
  t2hat <- t+shat
  ## derive covariance function Paper 4.1
  ### 1. same group covariance
  Ahat = (ahat^2)*delta+1
  D <-  plgp::distance(c(t,t2hat))/Ahat
  # eps <- sqrt(.Machine$double.eps)
  # D2 <- dist(X, diag = T,upper = T)^2 %>% as.matrix() ## althernative way
  Chat <- exp(-bhat*D)
  Vk = (1/sqrt(Ahat)*Chat) %>% as.matrix()+diag(ghat,nrow(Y.shift)) ## covariance matrix
  Vi <- pinv(Vk) ## psduo-inverse covariance matrix
  group.size$values <- c(0,muhat)
  mu.expand <- inverse.rle(group.size)
  Y.shifthat <- Y.shift[,i] - mu.expand

  tau2hat2 <- drop(t(Y.shifthat) %*% Vi %*% Y.shifthat / length(tinput))
  return(c(bhat,ahat,shat,muhat,ghat,tau2hat2))
},cl=18)

# knitr::kable(
#   res %>% `rownames<-`(c("b","a","s","mu","tau2")),
#   col.names = colnames(Y)
# )
write_csv(res2 %>% as.data.frame() %>% `colnames<-`(colnames(dat_filter)),file="LIMA_binnedEPCounts_nonstatic.rbf.estimated.param.csv")

# load("LIMA_binnedEPCounts_nonstatic_res.rda")
# library(ggplot2)
# pdf("../plots/Gene-epCounts_filter-posterior.pdf",         # File name
#     width = 8, height = 5, # Width and height in inches
#     colormodel = "cmyk")
# par(mfrow = c(2, 3))
ci <-pbsapply(1: ncol(dat_filter),function(i){
  params <- res2[,i]
  cor <- interpolation.rbf(t,Y_shift=Y.shift[,i],delta,group.size,params,ntest=100, plot = FALSE)
  return(c(cor[[1]],quantile(cor[[2]],c(0.05,0.95))))
},cl=18)

# plotPair(Y.shift,2,t, s=c(0,shat), mu=muhat)
write_csv(ci %>% as.data.frame() %>% `colnames<-`(colnames(dat_filter)),file="LIMA_binnedEPCounts_nonstatic.rbf.ci.csv")


