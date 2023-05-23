library(plgp)
library(mvtnorm)
library(tidyverse)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(pbapply)
library(GenomicRanges)
library(nullranges)
res2<-read.csv(file="lima/data/LIMA_binnedEPCounts_nonstatic.rbf.estimated.param.csv")
ci <- read.csv(file="lima/data/LIMA_binnedEPCounts_nonstatic.rbf.ci.csv")
load(file = "lima/data/LIMA_binnedEPCounts_nonstatic.rda")

t=c(0,0.5,1,1.5,2,4,6)
n=length(t)
tinput <- c(t,t) %>% as.matrix()
peak<-epCounts_filter %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(anchor1.K27_m0000_VST:anchor1.K27_m0360_VST) %>% t() %>% `colnames<-`(paste0(seqnames(anchors(epCounts_filter, type="first"))," ",start(anchors(epCounts_filter, type="first")),"-", end(anchors(epCounts_filter, type="first"))))

gene<-epCounts_filter %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(anchor2.RNA_m0000_VST:anchor2.RNA_m0360_VST) %>% t() %>% `colnames<-`(epCounts_filter$anchor2.RNA_hgnc)

dat_filter <- rbind(gene,peak) %>% `colnames<-`(paste(colnames(peak),colnames(gene),sep = "-"))

meanGene = colMeans(dat_filter[1:n,])
meanGeneMat <- matrix(rep(meanGene, each = 2 * n), nrow = 2 * n)
Y.shift <- dat_filter - meanGeneMat
write.csv(Y.shift, file = "lima/data/Yshift.csv", row.names = FALSE)

### calucalte ccf
negative_index = floor(10*log10(n/2))+1
shiftcorr <- pbsapply(1:ncol(Y.shift), function(i){
  cc <- ccf(Y.shift[1:n,i],Y.shift[(n+1):(2*n),i],plot=FALSE)
  cor = max(cc$acf[negative_index:(2*negative_index-2)])
  s = t[which.max(cc$acf[negative_index:(2*negative_index-2)])]
  return(c(cor,s))
},cl=18)

### calucalte dtw
library(dtw)
dtwdistance<- pbsapply(1:ncol(Y.shift), function(i){
  dtw_result <- dtw(Y.shift[1:n,i],Y.shift[(n+1):(2*n),i])
  dtw_distance <- dtw_result$distance
  aligned_ts1 <- dtw_result$index1
  aligned_ts2 <- dtw_result$index2
  correlation <- cor(aligned_ts1, aligned_ts2)
  return(c(dtw_distance, correlation))
},cl=18)

### calucalte softdtw
library(dtwclust)
softdtwdistance<- pbsapply(1:ncol(Y.shift), function(i){
  softdtw_distance <- sdtw(Y.shift[1:n,i],Y.shift[(n+1):(2*n),i])
  return(softdtw_distance)
},cl=18)

### calucalte softdtwdiv in Python
softdtwdiv <- read.csv(file = "lima/data/epcount_softdtwdiv.csv")

## validate with HiC data
epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic <- 1 + (epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic - epCounts_filter$expectedCounts) / (epCounts_filter$expectedCounts + 1)
hist(epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic)


## Use observed / expected
hic = mcols(epCounts_filter)[,183:189] %>% as.matrix()

## filter not converged pairs; when a is small, all cor variance come from noise
# when a is too small or noise is large, cor will either =1 or cor not fall in corrlow and corrhigh
index1 <- (ci[2,]<= ci[1,] & ci[3,] >= ci[1,] & res2[5,] <1e-2) %>% as.vector()
# index1 <-  (res2[5,] <1e-4 & res2[2,] > 1e-4)  %>% as.vector()
table(index1)
plotPair(Y.shift,1,t, s=c(0,res2[3,i]), mu=res2[4,i])

## second filter to high dynamic change
peak_max = apply(peak,2,max)
peak_min = apply(peak,2,min)
peak_diff = peak_max-peak_min
peak_index = peak_diff > 1
gene_max = apply(gene,2,max)
gene_min = apply(gene,2,min)
gene_diff = gene_max-gene_min
gene_index = gene_diff > 1

hic_index = epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic>100

filter2 = peak_index & gene_index & index1 & hic_index
# filter2 = peak_index & gene_index
epCounts_filter2 = epCounts_filter[filter2] # 4776

Y.shift2 <- Y.shift[,filter2]

# summarize results from all methods
data = data.frame(
  corr = ci[1,]  |>  unlist(),
  corrhigh = ci[3,]  |>  unlist(),
  corrlow = ci[2,]  |>  unlist(),
  hic = epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic,
  hic_max = epCounts_filter$hic_max,
  obs_hic = epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic,
  obs_hic_max = rowMaxs(hic),
  pvalue_adj = epCounts_filter$HIC_adjusted_pval,
  loop = !is.na(epCounts_filter$HIC_value),
  similarity = res2[2,] %>% unlist(),
  noise = res2[5,] %>% unlist(),
  bhat = res2[1,] %>% unlist(),
  binit = theta,
  time_shift = res2[3,] %>% unlist(),
  spatial_variance = res2[6,] %>% unlist(),
  y_shift = res2[4,] %>% unlist(),
  ccf = shiftcorr[1,],
  dtw = dtwdistance[2,],
  sdtw = softdtwdistance,
  sdtwdiv = softdtwdiv[,1]
)
# rownames(data) <- colnames(dat_filter)
# filter pairs
data2 <- data[filter2,]

# create candidate pairs group, each methods has its own threshold
epCounts_filter2$candidate1 = ifelse(ci[2,filter2]>0.8, TRUE, FALSE) #20%
threshold1 <- quantile(data2$similarity, 0.2)
threshold2 <- quantile(data2$ccf, 0.8)
threshold3 <- quantile(data2$dtw, 0.8)
threshold4 <- quantile(data2$sdtw, 0.2)
threshold5 <- quantile(data2$sdtwdiv, 0.2)
# Check if each value is in the smallest 15% or not
epCounts_filter2$candidate1 <- data2$similarity < threshold1
epCounts_filter2$candidate2 <- data2$ccf > threshold2
epCounts_filter2$candidate3 <- data2$dtw > threshold3
epCounts_filter2$candidate4 <- data2$sdtw < threshold4
epCounts_filter2$candidate5 <- data2$sdtwdiv < threshold5
epCounts_filter2$loop = data2$loop
epCounts_filter2$corr = data2$corrlow
epCounts_filter2$similarity = data2$similarity
set.seed(123)
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate1],
                   pool = epCounts_filter2[!epCounts_filter2$candidate1],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate2],
                   pool = epCounts_filter2[!epCounts_filter2$candidate2],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate3],
                   pool = epCounts_filter2[!epCounts_filter2$candidate3],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate4],
                   pool = epCounts_filter2[!epCounts_filter2$candidate4],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate5],
                   pool = epCounts_filter2[!epCounts_filter2$candidate5],
                   covar = ~ distance,
                   method = 'stratified')
# wilcox test for candidate pairs vs other pairs
wilcox.test(focal(mgi)$obs_exp_CMB_omegaMap_inter_withNorms.hic,mgi$obs_exp_CMB_omegaMap_inter_withNorms.hic)


## Enrichment with loop # no identify
# epCounts_filter2$candidate2 = ifelse(res2[2,filter2] %>% as.numeric()<1,TRUE, FALSE)

plotPair(Y.shift2,2,t, s=c(0,(res2[3,filter2])[2] %>% as.numeric()), mu=(res2[4,filter2])[2] %>% as.numeric())

all <- c(focal(mgi), mgi)
table <- table(all$candidate1, all$loop) #
table <- table(all$candidate2, all$loop)
table <- table(all$candidate3, all$loop)
table <- table(all$candidate4, all$loop)
table <- table(all$candidate5, all$loop)
## p-value = 0.003567 when cor > 0.9 when filter non-converged pairs;
chisq.test(table)
fisher.test(table)


jpeg(file="Top-similarity-pairs.jpg",width = 15, height = 3,units = "in",res=450)
par(mfrow = c(1, 4))
data3 = data2[order(data2$similarity),]
pair.names <- rownames(data3[which(data3$similarity < 0.1 & data3$time_shift != 2 & data3$y_shift != -7),])
# pair.names2 <- str_replace_all(pair.names,"\\.","-")
for(j in c(2,3,5,6)){
  i <- which(rownames(data) == pair.names[j])
  bound <-quantile(c(Y.shift[,i],Y.shift[(n+1):(2*n),i]-res2[4,i]),probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1,xlab="t(h)", ylab="Expression", ylim = bound,lwd=3)

  title(main=colnames(Y.shift)[i],
        sub = paste0("GPlag: a =", round(res2[2,i],3), ";  ", "s = ", round(res2[3,i],3), ";  ", " TLCC: s = ", round(shiftcorr[2,i],3), "; "
                      # "hic = ", round(epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic[i],3), "obs hic = ", round(epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic[i],3)
        ), cex.main = 2, cex.sub = 1.5)

  points(t,Y.shift[1:n,i], pch=20, cex=2)

  lines(t,Y.shift[(n+1):(2*n),i],lty=1,col="blue",lwd = 3)
  points(t,Y.shift[(n+1):(2*n),i], pch=17, cex=2,col="blue")
  # lines(t,hic[,i]-mean(hic[,i]),lty=2,col="green")
  # points(t,hic[,i]-mean(hic[,i]), pch=20, cex=1,col="green")
  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red",lwd=3)
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=21, cex=2,col="red")
  lines(t+shiftprior[i],Y.shift[(n+1):(2*n),i],lty=2,col="orange3",lwd=3)
  points(t+shiftprior[i],Y.shift[(n+1):(2*n),i], pch=24, cex=2,col="orange3")

  legend("bottomright", legend=c("Gene", "Enhancer","GPlag-shifted Enhancer", "TLCC-shifted Enhancer"),col=c("black","blue","red","orange3"), lty=c(1,1,2,2),pch=c(20,17,21,24), cex=0.8,text.font=2)

}
dev.off()


