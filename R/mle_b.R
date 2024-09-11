library(plgp)
library(mvtnorm)
library(dplyr)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(kernlab)
library(MagmaClustR)
library(splines)
library(wesanderson)

library(ggplot2)
library(cowplot)
library(egg)
library(data.table)
library(ggsci)

source("R/GPlag.R")

# For Biometrics reviewer 2
### larger range of b to show its impact on a and s estimation

timestamp()
# n is number of timepoints
for(n in c(50,100,200,300)){
  set.seed(781)
  a = 1
  s <- c(0,2) # simulate 2 time series with shift 0 and 2
  b = 10
  Y_exp100_b10 = simulateData(b,a,s,tau2, n, kernel = "exp", n_replicate = 20)
  b = 0.3
  Y_exp100_b03 = simulateData(b,a,s,tau2, n, kernel = "exp", n_replicate = 20)
  b = 5
  Y_exp100_b5 = simulateData(b,a,s,tau2, n, kernel = "exp", n_replicate = 20)
  res100_b <- list()
  res100_b[["b0.3"]] <- deriveEstimation(Y_exp100_b03, kernel ="exp",b0 = 0.3, s0=2.0)
  res100_b[["b5"]] <- deriveEstimation(Y_exp100_b5, kernel ="exp",b0 = 3, s0=2.0,bu=10)
  res100_b[["b10"]] <- deriveEstimation(Y_exp100_b10, kernel ="exp", b0 = 10,s0=2.0)
  res_filter100_b <- lapply(1:3,function(i){
    res100_b[[i]][res100_b[[i]][,1] != softplus(-10) & res100_b[[i]][,1] != softplus(20) & res100_b[[i]][,2] < 10,]
  })
  df100_b <- do.call("rbind", res_filter100_b) %>% as.data.frame() %>%
    `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
    mutate(n = 100,
           kernel = rep(c("0.3","10","5"), unname(unlist(lapply(res_filter100_b,nrow)))))
  df100_b %>% group_by(kernel) %>%
    mutate(outlier.p = is.outlier(ahat),outlier.s = is.outlier(shat)) %>%
    ungroup() -> df100_b.out
  df100_b.out$kernel <- factor(df100_b.out$kernel, levels = c(0.3, 5, 10))
  
  fwrite(df100_b.out,paste0('n',n,'.txt'))
}

n50 = fread('n50.txt')
n100 = fread('n100.txt')
n200 = fread('n200.txt')
n300 = fread('n300.txt')

result_all = rbind(n50,n100,n200,n300)

result_all$n <- factor(result_all$n, levels = c(50,100,200,300))

result_all$kernel = paste0('b=',result_all$kernel)
result_all$kernel <- factor(result_all$kernel, levels = c('b=0.3','b=5','b=10'))



## plot various b influence on a and s
k1 = 
  ggplot(result_all,aes(y = ahat, x= n, fill = n))+
  facet_wrap(~kernel,nrow=1)+
  geom_boxplot(outlier.size = 0.5)+
  #geom_point(data = df100_b.out, aes(col = kernel))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Moonrise3"))+
  scale_colour_manual(values = wes_palette("Moonrise3"))+
  #scale_fill_bmj()+
  #scale_color_bmj()+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Sample size")+
  ylab("Fitted value of a")

k2 = 
  ggplot(result_all,aes(y = shat, x= n, fill = n))+
  facet_wrap(~kernel)+
  geom_boxplot(outlier.size = 0.5)+
  #geom_point(data = df100_b.out, aes(col = kernel))+
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Moonrise3"))+
  scale_colour_manual(values = wes_palette("Moonrise3"))+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Sample size")+
  ylab("Fitted value of s") +
  coord_cartesian(ylim=c(1.3,2.7))
#ylim(-1,3)

k=plot_grid(k1,k2)
ggsave('MLE_b.png',k,bg='white',height = 3,width = 10)
