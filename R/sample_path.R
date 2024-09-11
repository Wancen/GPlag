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

is.outlier <- function (x) {
  x < quantile(x, .25) - 1.5 * IQR(x) |
    x > quantile(x, .75) + 1.5 * IQR(x)
}



### sample path

plot_samplepath <- function(b,a,s,tau2=4,n=100,kernel,timemax=100){
  Y = simulateData(b = b,a = a,s = s, tau2, n, kernel, timemax)
  n = nrow(Y)/2
  cumsum.index <- cumsum(c(n,n))
  y1 = Y[1:cumsum.index[1],i]
  y2 = Y[(cumsum.index[1]+1):(cumsum.index[2]),i]
  t1 = rownames(Y)[1:n] %>% as.numeric()
  kernel_nice = ''
  if(kernel=='exp'){kernel_nice='LExp';kernel_color=wes_palette("Rushmore1")[1]}
  else{
    if(kernel=='rbf'){kernel_nice='LRBF';kernel_color=wes_palette("Rushmore1")[2]}
    else{
      if(kernel=='matern'){kernel_nice='LMat3/2';kernel_color=wes_palette("Rushmore1")[3]}
      else{print('not implemented')}}
    }
  title = paste0(kernel_nice,' b=',b," a=",a, ' s=', s[2])
  ggplot()+geom_point(aes(x=t1+s[1],y=y1),col="gray")+geom_line(aes(x=t1+s[1],y=y1),col="gray",lty=1)+
    geom_point(aes(x=t1+s[2],y=y2),col="black")+geom_line(aes(x=t1+s[2],y=y2),col="black",lty=2)+
    ylim(-10,10)+theme_article()+ggtitle(title)+xlab('t')+ylab('y')+
    theme(plot.title = element_text(color = 'black'))
  #matplot(t1+s[1], Y[1:cumsum.index[1],i], pch = 16, cex = pointSize, col="gray", lty=1, xlab="x", ylab="y", ylim = c(-10,10), cex.axis=textSize, cex.lab=textSize)
  #points(t1+s[2], Y[(cumsum.index[1]+1):(cumsum.index[2]),i], lty=2, col = "black", cex = pointSize)
}

k=plot_grid(
plot_samplepath(b = 0.1,a = 0.1,s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 1,a = 0.1,  s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 10,a = 0.1, s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 0.1,a = 1,  s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 1. ,a = 1,  s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 10 ,a = 1,  s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 0.1,a = 10, s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 1. ,a = 10, s = c(0,5),kernel = "rbf"),
plot_samplepath(b = 10 ,a = 10, s = c(0,5),kernel = "rbf"),ncol=3)

ggsave('rbf_samplepath.png',k,bg='white',height = 15,width=15)

k=plot_grid(
plot_samplepath(b = 0.1,a = 0.1,s = c(0,5),kernel = "exp"),
plot_samplepath(b = 1,a = 0.1,  s = c(0,5),kernel = "exp"),
plot_samplepath(b = 10,a = 0.1, s = c(0,5),kernel = "exp"),
plot_samplepath(b = 0.1,a = 1,  s = c(0,5),kernel = "exp"),
plot_samplepath(b = 1. ,a = 1,  s = c(0,5),kernel = "exp"),
plot_samplepath(b = 10 ,a = 1,  s = c(0,5),kernel = "exp"),
plot_samplepath(b = 0.1,a = 10, s = c(0,5),kernel = "exp"),
plot_samplepath(b = 1. ,a = 10, s = c(0,5),kernel = "exp"),
plot_samplepath(b = 10 ,a = 10, s = c(0,5),kernel = "exp"),ncol=3)

ggsave('exp_samplepath.png',k,bg='white',height = 15,width=15)

k= plot_grid(
plot_samplepath(b = 0.1,a = 0.1,s = c(0,5),kernel = "matern"),
plot_samplepath(b = 1,a = 0.1,  s = c(0,5),kernel = "matern"),
plot_samplepath(b = 10,a = 0.1, s = c(0,5),kernel = "matern"),
plot_samplepath(b = 0.1,a = 1,  s = c(0,5),kernel = "matern"),
plot_samplepath(b = 1. ,a = 1,  s = c(0,5),kernel = "matern"),
plot_samplepath(b = 10 ,a = 1,  s = c(0,5),kernel = "matern"),
plot_samplepath(b = 0.1,a = 10, s = c(0,5),kernel = "matern"),
plot_samplepath(b = 1. ,a = 10, s = c(0,5),kernel = "matern"),
plot_samplepath(b = 10 ,a = 10, s = c(0,5),kernel = "matern"),ncol=3)

ggsave('matern1.5_samplepath.png',k,bg='white',height = 15,width=15)

