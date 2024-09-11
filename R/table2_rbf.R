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

source("GPlag_server.R")

# Table 2 left

## Simulation 3 -clustering/ranking
n <- 100
X <- matrix(seq(-2, 2, length=n), ncol=1)
k=0.01
y1 <- atan(k*X)/atan(k)
plot(X,y1)
k=1
y2 <- atan(k*sin(X))/atan(k)
k=0.01
y2 <- atan(k*sin(X))/atan(k)

dat <- data.frame(x= X, y1= y1, y2 = y2)
ggplot(dat, aes(x = X))+
  geom_point(aes(y=y1,col = "blue"))+
  geom_point(aes(y=y2,col = "red"))+
  geom_line(aes(y=y1,col = "blue"))+
  geom_line(aes(y=y2,col = "red"))
# geom_line(aes(y=y3,col = "black"))
delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
group.index <- rep(c(1,2),each=n)
group.size <- rle(group.index)

k_across = rep(c(0.01,1,10),3)
s_across = rep(c(0,0.5,1),each=3)

yother <- sapply(1:length(s_across), function(i){
  k = k_across[i]
  s = s_across[i]
  # y2 <- atan(k*sin(X+s))/atan(k)
  atan(k*(X+s))/atan(k)
})


dat1 <- cbind(X, yother) %>% as.data.frame() %>% `colnames<-`(c("t",paste0("ts",1:9))) %>%
  tidyr::pivot_longer(!t, names_to = "groups") %>%
  mutate(k = rep(k_across %>% as.factor(), n),
         s = rep(s_across %>% as.factor(), n))
dat1$label <- sprintf(
  "k = %s, s = %s",
  dat1$k,
  dat1$s
)
mypal = pal_npg("nrc", alpha = 0.8)(9)
mypal = pal_d3("category10")(9)

ggplot() +
  geom_line(data = dat1,aes(x = t, y = value,col = groups),linewidth = 2)+
  geom_line(data = data.frame(t = X, target = y1), aes(x = t,y = target), linewidth = 2)+
  facet_wrap(~groups) +
  scale_color_manual(name = "Groups", values = mypal) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none")+
  geom_text(
    size    = 4,
    data    = dat1,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.5,
    vjust   = 1.5
  )
# xlab("Kernel from which data is generated")+
# ylab("Fitted value of a")

Xtest = matrix(seq(-2, 2, length=200), ncol=1)

# rbf

a_across_rbf = matrix(rep(0,length(s_across)*2),nrow=2)

##########################################################################################
i=1
y2 <- yother[,i]
set.seed(1)
out <- optim(c(-5,-1,0), nl_rbf, method="L-BFGS-B",lower=c(-20,-10,-0.05),
             upper=c(20,20,0.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par))
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]
##########################################################################################

i=4
set.seed(1)
y2 <- yother[,i]
out <- optim(c(0,0,0.5), nl_rbf, method="L-BFGS-B",lower=c(-20,-10,0.45),
             upper=c(20,20,0.55),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]
##########################################################################################

i=7
set.seed(1)
y2 <- yother[,i]
out <- optim(c(-5,-10,1), nl_rbf, method="L-BFGS-B",lower=c(-20,-20,0.95),
             upper=c(20,20,1.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]
##########################################################################################

i=2
y2 <- yother[,i]
set.seed(1)
out <- optim(c(0.5,-0.5,0), nl_rbf, method="L-BFGS-B",lower=c(-20,-10,-0.05),
             upper=c(20,20,0.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par))
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]

##########################################################################################

i=5
set.seed(1)
y2 <- yother[,i]
out <- optim(c(0.5,-0.5,0.5), nl_rbf, method="L-BFGS-B",lower=c(-20,-10,0.45),
             upper=c(20,20,0.55),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]

##########################################################################################

i=8
set.seed(1)
y2 <- yother[,i]
out <- optim(c(0,0.5,1), nl_rbf, method="L-BFGS-B",lower=c(-20,-20,0.95),
             upper=c(20,20,1.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]
##########################################################################################

i=3
y2 <- yother[,i]
set.seed(1)
out <- optim(c(-3,2,0), nl_rbf, method="L-BFGS-B",lower=c(-20,-20,-0.05),
             upper=c(20,20,0.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$par)

a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]

##########################################################################################

i=6
set.seed(1)
y2 <- yother[,i]
out <- optim(c(-2,2,0.5), nl_rbf, method="L-BFGS-B",lower=c(-20,-10,0.45),
             upper=c(20,20,0.55),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]

##########################################################################################

i=9
set.seed(1)
y2 <- yother[,i]
out <- optim(c(0,4,1), nl_rbf, method="L-BFGS-B",lower=c(-20,-20,0.95),
             upper=c(20,20,1.05),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size,control=list(trace=3))
print(out$message)
print(softplus(out$par)[1:2]);print(out$par[3])
a_across_rbf[1,i] = softplus(out$par[2])
a_across_rbf[2,i] = out$par[3]
##########################################################################################
mypal = pal_d3("category10")(9)
dat <- data.frame(k= k_across %>% as.factor(), s =s_across %>% as.factor(), a = a_across_rbf[1,], shat = a_across_rbf[2,],
                  groups = paste0("ts",1:9))
fwrite(dat,'table2_rbf.csv')
p_rbf = ggplot(dat, aes(x= k,y=a, col=groups, shape=s))+
  geom_point(size = 5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = c(.05, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  scale_color_manual(name = "Groups", values = mypal, guide = "none")

# calculate ARI and accuracy for cluster performance
set.seed(2)
km <- kmeans(a_across_rbf[1,],3)

library(aricode)
ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)

