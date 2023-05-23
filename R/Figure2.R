library(plgp)
library(mvtnorm)
library(dplyr)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(ggplot2)
library(ggsci)
library(patchwork)
library(wesanderson)
source("R/GPlag.R")

## Simulation 2 - evaluate a is close to simularity
n = 50
tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4

## data generated from mtsgp kernel
Y_exp = simulateData(b = 1,a = 0.3,s = c(0,2), tau2, n, kernel = "exp", timemax = n)
plotSimulate(Y_exp,i=1,s= c(0,2))

## Prediction tasks on irregualr data
set.seed(666)
train.index = sample(1:50,25)
Y.mtsgp.train <- Y_exp[c(train.index,train.index+50),]
mtsgp.mtsgp.train <- deriveEstimation(Y.mtsgp.train, "exp", sl = -1, su = 5) #perform MTSE - return b,a,s, tau2, ll in each column
mtsgp.sgp.train <- deriveEstimation(Y.mtsgp.train, "sep.exp") # perform SGP

# derive prediction
Y.mtsgp.test <- Y_exp[-c(train.index,train.index+50),]
t <- (rownames(Y.mtsgp.train) %>% as.numeric())
ttest = (rownames(Y.mtsgp.test) %>% as.numeric())
ntest = length(ttest)/2

mse.mtsgp <- sapply(1:100, function(i){
  Y = Y.mtsgp.train[,i]
  params = mtsgp.mtsgp.train[i,]
  yhat <- interpolation_sim(t,Y,params,ttest)
  res.mtsgp <- mean((Y.mtsgp.test[,i] - yhat)^2)

  # fit splines
  n = length(Y)/2
  t1 <- t[1:n]
  t2 <- t[(n+1):(2*n)]
  fit = lm(Y[1:n] ~ ns(t1, knots=NULL))
  ntest = length(ttest)/2
  y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
  fit2 = lm(Y[(n+1):(2*n)] ~ ns(t2, knots=NULL))
  y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
  res.splines = mean((Y.mtsgp.test[,i] - c(y1hat,y2hat))^2)

  # fit sep GP
  params = mtsgp.sgp.train[i,]
  yhat <- interpolation_sim(t,Y,params,ttest, kernel = "sep.matern")
  res.sep <- mean((Y.mtsgp.test[,i] - yhat)^2)

  # fit MAGMA
  magma_train <- tibble(ID = rep(c("1","2"), each = n), Input = t, Output = Y)
  magma_pred1 <- tibble(ID = rep("1", n), Input = ttest[1:n], Output = unname(Y.mtsgp.test[1:n,i]))
  magma_pred2 <- tibble(ID = rep("2", n), Input = ttest[(n+1):(2*n)], Output = unname(Y.mtsgp.test[(n+1):(2*n),i]))
  model <- train_magma(data = magma_train, common_hp = T)
  pred1  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[1:n])
  pred1 <- pred1[order(pred1$Input),]
  y1hat = pred1$Mean
  pred2  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[(n+1):(2*n)])
  pred2 <- pred2[order(pred2$Input),]
  y2hat = pred2$Mean
  res.magma = mean((Y.mtsgp.test[,i] - c(y1hat,y2hat))^2)

  return(c(res.mtsgp,res.splines,res.sep, res.magma))
})

dat3.1 <- t(mse.mtsgp) %>% as.data.frame() %>% `colnames<-`(c("GPlag","Splines","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>% mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

p6 <-ggplot(dat3.1,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  # geom_point(data = dat3.1[dat3.1$outlier.p,],aes(  col = model))+
  # facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(2,3,4,5)], labels = c("GPlag","Seperate GPs", "Splines","MAGMA"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(2,3,4,5)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        legend.position = "none"
  ) +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")
p6

## Simulation 3 -clustering/ranking
n <- 50
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
out <- optim(c(1,1), nl_mggp, method="L-BFGS-B",lower=c(-10,-10),
             upper=c(20,Inf),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
bhat = softplus(out$par[1])
ahat = softplus(out$par[2])

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
p1 <- ggplot() +
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
    hjust   = 1.0,
    vjust   = 1.5
  )
# xlab("Kernel from which data is generated")+
# ylab("Fitted value of a")

Xtest = matrix(seq(-2, 2, length=200), ncol=1)
a_across<- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  out <- optim(c(1,1,0), nl_exp, method="L-BFGS-B",lower=c(-10,-10,-1),
               upper=c(20,Inf,4),
               t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
  return(c(softplus(out$par[2]),out$par[3]))
})

dat <- data.frame(k= k_across %>% as.factor(), s =s_across %>% as.factor(), a = a_across[1,], shat = a_across[2,],
                  groups = paste0("ts",1:9))

p2 <- ggplot(dat, aes(x= k,y=a, col=groups, shape=s))+
  geom_point(size = 5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = c(.05, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  scale_color_manual(name = "Groups", values = mypal, guide = "none")
p1+p2

jpeg(file="plots/simulation3.jpg",width = 15, height = 5.5,units = "in",res=450)
p6+ p1 + p2 + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 18),
        plot.tag = element_text(size = 14)
  )
dev.off()

write.csv(yother, file = "data/cluster.csv",row.names = F)
write.csv(y1, file = "data/target.csv",row.names = F)

# calculate ARI and accuracy for cluster performance
km <- kmeans(a_across[1,],3)

library(aricode)
ari_gplag <- adjustedRandIndex(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform TLCC
tlcc <- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  cc <- ccf(y1[,1],y2,plot=FALSE)
  return(max(cc$acf))
})
km <- kmeans(tlcc,3)
ari_tlcc <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform magmaclust
magmaclust_train <- yother %>% as.data.frame() %>% pivot_longer(everything(), names_to = "ID", values_to = "Output") %>%
  cbind(Input = rep(X,each = 9))
model_clust <- train_magmaclust(data = magmaclust_train)
data_train_with_clust = data_allocate_cluster(model_clust)
cl = data_train_with_clust %>% group_by(ID) %>% summarise(cluster = Cluster %>% unique())
ari_magmaclust <- ARI(c(1,2,3,1,2,3,1,2,3), cl$cluster)
# Calculate accuracy
NMI(c(1,2,3,1,2,3,1,2,3), cl$cluster)

# perform dtw
library(dtw)
dtwdistance<- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  dtw_result <- dtw(y1, y2)
  dtw_distance <- dtw_result$distance
  return(dtw_distance)
})

km <- kmeans(dtwdistance,3)
ari_dtw <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform soft-dtw
library(dtwclust)
# softdtw_result <- tsclust(t(yother), k=3, distance = "sdtw", centroid = "dba")
# plot(softdtw_result, type = "sc")
set.seed(101)
softdtwdistance<- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  softdtw_distance <- sdtw(y1, y2)
  return(softdtw_distance)
})

softdtwdivdistance<- read.csv(file = "data/softdtwdiv.csv")
km <- kmeans(softdtwdivdistance$softdtwdiv,3)
ari_softdtwdiv <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)


