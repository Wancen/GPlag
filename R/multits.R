library(plgp)
library(mvtnorm)
library(tidyverse)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(ggplot2)
library(ggsci)
library(patchwork)
library(wesanderson)
library(splines)
library(MagmaClustR)
source("R/GPlag.R")


n = 50
tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4
b = 0.3
a = c(1,1,1)
s <- c(0,2,4) # simulate 2 time series with shift 0 and 2

## data generated from mtsgp kernel
Ymultiexp = simulateMultiData(b = b,a = a,s = s, tau2, n, kernel = "exp", nrep = 30)
Ymultirbf = simulateMultiData(b,a,s,tau2, n, kernel = "rbf", nrep = 30)
resMultiexp <- deriveMultiEstimation(Ymultiexp, kernel = "exp")
resMultirbf <- deriveMultiEstimation(Ymultirbf, kernel = "rbf")
colnames(resMultirbf) <- c("b","a12","a13","a23","s2","s3","tau2","loglik")
datMultirbf <- resMultirbf %>% as.data.frame() %>%
  # filter(kernel == "exp") %>%
  dplyr::select(a12:s3) %>%
  tidyr::pivot_longer(a12:s3,names_to ="paramater")%>% group_by(paramater) %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

colnames(resMultiexp) <- c("b","a12","a13","a23","s2","s3","tau2","loglik")
datMultiexp <- resMultiexp %>% as.data.frame() %>%
  # filter(kernel == "exp") %>%
  dplyr::select(a12:s2) %>%
  tidyr::pivot_longer(a12:s2,names_to ="paramater")%>% group_by(paramater) %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

pltMultiparam <- ggplot(datMultiexp,aes(y = value, x= paramater, fill = paramater))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = dat_exp[dat_exp$outlier.p,], aes(fill = n, col = n))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 4, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  scale_fill_startrek() +
  scale_color_startrek() +
  scale_x_discrete(labels = c("a12","a13","a23","s1","s2"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey10"),
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("Parameter")+
  ylab("Fitted value")+
  ylim(-1,5)

pltMultiparam <-ggplot(datMultirbf,aes(y = value, x= paramater, fill = paramater))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = dat_exp[dat_exp$outlier.p,], aes(fill = n, col = n))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 4, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  scale_fill_startrek() +
  scale_color_startrek() +
  scale_x_discrete(labels = c("a12","a13","a23","s2","s3"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey10"),
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("Parameter")+
  ylab("Fitted value")+
  ylim(-1,5)


set.seed(666)
train.index = sample(1:n,n/2)
train.index <- train.index[order(train.index)]
Ymultiexp.train <- Ymultiexp[c(train.index,train.index+n,train.index+2*n),]
plotMultiSimulate(Ymultiexp.train,i=2,s= c(0,2,4))
gplag.train <- deriveMultiEstimation(Ymultiexp.train, "exp") #perform MTSE - return b,a,s, tau2, ll in each column
gplag.sgp.train <- deriveMultiEstimation(Ymultiexp.train, "sep.exp") # perform SGP

# derive prediction
Ymultiexp.test <- Ymultiexp[-c(train.index,train.index+n,train.index+2*n),]
t <- (rownames(Ymultiexp.train) %>% as.numeric())
ttest = (rownames(Ymultiexp.test) %>% as.numeric())
ntest = length(ttest)/3



mse.gplag <- sapply(1:30, function(i){
  Y = Ymultiexp.train[,i]
  params = gplag.train[i,]
  yres <- interpolation_Multisim(t,Y,Ymultiexp.test[,i],params,ttest, kernel = "lexp")
  yhat <- yres[1:length(ttest)]
  # yhat_capped <- pmin(pmax(yhat, min(Ymultiexp.test[,i])), max(Ymultiexp.test[,i]))
  res.mtsgp <- mean((Ymultiexp.test[,i] - yhat)^2)
  gplag.logL <- yres[length(ttest)+1]

  # fit splines
  n = length(Y)/3
  t1 <- t[1:n]
  t2 <- t[(n+1):(2*n)]
  t3 <- t[(2*n+1):(3*n)]
  fit = lm(Y[1:n] ~ ns(t1, knots=NULL))
  ntest = length(ttest)/3
  y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
  fit2 = lm(Y[(n+1):(2*n)] ~ ns(t2, knots=NULL))
  y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
  fit3 = lm(Y[(2*n+1):(3*n)] ~ ns(t3, knots=NULL))
  y3hat = predict(fit3,data.frame(t3=ttest[(2*ntest+1):(3*ntest)]))
  res.splines = mean((Ymultiexp.test[,i] - c(y1hat,y2hat,y3hat))^2)

  # fit sep GP
  sigma <- 1 # Scale parameter
  kernel <- besseldot(sigma = sigma, order = 0, degree = 1)
  gp_model <- gausspr(as.matrix(t1), Y[1:n], kernel = kernel)
  predictions1 <- predict(gp_model, as.matrix(ttest[1:ntest]))
  gp_model <- gausspr(as.matrix(t2), Y[(n+1):(2*n)], kernel = kernel)
  predictions2 <- predict(gp_model, as.matrix(ttest[(ntest+1):(2*ntest)]))
  gp_model <- gausspr(as.matrix(t3), Y[(2*n+1):(3*n)], kernel = kernel)
  predictions3 <- predict(gp_model, as.matrix(ttest[(2*ntest+1):(3*ntest)]))

  res.sgp <- mean((Ymultiexp.test[,i] - c(predictions1, predictions2, predictions3))^2)

  # fit MAGMA
  magma_train <- tibble(ID = rep(c("1","2","3"), each = n), Input = t, Output = Y)
  predict.index = sample(n,size = n/2)
  predict.index <- predict.index[order(predict.index)]
  magma_tr <- magma_train[c(predict.index, predict.index+n),]

  magma_pred1 <- magma_train[setdiff(1:(n),predict.index),]
  magma_pred2 <- magma_train[setdiff((n+1):(2*n),predict.index+n),]
  magma_pred3 <- magma_train[setdiff((2*n+1):(3*n),predict.index+2*n),]

  model <- train_magma(data = magma_tr, common_hp = T)

  pred1  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[1:n])
  pred1 <- pred1[order(pred1$Input),]
  y1hat = pred1$Mean
  pred2  <- pred_magma(data = magma_pred2,
                       trained_model = model,
                       grid_inputs = ttest[(n+1):(2*n)])
  pred2 <- pred2[order(pred2$Input),]
  y2hat = pred2$Mean
  pred3  <- pred_magma(data = magma_pred3,
                       trained_model = model,
                       grid_inputs = ttest[(2*n+1):(3*n)])
  pred3 <- pred3[order(pred3$Input),]
  y3hat = pred3$Mean
  res.magma = mean((Ymultiexp.test[,i] - c(y1hat,y2hat,y3hat))^2)

  return(c(res.mtsgp,res.splines,res.sgp, res.magma))
})

datMulti <- t(mse.gplag) %>% as.data.frame() %>%  `colnames<-`(c("GPlag","Splines","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

pltMulti <-ggplot(datMulti,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  # geom_point(data = dat3.1[dat3.1$outlier.p,],aes(  col = model))+
  # facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(2,3,4,5)], labels = c("GPlag","Seperate GPs", "Splines","MAGMA"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(2,3,4,5)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        legend.position = "none"
  ) +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/Rev4.jpg",width = 15, height = 8,units = "in",res=450)
pltMultiparam + pltMulti + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 18),
        plot.tag = element_text(size = 14)
  )
dev.off()
