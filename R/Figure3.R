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

## Simulation 2 - evaluate a is close to simularity
n = 50
tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4

## data generated from mtsgp kernel
Y_exp = simulateData(b = 1,a = 0.3,s = c(0,2), tau2, n, kernel = "exp", timemax = n)
Y_exp.smooth = simulateData(b = 0.1,a = 0.3,s = c(0,2), tau2, n, kernel = "exp", timemax = n)
#jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer1_b1viz.jpg",width = 15, height = 10,units = "in",res=450)
plotSimulate(Y_exp,i=1,s= c(0,2), pointSize = 4, textSize = 5)
#dev.off()

#jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer1_b0.1viz.jpg",width = 15, height = 10,units = "in",res=450)
plotSimulate(Y_exp.smooth,i=4,s= c(0,2), pointSize = 4, textSize = 5)
#dev.off()

## Prediction tasks on irregualr data
set.seed(666)
train.index = sample(1:50,25)
train.index <- train.index[order(train.index)]
Y.mtsgp.train <- Y_exp[c(train.index,train.index+50),]
plotSimulate(Y.mtsgp.train,i=2,s= c(0,2))
gplag.time <- system.time({mtsgp.mtsgp.train <- deriveEstimation(Y.mtsgp.train, "exp", sl = -1, su = 5)}) #perform MTSE - return b,a,s, tau2, ll in each column
sgp.time <- system.time({mtsgp.sgp.train <- deriveEstimation(Y.mtsgp.train, "sep.exp")}) # perform SGP

# derive prediction
Y.mtsgp.test <- Y_exp[-c(train.index,train.index+50),]
t <- (rownames(Y.mtsgp.train) %>% as.numeric())
ttest = (rownames(Y.mtsgp.test) %>% as.numeric())
ntest = length(ttest)/2

# save date to compare amtsgp and gplva
write.table(t, file = "/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/train_t_wiggle.csv",row.names = F, col.names = F, sep = ",")
write.table(ttest, file = "/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/test_t_wiggle.csv",row.names = F, col.names = F, sep = ",")
write.table(Y.mtsgp.train, file = "/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/train_y_wiggle.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(Y.mtsgp.test, file = "/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/test_y_wiggle.csv",row.names = F, col.names = F, sep = ",")




# helper function to calulate loglikelihood from magma
logL_GP_mod <- function(hp, db, mean, kern, pen_diag) {
  if (length(mean) == 1) {
    mean <- rep(mean, nrow(db))
  }
  ## mean is equal for all timestamps
  
  ## Extract the input variables (reference Input + Covariates)
  inputs <- db %>% dplyr::select(-.data$Output)
  ## Compute the inverse of the covariance matrix
  inv <- kern_to_inv(inputs, kern, hp, pen_diag)
  
  ## Classical Gaussian log-likelihood
  LL_norm <- -dmnorm(db$Output, mean, inv, log = T)
  ## Correction trace term (- 1/2 * Trace(inv %*% post_cov))
  cor_term <- 0.5 * sum(inv * post_cov)
  
  return(LL_norm)
}
logL_monitoring <- function(hp_0,
                            hp_i,
                            db,
                            m_0,
                            kern_0,
                            kern_i,
                            post_mean,
                            pen_diag) {
  ## Compute the modified logL for the mean process
  ll_0 <- logL_GP_mod(
    hp = hp_0,
    db = post_mean,
    mean = m_0,
    kern = kern_0,
    pen_diag = pen_diag
  )
  
  ## Sum over the individuals
  funloop <- function(i) {
    ## Extract the i-th specific reference inputs
    input_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::pull(.data$Reference)
    ## Extract the i-th specific hyper-parameters
    hp_i_i <- hp_i %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## Extract the i-th specific Inputs and Output
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## Extract the mean values associated with the i-th specific inputs
    post_mean_i <- post_mean %>%
      dplyr::filter(.data$Reference %in% input_i) %>%
      dplyr::pull(.data$Output)
    
    ## Compute the modified logL for the individual processes
    logL_GP_mod(hp_i_i,
                db_i,
                post_mean_i,
                kern_i,
                pen_diag) %>%
      return()
  }
  sum_ll_i <- sapply(unique(db$ID), funloop) %>% sum()
  
  return(- sum_ll_i )
}

mse.mtsgp.smooth <- sapply(1:20, function(i){
  Y = Y.mtsgp.train[,i]
  params = mtsgp.mtsgp.train[i,]
  yres <- interpolation_sim(t,Y,Y.mtsgp.test[,i],params,ttest)
  yhat <- yres[1:length(ttest)]
  res.mtsgp <- mean((Y.mtsgp.test[,i] - yhat)^2)
  gplag.logL <- yres[length(ttest)+1]
  
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
  yres <- interpolation_sim(t,Y,Y.mtsgp.test[,i],params,ttest, kernel = "sep.matern")
  yhat <- yres[1:length(ttest)]
  res.sep <- mean((Y.mtsgp.test[,i] - yhat)^2)
  sgp.logL <- yres[length(ttest)+1]
  
  # fit MAGMA
  magma_train <- tibble(ID = rep(c("1","2"), each = n), Input = t, Output = Y)
  predict.index = sample(n,size = n/2)
  predict.index <- predict.index[order(predict.index)]
  magma_tr <- magma_train[c(predict.index, predict.index+n),]
  
  magma_pred1 <- magma_train[setdiff(1:(n),predict.index),]
  magma_pred2 <- magma_train[setdiff((n+1):(2*n),predict.index+n),]
  
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
  res.magma = mean((Y.mtsgp.test[,i] - c(y1hat,y2hat))^2)
  
  # modify the dataset to calculate likelihood from magma to be compariable with gplag
  magma_test <- tibble(ID = rep(c("1","2"), each = ntest), Input = ttest, Output = c(y1hat,y2hat))
  data = magma_test
  if (!("Reference" %in% (data %>% names()))) {
    names_col <- data %>%
      dplyr::select(- c(.data$ID, .data$Output)) %>%
      names()
  } else {
    names_col <- data %>%
      dplyr::select(- c(.data$ID, .data$Output, .data$Reference)) %>%
      names()
  }
  data <- data %>%
    purrr::modify_at(tidyselect::all_of(names_col), signif) %>%
    tidyr::unite("Reference",
                 tidyselect::all_of(names_col),
                 sep = ":",
                 remove = FALSE) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::arrange(.data$Reference, .by_group = TRUE) %>%
    dplyr::ungroup()
  all_inputs <- data %>%
    dplyr::select(-c(.data$ID, .data$Output)) %>%
    unique()
  all_input <- all_inputs %>% dplyr::pull(.data$Reference)
  m_0 <- rep(0, length(all_input))
  magma_logL <- logL_monitoring(
    hp_0 = model$hp_0,
    hp_i = model$hp_i,
    db = data,
    m_0 = m_0,
    kern_0 = model$ini_args$kern_0,
    kern_i = model$ini_args$kern_i,
    post_mean = data[1:n,2:4],
    pen_diag = 1e-10
  )
  return(c(res.mtsgp,res.splines,res.sep, res.magma, gplag.logL, sgp.logL, magma_logL))
})



# change outliers color in boxplot
is.outlier <- function (x) {
  x < quantile(x, .25) - 1.5 * IQR(x) |
    x > quantile(x, .75) + 1.5 * IQR(x)
}

dat3.1 <- t(mse.mtsgp[1:4,]) %>% as.data.frame() %>% `colnames<-`(c("GPlag","Splines","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  mutate(outlier.p = is.outlier(value)) %>%
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

## plot for reviewer 1

# read mse of data b=1(less smooth)
mse.amtsgp.unsmooth <- read.csv("/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/test_mse_b1.csv")
# mse.amtsgp <- rbind(mse.amtsgp,data.frame(MSE = rep(NA,7)))

dat3.1.unsmooth <- cbind(t(mse.mtsgp.unsmooth[1:4,]) %>% as.data.frame() , mse.amtsgp)%>% `colnames<-`(c("GPlag","Splines","SGP", "MAGMA","AMTGP")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
p6reviewer1.unsmooth <-ggplot(dat3.1.unsmooth,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  # geom_point(data = dat3.1[dat3.1$outlier.p,],aes(  col = model))+
  # facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(1,2,3,4,5)], labels = c("GPlag","Seperate GPs", "Splines","MAGMA","AMTGP"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(1,2,3,4,5,6)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        legend.position = "none"
  ) +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")

# read prediction of data b=0.1(smooth)
mse.amtsgp.smooth <- read.csv("/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/data/test_mse_b0.1.csv")

dat3.1.smooth <- cbind(t(mse.mtsgp.smooth[1:4,1:20]) %>% as.data.frame() , mse.amtsgp.smooth)%>% `colnames<-`(c("GPlag","Splines","SGP", "MAGMA","AMTGP")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

p6reviewer1.smooth <-ggplot(dat3.1.smooth,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  # geom_point(data = dat3.1[dat3.1$outlier.p,],aes(  col = model))+
  # facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(1,2,3,4,5)], labels = c("GPlag","Seperate GPs", "Splines","MAGMA","AMTGP"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(1,2,3,4,5,6)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 35),axis.title.x = element_text(size = 35),
        legend.position = "none"
  ) +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")
p6reviewer1.smooth


jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer1.smooth.jpg",width = 15, height = 10,units = "in",res=450)
p6reviewer1.smooth
dev.off()

dat3.1logl <- t(mse.mtsgp[5:7,]) %>% as.data.frame() %>% `colnames<-`(c("GPlag","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>% filter(value>(-130)&value<(-10)) %>%
  ungroup()

dat3.1logl %>% group_by(model) %>% summarise(mean = mean(value), se = std_err(value))
# A tibble: 3 × 3
# model  mean    se
# <chr> <dbl> <dbl>
#   1 GPlag -45.9 4.02
# 2 MAGMA -71.1 0.597
# 3 SGP   -99.8 0.763
p6supp <- ggplot(dat3.1logl ,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  scale_alpha(0.3)+
  ylim(-150,0)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(2,3,5)], labels = c("GPlag","Seperate GPs", "MAGMA"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(2,3,5)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        legend.position = "none"
  ) +
  # xlab("Model from which data is generated") +
  ylab("Log-likelihood in test set")

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/loglikelihood.jpg",width = 10, height = 10,units = "in",res=450)
p6supp
dev.off()


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
gplagtime <- system.time({
  a_across<- sapply(1:length(s_across), function(i){
    y2 <- yother[,i]
    out <- optim(c(1,1,0), nl_exp, method="L-BFGS-B",lower=c(-10,-10,-1),
                 upper=c(20,Inf,4),
                 t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
    return(c(softplus(out$par[2]),out$par[3]))
  })
})

a_across_sep<- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  out <- optim(c(1,0,0), nl_expsep, method="L-BFGS-B",lower=c(-10,-10,-1),
               upper=c(20,Inf,4),
               t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
  return(c(softplus(out$par[2]),out$par[3]))
})

dat <- data.frame(k= k_across %>% as.factor(), s =s_across %>% as.factor(), a = a_across[1,], shat = a_across[2,],
                  groups = paste0("ts",1:9))
dat <- data.frame(k= k_across %>% as.factor(), s =s_across %>% as.factor(), a = a_across_sep[1,], shat = a_across_sep[2,],
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
write.csv(X, file = "data/time.csv",row.names = F)

# calculate ARI and accuracy for cluster performance
km <- kmeans(a_across[1,],3)

library(aricode)
ari_gplag <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform TLCC
tlcctime <- system.time({
  tlcc <- sapply(1:length(s_across), function(i){
    y2 <- yother[,i]
    cc <- ccf(y1[,1],y2,plot=FALSE)
    return(max(cc$acf))
  })
})
km <- kmeans(tlcc,3)
ari_tlcc <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform magmaclust
magmaclust_train <- yother %>% as.data.frame() %>% pivot_longer(everything(), names_to = "ID", values_to = "Output") %>%
  cbind(Input = rep(X,each = 9))
magmaclusttime <- system.time({
  model_clust <- train_magmaclust(data = magmaclust_train)
})
data_train_with_clust = data_allocate_cluster(model_clust)
cl = data_train_with_clust %>% group_by(ID) %>% summarise(cluster = Cluster %>% unique())
ari_magmaclust <- ARI(c(1,2,3,1,2,3,1,2,3), cl$cluster)
# Calculate accuracy
NMI(c(1,2,3,1,2,3,1,2,3), cl$cluster)

# perform dtw
library(dtw)

dtwtime <- system.time({
  dtwdistance <- sapply(1:length(s_across), function(i){
    y2 <- yother[,i]
    dtw_result <- dtw(y1, y2)
    dtw_distance <- dtw_result$distance
    return(dtw_distance)
  })
})

print(time_taken) # This will print the time taken for execution

km <- kmeans(dtwdistance,3)
ari_dtw <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)

# perform soft-dtw
library(dtwclust)
# softdtw_result <- tsclust(t(yother), k=3, distance = "sdtw", centroid = "dba")
# plot(softdtw_result, type = "sc")
set.seed(101)
softdtwtime <- system.time({
  softdtwdistance<- sapply(1:length(s_across), function(i){
    y2 <- yother[,i]
    softdtw_distance <- sdtw(y1, y2)
    return(softdtw_distance)
  })
})
# softdtwdivdistance<- read.csv(file = "data/softdtwdiv.csv")
km <- kmeans(softdtwdistance,3)
ari_softdtwdiv <- ARI(c(1,2,3,1,2,3,1,2,3), km$cluster)
ari_softdtwdiv
NMI(c(1,2,3,1,2,3,1,2,3), km$cluster)


