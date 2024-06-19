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

tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4
b = 0.3
a = 1
s <- c(0,2) # simulate 2 time series with shift 0 and 2
kernel = c("LRBF","LExp","LMat3/2")
## sample size 20 in Fig.1c
n = 20
set.seed(671)
Y_rbf20 = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_exp20 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern20 = simulateData(b,a,s,tau2, n, kernel = "matern")
res20 <- list()
res20[["rbf_rbf"]] <- deriveEstimation(Y_rbf20, kernel = "rbf")
res20[["exp_exp"]] <- deriveEstimation(Y_exp20, kernel ="exp")
res20[["matern_matern"]] <- deriveEstimation(Y_matern20, kernel = "matern")
## filter replicates where parameter estimation hit boundary
res_filter20 <- lapply(1:3,function(i){
  res20[[i]][res20[[i]][,1] != softplus(-10) & res20[[i]][,2] < 10 & res20[[i]][,2] != softplus(-10),]
})

## sample size 50 in Fig.1c
n = 50
set.seed(2024)
Y_rbf50 = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_exp50 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern50 = simulateData(b,a,s,tau2, n, kernel = "matern")
res50 <- list()
res50[["rbf_rbf"]] <- deriveEstimation(Y_rbf50, kernel = "rbf")
res50[["exp_exp"]] <- deriveEstimation(Y_exp50, kernel = "exp")
res50[["matern_matern"]] <- deriveEstimation(Y_matern50, kernel = "matern")
res_filter50 <- lapply(1:3,function(i){
  res50[[i]][res50[[i]][,1] != softplus(-10) & res50[[i]][,2] < 10,]
})

## sample size 100 in Fig.1 all panels
n = 100
set.seed(2024)
Y_rbf100 = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_exp100 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern100 = simulateData(b,a,s,tau2, n, kernel = "matern")
res100 <- list()
res100[["rbf_rbf"]] <- deriveEstimation(Y_rbf100, kernel = "rbf")
res100[["exp_exp"]] <- deriveEstimation(Y_exp100, kernel ="exp")
res100[["matern_matern"]] <- deriveEstimation(Y_matern100, kernel = "matern")
res_filter100 <- lapply(1:3,function(i){
  res100[[i]][res100[[i]][,1] != softplus(-10) & res100[[i]][,2] < 10,]
})


## sample size 200 in Fig.1c
n = 200
Y_exp200 = simulateData(b,a,s,tau2, n, kernel = "exp")
res200 <- list()
res200[["exp_exp"]] <- deriveEstimation(Y_exp200, "exp")
res_filter200 <- lapply(1:1,function(i){
  res200[[i]][res200[[i]][,1] != softplus(-10) & res200[[i]][,2] < 10,]
})

save(res20,res50,res100,res200, res_filter20,res_filter50,res_filter100,res_filter200, file = "data/simulation.rda")
load("data/simulation.rda")

# concatenate various sample size data together
df20 <- do.call("rbind", res_filter20) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 20,
         kernel = rep(kernel, unname(unlist(lapply(res_filter20,nrow)))))

df50 <- do.call("rbind", res_filter50) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 50,
         kernel = rep(kernel, unname(unlist(lapply(res_filter50,nrow)))))

df100 <- do.call("rbind", res_filter100) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 100,
         kernel = rep(kernel, unname(unlist(lapply(res_filter100,nrow)))))

df200 <- do.call("rbind", res_filter200) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 200,
         kernel = rep(kernel[2], unname(unlist(lapply(res_filter200,nrow)))))

dat <- rbind(df20,df50, df100,df200)
dat$n <- as.factor(dat$n)

# change outliers color in boxplot
is.outlier <- function (x) {
  x < quantile(x, .25) - 1.5 * IQR(x) |
    x > quantile(x, .75) + 1.5 * IQR(x)
}
df100 %>% group_by(kernel) %>%
  mutate(outlier.p = is.outlier(ahat),outlier.s = is.outlier(shat)) %>%
  ungroup() -> df100.out

## plot Figure 1a
p1 = ggplot(df100.out,aes(y = ahat, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df100.out[df100.out$outlier.p,], aes(col = kernel))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("GPlag kernel (n=100)")+
  ylab("Fitted value of a") +
  ylim(-1,3)
p1
# scale_fill_simpsons()

# read in Lead-Lag result and plot Figure 1b
leadlag.rbf = read.csv("data/lead-lag-rbf.csv")
leadlag.matern = read.csv("data/lead-lag-matern.csv")
leadlag.exp = read.csv("data/lead-lag-exp.csv")
df100.out$model <- "GPlag"
df100.out %>% dplyr::select(shat, kernel,model) %>%
  rbind(data.frame(shat = c(-leadlag.rbf$lead.lag,-leadlag.exp$lead.lag,-leadlag.matern$lead.lag),
                   kernel = rep(c("LRBF","LExp","LMat3/2"),each=100),
                   model = rep("Lead-Lag",300)))-> p2data
p2 = ggplot(p2data,aes(y = shat, x= model, fill = model))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = df100.out[df100.out$outlier.s,], aes(col = model))+
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col = wes_palette("Darjeeling1")[5]) +
  facet_wrap(~kernel)+
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  scale_colour_manual(values = wes_palette("Darjeeling1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab("Method (n=100)")+
  ylab("Fitted value of s")

# select LExp results in various sample size and plot figure 1c
dat_exp <- dat %>% filter(kernel == "LExp") %>%
  dplyr::select(ahat,shat,n) %>%
  tidyr::pivot_longer(ahat:shat,names_to ="paramater")%>% group_by(paramater,n) %>%
  mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
p4 <- ggplot(dat_exp,aes(y = value, x= paramater, fill = n))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = dat_exp[dat_exp$outlier.p,], aes(fill = n, col = n))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  scale_fill_startrek() +
  scale_color_startrek() +
  scale_x_discrete(labels = c("a", "s"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey10"),
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("Parameter")+
  ylab("Fitted value")+
  ylim(-1,4)
p4
# scale_fill_simpsons()

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/simulation_across_kernels.jpg",width = 15, height = 6,units = "in",res=450)
(p1 + p2 +p4) + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 22),
        plot.tag = element_text(size = 16)
  )
dev.off()

# For Biometrics reviewer 2
### MLE of sigma2 and b
dat_exp <- dat %>% filter(kernel == "LExp") %>%
  dplyr::select(bhat, tau2, n) %>%
  dplyr::mutate(bsigma2 = bhat*tau2) %>% 
  tidyr::pivot_longer(c(bhat,tau2,bsigma2),names_to ="paramater")%>% group_by(paramater,n) %>%
  mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
p4.supp <- ggplot(dat_exp,aes(y = value, x= paramater, fill = n))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = dat_exp[dat_exp$outlier.p,], aes(fill = n, col = n))+
  geom_hline(yintercept = 1.2, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 0.3, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 4, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  scale_fill_startrek() +
  scale_color_startrek() +
  scale_x_discrete(labels = c( "b",  expression(sigma^2*b), expression(sigma^2)))+
  theme(panel.background = element_rect(fill = "white", colour = "grey10"),
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("Parameter")+
  ylab("Fitted value")+
  ylim(-1,4)
p4.supp

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/simulation1_bsigma2.jpg",width = 15, height = 6,units = "in",res=450)
(p4.supp) + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 22),
        plot.tag = element_text(size = 16)
  )
dev.off()

### larger range of b to show its impact on a and s estimation
n = 100
b = 10
set.seed(781)
Y_exp100_b10 = simulateData(b,a,s,tau2, n, kernel = "exp")
b = 0.3
Y_exp100_b03 = simulateData(b,a,s,tau2, n, kernel = "exp")
b = 5
Y_exp100_b5 = simulateData(b,a,s,tau2, n, kernel = "exp")
b = 20
Y_exp100_b20 = simulateData(b,a,s,tau2, n, kernel = "exp")
res100_b <- list()
res100_b[["b0.3"]] <- res100[["exp_exp"]]
res100_b[["b5"]] <- deriveEstimation(Y_exp100_b5, kernel ="exp",b0 = 3, s0=1.5)
res100_b[["b10"]] <- deriveEstimation(Y_exp100_b10, kernel ="exp", b0 = 10,s0=1.5)
res100_b[["b20"]] <- deriveEstimation(Y_exp100_b20, kernel = "exp", b0 = 20,s0=1.5, bu = 40)
res_filter100_b <- lapply(1:4,function(i){
  res100_b[[i]][res100_b[[i]][,1] != softplus(-10) & res100_b[[i]][,1] != softplus(20) & res100_b[[i]][,2] < 10,]
})
df100_b <- do.call("rbind", res_filter100_b) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 100,
         kernel = rep(c("0.3","10","20","5"), unname(unlist(lapply(res_filter100_b,nrow)))))
df100_b %>% group_by(kernel) %>%
  mutate(outlier.p = is.outlier(ahat),outlier.s = is.outlier(shat)) %>%
  ungroup() -> df100_b.out
df100_b.out$kernel <- factor(df100_b.out$kernel, levels = c(0.3, 5, 10, 20))
## plot various b influence on a and s
p_R2Q6a = ggplot(df100_b.out,aes(y = ahat, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df100_b.out[df100_b.out$outlier.p,], aes(col = kernel))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("lengthscale b (n=100)")+
  ylab("Fitted value of a") +
  ylim(-1,3)

p_R2Q6b = ggplot(df100_b.out,aes(y = shat, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df100_b.out[df100_b.out$outlier.p,], aes(col = kernel))+
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("lengthscale b (n=100)")+
  ylab("Fitted value of s") +
  ylim(-1,3)
jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/simulation1_various_b.jpg",width = 15, height = 6,units = "in",res=450)
(p_R2Q6a + p_R2Q6b) + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 22),
        plot.tag = element_text(size = 16)
  )
dev.off()
