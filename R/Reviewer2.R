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
source("R/GPlag.R")
n=100
x <- seq(from = 1, to = 100, by = 1)
t <- rep(x,2)

# Generate y using the linear relationship and add some Gaussian noise
set.seed(42) # For reproducibility
y <- arima.sim(list(ar = 0.7), n = 200)
y1 <- y[1:n]
# Generate the second time series (ts2) with a time shift
x_shift <- 10  # Number of time steps to shift
y2 <- y[(1+x_shift) : (n+x_shift)]


# Function to perform the simulation and evaluation
simulate_and_evaluate <- function(noise_sd) {
  mse.noise <- sapply(1:100, function(i){
    # t_noise1 <- rnorm(n, sd = noise_sd)
    # t_noise2 <- rnorm(n, sd = noise_sd)
    # Generate the time series using arima.sim with custom noise
    custom_noise <- rnorm(n+x_shift, mean = 0, sd = noise_sd)  # For example, with a standard deviation of 2
    y <- arima.sim(list(ar = 0.5), n = n+x_shift,  innov = custom_noise)
    y1_noisy1 <- y[1:n]
    y2_noisy2 <- y[(1+x_shift) : (n+x_shift)]
    # Add the noise to y
    # y_noisy1 <- y1 + t_noise1
    # y_noisy2 <- y2 + t_noise2
    
    Y_noise <- data.frame(y=c(y_noisy1,y_noisy2)) %>% as.matrix()
    rownames(Y_noise) <- t
  
    train.index = sample(1:n,n/2)
    train.index <- train.index[order(train.index)]
    Y_noise.train <- Y_noise[c(train.index,train.index+n),] %>% as.matrix()
    rownames(Y_noise.train) <- t[c(train.index,train.index+n)]
  
    Y_noise.test <- Y_noise[-c(train.index,train.index+n),] %>% as.matrix()
    rownames(Y_noise.test) <- t[-c(train.index,train.index+n)]
    ttrain <- (rownames(Y_noise.train) %>% as.numeric())
    ttest = (rownames(Y_noise.test) %>% as.numeric())
    ytruth = c(y1,y2)[-c(train.index,train.index+n)]
    ## prepare data running GPlag
    Y_noise_param <- deriveEstimation(Y_noise.train, "rbf", b0=-1, a0 = -1, s0 = 10, sl = 5, su = 15) #perform MTSE - return b,a,s, tau2, ll in each column
    ## return value first n are prediction and the last value is -loglikelihood
    yres <- interpolation_sim(ttrain,Y_noise.train,ytruth,Y_noise_param,ttest, kernel = "lrbf")
    res.gplag.lrbf <- mean((ytruth - yres[1:n])^2)
    
    Y_noise_param <- deriveEstimation(Y_noise.train, "exp", b0=-1, a0 = -1, s0 = 10, sl = 5, su = 15) #perform MTSE - return b,a,s, tau2, ll in each column
    ## return value first n are prediction and the last value is -loglikelihood
    yres.exp <- interpolation_sim(ttrain,Y_noise.train,ytruth,Y_noise_param,ttest, kernel = "lexp")
    res.gplag.lexp <- mean((ytruth - yres.exp[1:n])^2)
    
    Y_noise_param <- deriveEstimation(Y_noise.train, "matern", b0=-1, bu=2, a0 = -1, s0 = 10, sl = 5, su = 15) #perform MTSE - return b,a,s, tau2, ll in each column
    ## return value first n are prediction and the last value is -loglikelihood
    yres.mat <- interpolation_sim(ttrain,Y_noise.train,ytruth,Y_noise_param,ttest, kernel = "lmatern")
    res.gplag.lmatern <- mean((ytruth - yres.mat[1:n])^2)
  
    # fit seperate rbf GP
    rbf_kernel <- rbfdot(sigma = 1)
    gp_model <- gausspr(as.matrix(ttrain[1:(n/2)]), Y_noise.train[1:(n/2)], kernel = rbf_kernel)
    predictions1 <- predict(gp_model, as.matrix(ttest[1:(n/2)]))
    gp_model <- gausspr(as.matrix(ttrain[(n/2+1):(n)]), Y_noise.train[(n/2+1):(n)], kernel = rbf_kernel)
    predictions2 <- predict(gp_model, as.matrix(ttest[(n/2+1):(n)]))
  
    res.sgp <- mean((ytruth - c(predictions1, predictions2))^2)
  
    ## prepare data running magma
    magma_train <- tibble(ID = rep(c("1","2"), each = n/2), Input = ttrain, Output = Y_noise.train[,1])
    predict.index = sample(n/2,size = n/4)
    predict.index <- predict.index[order(predict.index)]
    magma_tr <- magma_train[c(predict.index, predict.index+n/2),]
  
    magma_pred1 <- magma_train[setdiff(1:(n/2),predict.index),]
    magma_pred2 <- magma_train[setdiff((n/2+1):(n),predict.index+n/2),]
    model <- train_magma(data = magma_tr, common_hp = T)
  
    pred1  <- pred_magma(data = magma_pred1,
                         trained_model = model,
                         grid_inputs = ttest[1:(n/2)]+runif(n/2,-1/4,1/4), plot = T)
    pred1 <- pred1[order(pred1$Input),]
    y1hatm = pred1$Mean
    pred2  <- pred_magma(data =  magma_pred2,
                         trained_model = model,
                         grid_inputs = ttest[(n/2+1):(n)], plot = T)
    pred2 <- pred2[order(pred2$Input),]
    y2hatm = pred2$Mean
    res.magma = mean((ytruth - c(y1hatm,y2hatm))^2)
  
    ## prepare data running splines
    fit1 = lm(Y_noise.train[1:(n/2)] ~ ns(ttrain[1:(n/2)], knots=NULL))
    y1hat = predict(fit1,data.frame(t1=ttest[1:(n/2)]))
  
    fit2 = lm(Y_noise.train[(n/2+1):(n)] ~ ns(ttrain[(n/2+1):(n)], knots=NULL))
    y2hat = predict(fit2,data.frame(t2=ttest[(n/2+1):(n)]))
    res.splines = mean((ytruth - c(y1hat,y2hat))^2)
  
    return(c(res.sgp, res.magma, res.gplag.lrbf, res.gplag.lexp, res.gplag.lmatern, res.splines))
  })
return(mse.noise)
}

noise_levels <- c(0.4, 2, 10)
# try with both ar - 0.5 and 0.998
mse_results_all0.5 <- lapply(noise_levels, function(noise_sd) {
  mse_results <- simulate_and_evaluate(noise_sd)
  return(mse_results)
})

# Combine the results into a single data frame with a column for the noise level
mse_results_df0.5 <- map2_df(mse_results_all0.5, noise_levels, function(result, noise_sd) {
  t(result) %>% 
    as.data.frame() %>% 
    setNames(c("SGP", "MAGMA", "GPlag(LRBF)", "GPlag(LExp)", "GPlag(LMat3/2)", "Splines")) %>%
    mutate(noise_sd = noise_sd)
})

# Reshape the data frame to long format
datr <- mse_results_df0.5 %>%
  pivot_longer(-noise_sd, names_to = "model") %>%
  ungroup()

datr2 <- t(mse.noise) %>% as.data.frame() %>%  `colnames<-`(c("SGP", "MAGMA","GPlag(LRBF)", "GPlag(LExp)", "GPlag(LMat3/2)","Splines")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

pltr2 <-ggplot(datr %>% filter(value<100),aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  # geom_point(data = dat3.1[dat3.1$outlier.p,],aes(  col = model))+
  facet_wrap(~noise_sd)+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("IsleofDogs1"), labels = c("GPlag(LExp)", "GPlag(LMat3/2)","GPlag(LRBF)","Seperate GPs", "Splines","MAGMA"))+
  scale_colour_manual(values = wes_palette("IsleofDogs1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        legend.position = "none"
  ) +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")
pltr2

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer2.jpg",width = 10, height = 10,units = "in",res=450)
pltr2
dev.off()



# Set up the jpeg output
jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer2viz.jpg", width = 10, height = 10, units = "in", res=450)
par(mfrow = c(2, 2))

yl = -5
yu = 5
# Plot 1
plot(ttest, ytruth, main = "Fitting GPlag LRBF", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], yres[1:(n/2)], col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], yres[(n/2+1):(n)], col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 2
plot(ttest, ytruth, main = "Fitting GPlag LExp", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], yres.exp[1:(n/2)], col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], yres.exp[(n/2+1):(n)], col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 3
plot(ttest, ytruth, main = "Fitting GPlag LMat3/2", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], yres.mat[1:(n/2)], col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], yres.mat[(n/2+1):(n)], col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 4
plot(ttest, ytruth, main = "Fitting MAGMA", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], y1hatm, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], y2hatm, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 5
plot(ttest, ytruth, main = "Fitting SGP", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], predictions1, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], predictions2, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 6
plot(ttest, ytruth, main = "Fiting splines", xlab = "Input", ylab = "Output", ylim= c(yl,yu))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], y1hat, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], y2hat, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts2"), col = c("black","gray", "blue", "green"), lty = 1)


# Close the jpeg device
dev.off()

