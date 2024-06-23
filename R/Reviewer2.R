library(tidyverse)
library(Matrix)
source("R/GPlag.R")
n=100
x <- seq(from = 1, to = 100, by = 1)
t <- rep(x,2)
# Define the true linear relationship (e.g., slope = 2, intercept = 3)
slope <- 2
intercept <- 3

# Generate y using the linear relationship and add some Gaussian noise
set.seed(42) # For reproducibility
y1 <- slope * x + intercept
y2 <- slope * (x-20) + intercept


# Generate noise from the t-distribution
# Degrees of freedom for the t-distribution
df <- 5
# Standard deviation of the noise
noise_sd <- 5

mse.noise <- sapply(1:100, function(i){
  t_noise1 <- rt(n, df) * noise_sd
  t_noise2 <- rt(n, df) * noise_sd

  # Add the noise to y
  y_noisy1 <- y1 + t_noise1
  y_noisy2 <- y2 + t_noise2
  Y_noise <- data.frame(y=c(y_noisy1,y_noisy2)) %>% as.matrix()
  rownames(Y_noise) <- t
  # Plot the original and noisy data
  # plot(x, y1, main = "Original and Noisy Data", xlab = "Input", ylab = "Output", ylim= c(-40,200))
  # points(x, y2, col = "gray")
  # points(x, y_noisy1, col = "red")
  # points(x, y_noisy2, col = "green")
  # legend("topright", legend = c("Original","Noisy1", "Noisy2"), col = c("blue", "red", "green"), lty = 1)
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
  Y_noise_param <- deriveEstimation(Y_noise.train, "rbf", sl = 15, su = 25) #perform MTSE - return b,a,s, tau2, ll in each column
  yres <- interpolation_sim(ttrain,Y_noise.train,ytruth,Y_noise_param,ttest, kernel = "lrbf")

  # plot(ttest, ytruth, main = "Fitting GPlag", xlab = "Input", ylab = "Output", ylim= c(-40,200))
  # points(ttrain,Y_noise.train, col="gray")
  # lines(ttest[1:(n/2)], yres[1:(n/2)], col = "blue", lwd = 2)
  # lines(ttest[(n/2+1):(n)], yres[(n/2+1):(n)], col = "green", lwd = 2)
  # legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","blue", "red", "green"), lty = 1)

  res.gplag <- mean((ytruth - yres)^2)

  # fit seperate rbf GP
  rbf_kernel <- rbfdot(sigma = 1)
  gp_model <- gausspr(as.matrix(ttrain[1:(n/2)]), Y_noise.train[1:(n/2)], kernel = rbf_kernel)
  predictions1 <- predict(gp_model, as.matrix(ttest[1:(n/2)]))
  gp_model <- gausspr(as.matrix(ttrain[(n/2+1):(n)]), Y_noise.train[(n/2+1):(n)], kernel = rbf_kernel)
  predictions2 <- predict(gp_model, as.matrix(ttest[(n/2+1):(n)]))

  # plot(ttest, ytruth, main = "Fitting SGP", xlab = "Input", ylab = "Output", ylim= c(-40,200))
  # points(ttrain,Y_noise.train, col="gray")
  # lines(ttest[1:(n/2)], predictions1, col = "blue", lwd = 2)
  # lines(ttest[(n/2+1):(n)], predictions2, col = "green", lwd = 2)
  # legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","blue", "red", "green"), lty = 1)

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

  # plot(ttest, ytruth, main = "Fiting splines", xlab = "Input", ylab = "Output", ylim= c(-40,200))
  # points(ttrain,Y_noise.train, col="gray")
  # lines(ttest[1:(n/2)], y1hat, col = "blue", lwd = 2)
  # lines(ttest[(n/2+1):(n)], y2hat, col = "green", lwd = 2)
  # legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","blue", "red", "green"), lty = 1)

  return(c(res.sgp, res.magma, res.gplag, res.splines))
})

datr2 <- t(mse.noise) %>% as.data.frame() %>%  `colnames<-`(c("SGP", "MAGMA","GPlag","Splines")) %>%
  tidyr::pivot_longer(everything(),names_to = "model") %>%
  # mutate(outlier.p = is.outlier(value)) %>%
  ungroup()

pltr2 <-ggplot(datr2 %>% filter(value<100),aes(x=model,y=value, fill = model))+
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
pltr2

jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer2.jpg",width = 10, height = 10,units = "in",res=450)
pltr2
dev.off()



# Set up the jpeg output
jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/reviewer2viz.jpg", width = 10, height = 10, units = "in", res=450)
par(mfrow = c(2, 2))

# Plot 1
plot(ttest, ytruth, main = "Fitting GPlag", xlab = "Input", ylab = "Output", ylim= c(-40,200))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], yres[1:(n/2)], col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], yres[(n/2+1):(n)], col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 3
plot(ttest, ytruth, main = "Fitting MAGMA", xlab = "Input", ylab = "Output", ylim= c(-40,200))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], y1hatm, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], y2hatm, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 3
plot(ttest, ytruth, main = "Fitting SGP", xlab = "Input", ylab = "Output", ylim= c(-40,200))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], predictions1, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], predictions2, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","gray", "blue", "green"), lty = 1)

# Plot 2
plot(ttest, ytruth, main = "Fiting splines", xlab = "Input", ylab = "Output", ylim= c(-40,200))
points(ttrain,Y_noise.train, col="gray")
lines(ttest[1:(n/2)], y1hat, col = "blue", lwd = 2)
lines(ttest[(n/2+1):(n)], y2hat, col = "green", lwd = 2)
legend("topleft", legend = c("test truth","observed train","predicted ts1", "predicted ts1"), col = c("black","gray", "blue", "green"), lty = 1)


# Close the jpeg device
dev.off()

