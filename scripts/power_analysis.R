# making sure that the residual method works

library(lmvar)
library(lmtest)
library(pbapply)

mean.sigma.test <- function(doy, explan, plot=FALSE, counter=FALSE, lmvar_method=FALSE, sqrt=FALSE, gls=FALSE, quant_reg=FALSE, ...){
  if(counter) global.counter()
  if(length(doy) != length(explan)) return(list(mean_slope = -8888, mean_se = -8888, mean_p = -8888, sigma_slope = -8888, sigma_se = -9999, sigma_p = -8888, bp_p = -8888, n=length(doy))) # doy and explan differ in lengths - this should never happen
  if(sum(!is.na(doy)) < 10 | sum(!is.na(explan)) < 10) return(list(mean_slope = -9999, mean_se = -9999, mean_p = -9999, sigma_slope = -9999, sigma_se = -9999, sigma_p = -9999, bp_p = -9999, n=length(doy))) # output for data.table all needs to be of the same type, so I can't output NA
  # if(sum(!is.na(doy)) != sum(!is.na(explan))) return(list(mean_slope = -7777, mean_p = -7777, sigma_slope = -7777, sigma_p = -7777, bp_p = -7777, n=length(doy))) # there in an NA in the doy or explan timeseries, but not both ... this was causing problems with manomet because some doys are NAs because phest couldn't estimate onset properly. I'm changing this to remove doys with NAs
  if(sum(!is.na(doy)) != sum(!is.na(explan))){
    warning("Unequal number of NAs in doy vs. explan") # not a problem if there are supposed to be some NAs in doy
    doy_explan <- data.frame(doy, explan)
    doy_explan <- doy_explan[which(complete.cases(doy_explan)),] # getting rid of years with either missing doy or explan data
    doy <- doy_explan$doy
    explan <- doy_explan$explan
  }
  
  mean_model <- lm(doy ~ explan)
  mean_summ <- summary(mean_model)
  
  #res <- abs(residuals(mean_model))
  if(quant_reg){
    res <- abs(residuals(mean_model))
    res_model <- rq(res ~ explan, tau=0.6827)
    res_summ <- summary(res_model, se="boot")
  } else{
    res <- abs(residuals(mean_model)) * sqrt(2)
    res_model <- lm(res ~ explan)
    res_summ <- summary(res_model)
  }
  
  if(sqrt){
    sign <- sign(res_summ$coefficients[2])
    sqrt_coef <- sqrt(abs(res_summ$coefficients[2])) * sign
    res_summ$coefficients[2] <- sqrt_coef
  } 
  bp_test <- bptest(mean_model)
  if(gls){
    #plot(doy ~ explan)
    #print(doy)
    #print(explan)
    sigma_model <- gls(doy ~ explan, weights = ~ explan)
    sigma_summ <- summary(sigma_model)
    res_summ$coefficients[2] <- sigma_summ$sigma
  }
  
  if(lmvar_method){
    X = model.matrix(~ explan - 1)
    model = lmvar(doy, X_mu = X, X_sigma = X)
    summ <- summary(model)
    if(plot){
      plot(doy ~ explan)
      fitted <- fitted(model)
      lines(fitted[,1] ~ explan)
      lines(fitted[,1] + fitted[,2] ~ explan, lty=2, col="red")
      lines(fitted[,1] - fitted[,2] ~ explan, lty=2, col="red")
    }
    results <- list(mean_slope = as.numeric(summ$coefficients[2,1]),
                    mean_p = as.numeric(summ$coefficients[2,4]),
                    mean_se = as.numeric(summ$coefficients[2,2]),
                    sigma_slope = as.numeric(summ$coefficients[4,1]),
                    sigma_se = as.numeric(summ$coefficients[4,2]),
                    sigma_p = as.numeric(summ$coefficients[4,4]),
                    bp_p = as.numeric(bp_test$p.value),
                    n = length(doy))
    return(results)
  }
  
  if(plot & !lmvar_method){
    par(mfrow=c(2,1), mar=c(2.5, 4.1, 1, 2.1))
    visreg(mean_model, ylab="Phenophase DOY", ...)
    par(mar=c(4.1, 4.1, 0, 2.1))
    visreg(res_model, ylab="Absolute value of residuals", ...)
    par(mfrow=c(1,1), mar=c(4.1, 4.1, 2.1, 2.1))
  }
  
  results <- list(mean_slope = as.numeric(mean_summ$coefficients[2]),
                  mean_p = as.numeric(mean_summ$coefficients[8]),
                  mean_se = as.numeric(mean_summ$coefficients[4]),
                  sigma_slope = as.numeric(res_summ$coefficients[2]),
                  sigma_se = as.numeric(res_summ$coefficients[4]),
                  sigma_p = as.numeric(res_summ$coefficients[8]),
                  bp_p = as.numeric(bp_test$p.value),
                  n = length(doy))
  
  return(results)
}

simulate.data <- function(n, mean_baseline, mean_change, sigma_baseline, sigma_change, plot=FALSE){
  mean_seq <- seq(mean_baseline, mean_baseline+mean_change, length.out = n)
  error_seq <- seq(sigma_baseline, sigma_baseline+sigma_change, length.out = n)
  error <- sapply(error_seq, function(x) rnorm(1,0,x))
  #error_seq <- seq(sigma_baseline, sigma_baseline+sigma_change, length.out = n)
  #error <- sapply(error_seq, function(x) rnorm(1,0,sqrt(x)))
  #error <- mean_seq + rnorm(n, 0, error_seq*sqrt(1:n))
  data <- mean_seq + error
  if(plot) plot(data)
  return(data)
}

power.analysis <- function(n, mean_baseline, mean_change, sigma_baseline, sigma_change, lmvar_method=F, sqrt=F, gls=F, quant_reg=F, plot=F){
  data <- simulate.data(n, mean_baseline, mean_change, sigma_baseline, sigma_change)
  years <- 1:n
  #years <- 0:(n-1)
  test <- mean.sigma.test(data, years, lmvar_method = lmvar_method, sqrt=sqrt, gls=gls, quant_reg=quant_reg, plot=plot)
  
  actual_mean_slope <- mean_change/n
  #actual_sigma_slope <- (sqrt(sigma_baseline + sigma_change)-sqrt(sigma_baseline))/n
  #actual_sigma_slope <- (sqrt(sigma_baseline + sigma_change)-sqrt(sigma_baseline))/n
  actual_sigma_slope <- sigma_change/n
  
  # are the estimated slopes within 2*SE of estimate?
  mean_slope_within <- actual_mean_slope < test$mean_slope+test$mean_se*2 & actual_mean_slope > test$mean_slope-test$mean_se*2
  sigma_slope_within <- actual_sigma_slope < test$sigma_slope+test$sigma_se*2 & actual_sigma_slope > test$sigma_slope-test$sigma_se*2
  
  return(c(n = n,
           actual_mean_slope = actual_mean_slope,
           actual_sigma_slope = actual_sigma_slope,
           predicted_mean_slope = test$mean_slope,
           predicted_sigma_slope = test$sigma_slope,
           mean_correct = mean_slope_within,
           sigma_correct = sigma_slope_within))
}

# n <- 30
# mean_baseline <- 100
# mean_change <- 20
# sigma_baseline <- 1
# sigma_change <- 20
# data <- simulate.data(n, mean_baseline, mean_change, sigma_baseline, sigma_change)
# years <- 1:n
# mean.sigma.test(data, years, lmvar_method=F, plot=T)$sigma_slope
# mean.sigma.test(data, years, lmvar_method=F, plot=T, sqrt=T)$sigma_slope
# #mean.sigma.test(data, years, lmvar_method=T, plot=T)$sigma_slope
# power.analysis(32, 100, 20, 1, 100, lmvar_method = F, sqrt=F, plot=T)
# power.analysis(32, 100, 20, 10, -10, lmvar_method = F, sqrt=F, plot=T)
# #power.analysis(32, 100, 20, 1, 20, gls=T, plot=T)
# #mean.sigma.test(data, years, gls = T, plot=T)
# power.analysis(32, 100, 20, 1, 100, quant_reg = T, plot=F)
# 
# n_vals <- round(seq(10,60,length.out = 10))
# sigma_change_vals <- seq(-10,10,length.out = 10)
#sigma_change_vals <- seq(-50,100,length.out = 10)
#sigma_change_vals <- c(-15,  -10, -5, 0, 5, 10, 15, 20, 30)

iters <- 20
iters <- 100
n_iters <- rep(n_vals, each = iters*length(sigma_change_vals))
sigma_change_iters <- rep(sigma_change_vals, times = iters*length(n_vals))

total_runs <- length(n_iters)

run_power <- pbmapply(power.analysis, n=32, sigma_change=sigma_change_iters,
                    mean_baseline = 100, mean_change = runif(total_runs,-20,20),
                    sigma_baseline = 12, sqrt=F, quant_reg=T)
# run_power <- pbmapply(power.analysis, n=32, sigma_change=sigma_change_iters,
#                       mean_baseline = 100, mean_change = runif(total_runs,-20,20),
#                       sigma_baseline = 100, gls=T)
power <- as.data.table(t(run_power))
percentages <- power[, (sum(sigma_correct)/.N)*100, by = actual_sigma_slope]

mean_power_model <- lm(predicted_mean_slope ~ actual_mean_slope, data=power[mean_correct == 1,])
summary(mean_power_model)
sigma_power_model <- lm(predicted_sigma_slope ~ actual_sigma_slope, data=power)
summary(sigma_power_model)

png("Figures/power_2020_04_15_quantreg.png", width=600, height=400)
plot(predicted_sigma_slope ~ jitter(actual_sigma_slope), data=power, col=adjustcolor(-sigma_correct+2,alpha.f = 0.2),
     xlab = "Actual variance change", ylab = "Predicted variance change", pch=20, cex=0.7,
     xaxt="n")
axis(1,at=percentages$actual_sigma_slope, labels= round(percentages$actual_sigma_slope,2))
abline(0,1)
abline(v=0, col="red")
abline(sigma_power_model, lty=2)
text(seq(par("usr")[1]+0.03, par("usr")[2]-0.03, length.out = 10), y = seq(0.1,1.15, length.out = 10),
     labels = paste0(round(percentages$V1), "%"))
legend("topleft", inset=0.05,
       legend=c("1:1 line", "Method prediction", "Baseline variance", "Estimate within CI", "Estimate outside CI"),
       lty=c(1,2,1, NA, NA), col=c(1,1,2,1,2), pch=c(NA,NA,NA,20,20), cex=0.8)
dev.off()

# multiplying residuials by sqrt(2) gets rid of the bias

plot(predicted_mean_slope ~ actual_mean_slope, data=power)
write.csv(power, "Edited Data/power_analysis/power.csv")

mean(power$sigma_correct)
mean(power$mean_correct)

### troubleshooting residual method
x <- 1:32
var_slope <- -0.3
var_intercept <- 10
mean_slope <- 0
mean_intercept <- 100
error <- rnorm(length(x), sd = var_intercept + var_slope*x)
y <- mean_intercept + x*mean_slope + error
plot(y~x)

mean_model <- lm(y~x)
residuals <- abs(residuals(mean_model)) * sqrt(2)
res_model <- lm(residuals ~ x)
summary(mean_model)
summ <- summary(res_model); summ
plot(abs(error) ~ x, col="red", pch=20, ylim=c(0,20), ylab="Absolute error or residuals")
points(residuals ~ x)
legend("topright", inset=0.05, legend=c("Abs. Error", "Abs. Residuals * sqrt(2)", "Error estimate", "Residual method estimate", "Actual error simulation"),
       col=c("red", "black","red","black","blue"),
       pch=c(20,1,NA,NA,NA),
       lty=c(NA,NA,1,1,1),
       lwd=c(NA,NA,1,1,2))
abline(res_model, col="black") # estimate from residual method

error_model <- lm(abs(error) ~ x)
abline(error_model, col="red") # estimate if error is known

abline(var_intercept, var_slope, lwd=2, col="blue") # actual variance model
paste("actual:", var_slope, "residual estimated:", round(summ$coefficients[2],3), "error actual:", round(summary(error_model)$coefficients[2],3))

quant_model <- rq(abs(residuals(mean_model)) ~ x, tau=0.6827)
quant_summ <- summary(quant_model, se="boot")
quant_summ

paste("actual:", var_slope, "residual estimate:", round(summ$coefficients[2],3), "quantile estimate", round(quant_summ$coefficients[2,1], 3))
abline(quant_model, col="green")






### figured it out

set.seed(0)
error <- rnorm(1000, 0, 10)
abs_error <- abs(error)
sd(error)
mean(abs_error)

hist(error, seq(-40,40,by=1), ylim=c(0,100), col=adjustcolor("red", 0.5))
abline(v=c(-sd(error),sd(error)), lwd=2, col="red")

hist(abs_error, seq(-40,40,by=1), add=T, col=adjustcolor("blue", 0.5))
abline(v=mean(abs_error), lwd=2, col="blue", lty=1)
#abline(v=mean(abs_error*sqrt(2)), lwd=2, col="blue", lty=2)
#abline(v=quantile(abs_test,0.6827), lwd=2, col="blue", lty=3)

tau <- 0.682689492
quant <- rq(abs_error ~ 1, tau=tau)
quant_summ <- summary(quant, se="boot")
quant_summ
abline(v=quant_summ$coefficients[1,1], lwd=2, col="green", lty=3)

tau/0.5
sqrt(2)

# THE SOLUTION
res_model <- rq(abs_residuals ~ year, tau=0.682689492)































### trying dglm
# library(dglm)
# n <- 30
# mean_baseline <- 100
# mean_change <- 20
# sigma_baseline <- 1
# sigma_change <- 10
# data <- simulate.data(n, mean_baseline, mean_change, sigma_baseline, sigma_change)
# years <- 1:n
# 
# plot(data ~ years)
# 
# 
# 
# model <- dglm(data ~ years, ~ years, gaussian,ykeep=TRUE,xkeep=TRUE,zkeep=TRUE)
# summary(model)
# summ <- summary(model)
# summary(model$dispersion.fit)
# 
# plot(fitted(model$dispersion.fit),residuals(model$dispersion.fit))
# 
# paste("actual mean change:", mean_change/n, "predicted mean change:", summ$coefficients[2,1])
# #paste("actual sigma change:",sigma_change/n,  "predicted mean change:", sqrt(summ$dispersion.summary$coefficients[2,1]))
# paste("actual sigma change:",sigma_change/n,  "predicted mean change:", summ$dispersion.summary$coefficients[2,1])
# 
# plot(simulate(model)$sim_1 ~ years)
# ?simulate



### trying gls and weights variance function
# n <- 100
# mean_baseline <- 100
# mean_change <- 20
# sigma_baseline <- 100
# sigma_change <- 200
# data <- simulate.data(n, mean_baseline, mean_change, sigma_baseline, sigma_change)
# years <- 1:n
# #years <- years/100
# plot(data ~ years)

# dat <- data.frame(data, years)
# brm_model <- brm(bf(data ~ years, sigma ~ years), data=dat)
# summary(brm_model)
# plot(brm_model)
# plot(conditional_effects(brm_model), points = TRUE)

# model <- gls(data ~ years)
# model <- gls(data ~ years, weights = ~ years)
# # model <- gls(data ~ years, weights = varFixed(~ years)) # identical to above
# model <- gls(data ~ years, weights = varConstPower(form = ~ years)) 
# #model <- gls(data ~ years, weights = varConstPower(form = ~ years, fixed=list(power=1, const=1))) 
# summary(model)
# summ <- summary(model)
# paste("actual mean change:", mean_change/n*100, "predicted mean change:", summ$coefficients[2])
# paste("actual sigma change:",sqrt(sigma_change)/n*100,  "predicted sigma change:", summ$sigma)
# 
# plot((1/attr(model$model$varStruct, "weights")))
# plot((1/attr(model$model$varStruct, "weights"))^2)
# (1/attr(model$model$varStruct, "weights"))^2
# 
# attr(model$model$varStruct, "formula")
# str(model$model$varStruct)$varStruct
# 
# str(model)
# 
# model$sigma
# 
# model$varBeta
# 
# intervals(model)
# 
# x <- seq(1,10, length.out = 100)
# noise <- rnorm(n = 100, mean = 0, sd = 2 * sqrt(x)) # square root here because SD = sqrt(variance)
# y <- 1.2 + 2.1 * x + noise
# plot(x, y)




### doing a sanity test with two x values
# y1 <- rnorm(1000,100,sqrt(10))
# y2 <- rnorm(1000,100,sqrt(20))
# y <- c(y1, y2)
# x <- rep(c(0,1), each=1000)
# plot(y~x)
# # should be variance change of 10
# test <- mean.sigma.test(y, x, plot=T)
# model <- lm(y~x)
# bp_test <- ols_test_breusch_pagan(model)
# bp_test$
# sd(y1)
# sd(y2)
# sqrt(var(y1))
# sqrt(var(y2))
# paste("SD change:", (sd(y2)-sd(y1))) 
# paste("variance change:", (var(y2)-var(y1))/10, "estimated:", test$sigma_slope) 
# paste("bias: "(var(y2)-var(y1))/10)










### unrelated stuff for Claire

# library(measurements)
# 
# data <- data.frame(weight = c("150", "140", "200kg", "250KG", "300 Kg", "150", "150", "150", "150"),
#                    height = c("4'9", "5'2'", "166cm", "168 cm", "5 feet 4 inches", "5 feet 4 in", "66 inches", "6 foot", "7'"),
#                    stringsAsFactors = FALSE)
# 
# # testing whether height was written as two separate numbers - not necessarily in metric because some people wrote 66 inches
# data$height_is_standard <- grepl("\\d*\\D\\d", data$height)
# 
# # one number in american standard
# data[which(!data$height_is_standard & grepl("in", data$height)), "height"] <- gsub("\\D", "", data[which(!data$height_is_standard & grepl("in", data$height)), "height"])
# # one number in metric
# data[which(!data$height_is_standard & grepl("cm", data$height)), "height"] <- conv_unit(as.numeric(gsub("\\D", "", data[which(!data$height_is_standard & grepl("cm", data$height)), "height"])),
#                                                                                         "cm", "inch")
# # one number in feet
# data[which(!data$height_is_standard & grepl("ft|foot|feet|\\'", data$height)), "height"] <- conv_unit(as.numeric(gsub("\\D", "", data[which(!data$height_is_standard & grepl("ft|foot|feet|\\'", data$height)), "height"])),
#                                                                                                   "ft", "inch")
# 
# # two numbers - this must be american standard
# data[which(data$height_is_standard), "height"] <- 12 * as.numeric(gsub("^(\\d)\\D.*","\\1",data[which(data$height_is_standard), "height"])) + as.numeric(gsub("^\\d+\\D*(\\d*)\\D*","\\1",data[which(data$height_is_standard), "height"]))
# 
# 
# data[which(grepl("kg|Kg|KG", data$weight)), "weight"] <- as.character(2.205 * as.numeric(gsub("([0-9]+).*$", "\\1",data[which(grepl("kg|Kg|KG", data$weight)), "weight"])))
# data$weight <- as.numeric(data$weight)
# 
# 
# 
# data[which(data$unit %in% c("kg", "Kg", "KG")),"weight"] <- 2.205*data[which(data$unit %in% c("kg", "Kg", "KG")),"weight"]





