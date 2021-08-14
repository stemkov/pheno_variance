# exploring whether shift in temperature to be more important variable changes observed variance

run.sim <- function(precip_mean, precip_sd, temp_mean, temp_sd, dev_mean, dev_sd){
  precip <- rnorm(1,precip_mean, precip_sd)
  temp <- rnorm(1,temp_mean, temp_sd)
  dev <- rnorm(1,dev_mean, dev_sd)
  
  pheno <- min(precip, temp, dev)
  cue <- c("precip", "temp", "dev")[which(c(precip, temp, dev) == min(c(precip, temp, dev)))]
  return(c(pheno = pheno, cue = cue))
}

run.sim(100, 10, 100, 10, 100, 10)

# historic phenology with ques equal
precip_means <- rep(100, 1000)
precip_sds <- rep(10, 1000)
temp_means <- rep(100, 1000)
temp_sds <- rep(5, 1000)
dev_means <- rep(105, 1000)
dev_sds <- rep(5, 1000)
historic_sim <- as.data.table(t(
  mapply(run.sim, precip_means, precip_sds, temp_means, temp_sds, dev_means, dev_sds, SIMPLIFY = T)))
historic_sim$pheno <- as.numeric(historic_sim$pheno)

plot(historic_sim$pheno)
plot(table(historic_sim$cue))
paste("mean", mean(historic_sim$pheno), "sd", sd(historic_sim$pheno))

# modern phenology with advanced, more important temp
precip_means <- rep(100, 1000)
precip_sds <- rep(10, 1000)
temp_means <- rep(90, 1000)
temp_sds <- rep(5, 1000)
dev_means <- rep(105, 1000)
dev_sds <- rep(5, 1000)
modern_sim <- as.data.table(t(
  mapply(run.sim, precip_means, precip_sds, temp_means, temp_sds, dev_means, dev_sds, SIMPLIFY = T)))
modern_sim$pheno <- as.numeric(modern_sim$pheno)

plot(modern_sim$pheno)
plot(table(modern_sim$cue))
paste("mean", mean(modern_sim$pheno), "sd", sd(modern_sim$pheno))

### population size affecting measured variance
pop.size.sim <- function(n){
  pop <- rnorm(n,100,10)
  pheno <- min(pop)
  return(pheno)
}

n_vals <- rep(c(10,1000,100000), each=1000)
x_pos <- rep(c(1,2,3), each=1000)
sim <- sapply(n_vals, pop.size.sim)

data <- data.table(n_vals, sim, x_pos)

plot(sim ~ jitter(x_pos), data=data, xlim=c(0.5,3.5), xaxt="n", ylab="Day Of Year Observation", xlab="Population Size")
axis(1, at=c(1,2,3), labels=c("n=10","n=1,000","n=100,000"))


par(mfrow=c(1,3))
plot(data[n_vals==10, sim], ylim=c(40,100), main="n = 10", ylab="Day Of Year")
plot(data[n_vals==1000, sim], ylim=c(40,100), main="n = 1,000")
plot(data[n_vals==100000, sim], ylim=c(40,100), main="n = 100,000")
par(mfrow=c(1,1))


### phenology coming up against an early season "wall"
library(shape)

time <- seq(50,150,0.1)
historic_pheno <- dnorm(time, mean=100, sd=10)
advanced_pheno <- dnorm(time, mean=80, sd=10)
early_season_filter <- pnorm(time, mean=70, sd=4)

plot(historic_pheno*25 ~ time, type="l", lwd=2, col="blue",
     xlab="Day of year", ylab="Abundance", yaxt="n")
polygon(c(time, max(time), min(time)),
        c(-early_season_filter + 1, 0, 0),
        col=adjustcolor("red",0.3), border=NA)
lines(-early_season_filter + 1 ~ time, lwd=2, col="red")
lines(historic_pheno*25 ~ time, lwd=2, col="blue")
lines(advanced_pheno*25 ~ time, lwd=2, col="#38c0ff")


historic_realized <- historic_pheno * early_season_filter
advanced_realized <- advanced_pheno * early_season_filter

#lines(historic_realized*25, type="l")
lines(advanced_realized*25 ~ time, col="#38c0ff", lwd=2, lty=2)

legend("topright", inset=0.03, cex=0.75,
       legend=c("Environmental filter", "Historic phenology", "Advanced phenology (underlying)", "Advanced phenology (filterd)"),
       lty=c(1,1,1,2), lwd=2, col=c("red", "blue", "#38c0ff", "#38c0ff"))

Arrows(100, 0.9, 80, 0.9,
       arr.type = "triangle", arr.adj = 1)
Arrows(56, 0.05, 65, 0.05,
       arr.type = "triangle", arr.adj = 1,
       col="red")

text(90, 0.85, "Phenological advance", cex=0.75)
text(60, 0.13, "Early season \nenvironmental constraint", cex=0.75)

### residual method demo/ 4 pheno metrics description

mean.sigma.test <- function(doy, explan, plot=FALSE, counter=FALSE, ...){
  
  if(FALSE){
    test <- all_data[which(all_data$species == "Spinus tristis"),]
    doy <- test$doy
    explan <- test$year
    
    test <- all_data[species == "Pieris rapae" & site == 100 & phenophase == "first_occurence",]
    doy <- test$doy
    explan <- test$year
    
    test <- all_data[species == "Hirundo rustica" & site == 100 & phenophase == "first_occurence",]
    doy <- test$doy
    explan <- test$year
  }
  
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
  res <- abs(residuals(mean_model)) * sqrt(2) # it took so long to figure out to multiply by sqrt(2) >:/
  res_model <- lm(res ~ explan)
  res_summ <- summary(res_model)
  
  bp_test <- bptest(mean_model)
  
  if(plot){
    par(mfrow=c(2,1), mar=c(2.5, 4.1, 1, 2.1))
    visreg(mean_model, ylab="Phenophase DOY", ...)
    par(mar=c(4.1, 4.1, 0, 2.1))
    visreg(res_model, ylab="Absolute residuals * sqrt(2)", ...)
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
  
  # if(may_plot){
  #   if(results$mean_slope > slope_lower & results$mean_slope < slope_upper){
  #     if(!(any(doy > out_upper) | any(doy < out_lower)))
  #       if(may_plot_frame == "mean") points(doy ~ explan, col=rgb(0.2,0,0.2,0.1), pch=19)
  #       if(may_plot_frame == "sigma") points(res ~ explan, col=rgb(0,0.2,0.2,0.1), pch=19)
  #   }
  # }
  
  return(results)
}

set.seed(33)
n <- 35
n <- 10
years <- 1:n
underlying_pheno <- seq(120,110, length.out = n)
underlying_temp <- seq(110,120, length.out = n)
underlying_temp_sd <- seq(0,3,length.out = n)
underlying_pheno_sd <- seq(6,0,length.out = n)
temp_data <- underlying_temp + rnorm(n,0,underlying_temp_sd)
pheno_data <- underlying_pheno + rnorm(n,0,underlying_pheno_sd)
temp_data <- -pheno_data + 230 + rnorm(n,0,1)

pheno_model <- lm(pheno_data ~ years)
pheno_res <- abs(residuals(pheno_model)) * sqrt(2) 
#pheno_res <- lm(res ~ years)
temp_model <- lm(pheno_data ~ temp_data)
temp_res <- abs(residuals(temp_model)) * sqrt(2) 
#temp_res <- lm(res ~ years)


# layout(matrix(c(1,1,1,1,2,3,4,5), 2,4))
# plot(temp_data, type="l", col="red")
# points(pheno_data)
# plot(pheno_data ~ years)
# plot(pheno_data ~ temp_data)
# plot(pheno_res ~ years)
# plot(temp_res ~ temp_data)

png("Figures/res_method_demo.png", width=625, height=600)
set.seed(30)
n <- 35
pheno <- seq(125, 100, length.out = n) + rnorm(n, 0, 4)
#years <- 1986:2020
years <- round(seq(1986,2020, length.out = n))
pheno_temp <- seq(125, 100, length.out = n) + rnorm(n, 0, seq(0,15, length.out = n))
temp <- sort(runif(n, 5, 20))

layout(matrix(c(1,2,3,4), 2,2))
res_pal <- c("#ff2e1f", "#346cfa")
point_pal <- c("#24b6d6","#e39536")

par(cex = 0.9)
par(mar = c(0, 0, 0, 0), oma = c(5, 5, 1, 1))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

plot(pheno ~ years, xaxt="n", las=1, col="white")
pheno_mod <- lm(pheno ~ years); pheno_res <- abs(residuals(pheno_mod)) * sqrt(2) 
abline(pheno_mod)
model_pred <- predict(pheno_mod)
above_below <- sapply(pheno - model_pred, function(x) ifelse(x >= 0, 1, 2)) # recording whether residual is above or below prediction
arrows(years, model_pred,
       years, pheno, length=0,
       col=res_pal[above_below], lwd=2)
points(pheno ~ years, col=point_pal[1], pch=19, cex=1)

mtext("Phenology (day of year)", side=2, line=3)
mtext("(-) Mean shift  ", side=3, line=-1.5, adj=1)

plot(pheno_res ~ years, las=2, col="white")
pheno_res_mod <- lm(pheno_res ~ years)
abline(pheno_res_mod)
arrows(years, 0,
       years, pheno_res, length=0,
       col=res_pal[above_below], lty=3, lwd=2)
points(pheno_res ~ years, col=point_pal[1], pch=19, cex=1)

#mtext(paste("|R|",expression(sqrt(2)), "(days)"), side=2, line=3)
mtext(bquote(abs(R)*" × "*sqrt(2)*"  (days)"), side=2, line=3)
mtext("Year", side=1, line=3)
mtext("(∅) Variance shift  ", side=3, line=-1.5, adj=1)

plot(pheno_temp ~ temp, xaxt="n", yaxt="n", col="white")
temp_mod <- lm(pheno_temp ~ temp); temp_res <- abs(residuals(temp_mod)) * sqrt(2) 
abline(temp_mod)
temp_pred <- predict(temp_mod)
above_below_temp <- sapply(pheno_temp - temp_pred, function(x) ifelse(x >= 0, 1, 2)) # recording whether residual is above or below prediction
arrows(temp, temp_pred,
       temp, pheno_temp, length=0,
       col=res_pal[above_below_temp], lwd=2)
points(pheno_temp ~ temp, col=point_pal[2], pch=19, cex=1)
mtext("(-) Mean  \n sensitivity  ", side=3, line=-2.5, adj=1)

plot(temp_res ~ temp, yaxt="n", col="white")
temp_res_mod <- lm(temp_res ~ temp)
abline(temp_res_mod)
arrows(temp, 0,
       temp, temp_res, length=0,
       col=res_pal[above_below_temp], lty=3, lwd=2)
points(temp_res ~ temp, col=point_pal[2], pch=19, cex=1)

mtext("Temperature (°C)", side=1, line=3)
mtext("(+) Variance  \n sensitivity  ", side=3, line=-2.5, adj=1)

par(fig=c(0.11,0.29,0.42,0.6), new=T)
plot((model_pred + abs(residuals(pheno_mod))) ~ years, axes=F,
     xlim=range(years)+c(0,0), ylim=range(pheno)+c(0,5))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
Arrows(par("usr")[2]+3.5, par("usr")[4]-4,
       par("usr")[2]+3.5, par("usr")[3]+4,
       arr.adj=1, arr.type="triangle", lwd=3, col="gray", xpd=NA)
#points((model_pred + abs(residuals(pheno_mod))) ~ years)
abline(pheno_mod)
arrows(years, model_pred,
       years, (model_pred + abs(residuals(pheno_mod))),
       length=0, col=res_pal[above_below], lwd=1)


par(fig=c(0.56,0.74,0.42,0.6), new=T)
plot((temp_pred + abs(residuals(temp_mod))) ~ temp, axes=F,
     xlim=range(temp)+c(0,0), ylim=range(pheno_temp)+c(0,5))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
Arrows(par("usr")[2]+1.5, par("usr")[4]-5,
       par("usr")[2]+1.5, par("usr")[3]+5,
       arr.adj=1, arr.type="triangle", lwd=3, col="gray", xpd=NA)
#points((temp_pred + abs(residuals(temp_mod))) ~ temp)
abline(temp_mod)
arrows(temp, temp_pred,
       temp, (temp_pred + abs(residuals(temp_mod))),
       length=0, col=res_pal[above_below_temp], lwd=1)

dev.off()



### same as above but defining data manually and less points
png("Figures/res_method_demo_manual.png", width=500, height=475)
n <- 10
pheno <- seq(125, 100, length.out = n) +
  c(3,-3,3,-3,3,-3,3,-3,3,-3)
years <- c(2001:2010)
#years <- round(seq(1986,2020, length.out = n))
pheno_temp <- seq(125, 100, length.out = n) +
  c(0.6,-0.2,0.3,-1.5,2,-4.2,4.5,-5,8,-7)
temp <- c(5,8,10,11,12,14.5,15,17,19,20)


layout(matrix(c(1,2,3,4), 2,2))
res_pal <- c("#ff2e1f", "#346cfa")
point_pal <- c("#24b6d6","#e39536")

par(cex = 0.9)
par(mar = c(0, 0, 0, 0), oma = c(5, 5, 1, 1))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

plot(pheno ~ years, xaxt="n", yaxt="n", las=1, col="white")
pheno_mod <- lm(pheno ~ years); pheno_res <- abs(residuals(pheno_mod)) * sqrt(2) 
abline(pheno_mod)
model_pred <- predict(pheno_mod)
above_below <- sapply(pheno - model_pred, function(x) ifelse(x >= 0, 1, 2)) # recording whether residual is above or below prediction
arrows(years, model_pred,
       years, pheno, length=0,
       col=res_pal[above_below], lwd=3)
points(pheno ~ years, col=point_pal[1], pch=19, cex=1.5)

mtext("Phenology (day of year)", side=2, line=2)
mtext("(-) Mean shift  ", side=3, line=-1.5, adj=1)

plot(pheno_res ~ years, yaxt="n", las=2, col="white", ylim=c(0,10))
pheno_res_mod <- lm(pheno_res ~ years)
pheno_res_mod <- rq(pheno_res ~ years, tau=0.68)
abline(pheno_res_mod)
arrows(years, 0,
       years, pheno_res, length=0,
       col=res_pal[above_below], lty=3, lwd=3)
points(pheno_res ~ years, col=point_pal[1], pch=19, cex=1.5)

#mtext(paste("|R|",expression(sqrt(2)), "(days)"), side=2, line=3)
#mtext(bquote(abs(R)*" × "*sqrt(2)*"  (days)"), side=2, line=2)
mtext("Absolute residuals (days)", side=2, line=2)
mtext("Year", side=1, line=3)
mtext("(∅) Variance shift  ", side=3, line=-1.5, adj=1)

plot(pheno_temp ~ temp, xaxt="n", yaxt="n", col="white", ylim=c(93,128))
temp_mod <- lm(pheno_temp ~ temp); temp_res <- abs(residuals(temp_mod)) * sqrt(2) 
abline(temp_mod)
temp_pred <- predict(temp_mod)
above_below_temp <- sapply(pheno_temp - temp_pred, function(x) ifelse(x >= 0, 1, 2)) # recording whether residual is above or below prediction
arrows(temp, temp_pred,
       temp, pheno_temp, length=0,
       col=res_pal[above_below_temp], lwd=3)
points(pheno_temp ~ temp, col=point_pal[2], pch=19, cex=1.5)
mtext("(-) Mean  \n sensitivity  ", side=3, line=-2.5, adj=1)

plot(temp_res ~ temp, yaxt="n", col="white", ylim=c(0,20))
temp_res_mod <- lm(temp_res ~ temp)
temp_res_mod <- rq(temp_res ~ temp, tau=0.68)
abline(temp_res_mod)
arrows(temp, 0,
       temp, temp_res, length=0,
       col=res_pal[above_below_temp], lty=3, lwd=3)
points(temp_res ~ temp, col=point_pal[2], pch=19, cex=1.5)

mtext("Temperature (°C)", side=1, line=3)
mtext("(+) Variance  \n sensitivity  ", side=3, line=-2.5, adj=1)

par(fig=c(0.04,0.22,0.42,0.6), new=T)
plot((model_pred + abs(residuals(pheno_mod))) ~ years, axes=F,
     xlim=range(years)+c(0,0), ylim=range(pheno)+c(0,5))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
Arrows(par("usr")[2]+1, par("usr")[4]-4,
       par("usr")[2]+1, par("usr")[3]+4,
       arr.adj=1, arr.type="triangle", lwd=3, col="gray", xpd=NA)
#points((model_pred + abs(residuals(pheno_mod))) ~ years)
abline(pheno_mod)
arrows(years, model_pred,
       years, (model_pred + abs(residuals(pheno_mod))),
       length=0, col=res_pal[above_below], lwd=3)


par(fig=c(0.56,0.74,0.42,0.6), new=T)
plot((temp_pred + abs(residuals(temp_mod))) ~ temp, axes=F,
     xlim=range(temp)+c(0,0), ylim=range(pheno_temp)+c(0,5))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
Arrows(par("usr")[2]+1.5, par("usr")[4]-5,
       par("usr")[2]+1.5, par("usr")[3]+5,
       arr.adj=1, arr.type="triangle", lwd=3, col="gray", xpd=NA)
#points((temp_pred + abs(residuals(temp_mod))) ~ temp)
abline(temp_mod)
arrows(temp, temp_pred,
       temp, (temp_pred + abs(residuals(temp_mod))),
       length=0, col=res_pal[above_below_temp], lwd=3)
dev.off()


