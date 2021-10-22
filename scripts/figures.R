# making figures to illustrate the analysis
# need to run data_analysis.R first

### temp over time - how is climate at the sites changing?
png("figures/temp_overtime.png", width=800, height=400)
par(mfrow=c(1,2))
hist(clim_trends$mean_slope, 100,
     xlab = "Mean temperature over time (degree C/year)",
     main="") # temps are increasing
abline(v=0, col="red", lwd=2)
hist(clim_trends$sigma_slope, 100,
     xlab = "Temperature variance over time (degree C/year)",
     main="") # almost no change in variance, but a small significant reduction
abline(v=0, col="red", lwd=2)
par(mfrow=c(1,1))
clim_mean <- lm(mean_slope ~ 1, data=clim_trends)
summary(clim_mean)
clim_sigma <- lm(sigma_slope ~ 1, data=clim_trends)
summary(clim_sigma)
dev.off()

### overall distributions
mean_sd_range <- c(-5*sd(trends$mean_slope), 5*sd(trends$mean_slope))
sigma_sd_range <- c(-5*sd(trends$sigma_slope), 5*sd(trends$sigma_slope))
mean_sd_range_temp <- c(-5*sd(trends$mean_slope_temp), 5*sd(trends$mean_slope_temp))
sigma_sd_range_temp <- c(-5*sd(trends$sigma_slope_temp), 5*sd(trends$sigma_slope_temp))

mean_sd_range <- range(trends$mean_slope)
sigma_sd_range <- range(trends$sigma_slope)
mean_sd_range_temp <- range(trends$mean_slope_temp)
sigma_sd_range_temp <- range(trends$sigma_slope_temp)

png("figures/overall_dists.png", width=700, height=500)
layout(matrix(1:4,nrow=2,byrow=F))
par(mar=c(4,4.1,2,0))
hist(trends$mean_slope, breaks=200, col=rgb(0, 0, 1, 0.7),
     main="Original", xlab="Mean shift", xlim=mean_sd_range)
abline(v=0, lty=2, lwd=2)
report.outliers(trends$mean_slope, "x")

hist(trends$mean_slope_temp, breaks=150, col=rgb(1, 0, 0, 0.7),
     main="", xlab="Mean sensitivity", xlim=mean_sd_range_temp)
abline(v=0, lty=2, lwd=2)
report.outliers(trends$mean_slope_temp, "x")

par(mar=c(4,4.1,2,2))
hist(trends$sigma_slope, breaks=300, col=rgb(0, 0.4, 1, 0.7),
     main="", ylab="", xlab="Variance shift", xlim=sigma_sd_range)
abline(v=0, lty=2, lwd=2)
report.outliers(trends$sigma_slope, "x")

hist(trends$sigma_slope_temp, breaks=200, col=rgb(1, 0.4, 0, 0.7),
     main="", ylab="", xlab="Variance sensitivity", xlim=sigma_sd_range_temp)
abline(v=0, lty=2, lwd=2)
report.outliers(trends$sigma_slope_temp, "x")
par(mfrow= c(1,1))
dev.off()


### summary of raw data figure (GIF)
palette <- brewer.pal(length(levels(as.factor(trends$dataset))), "Set1") # need to fix the rothamsted shit stain...
palette[7] <- "#6b6b6b"
palette[6] <- "#59f7ff"

saveGIF({
  ani.options(interval = 2, nmax = 10)
  for(i in unique(all_data$dataset)){
    color <- palette[which(levels(as.factor(all_data$dataset)) == i)]
    plot(doy ~ tmax, data=all_data[dataset == i,], pch=20, cex=0.7, col=adjustcolor(color,alpha.f = 0.2),
         xlab="Temperature in the typical month (C)", ylab="Date of phenophase (day-of-year)",
         xlim=c(-20,40), ylim=c(0,365), main=i, col.main=color, cex.main=2)
  }
}, movie.name = "raw_data.gif",ani.width = 800, ani.height = 600)
file.copy(from = "raw_data.gif", to = "figures/raw_data.gif")
file.remove("raw_data.gif") # doing this because saveGIF isn't letting me put the file into another directory

### exploratory climate metrics
hist(trends$seasonality, 100)
hist(trends$mean_temp, 100)
plot(seasonality ~ lat, data=sites, col=palette[as.factor(dataset)], pch=20,
     ylab="Seasonality (max temp - min temp)", xlab="Latitude")
legend("topleft", inset = 0.05, legend=levels(as.factor(trends$dataset)),
       col = palette, pch=20, title="Dataset", cex=0.8)

plot(mean_temp ~ lat, data=sites, col=palette[as.factor(dataset)], pch=20,
     ylab="Yearly mean temperature", xlab="Latitude")
legend("topright", inset = 0.05, legend=levels(as.factor(all_data$dataset)),
       col = palette, pch=20, title="Dataset")

plot(seasonality ~ mean_temp, data=sites, col=palette[as.factor(dataset)], pch=20,
     ylab="Seasonality (max temp - min temp)", xlab="Mean temperature")
legend("topright", inset = 0.05, legend=levels(as.factor(trends$dataset)),
       col = palette, pch=20, title="Dataset", cex=0.8)

### world map
png("figures/world_map.png", width=700, height=500)  
# svg("figures/world_map.svg", width=7, height=5) 
main <- ggplot() + geom_sf(data = world) +
  geom_sf(data = sites_sf, aes(col = seasonality), show.legend = "T") + # seasonality in aes()
  coord_sf(xlim = c(-20, 173), ylim = c(10, 73), expand = FALSE) +
  scale_color_gradient(low = "orange", high = "#85039c") + # let ggplot make the scale for you
  theme_void(base_rect_size = 1) +
  theme(legend.position = c(0.93, 0.27),
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) # move the legend around here
inset <- ggplot() + geom_sf(data = world) +
  geom_sf(data = sites_sf, aes(col = seasonality), show.legend = "T") +
  coord_sf(xlim = c(-120, -60), ylim = c(15, 50), expand = FALSE) +
  theme_void(base_rect_size = 1) + scale_color_gradient(low = "orange", high = "#85039c") +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none",
        panel.grid.major = element_line(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.background = element_rect(fill="white")) # probably want to remove inset legend?
full = ggdraw() + draw_plot(main) + draw_plot(inset, x = 0, y = 0.191, width = 0.25, height = 0.25)
full # I hate ggplot
dev.off()


### timeseries figure to go with map
pheno_group_pallete <- brewer.pal(length(levels(as.factor(trends_focus$pheno_group))), "Dark2")
all_data_plot <- all_data
all_data_plot$col <- pheno_group_pallete[as.factor(all_data_plot$phenophase)]
pheno_group_pallete <- brewer.pal(length(levels(as.factor(trends_focus$pheno_group))), "Dark2")

# long version
png("figures/trends_summary_long.png", width=4, height=8, units="in", res=100)
# svg("figures/trends_summary_long.svg", width=4, height=5)
plot(1,type = "n", yaxt = "n", xaxt="n",
     ylim = c(1955,2020), xlim = c(0,max(trends_summary$n_species)*1.1),
     ylab = "", xlab = "Number of species",
     bty="l")

axis(2, seq(1960,2020,by=10), TRUE, las=2)
axis(1, seq(0,max(trends_summary$n_species)*1.1, by=200), seq(0,max(trends_summary$n_species)*1.1, by=200),
     las=2)
trends_summary[, arrows(n_species, mean_min, n_species, mean_max,
                        code=3, angle=90, length=0.07, lwd=4,
                        col=col)]
trends_summary[, segments(rep(n_species,2), c(min25, max25),
                          rep(n_species,2), c(min75, max75),
                          lwd=16, col=adjustcolor(col, 0.5))]
trends_summary[, text(n_species-max(trends_summary$n_species)*0.06, mean(c(mean_min, mean_max))+4,
                      group_name, col=col, srt=90)]
box(bty="l",lwd=2)
summary_plot <- recordPlot()
dev.off()

### map and data summary figures together
png("figures/map_and_summary.png", width=900, height=400)
plot_grid(full, summary_plot,
          ncol=2, rel_widths = c(7,3))
dev.off()


### Phenogroup histograms
# shading not-significant trends
#png("figures/pheno_group_hist_sigma_year_signif.png", width=600, height=400)
png("figures/pheno_group_hist_sigma_year_signif.png", width=900, height=600, res=100)
shading_opacity <- 0.7
shading_density <- 30

main_hist <- quote(mean_slope_temp) # you can change which metric gets plotted here
sub_hist <- quote(sigma_slope_temp)
main_hist_p <- quote(mean_p_temp)
sub_hist_p <- quote(sigma_p_temp)
#sub_xlim <- c(-1,1) # if sub_hist is sigma_slope
sub_xlim <- c(-4,4) # if sub_hist is sigma_slope_temp
sub_breaks1 <- 150 # set to 100 for sigma_slope
sub_breaks2 <- 100
sub_breaks3 <- 50

par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
hist(trends[pheno_group == "flower", eval(main_hist)], 100, col=pheno_group_pallete[2], border=pheno_group_pallete[2], xlim=c(-10,12), main="", xlab="Mean sensitivity (days/degree C)", ylab="Number of time-series")
hist(trends[pheno_group == "flower" & eval(main_hist_p) > 0.05, eval(main_hist)], 100, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "leaf", eval(main_hist)], 100, col=pheno_group_pallete[1], border=pheno_group_pallete[1], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "leaf" & eval(main_hist_p) > 0.05, eval(main_hist)], 100, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "bird", eval(main_hist)], 100, col=pheno_group_pallete[3], border=pheno_group_pallete[3], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "bird" & eval(main_hist_p) > 0.05, eval(main_hist)], 100, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "insect", eval(main_hist)], 100, col=pheno_group_pallete[4], border=pheno_group_pallete[4], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "insect" & eval(main_hist_p) > 0.05, eval(main_hist)], 100, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
abline(v=0, lwd=1, lty=2)
mtext("  a", side=3, line=-2, adj=0)
legend("left", inset=.02, c("Flowers","Leaves","Birds","Insects"), cex=0.8,
       fill=pheno_group_pallete[c(2,1,3,4)],
       bty="n", border="white")
#report.outliers(trends$mean_slope_temp, "x", zero=F)

par(fig = c(0.5,1, 0.3, 1), new = T, mgp=c(3,0.7,0)) 
hist(trends[pheno_group == "flower", eval(sub_hist)], sub_breaks1, col=pheno_group_pallete[2], border=pheno_group_pallete[2], main="", xlab="Change in variance (days/year)", xlim=sub_xlim, ylab="", ann=F)
hist(trends[pheno_group == "flower" & eval(sub_hist_p) > 0.05, eval(sub_hist)], sub_breaks2, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "leaf", eval(sub_hist)], sub_breaks2, col=pheno_group_pallete[1], border=pheno_group_pallete[1], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "leaf" & eval(sub_hist_p) > 0.05, eval(sub_hist)], sub_breaks2, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "bird", eval(sub_hist)], sub_breaks2, col=pheno_group_pallete[3], border=pheno_group_pallete[3], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "bird" & eval(sub_hist_p) > 0.05, eval(sub_hist)], sub_breaks3, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
hist(trends[pheno_group == "insect", eval(sub_hist)], sub_breaks2, col=pheno_group_pallete[4], border=pheno_group_pallete[4], add=T, xlim=c(-10,5))
hist(trends[pheno_group == "insect" & eval(sub_hist_p) > 0.05, eval(sub_hist)], sub_breaks2, border=F, add=T,
     density = shading_density, angle=45, col=adjustcolor("white",shading_opacity))
mtext(side = 1, text = "Variance sensitivity (days/degree C)", line = 2)
mtext("  b", side=3, line=-1.5, adj=0)
abline(v=0, lwd=1, lty=2)
#report.outliers(trends[, eval(sub_hist)], "x", zero=F)
dev.off()


### year vs. temp shift - colored by pheno_position w/ contours and/or hex bins

png("figures/mean_scatter_pheno_position_contours.png", width=600, height=500)
#svg("figures/year_vs_temp_pheno_position_contours.svg", width=7, height=6)

x_var <- quote(mean_slope_temp)
y_var <- quote(mean_slope)
xlim <- c(-10,7); ylim <- c(-2,1.8) # for mean_slope~mean_slope_temp
#xlim <- c(-5,5); ylim <- c(-1.5,1.5)

blue <- "#0044f2"; red <- "#e6073b"; purple <- "#9333d4"  # original colors
blue <- "#0031b0"; red <- "#ff2b2b"; purple <- "#7e2bd6" # higher contrast colors
color_reps <- 6
color_breaks <- 200
pal <- colorRampPalette(c(rep(blue, color_reps), adjustcolor(purple, alpha.f = 0.7), rep(red, color_reps)),
                        alpha=T)
c_pal <- colorRampPalette(c(rep(blue, color_reps), adjustcolor(purple, alpha.f = 1), rep(red, color_reps)),
                          alpha=T)

# pal <- colorRampPalette(adjustcolor(c(rep(brewer.pal(11, "RdYlBu")[11], color_reps),
#                                       brewer.pal(11, "RdYlBu")[11:1],
#                                       rep(brewer.pal(11, "RdYlBu")[1], color_reps)), 0.7))
# c_pal <- colorRampPalette(adjustcolor(c(rep(brewer.pal(11, "RdYlBu")[11], color_reps),
#                                       brewer.pal(11, "RdYlBu")[11:1],
#                                       rep(brewer.pal(11, "RdYlBu")[1], color_reps)), 1))


trends_plot <- trends_focus
breaks <- seq(-max(abs(trends_plot[, pheno_position])), max(abs(trends_plot[, pheno_position])), length.out = color_breaks) # so that color ramp is centered at zero
trends_plot$col <- adjustcolor(pal(color_breaks)[as.numeric(cut(trends_plot[, pheno_position],breaks = breaks))],
                               alpha.f = 0.6)
layout(matrix(1:2, nrow=1), widths = c(0.8,0.2))
old_mar <- par("mar")
par(mar = c(5,4,3,0), xpd=FALSE)

plot(eval(y_var) ~ eval(x_var), data=trends_focus,
     pch=16, col=trends_plot[, col], cex=0.35,
     xlab="Variance sensitivity (days/degree C)", ylab="Variance shift (days/year)",
     xlim=xlim, ylim=ylim, axes=F)
axis(1, at=seq(-10,7,by=2), labels=seq(-10,7,by=2))
axis(2, at=seq(-2,1.8,by=0.5), labels=seq(-2,1.8,by=0.5))
box(which = "plot", lty = "solid", bty="l", lwd=1.5)

kd <- kde2d(trends_focus[, eval(x_var)], trends_focus[, eval(y_var)],
            h = c(1.5,0.5),n=100) #kernel density for contours
c_lines <- contourLines(kd, levels = c(0.003,0.01,0.03,0.1,0.25))

# find average pheno position color in a range around a give point - used to define color of contour segment
find.avg.col <- function(x_val, y_val, x_width=0.5, y_width=0.25){
  rel_data <- trends_plot[eval(x_var) > x_val-x_width & eval(x_var) < x_val+x_width & eval(y_var) > y_val-y_width & eval(y_var) < y_val+y_width,
                          pheno_position]
  rel_val <- mean(rel_data)
  c_col <- c_pal(color_breaks)[which.min(abs(breaks - rel_val))]
  return(c_col)
}

# plot contours colored according to the average pheno_position around them
plot.colored.contours <- function(c_list, c_cols){
  segments(c_list$x[-length(c_list$x)],
           c_list$y[-length(c_list$y)],
           c_list$x[-1L],
           c_list$y[-1L],
           col=c_cols, lwd=4)
}

c_cols <- lapply(c_lines, function(c_line) mapply(find.avg.col, x_val=c_line$x, y_val=c_line$y))
lapply(c_lines, lines, col="white", lwd=6)
mapply(plot.colored.contours, c_list=c_lines, c_cols=c_cols)

report.outliers(trends_focus[, eval(y_var)], "y", just_report=T)
report.outliers(trends_focus[, eval(x_var)], "x", just_report=T)
model <- lm(eval(y_var) ~ eval(x_var), data=trends_focus)
summary(model)
abline(v=0, lty=2)
abline(h=0, lty=2)
abline(model, lwd=2)

mtext(paste(" ",nrow(trends_focus[eval(y_var) < 0 & eval(x_var) < 0,])),
      1, line=-1.2, adj=0)
mtext(paste(" ",nrow(trends_focus[eval(y_var) > 0 & eval(x_var) < 0,])),
      3, line=-1.2, adj=0)
mtext(paste(" ",nrow(trends_focus[eval(y_var) < 0 & eval(x_var) > 0,]), " "),
      1, line=-1.2, adj=1)
mtext(paste(nrow(trends_focus[eval(y_var) > 0 & eval(x_var) > 0,]), " "),
      3, line=-1.2, adj=1)
par(mar = c(7,3.5,7,3), xpd=NA)

legend_range <- c(-80,80)
legend_pal_range <- Closest(breaks, legend_range, which=TRUE)
legend_labels <- seq(legend_range[1],legend_range[2],by=20); legend_labels[1] <- -150; legend_labels[length(legend_labels)] <- 150
legend_axis <- list(at = seq(legend_range[1],legend_range[2],by=20),
                    labels = legend_labels,
                    las = 2, cex.axis=0.8)
legend.scale(legend_range, col = pal(color_breaks)[c(legend_pal_range[1]:legend_pal_range[2])],
             horizontal=F, axis.args = legend_axis)
axis(2, at=0, labels="Phenological position (days)", tick=F, cex.axis=0.8)
scale_break <- readPNG("figures/scale_break.png")
rasterImage(scale_break, -0.2,67,1.2,73)
rasterImage(scale_break, -0.2,-73,1.2,-67)
par(mar = old_mar)
dev.off()


### Main models - trend viz panels
png("figures/year_mean_panels.png", width=1000, height = 800)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
plot <- visreg(year_mean_model, "pheno_group", ylab="Predicted shift w.r.t. time", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(year_mean_model, "seasonality", ylab="Predicted shift w.r.t. time", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mar=c(4,1,1,2))
plot <- visreg(year_mean_model, "mean_temp", ylab="Predicted rate of shift", points.par=list(col="#9cbbff"), yaxt="n")
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(year_mean_model, "pheno_position", ylab="Predicted rate of shift", points.par=list(col="#9cbbff"), yaxt="n")
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()

png("figures/temp_mean_panels.png", width=1000, height = 800)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
plot <- visreg(temp_mean_model, "pheno_group", ylab="Predicted shift w.r.t. temperature", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(temp_mean_model, "seasonality", ylab="Predicted shift w.r.t. temperature", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mar=c(4,1,1,2))
plot <- visreg(temp_mean_model, "mean_temp", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(temp_mean_model, "pheno_position", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()

png("figures/year_sigma_panels.png", width=1000, height = 800)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
plot <- visreg(year_sigma_model, "pheno_group", ylab="Predicted variance change w.r.t. time", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(year_sigma_model, "seasonality", ylab="Predicted variance change w.r.t. time", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mar=c(4,1,1,2))
plot <- visreg(year_sigma_model, "mean_temp", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(year_sigma_model, "pheno_position", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#9cbbff"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()


png("figures/temp_sigma_panels.png", width=1000, height = 800)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
plot <- visreg(temp_sigma_model, "pheno_group", ylab="Predicted variance change w.r.t. temperature", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(temp_sigma_model, "seasonality", ylab="Predicted variance change w.r.t. temperature", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mar=c(4,1,1,2))
plot <- visreg(temp_sigma_model, "mean_temp", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
plot <- visreg(temp_sigma_model, "pheno_position", ylab="Predicted rate of shift", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()


### standardized effects size figure
arrows.label <- function(label="",col="black",angle=45,length=0.1,i=0,letter=FALSE, cex=0.75){
  extremes = par("usr")
  xx = (extremes[2] - extremes[1])/16
  yy = (extremes[4] - extremes[3])/16
  
  if(!isFALSE(letter)){
    length=0
    text(x = extremes[1]+xx, extremes[4]-2*yy-i, labels = substitute(italic(bold(letter))), col=col)
  }
  
  arrows(x0 = extremes[1]+xx, x1 = extremes[1]+2*xx,
         y0 = extremes[4]-2*yy-i, y1 = extremes[4]-2*yy-i,
         col = adjustcolor(col,0.8), angle=angle, length=length, code = 1, lwd=2)
  
  text(x = extremes[1]+2*xx,
       y = extremes[4]-2*yy-i,
       labels = label, pos=4, cex=cex)
}

coefs <- data.table(type = rep(c("cat", "cont"),times=c(16,12)),
                    variable = rep(c("leaf", "flower", "insect", "bird", "seasonality", "mean_temp", "pheno_position"), each=4),
                    term = c("mean", "mean", "sigma", "sigma"),
                    explan = c("year", "temp"),
                    coef = rep(999, 28),
                    se = rep(999, 28),
                    p = rep(999, 28))
coefs[term=="mean" & explan=="year", "coef"] <- as.numeric(summary(year_mean_model)$coefficients[,1])
coefs[term=="mean" & explan=="temp", "coef"] <- as.numeric(summary(temp_mean_model)$coefficients[,1])
coefs[term=="sigma" & explan=="year", "coef"] <- as.numeric(summary(year_sigma_model)$coefficients[,1])
coefs[term=="sigma" & explan=="temp", "coef"] <- as.numeric(summary(temp_sigma_model)$coefficients[,1])
coefs[term=="mean" & explan=="year", "se"] <- as.numeric(summary(year_mean_model)$coefficients[,2])
coefs[term=="mean" & explan=="temp", "se"] <- as.numeric(summary(temp_mean_model)$coefficients[,2])
coefs[term=="sigma" & explan=="year", "se"] <- as.numeric(summary(year_sigma_model)$coefficients[,2])
coefs[term=="sigma" & explan=="temp", "se"] <- as.numeric(summary(temp_sigma_model)$coefficients[,2])
coefs[term=="mean" & explan=="year", "p"] <- as.numeric(summary(year_mean_model)$coefficients[,5])
coefs[term=="mean" & explan=="temp", "p"] <- as.numeric(summary(temp_mean_model)$coefficients[,5])
coefs[term=="sigma" & explan=="year", "p"] <- as.numeric(summary(year_sigma_model)$coefficients[,5])
coefs[term=="sigma" & explan=="temp", "p"] <- as.numeric(summary(temp_sigma_model)$coefficients[,5])
coefs$p <- round(coefs$p, 3)
padding <- 2
label_padding <- 0.11
x_lim <- c(min(coefs$coef - 2*coefs$se),max(coefs$coef + 2*coefs$se))
y_lim_top <- c(0 - padding, nrow(coefs[type == "cont"]) + padding)
y_lim_bottom <- c(0 - padding, nrow(coefs[type == "cat"]) + padding+1)

cont_coefs <- coefs[type == "cont",]
unit <- nrow(cont_coefs)*0.075
cont_coefs[term=="mean", "ypos"] <- rep(seq(0, nrow(cont_coefs), length.out = 3), each=2) + rep(c(0*unit,2*unit), 3)
cont_coefs[term=="sigma", "ypos"] <- rep(seq(0, nrow(cont_coefs), length.out = 3), each=2) + rep(c(-1*unit,1*unit), 3)

cat_coefs <- coefs[type == "cat",]
unit <- nrow(cat_coefs)*0.05
cat_coefs[term=="mean", "ypos"] <- rep(seq(0, nrow(cat_coefs), length.out = 4), each=2) + rep(c(0*unit,2*unit), 4)
cat_coefs[term=="sigma", "ypos"] <- rep(seq(0, nrow(cat_coefs), length.out = 4), each=2) + rep(c(-1*unit,1*unit), 4)


png("figures/standardized_effects.png", width=400, height=600)
#svg("figures/standardized_effects.svg", width=5, height=8)
# continuous
layout(matrix(c(1,2), nrow=2), heights=c(y_lim_top[2], y_lim_bottom[2]))
par(mar=c(0, 4.1, 2.1, 2.1))
plot(1,type = "n", yaxt = "n", xaxt="n",
     xlim = x_lim, ylim = y_lim_top,
     xlab = "", ylab = "",
     bty="l")
text(x = par("usr")[1]-label_padding, y=seq(2,nrow(cont_coefs)+3, length.out=3),
     labels=c("Seasonality", "Regional\ntemperature  ", "Phenological\nposition  "),
     xpd=NA, srt=60, adj=1)
mtext(" Continuous", 1, line=-1.2, adj=0)
segments(x0 = cont_coefs$coef - 2*cont_coefs$se, x1 = cont_coefs$coef + 2*cont_coefs$se,
         y0 = cont_coefs$ypos,
         lwd = 10, lend=2, col=adjustcolor(c("blue", "orange"), alpha.f = 0.4))
arrows(x0 = rep(0,nrow(cont_coefs[term=="mean",])), y0 = cont_coefs[term=="mean",ypos],
       x1 = cont_coefs[term=="mean",coef], y1 = cont_coefs[term=="mean",ypos],
       length=0, col=c("#4054c7", "orange"), angle=45, lwd=2) # mean shifts
text(x = cont_coefs[term=="mean",coef], cont_coefs[term=="mean",ypos],
     labels = substitute(italic(bold("μ"))), col=c("blue", "#e69100"))
arrows(x0 = rep(0,nrow(cont_coefs[term=="sigma",])), y0 = cont_coefs[term=="sigma",ypos],
       x1 = cont_coefs[term=="sigma",coef], y1 = cont_coefs[term=="sigma",ypos],
       length=0, col=c("#4054c7", "orange"), angle=90, lwd=2) # sigma shifts
text(x = cont_coefs[term=="sigma",coef], cont_coefs[term=="sigma",ypos],
     labels = substitute(italic(bold("σ"))), col=c("blue", "#e69100"))
text(x = c(cont_coefs[term=="mean",coef], cont_coefs[term=="sigma",coef])+label_padding,
     y = c(cont_coefs[term=="mean",ypos], cont_coefs[term=="sigma",ypos])+label_padding,
     labels= sapply(c(cont_coefs[term=="mean",p], cont_coefs[term=="sigma",p]), function(x) ifelse(x < 0.01, "*", "")),
     cex=1.5) # significants at a=0.01
abline(v=0, lty=2)
box(which = "plot", lty = "solid", bty="l", lwd=1.5)

#legend
arrows.label("Mean sensitivity", "orange", 45, 0.1, 0, letter="μ", cex=0.85)
arrows.label("Variance sensitivity", "orange", 90, 0.07, 1.5, letter="σ", cex=0.85)
arrows.label("Mean shift", "blue", 45, 0.1, 3, letter="μ", cex=0.85)
arrows.label("Variance shift", "blue", 90, 0.07, 4.5, letter="σ", cex=0.85)


#categorical
par(mar=c(4.1, 4.1, 0, 2.1))
plot(1,type = "n", yaxt = "n",
     xlim = x_lim, ylim = y_lim_bottom,
     xlab = "Standardized effect size", ylab = "",
     bty="l")
text(x = par("usr")[1]-label_padding, y=seq(2,nrow(cat_coefs)+1, length.out=4),
     labels=c("Leaves", "Flowers", "Insects", "Birds"),
     xpd=NA, srt=60, adj=1)
mtext(" Categorical", 1, line=-1.2, adj=0)
segments(x0 = cat_coefs$coef - 2*cat_coefs$se, x1 = cat_coefs$coef + 2*cat_coefs$se,
         y0 = cat_coefs$ypos,
         lwd = 10, lend=2, col=adjustcolor(c("blue", "orange"), alpha.f = 0.4))
arrows(x0 = rep(0,nrow(cat_coefs[term=="mean",])), y0 = cat_coefs[term=="mean",ypos],
       x1 = cat_coefs[term=="mean",coef], y1 = cat_coefs[term=="mean",ypos],
       length=0, col=c("#4054c7", "orange"), angle=45, lwd=2) # mean shifts
text(x = cat_coefs[term=="mean",coef], cat_coefs[term=="mean",ypos],
     labels = substitute(italic(bold("μ"))), col=c("blue", "#e69100"))
arrows(x0 = rep(0,nrow(cat_coefs[term=="sigma",])), y0 = cat_coefs[term=="sigma",ypos],
       x1 = cat_coefs[term=="sigma",coef], y1 = cat_coefs[term=="sigma",ypos],
       length=0, col=c("#4054c7", "orange"), angle=90, lwd=2) # sigma shifts
text(x = cat_coefs[term=="sigma",coef], cat_coefs[term=="sigma",ypos],
     labels = substitute(italic(bold("σ"))), col=c("blue", "#e69100"))
text(x = c(cat_coefs[term=="mean",coef], cat_coefs[term=="sigma",coef])-label_padding,
     y = c(cat_coefs[term=="mean",ypos], cat_coefs[term=="sigma",ypos])-label_padding,
     labels= sapply(c(cat_coefs[term=="mean",p], cat_coefs[term=="sigma",p]), function(x) ifelse(x < 0.01, "*", "")),
     cex=1.5) # significants at a=0.01
abline(v=0, lty=2)
box(which = "plot", lty = "solid", bty="l", lwd=1.5)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()


### Plant traits

png("figures/flowers_year_mean.png", width=700, height=500)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
visreg(flower_traits_year_mean, "growth_form", overlay=FALSE,
       ylab="Shift in mean over time", main="")
abline(h=0, lty=2)

visreg(flower_traits_year_mean, "height", overlay=FALSE,
       ylab="Shift in mean over time", main="", xlab="log(height)")
abline(h=0, lty=2)
par(mar=c(4,1,1,2))
visreg(flower_traits_year_mean, "seed_mass", overlay=FALSE,
       ylab="Shift in mean over time", main="", xlab="log(seed_mass)", yaxt="n")
abline(h=0, lty=2)
plot <- visreg(flower_traits_year_mean, "leaf_area_per_mass", overlay=FALSE,
               ylab="Shift in mean over time", main="", xlab="log(leaf_area_per_mass)", yaxt="n")
report.outliers(plot$res$visregRes, "y", zero=T)
#mtext(" *", 3, line=-1.8, adj=0, cex=1.5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()

png("figures/flowers_temp_mean.png", width=700, height=500)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
visreg(flower_traits_temp_mean, "growth_form", overlay=FALSE,
       ylab="Shift in mean over temp", main="", points.par=list(col="#f0a21d"))
mtext(" *", 3, line=-1.8, adj=0, cex=1.5)
abline(h=0, lty=2)

visreg(flower_traits_temp_mean, "height", overlay=FALSE,
       ylab="Shift in mean over temp", main="", xlab="log(height)", points.par=list(col="#f0a21d"))
abline(h=0, lty=2)
par(mar=c(4,1,1,2))
visreg(flower_traits_temp_mean, "seed_mass", overlay=FALSE,
       ylab="Shift in mean over temp", main="", xlab="log(seed_mass)", yaxt="n", points.par=list(col="#f0a21d"))
abline(h=0, lty=2)
plot <- visreg(flower_traits_temp_mean, "leaf_area_per_mass", overlay=FALSE,
               ylab="Shift in mean over temp", main="", xlab="log(leaf_area_per_mass)", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()

png("figures/flowers_year_sigma.png", width=700, height=500)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
visreg(flower_traits_year_sigma, "growth_form", overlay=FALSE,
       ylab="Shift in variance over time", main="")
abline(h=0, lty=2)

visreg(flower_traits_year_sigma, "height", overlay=FALSE,
       ylab="Shift in variance over time", main="", xlab="log(height)")
abline(h=0, lty=2)
par(mar=c(4,1,1,2))
visreg(flower_traits_year_sigma, "seed_mass", overlay=FALSE,
       ylab="Shift in variance over time", main="", xlab="log(seed_mass)", yaxt="n")
abline(h=0, lty=2)
plot <- visreg(flower_traits_year_sigma, "leaf_area_per_mass", overlay=FALSE,
               ylab="Shift in variance over time", main="", xlab="log(leaf_area_per_mass)", yaxt="n")
report.outliers(plot$res$visregRes, "y", zero=T)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()

png("figures/flowers_temp_sigma.png", width=700, height=500)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,1,0))
visreg(flower_traits_temp_sigma, "growth_form", overlay=FALSE,
       ylab="Shift in variance over temp", main="", points.par=list(col="#f0a21d"))
abline(h=0, lty=2)

visreg(flower_traits_temp_sigma, "height", overlay=FALSE,
       ylab="Shift in variance over temp", main="", xlab="log(height)", points.par=list(col="#f0a21d"))
abline(h=0, lty=2)
par(mar=c(4,1,1,2))
visreg(flower_traits_temp_sigma, "seed_mass", overlay=FALSE,
       ylab="Shift in variance over temp", main="", xlab="log(seed_mass)", yaxt="n", points.par=list(col="#f0a21d"))
abline(h=0, lty=2)
plot <- visreg(flower_traits_temp_sigma, "leaf_area_per_mass", overlay=FALSE,
               ylab="Shift in variance over temp", main="", xlab="log(leaf_area_per_mass)", yaxt="n", points.par=list(col="#f0a21d"))
report.outliers(plot$res$visregRes, "y", zero=T)
# mtext(" *", 3, line=-1.8, adj=0, cex=1.5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()


### Plant phenophase and growth form

png("figures/leaves_flowers_growth_form.png", width=800, height=600)
layout(matrix(1:4, nrow=2, byrow = F))
par(mar=c(4,4.1,2,1))
visreg(plant_phase_year_mean, "phenophase", by="growth_form", overlay=TRUE, points.par=list(col=c("red","green",adjustcolor("blue", alpha.f = 0.3))),
       ylab="Predicted mean shift over time")
visreg(plant_phase_temp_mean, "phenophase", by="growth_form", overlay=TRUE, points.par=list(col=c("red","green",adjustcolor("blue", alpha.f = 0.3))),
       ylab="Predicted mean shift over temp")
visreg(plant_phase_year_sigma, "phenophase", by="growth_form", overlay=TRUE, points.par=list(col=c("red","green",adjustcolor("blue", alpha.f = 0.3))),
       ylab="Predicted variance shift over time")
visreg(plant_phase_temp_sigma, "phenophase", by="growth_form", overlay=TRUE, points.par=list(col=c("red","green",adjustcolor("blue", alpha.f = 0.3))),
       ylab="Predicted variance shift over temp")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 3.1, 2.1))
dev.off()


### Bird traits

png("figures/bird_traits.png", width=800, height=600)
par(mfrow=c(2,2))
visreg(bird_traits_model_mean_time, "diet", ylab="Shift in first occurence mean over time", main="Mean vs. temp",
       points.par=list(col=rgb(0, 0, 1, 0.7))); abline(h=0, lty=2)
visreg(bird_traits_model_mean_temp, "diet", ylab="Shift in first occurence mean over temperature", main="Mean vs. temp",
       points.par=list(col=rgb(1, 0, 0, 0.7))); abline(h=0, lty=2)
visreg(bird_traits_model_sigma_time, "diet", ylab="Shift in first occurence variance over time", main="Sigma vs. time",
       points.par=list(col=rgb(0, 0.4, 1, 0.7))); abline(h=0, lty=2)
visreg(bird_traits_model_sigma_temp, "diet", ylab="Shift in first occurence variance over temperature", main="Sigma vs. temp",
       points.par=list(col=rgb(1, 0.4, 0, 0.7))); abline(h=0, lty=2)
par(mfrow=c(1,1))
dev.off()






### MISC plots

# is Rothamsted causing the seasonality effect?
png("figures/rothamsted_seasonality.png", width=600, height=500)
plot(mean_slope_temp ~ seasonality, data=trends_focus,
     pch=20, col=adjustcolor(palette[as.factor(dataset)],alpha.f = 0.5), cex=0.5,
     xlab="Seasonality (Degrees C)", ylab="Shift in response to increased temperature")
abline(h=0, lty=2)
model_full <- lm(mean_slope_temp ~ seasonality, data=trends_focus)
abline(model_full, lwd=2, col=rgb(0.1,0.1,0.1), lty=3)
model_no_roth <- lm(mean_slope_temp ~ seasonality, data=trends_focus[dataset != "rothamsted",])
abline(model_no_roth, lwd=2, col=rgb(0,0.2,1))
legend("bottomright", inset = 0.05, legend=levels(as.factor(trends_focus$dataset)),
       col = palette, pch=20, title="Dataset", cex=0.8)
dev.off()

# seasonality trends by excluded dataset
png("figures/datasets_seasonality.png", width=600, height=600)
plot(mean_slope_temp ~ seasonality, data=trends_focus,
     pch=as.numeric(as.factor(dataset)), col=adjustcolor(palette[as.factor(dataset)],alpha.f = 0.5), cex=0.8,
     xlab="Seasonality (degrees C)", ylab="Mean sensitivity (days/degree C)")
abline(h=0, lty=2)
model_full <- lm(mean_slope_temp ~ seasonality, data=trends_focus)
counter <- 1
for(dset in sort(unique(trends_focus$dataset))){
  dataset_model <- lm(mean_slope_temp ~ seasonality, data=trends_focus[dataset != dset,])
  abline(dataset_model, lwd=2, col=palette[counter])
  counter <- counter + 1
}
legend("bottomright", inset = 0.05, legend=levels(as.factor(trends_focus$dataset)),
       col = palette, pch=1:8, title="Dataset", cex=0.8)
#abline(model_full, lwd=3, col="black", lty=1)
dev.off()

# seasonality over latitude by dataset
png("figures/seasonality_lat.png", width=600, height=400)
plot(seasonality ~ lat, data=sites, col=palette[as.factor(dataset)], pch=as.numeric(as.factor(dataset)),
     ylab="Seasonality (max temp - min temp)", xlab="Latitude")
legend("topleft", inset = 0.05, legend=levels(as.factor(trends$dataset)),
       col = palette, pch=1:8, title="Dataset", cex=0.8)
dev.off()

beep(4)

### cleanup

rm(all_data, datasets, bird_trait_data, countries, model, model_full, model_no_roth, clim_mean, clim_sigma)
gc()
