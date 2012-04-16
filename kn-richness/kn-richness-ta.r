######################
# SAA2012 
# 
# Analysis of trait richness in Ewens samples (K_n) from Wright-Fisher and conformist CT models
# under time-averaging.
#
# Mark E. Madsen (c) 2012 All rights reserved.  
# Software and Analysis for Mark E. Madsen SAA 2012 Papers is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# Based on a work at github.com.  See license.html for link to CC 3.0 license
#
######################

setwd("~/Dropbox/Research/Dissertation Project/analysis/saa2012/kn-richness")

# epicalc contains a version of "aggregate" that takes multiple analysis functions
#library("epicalc")
library("ggplot2")
library("plyr")
library("foreach")
library("xtable")

# read in data sets
#rich <- read.delim("../rawdata/trait-richness-wf-conf-5000.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer"))
rich.longtheta <- read.delim("../rawdata/wf-0.1-100-e100-combined-richness.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer"))
#rich.40ktheta <- read.delim("../rawdata/kn-richness-0.1-100-40k-100.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer"))
lifetime.longtheta <- read.delim("../rawdata/wf-0.1-100-e100-combined-lifetime.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","numeric","numeric","numeric","numeric","numeric","integer"))
#lifetime.40ktheta <- read.delim("../rawdata/lifetime-0.1-100-40k-100.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","numeric","numeric","numeric","numeric","numeric","integer"))


# determine mean lifetime for each theta value (windowsize doesn't figure into this)
lifetime.longtheta.stat <- ddply(lifetime.longtheta, .(Theta), summarise, mean_lifetime = mean(MeanLifetime), stdev_lifetime = sd(MeanLifetime))
#lifetime.40ktheta.stat <- ddply(lifetime.40ktheta, .(Theta), summarise, mean_lifetime = mean(MeanLifetime))

# create factors for plotting etc
#rich$WinFactor <- factor(rich$Windowsize)
#rich$ThetaFactor <- factor(rich$Theta)
#rich.40ktheta$ThetaFactor <- factor(rich.40ktheta$Theta)
#rich.40ktheta$WinFactor <- factor(rich.40ktheta$Windowsize)
#rich.wf <- subset(rich, Model == "WrightFisherDrift")
rich.longtheta$ThetaFactor <- factor(rich.longtheta$Theta)
rich.longtheta$WinFactor <- factor(rich.longtheta$Windowsize)

# Create data frame with mean value of k_n for each combination of windowsize and theta
#rich.wf.stat <- ddply(rich.wf, .(Theta, Windowsize), summarise, mean_kn = mean(Kn), sd_kn = sd(Kn), num = length(Kn))
rich.longtheta.stat <- ddply(rich.longtheta, .(Theta, Windowsize), summarise, mean_kn = mean(Kn), sd_kn = sd(Kn), num = length(Kn))
rich.longtheta.stat$ThetaFactor <- factor(rich.longtheta.stat$Theta)
#rich.40ktheta.stat <- ddply(rich.40ktheta, .(Theta, Windowsize), summarise, mean_kn = mean(Kn), sd_kn = sd(Kn), num = length(Kn))
#rich.40ktheta.stat$ThetaFactor <- factor(rich.40ktheta.stat$Theta)

#################### Plots of Raw Kn Data #########################
# Plot a line graph of mean Kn for a given theta level, versus log(windowsize) on the x-axis
# qplot(log(Windowsize), mean_kn, data=rich.wf.stat, main="Mean Kn vs Windowsize by Theta 2-40 for 50K run length", facets = Theta ~ ., geom=c("point", "smooth"))
# ggsave(file = "../images/mean-kn-vs-logwindowsize-bytheta-2-40-50k.pdf")

qplot(Windowsize, mean_kn, data=rich.longtheta.stat, main="Mean Kn vs Windowsize by Theta 0.1-100", facets = Theta ~ ., geom=c("point", "smooth")) + theme_bw() + scale_x_log10()
ggsave(file = "../images/kn-0.1-100-e100-combined-mean-kn-vs-logwindowsize-bytheta.pdf")

#qplot(log(Windowsize), mean_kn, data=rich.40ktheta.stat, main="Mean Kn vs Windowsize by Theta 0.1-100 for 40K run length", facets = Theta ~ ., geom=c("point", "smooth"))
#ggsave(file = "../images/mean-kn-vs-logwindowsize-bytheta-0.1-100-40k.pdf")


# Create a set of histogram stacks, one for each theta level, showing the distribution of Kn by window size
# Save each plot to PDF
subset_and_plot <- function(x, frame, filebase, titlebase, runlengthstring) {
  print(x)
  data.subset <- subset(frame, ThetaFactor == x)
  title <- paste(titlebase, x, runlengthstring)
  #qplot(Kn, ..density.., data = data.subset, main=title, facets = Windowsize ~., geom = "histogram", binwidth=1 )
  knplot <- qplot(Kn, ..density.., data = data.subset, main = title, geom="freqpoly", binwidth = 3, color = WinFactor) + theme_bw()
  #knplot + scale_colour_brewer(type="seq")
  fname <- paste("../images/", filebase, x, ".pdf", sep="")
  ggsave(filename = fname)          
}

#foreach(lvl=levels(rich.wf$ThetaFactor)) %do% subset_and_plot(lvl, rich.wf, "kn-2-40-freqpoly-bywindowsize-fortheta-", "Kn Distribution by Windowsize for Theta", "50K Run Length")
foreach(lvl=levels(rich.longtheta$ThetaFactor)) %do% subset_and_plot(lvl, rich.longtheta, "kn-0.1-100-e100-combined-freqpoly-bywindowsize-fortheta-5k-40k-", "Kn Distribution by Windowsize for Theta", "5K-40K Run Length")

# Plot the mean trait lifetime by theta value (this really ought to be the same for every single data set)
errorbars_lifetime <- aes(ymax = mean_lifetime + (stdev_lifetime), ymin = mean_lifetime - (stdev_lifetime) )
qplot(Theta, mean_lifetime, data = lifetime.longtheta.stat) + scale_y_continuous(limits = c(0,40)) + scale_x_continuous(limits=c(0,100)) + geom_point(colour = alpha("#0000FF", 1/3)) + geom_line(colour = alpha("black", 1/5)) + geom_smooth(errorbars_lifetime, stat="identity") + theme_bw()
ggsave(filename = "../images/kn-0.1-100-e100-combined-lifetime-by-theta.pdf")
#qplot(Theta, mean_lifetime, data = lifetime.40ktheta.stat, geom=c("point", "line"))


### need an unscaled graph to compare mean kn versus duration with all theta on the same plot

unscaled_title <- "Mean K_n versus Raw Assemblage Duration"
rich.longtheta.stat$ThetaFactor <- factor(rich.longtheta.stat$ThetaFactor, rev(levels(rich.longtheta.stat$ThetaFactor))) 
u <- ggplot(rich.longtheta.stat, aes(x = Windowsize, y = mean_kn, group=ThetaFactor)) + scale_colour_discrete(name="Theta") 
u + geom_line(aes(color = ThetaFactor)) + scale_x_continuous(name = "Assemblage Duration (Raw Simulation Steps)", trans="log10") + scale_y_continuous(name = "Mean K_n")
ggsave(filename = "../images/unscaled_kn_by_duration_stacked.pdf")

##################### Analyze Kn vs Time Averaging Window by Theta #################
# Add mean lifetime to each of the statistical summaries
#rich.wf.stat <- merge(rich.wf.stat, lifetime.wf.stat, by = "Theta")   ## need to generate this one, find the orig data?
rich.longtheta.stat <- merge(rich.longtheta.stat, lifetime.longtheta.stat, by = "Theta")
#rich.40ktheta.stat <- merge(rich.40ktheta.stat, lifetime.40ktheta.stat, by = "Theta")

# Add a column scaling windowsize by mean trait lifetime
rich.longtheta.stat$TLScaledTADur <- rich.longtheta.stat$Windowsize / rich.longtheta.stat$mean_lifetime
#rich.40ktheta.stat$TLScaledTADur <- rich.40ktheta.stat$Windowsize / rich.40ktheta.stat$mean_lifetime

# Add a column scaling the mean_kn value to the Kn value of the unaveraged process at each value of theta
# -- first, create a subset with just the mean_kn values for each theta value
unaveraged_kn <- subset(rich.longtheta.stat, Windowsize == 1, select = c("Theta", "mean_kn"))
# -- second, change the name of the Kn column so it doesn't conflict when we merge it back in
names(unaveraged_kn)[names(unaveraged_kn)=="mean_kn"] = "unaveraged_kn"
# -- next, merge it back in.  The "recycling" property of R vectors will ensure that the value is used for all rows matching a given theta
rich.longtheta.stat <- merge(rich.longtheta.stat, unaveraged_kn, by = "Theta")
# -- finally, calculate a new column which is the ratio of a specific Kn value for (theta, windowsize), and unaveraged Kn for that theta.  
#      any row with a windowsize of 1 will, of course, have 1.0000 for this ratio
rich.longtheta.stat$Kn_Ratio_to_Unaveraged <- rich.longtheta.stat$mean_kn / rich.longtheta.stat$unaveraged_kn



title2 <- "K_n Ratio by Log(Lifetime Scaled Duration) for Theta 0.1-100"
qplot(TLScaledTADur, Kn_Ratio_to_Unaveraged, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Time Averaged K_n / Non-time-averaged K_n", data = rich.longtheta.stat, main = title2) + geom_point(colour = alpha("#0000FF", 1/2)) + geom_line(colour = alpha("black", 1/10)) + facet_wrap(~ Theta, ncol = 3) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3)) + theme_bw()
ggsave(file = "../images/kn-0.1-100-e100-combined-scaled-richness-by-scaled-duration-vline.pdf")

minmaxmedianBreaks <- function(x){
  breaks <- c(min(x),median(x),max(x))
  #names(breaks) <- attr(breaks,"labels")
  breaks
}

errorbars <- aes(ymax = mean_kn + (sd_kn), ymin = mean_kn - (sd_kn) )
title3 <- "Mean K_n by Log(Lifetime Scaled Duration) for Theta 0.1-100"
qplot(TLScaledTADur, mean_kn, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Mean K_n", data = rich.longtheta.stat, main = title3) + scale_x_log10() + geom_point(colour = alpha("#0000FF", 1/2)) + geom_line(colour = alpha("black", 1/10)) + facet_wrap(~ Theta, ncol = 3, scales = "free_y") + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size = 1.5) + geom_smooth(errorbars, stat="identity") + theme_bw()
ggsave(file = "../images/kn-0.1-100-e100-combined-mean-richness-by-scaled-duration-errorbars-vline.pdf")



subset_and_plot_kn <- function(x, frame, filebase, titlebase, runlengthstring) {
  print(x)
  data.subset <- subset(frame, ThetaFactor == x)
  title <- paste(titlebase, x, runlengthstring)
  knplot <- qplot(log(TLScaledTADur), Kn_Ratio_to_Unaveraged, data = data.subset, main = title, geom = c("point", "line")) + theme_bw()
  fname <- paste("../images/", filebase, x, ".pdf", sep="")
  ggsave(filename = fname)          
}

#foreach(lvl=levels(rich.longtheta.stat$ThetaFactor)) %do% subset_and_plot_kn(lvl, rich.longtheta.stat, "kn-0.1-100-e100-combined-scaledkn-by-scaleduration-", "Scaled Kn by Scaled TA Duration by Theta ", " - 5K-40K Run Length")

# Create a faceted plot of histograms showing the distribution of trait richness for each combination 
# of theta and time averaging windowsize.  
qplot(Kn, ..density.., data = rich.longtheta, facets = ThetaFactor ~ WinFactor, geom = "histogram", binwidth=3 ) + theme_bw()
ggsave(file = "../images/kn-0.1-100-e100-combined-distribution-bytheta-bywindowsize.pdf")


# create tables useful for paper
sample_sizes <- ddply(rich.longtheta.stat, .(Windowsize), summarise,  min_sample = min(num), max_sample = max(num))
names(sample_sizes)[names(sample_sizes)=="Windowsize"] = "TA Duration"
names(sample_sizes)[names(sample_sizes)=="min_sample"] = "Min Sample Size"
names(sample_sizes)[names(sample_sizes)=="max_sample"] = "Max Sample Size"
write.csv(sample_sizes, file = "tables/wf-0.1-100-e100-combined-richness-sample-size-by-windowsize.csv")
ss.xtable <- xtable(sample_sizes,align="|c|c|c|c|", label="tab:sample-size-kn")
caption(ss.xtable) <- "Breakdown of sample sizes for analysis of trait richness, by size of time-averaging ``window.''  Some values of $\theta$ required larger numbers of simulation runs to achieve stable result, thus the difference between samples sizes at the same TA duration."
print(ss.xtable, file = "tables/wf-0.1-100-e100-combined-richness-sample-size-by-windowsize.tex", include.rownames=FALSE, latex.environment="ruledtabular")

mean_kn_subset_example <- subset(rich.longtheta.stat, Theta %in% c(0.5, 1.0, 20.0 ), select = c(Theta,Windowsize,mean_kn))
melted_subset <- melt(mean_kn_subset_example, id = c("Theta", "Windowsize"))
mean_kn_example_table <- cast(melted_subset, Windowsize ~ Theta | variable)
# I wasn't able to directly export this table to LaTeX with xtable, so copy and paste this result into \label{{tab:example-kn-ta-unscaled}
# in the paper! 

# save a table of mean trait lifetime
lifetime.xtable <- xtable(lifetime.longtheta.stat, align = "|c|c|c|", label="tab:mean-trait-lifetime")
names(lifetime.longtheta.stat)[names(lifetime.longtheta.stat)=="mean_lifetime"] = "Mean Trait Lifetime"
caption(lifetime.xtable) <- "Mean lifetime (in model generations) of traits, by $\theta$."
print(lifetime.xtable, file = "tables/wf-0.1-100-e100-combined-mean-lifetime.tex", include.rownames=FALSE, latex.environment="ruledtabular" )








######### end of code #################



