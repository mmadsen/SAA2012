######################
# SAA2012 
# 
# Analysis of Slatkin Exact test results from Wright-Fisher models
# under time-averaging.
#
# Mark E. Madsen (c) 2012 All rights reserved.  
# Software and Analysis for Mark E. Madsen SAA 2012 Papers is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# Based on a work at github.com.  See license.html for link to CC 3.0 license
#
######################
setwd("~/Dropbox/Research/Dissertation Project/analysis/saa2012/slatkin-tests")
library(ggplot2)
library(plyr)

#### function definitions ####

sigslat <- function(x) { return(ifelse(x < 0.05 | x > 0.95, 1, 0)) }


#### pull in data files ####

#slat <- read.delim("../rawdata/unified-slatkin-results-by-windowsize-and-gen-5000.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer","character","numeric"))
#slatbig <- read.delim("../rawdata/unified-wf-sweep-50k.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer","character","numeric"))
#slat100 <- read.delim("../rawdata/sweep-theta-2-40-ewens-100.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer","character","numeric"))
slatsweep <- read.delim("../rawdata/wf-0.1-100-e100-combined-slatkin.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","integer","numeric","integer","character","numeric", "numeric"))
lifetime.longtheta <- read.delim("../rawdata/wf-0.1-100-e100-combined-lifetime.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","numeric","numeric","numeric","numeric","numeric","integer"))

# determine mean lifetime for each theta value (windowsize doesn't figure into this)
lifetime.longtheta.stat <- ddply(lifetime.longtheta, .(Theta), summarise, mean_lifetime = mean(MeanLifetime))

# now let's post-process the data from sweep with theta 0.1-100, but ewens sample sizes of 100
slatsweep$sslat <- sigslat(slatsweep$SlatkinExactP)
slatsweep$WinFactor <- factor(slatsweep$Windowsize)
slatsweep$ThetaFactor <- factor(slatsweep$Theta)
slatsweep$sigfactor <- factor(slatsweep$sslat)



slatsweep.wf <- subset(slatsweep, Model == "WrightFisherDrift")
attach(slatsweep.wf)
table(slatsweep.wf$sigfactor,slatsweep.wf$WinFactor,slatsweep.wf$ThetaFactor)
slatsweep.passfail.table <- table(slatsweep.wf$sigfactor,slatsweep.wf$WinFactor,slatsweep.wf$ThetaFactor)
write.table(slatsweep.passfail.table, file="tables/slatkin-wf-passfail-byfactors-ewens100-5k-40k.txt", append = FALSE, quote=FALSE, sep="\t")
slatsweep.wf.extrafailures <- ddply(slatsweep.wf, .(Windowsize,Theta), summarise, FailuresOver = mean(sslat) - 0.10, StdevFail = sd(sslat), StdErrFail = sd(sslat)/sqrt(length(sslat)))
slatsweep.wf.thetaestimates <- ddply(slatsweep.wf, .(Theta, Windowsize), summarise, MeanEstTheta = mean(ThetaEst), StdevFail = sd(ThetaEst), StdErrEstTheta = sd(ThetaEst)/sqrt(length(ThetaEst)))

# Add mean lifetime to the summary of extra failures of the slatkin test

slatsweep.wf.extrafailures <- merge(slatsweep.wf.extrafailures, lifetime.longtheta.stat, by = "Theta")
slatsweep.wf.extrafailures$TLScaledTADur <- slatsweep.wf.extrafailures$Windowsize / slatsweep.wf.extrafailures$mean_lifetime
slatsweep.wf.extrafailures$LifetimeScaledDuration <- factor(slatsweep.wf.extrafailures$TLScaledTADur)
slatsweep.wf.extrafailures$ThetaFactor <- factor(slatsweep.wf.extrafailures$Theta)

slatsweep.wf.thetaestimates <- merge(slatsweep.wf.thetaestimates, lifetime.longtheta.stat, by = "Theta")
slatsweep.wf.thetaestimates$TLScaledTADur <- slatsweep.wf.thetaestimates$Windowsize / slatsweep.wf.thetaestimates$mean_lifetime
slatsweep.wf.thetaestimates$ThetaEstRatio <- slatsweep.wf.thetaestimates$MeanEstTheta / slatsweep.wf.thetaestimates$Theta
#qplot(log(TLScaledTADur), ThetaEstRatio, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Ratio of Slatkin Theta Estimate to Actual Theta", main = "Ratio of Slatkin Theta Estimates to Model Theta under Time-Averaging - 0.1-100", data = slatsweep.wf.thetaestimates) + geom_point(colour = alpha("#0000FF", 1/2)) + geom_line(colour = alpha("black", 1/5)) + facet_wrap(~ Theta, ncol = 3, scales = "free_y") + geom_vline(xintercept = log(1), colour = alpha("#FF0000", 1/3))
#ggsave(file = "../images/slatkin-theta-estimate-ratio-0.1-100-5k-40k.pdf")

pdf(file = "../images/slatkin-tests-theta-0.1-100-ewens100-5k-40k-failure-heatmap.pdf", width = 6, height=6)
slatsweep.heat.fail <- ggplot(slatsweep.wf.extrafailures, aes(Windowsize, Theta) , main="Average Slatkin Exact Excess Failure Rate - Theta 0.1-100") + theme_bw()
slatsweep.heat.fail + geom_tile(aes(fill = FailuresOver), colour="white")
dev.off()
#pdf(file = "../images/slatkin-tests-theta-0.1-100-ewens100-5k-40k-failure-scaled-heatmap.pdf", width = 6, height=6)
#slatsweep.heat.fail <- ggplot(slatsweep.wf.extrafailures, aes(LifetimeScaledDuration, ThetaFactor), main="Average Slatkin Exact Excess Failure Rate (Scaled) - Theta 0.1-100")
#slatsweep.heat.fail + geom_tile(aes(fill = FailuresOver), colour="white")
#dev.off()

errorbars <- aes(ymax = FailuresOver + (StdErrFail), ymin = FailuresOver - (StdErrFail) )
# for print
qplot(TLScaledTADur, FailuresOver, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Slatkin Exact Failure Rate beyond alpha = 10%", data = slatsweep.wf.extrafailures, main="Mean Failure Rate for Slatkin Exact Tests Beyond 10% Expected - Theta 0.1-100") + geom_point(colour = alpha("#0000FF", 2/3)) + facet_wrap(~ Theta, ncol = 3) + geom_line(colour = alpha("black", 1/4)) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size=1.5) + theme_bw() + geom_smooth(errorbars, stat="identity")
ggsave(file = "../images/slatkin-test-theta-0.1-100-100-5k-40k-failure-means-with-abline-labels.pdf")
# for presentation - NOT GOOD
#qplot(TLScaledTADur, FailuresOver, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Slatkin Exact Failure Rate beyond alpha = 10%", data = slatsweep.wf.extrafailures, main="Mean Failure Rate for Slatkin Exact Tests Beyond 10% Expected - Theta 0.1-100") + geom_point(colour = alpha("#0000FF", 2/3)) + facet_wrap(~ Theta, ncol = 2) + geom_line(colour = alpha("black", 1/4)) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size=1.5) + theme_bw() + geom_smooth(errorbars, stat="identity")
#ggsave(file = "../images/slatkin-test-theta-0.1-100-100-5k-40k-failure-means-with-abline-labels-slide.pdf")


############ end of code ##########

# slat$sslat <- sigslat(slat$SlatkinExactP)
# slat$WinFactor <- factor(slat$Windowsize)
# slat$ThetaFactor <- factor(slat$Theta)
# slat$sigfactor <- factor(slat$sslat)
# slat.wf <- subset(slat, Model == "WrightFisherDrift")
# attach(slat.wf)
# table(slat.wf$sigfactor,slat.wf$WinFactor,slat.wf$ThetaFactor)
# slat.passfail.table <- table(slat.wf$sigfactor,slat.wf$WinFactor,slat.wf$ThetaFactor)
# write.table(slat.passfail.table, file="tables/slatkin-wf-passfail-byfactors.txt", quote=FALSE, append = FALSE, sep="\t")
# slat.wf.extrafailures <- ddply(slat.wf, .(Windowsize,Theta), summarise, FailuresOver = mean(sslat) - 0.10)
# slat.wf.thetaestimates <- ddply(slat.wf, .(Theta, Windowsize), summarise, MeanEstTheta = mean(ThetaEst))
# qplot(log(Windowsize), FailuresOver, data = slat.wf.extrafailures, main="Average Slatkin Exact Excess Failure Rate - Theta 2-40 Ewens 30", geom=c("point", "line"), facets = Theta ~ .)
# ggsave(file="../images/slatkin-tests-theta-2-40-ewens30-failure-means.pdf")
# slat.wf.extrafailures$WindowFactor <- factor(slat.wf.extrafailures$Windowsize)
# slat.wf.extrafailures$ThetaFactor <- factor(slat.wf.extrafailures$Theta)
# pdf(file = "../images/slatkin-tests-theta-2-40-ewens30-failure-heatmap.pdf", width = 6, height=6)
# slat.heat.fail <- ggplot(slat.wf.extrafailures, aes(WindowFactor, ThetaFactor), main="Average Slatkin Exact Excess Failure Rate - Theta 2-40 Sample 30 - Runlength 5K")
# slat.heat.fail + geom_tile(aes(fill = FailuresOver), colour="white")
# dev.off()

# now let's analyze the data set with long 50K step runs
# slatbig$sslat <- sigslat(slatbig$SlatkinExactP)
# slatbig$WinFactor <- factor(slatbig$Windowsize)
# slatbig$ThetaFactor <- factor(slatbig$Theta)
# slatbig$sigfactor <- factor(slatbig$sslat)
# slatbig.wf <- subset(slatbig, Model == "WrightFisherDrift")
# attach(slatbig.wf)
# table(slatbig.wf$sigfactor,slatbig.wf$WinFactor,slatbig.wf$ThetaFactor)
# slatbig.passfail.table <- table(slatbig.wf$sigfactor,slatbig.wf$WinFactor,slatbig.wf$ThetaFactor)
# write.table(slatbig.passfail.table, file="tables/slatkin-wf-passfail-byfactors-50k-runs.txt", append = FALSE, quote=FALSE, sep="\t")
# slatbig.wf.extrafailures <- ddply(slatbig.wf, .(Windowsize,Theta), summarise, FailuresOver = mean(sslat) - 0.10)
# slatbig.wf.thetaestimates <- ddply(slatbig.wf, .(Theta, Windowsize), summarise, MeanEstTheta = mean(ThetaEst))
# qplot(log(Windowsize), FailuresOver, data = slatbig.wf.extrafailures, main="Average Slatkin Exact Excess Failure Rate - Theta 2-40 Ewens 30 Runlength 50K",geom=c("point", "line"),facets = Theta ~ .)
# ggsave(file="../images/slatkin-tests-theta-2-40-ewens30-50k-failure-means.pdf")
# slatbig.wf.extrafailures$WindowFactor <- factor(slatbig.wf.extrafailures$Windowsize)
# slatbig.wf.extrafailures$ThetaFactor <- factor(slatbig.wf.extrafailures$Theta)
# pdf(file = "../images/slatkin-tests-theta-2-40-ewens30-50k-failure-heatmap.pdf", width = 6, height=6)
# slatbig.heat.fail <- ggplot(slatbig.wf.extrafailures, aes(WindowFactor, ThetaFactor), main="Average Slatkin Exact Excess Failure Rate - Theta 20-40 Sample 30 - Runlength 50K")
# slatbig.heat.fail + geom_tile(aes(fill = FailuresOver), colour="white")
# dev.off()

# # now let's post-process the data from sweep with theta 2-40, but ewens sample sizes of 100
# slat100$sslat <- sigslat(slat100$SlatkinExactP)
# slat100$WinFactor <- factor(slat100$Windowsize)
# slat100$ThetaFactor <- factor(slat100$Theta)
# slat100$sigfactor <- factor(slat100$sslat)
# slat100.wf <- subset(slat100, Model == "WrightFisherDrift")
# attach(slat100.wf)
# table(slat100.wf$sigfactor,slat100.wf$WinFactor,slat100.wf$ThetaFactor)
# slat100.passfail.table <- table(slat100.wf$sigfactor,slat100.wf$WinFactor,slat100.wf$ThetaFactor)
# write.table(slat100.passfail.table, file="tables/slatkin-wf-passfail-byfactors-ewens100.txt", append = FALSE, quote=FALSE, sep="\t")
# slat100.wf.extrafailures <- ddply(slat100.wf, .(Windowsize,Theta), summarise, FailuresOver = mean(sslat) - 0.10)
# slat100.wf.thetaestimates <- ddply(slat100.wf, .(Theta, Windowsize), summarise, MeanEstTheta = mean(ThetaEst))
# qplot(log(Windowsize), FailuresOver, data = slat100.wf.extrafailures, main="Average Slatkin Exact Excess Failure Rate - Theta 2-40 Ewens 100",geom=c("point", "line"),facets = Theta ~ .)
# ggsave(file="../images/slatkin-tests-theta-2-40-ewens100-failure-means.pdf")
# slat100.wf.extrafailures$WindowFactor <- factor(slat100.wf.extrafailures$Windowsize)
# slat100.wf.extrafailures$ThetaFactor <- factor(slat100.wf.extrafailures$Theta)
# pdf(file = "../images/slatkin-tests-theta-2-40-ewens100-failure-heatmap.pdf", width = 6, height=6)
# slat100.heat.fail <- ggplot(slat100.wf.extrafailures, aes(WindowFactor, ThetaFactor), main="Average Slatkin Exact Excess Failure Rate - Theta 2-40 - Runlength 5K")
# slat100.heat.fail + geom_tile(aes(fill = FailuresOver), colour="white")
# dev.off()


############ end of code ##########
