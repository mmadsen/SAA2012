######################
# SAA2012 
# 
# Analysis of Tf diversity test results from Wright-Fisher models
# under time-averaging.
#
# Mark E. Madsen (c) 2012 All rights reserved.  
# Software and Analysis for Mark E. Madsen SAA 2012 Papers is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# Based on a work at github.com.  See license.html for link to CC 3.0 license
#
######################
setwd("~/Dropbox/Research/Dissertation Project/analysis/saa2012/tf-diversity")
library(ggplot2)
library(plyr)

#### function definitions ####

#sigslat <- function(x) { return(ifelse(x < 0.05 | x > 0.95, 1, 0)) }


#### pull in data files ####

tfdiverse <- read.delim("../rawdata/wf-0.1-100-e100-combined-tf.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","character","numeric", "numeric","numeric","numeric","numeric"))
lifetime.longtheta <- read.delim("../rawdata/wf-0.1-100-e100-combined-lifetime.txt", header=TRUE, colClasses = c("factor","numeric","integer","numeric","numeric","numeric","numeric","numeric","numeric","integer"))

# statistical summary for add mean, stdev, stderr
tfdiverse.stat <- ddply(tfdiverse, .(Theta,Windowsize), summarise, MeanTf = mean(TfEst), StdevTf = sd(TfEst), StdErrTf = sd(TfEst)/sqrt(length(TfEst)), MeanIQV = mean(IQV), StdevIQV = sd(IQV), mean_sw = mean(ShannonWeaver), stdev_sw = sd(ShannonWeaver))

# determine mean lifetime for each theta value (windowsize doesn't figure into this)
lifetime.longtheta.stat <- ddply(lifetime.longtheta, .(Theta), summarise, mean_lifetime = mean(MeanLifetime))

# add mean lifetime and calculate mean lifetime scaled TA duration
tfdiverse.stat <- merge(tfdiverse.stat, lifetime.longtheta.stat, by = "Theta") 
tfdiverse.stat$TLScaledTADur <- tfdiverse.stat$Windowsize / tfdiverse.stat$mean_lifetime

errorbars <- aes(ymax = MeanTf + (StdevTf), ymin = MeanTf - (StdevTf) )
qplot(TLScaledTADur, MeanTf, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "T_f Diversity Index", data = tfdiverse.stat, main="T_f Diversity Index for Theta 0.1-100") + geom_point(colour = alpha("#0000FF", 2/3)) + facet_wrap(~ Theta, ncol = 3, scales = "free_y") + geom_line(colour = alpha("black", 1/4)) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size=1.5) + geom_smooth(errorbars, stat="identity") + theme_bw()
ggsave(file = "../images/wf-0.1-100-100-tf-index-abline-labels.pdf")

errorbars2 <- aes(ymax = MeanIQV + (StdevIQV), ymin = MeanIQV - (StdevIQV) )
qplot(TLScaledTADur, MeanIQV, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "IQV Diversity Index", data = tfdiverse.stat, main="IQV Diversity Index for Theta 0.1-100") + geom_point(colour = alpha("#0000FF", 2/3)) + facet_wrap(~ Theta, ncol = 3, scales = "free_y") + geom_line(colour = alpha("black", 1/4)) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size=1.5) + geom_smooth(errorbars2, stat="identity") + theme_bw()
ggsave(file = "../images/wf-0.1-100-100-iqv-index-abline-labels.pdf")

errorbars3 <- aes(ymax = mean_sw + (stdev_sw), ymin = mean_sw - (stdev_sw) )
qplot(TLScaledTADur, mean_sw, xlab = "Log(TA Duration Scaled by Trait Lifetime)", ylab = "Shannon Weaver Diversity Index", data = tfdiverse.stat, main="Shannon Weaver Diversity Index for Theta 0.1-100") + geom_point(colour = alpha("#0000FF", 2/3)) + facet_wrap(~ Theta, ncol = 3, scales = "free_y") + geom_line(colour = alpha("black", 1/4)) + scale_x_log10() + geom_vline(xintercept = 1, colour = alpha("#FF0000", 1/3), size=1.5) + geom_smooth(errorbars3, stat="identity") + theme_bw()
ggsave(file = "../images/wf-0.1-100-100-sw-index-abline-labels.pdf")

###### end of code #######
