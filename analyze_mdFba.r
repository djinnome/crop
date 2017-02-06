##jmd
##10.22.12
##analyze_mdFba.r
#analyze test non-essentials from md-fba, since ended up w/ redundant growth rates for same KO that were inconsistent

setwd('C:/Documents and Settings/jdreyf/Desktop/farm_data/matlab_data')

teNoness <- read.csv('mdfba_teNeGrowth.csv', as.is=TRUE, na='NaN')
teNoness[is.na(teNoness)] <- 0

growth <- tapply(teNoness$Growth.rate, teNoness$Name, max)

# get 60% whether I use teNoness$Name or teNoness$No, & teNoness$Growth.rate>0.01 or teNoness$Rel..GR>1 (Rel..GR is in percent)
mean(growth>0.01) #60%
