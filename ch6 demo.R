# Table 6.1 dataset
x<-c(1089,1092,rep(1094,2),rep(1095,4),rep(1098,8),rep(1100,9),
     rep(1104,12),rep(1105,4),rep(1107,5),rep(1108,5),rep(1110,4),
     rep(1112,3),rep(1115,2))

# histogram visualization: Figure 6.1
breaks <- seq(1085,1120,5)
hist(x, breaks = breaks, freq=TRUE, right = FALSE, xlab="Temperature range",ylab="Number of measurements",main="Histogram")
# you can a distribution curve 
#curve(dnorm(x, mean=mean(x), sd=sd(x)), add=TRUE, col="darkblue", lwd=2)  

# calculate statistics
mean(x)
getmode(x)
median(x)
var(x)
sd(x)
sd(x)/mean(x) #coefficient of variation
sd(x)/length(x) #standard error

#user function: getmode()
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#uncertainty of the mean
qnorm(p=0.975)*sd(x)/length(x)

#---------------------------------------------------------------------------------------
#Introduction to R: 
# 1). help(functionName), ?functionName
# 2). install.packages() and library(), take for example "asbio" and "chemometrics" 
# 3). ctrl+L or cat('\f')
# 4). Online-resources: 
#     https://www.rdocumentation.org
#     http://www.r-tutor.com/
#version
#install.packages("installr")
#library(installr)
#updateR()
#install.packages("asbio")
#"asbio" has been built in since R3.4
library("asbio")
Mode(x) # "asbio" function

#---------------------------------------------------------------------------------------
# Data fitting: lm
data("Fbird") # data frame: vol and freq
attach(Fbird)
Fbird.lm <- lm(freq~vol, data=Fbird)
summary(Fbird.lm)

coeffs = Fbird.lm$coefficients
plot(vol,freq, xlab="Vol", ylab="Freq")
abline(coeffs[1],coeffs[2])

#---------------------------------------------------------------------------------------
#Multivaraite 
