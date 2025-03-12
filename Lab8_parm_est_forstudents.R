## Lab 8: parameter estimation
## serial interval and R0
## install the packages first
## install two pacakges: flexsurv and R0 first if you don't have them yet
library(flexsurv)
library(R0)

####################################################################################
## Part 1. Estimate the serial inteval for SARS using data from Lipsitch et al. 2003
####################################################################################
## read in data
getwd()  ## use "getwd()" to see your current directory and the path
setwd('YOUR DIRECTORY SAVING THE DATA')
# Go to Session -> Set Working Directory -> 


# alternative method using the file.choose() function
# da.sars=read.csv(file.choose())
da.sars=read.csv('SARS_serial_intervals.csv')

# Now take a look at the data
## each row is 1 case
## 1st column 'tm1' is the lower bound of the interval
## 2nd column 'tm2' is the upper bound of the interval
## they are the same for this dataset, as we only have a single value for each case
## 3rd column 'event' is the The status indicator, normally 
## 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
## here we have event=1 (showing symptoms)

## before doing the analysis, plot the dataset to see how it looks
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
hist(da.sars[,1],breaks=20,col='grey',ylim=c(0,25), main='',
     xlab='Serial Interval (days)',ylab='Number of cases')

## to use the survival function in the package (either flexsurv or survival)
## we have to first covert the data to a survival object
## to do so, run the command:
xsurv=Surv(da.sars$tm1,da.sars$tm2,da.sars$event,type='interval')

####################################################
## EXAMPLE CODE: Fit to the Weibull distribution
####################################################
## First try the survival function, with a Weibull distribution:
surv1=flexsurvreg(xsurv~1,dist='weibull')

surv1; # check the model output
surv1$res;  # model parm estimates are save in 'res'

## Calculate the mean based on the estimates
## for Weibull distribution:
## Check the Weibull distriubtion parameters here:  https://en.wikipedia.org/wiki/Weibull_distribution
surv1.mean=surv1$res['scale','est']*gamma(1+1/surv1$res['shape','est']); # the mean for weibull
surv1.sd=sqrt(surv1$res['scale','est']^2*(gamma(1+2/surv1$res['shape','est'])-(gamma(1+1/surv1$res['shape','est']))^2))
print(paste(round(surv1.mean,1),'+/-',round(surv1.sd,1)))

## Get the AIC:
surv1.AIC=surv1$AIC

## compute the model fit:
tm=seq(min(da.sars[,'tm1']),max(da.sars[,'tm2']),by=.2) # more time points than the observed, so we can fill the gap
# compute the number of cases per the survival model
# dweibull is the density function for the Weibull distribution
fit1=nrow(da.sars)*dweibull(tm,scale=surv1$res['scale','est'],
                            shape=surv1$res['shape','est']);

# Super-impose the model fit on the data for comparison
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
hist(da.sars[,1],breaks=20,col='grey',ylim=c(0,25), main='',
     xlab='Serial Interval (days)',ylab='Number of cases')
lines(tm,fit1,col='red',lwd=2)
legend('topright',cex=.9,seg.len = .8,
       legend=c('Observed','Fitted (Weibull)'),
       lty=c(0,1),pch=c(22,NA),lwd=c(NA,2),pt.bg = c('grey',NA),
       col=c('grey','red'),bty='n')


####################################################
## CODE YOURSELF: Fit to the Exponential distribution
####################################################
# [Q1] Fit the data to an exponential distribution, 
# and compute the mean and standard deviation (SD) of the serial interval, based on your model.  
# To visualize model fit, you can also plot the data and model fit from this model. 
# Hint: Use the table from wikipedia to for the exponential distribution to get relevant parameters 
# for computing the mean and SD of the serial interval: https://en.wikipedia.org/wiki/Exponential_distribution



####################################################
## CODE YOURSELF: Fit to the log-normal distribution
####################################################
# [Q2] Fit the data to a log-normal distribution, 
# and compute the mean and standard deviation (SD) of the serial interval, based on your model.  
# To visualize model fit, you can also plot the data and model fit from this model.
# Hint: Use the table from wikipedia to for the log-normal distribution to get relevant parameters 
# for computing the mean and SD of the serial interval: https://en.wikipedia.org/wiki/Log-normal_distribution





# [Q3] Based on the AIC of the three models (Weibull, exponential, and log-normal), which fits the data best here? (0.5 pt)




####################################################################################
## Part 2: Estimating R0 from the exponential growth phase of the epidemic
####################################################################################
## data: daily incidence during 1918 influenza pandemic in Germany (from 'R0' library)
da.flu=read.csv('data_1918pandemic_Germany.csv')
da.flu$date = as.Date(da.flu$date)  # convert to date

## Always plot and check the data first
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
plot(da.flu[,1],da.flu[,2], xlab='',ylab='Cases', pch = 20)



# [LQ4] Calculate the cumulative incidence using the daily incidence dataset (‘data_1918pandemic_Germany.csv’) and plot the log(cumulative incidence) v. time (0.5pt)
## Hint: cumpute the cumulative incidence using the cumsum function
cumI=cumsum(da.flu[,2]); 

# plot and see:



# [LQ5] Assume a generation time of 3 days, estimate R0 from the exponential growth phase of the epidemic
# (try the first 7, 14, and 21 days)
# Report your R0 estimates for these three time periods (i.e., the first 7, 14, and 21 days)  (1pt)

# EXAMPLE: FOR THE FIRST 7 DAYS
# Fit to the data during the first 7 days:
D=3; # set the generation time to 3 days
Ndays=7;  # ADJUST THE NUMBER OF DAYS INCLUDED IN THE FIT HERE
tm1=1:Ndays;
fit1=lm(log(cumI[1:Ndays])~tm1)
summary(fit1)

# compute R0 based on the model-fit
R=1+fit1$coefficients[2]*D   # slope: fit1$coefficients[2]





####################################################################################
## Part 3: R0 - Maximum Likelihood Estimation (MLE)
####################################################################################
data("Germany.1918") # Request the data (it's from the package)
Germany.1918 # print the data to see the structure

# First we need the distribution of generation time (i.e. serial interval)
mGT<-generation.time("gamma", c(2.45, 1.38))


#  [LQ7] Based on the flu data and the above generation time, 
# estimate R0 using the est.R0.ML function and incidence 
# during the first 7, 14, and 21 days of the Germany.1918 time series.  
# Report your R0 estimates for these three time periods.  (1pt)

# EXAMPLE: FOR THE FIRST 7 DAYS
# Maximum Likelihood Estimation using the est.R0.ML function
est.R0.ML(Germany.1918, # the data
          mGT, # the distribution of serial interval
          begin=1, # the start of the data
          end=7, # ADJUST THE NUMBER OF DAYS TO INCLUDE IN THE MODEL HERE
          range=c(0.01,50) # the range of possible values to test
          )



# [LQ8] Base on your estimated R0 over the first 7 days, 
# how long would it take for the number of cases to double (i.e. doubling time) 
# given the serial interval of 2.45 days used here? How is your estimate compared to the observations? (0.5pt)



