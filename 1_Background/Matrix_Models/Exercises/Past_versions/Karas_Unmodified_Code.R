#  Demographic analysis using the popbio library and some other fun stuff
#  A population viability type analysis for a rare herbacious perennial plant,
#  Penstemon albomarginatus, for its only remaining California population

#  Objectives:
#  1) Estimating the population growth rate - deterministic and stochastic methods
#  2) Sensitivity and elasticity analysis - which transitions/rate are most variable, sensitive to change
#  3) Projecting the population stochastically for different scenarios, observing variation in growth rate
#  4) Quasi extinction probability
#  5) Estimating vital rates from real data, which will never be as good/much the data
#     that you would like to have for a PVA - some therapy for me, and maybe you too


#  Load popbio and some other helpful libraries. I use the dropdown menu in RStudio
#  to install libraries, it seems to work more consistently for my mac than scripting the installs, sorry.  
#  Am I missing something basic about installing packages? 

library(popbio)
library(plyr)
library(reshape)


# I have already done a bit of magic on my data get here, most or all of which is suggested
# in the Morris and Doak 2005 book Quantitative population biology. That book was an excellent
# guide to this kind of analysis.
# If you want to go deep (and deeper) into this kind of analysis you should get borrow or steal
# Morris and Doak 2005 and Caswell 2002 from the start!
# They are both really excellent and straightforward and you'll be able to answer a lot of the 
# questions that you have along the way.
# Most of the functions in popbio are based on functions in these two books. Many of the popbio
# examples are on their data and generate the same graphs etc. 

# what you need for a class/stage structured demographic model:
  # - a bunch of individuals (maybe in different populations)
  # - annual survival rates
  # - annual class or stage transition rates
  # - annual fecundity, 
  #       ie. probability of contribution to the juvenile class in the next year
  # - these data for a lot of years, preferably 10 or more (but we do what we can with less)

# The crux for plants is that its challenging to 1) count all seeds produced annually and
# 2) know how many seeds really yeilds a juvenile in any year.  I am doing a bit of work
# on these issues from limited field data.
# If you have a well behaved penguin  that you can collar and track and you 
# know produces 2 live juveniles each season, your PVA might be a bit simpler.


# Our field data had surivival for each year, mean plant diameter, and inflorescence count
# classes.  From this I found the median inflorescence number for each class MEDIAN_INLF
# I assigned each plant to CLASS based on its xDIAM_cm.  I 
# convinced myself that this had biological meaning by looking at the relationship between
# size class and survival.  Whatever classes or stages you use, you should be confident 
# that they are meaningful for your study species. OR use an integral projection model (IPM)
# instead. These allow you to use continuous variables like size or age rather than classes.

# nice abbreviated dataset
andre <-read.csv ("karadat.csv")
str(andre)

# Let's look at the stages/classes
levels(andre$CLASS)  # note that you need "dead" as a class for the first year that an individual is
                     # dead. After that it can be omitted. 


#####################################################################################  
####  1.  ESTIMATING FECUNDITY
#####################################################################################    
# The big challenge for plants = estimating seeds/indiv --> juveniles produced the next year
# based on a few fruit and seed counts and a lot of inflorescence class data

# Part 1:  How to estimate seed production per plant?
# Multiply seeds/infl from field observations by inflorescence number.
# Ideally you would have seed counts for each plant, but in the absence of those nearly impossible
# data, I'm using a few seeds/fruit counts * fruits/infl counts from 2011 and 2012.

# Part 2:  How many seeds makes an emergent juvenile?

# I'll add some notes on how I made these estimates at the end of this script, but let's start today 
# with a dataset ready to go for popbio
# I'm sure that every PVA has a hurdle or five to leap before you have the vital rate data ready to go.

# I've added estimated seeds produced per individual per year as a column in this simple dataset
# more on how I generated that later if we are intersted
str(andre)

#####################################################################################  
####  2.  PICK A STARTING POPULATION VECTOR
#####################################################################################    
# this is the # of individuals in each class/stage at the start of your model. 
# you should play around with this to see how it effects the outcome.  My model is insensitive to realistic changes in this vector.

#  Here are the "options", the number in each class observed annually
n_options<-ddply (andre, c("YEAR"), function (df)
  return(table(df$CLASS)))
n_options  

# picked starting population vector from 1995, the first year with 9 observed populations
n95<-c(81,31,17,13,11)
n=n95  # this must be called "n" for popbio to work its magic

# you can also see that I only have J class individuals for a subset of the years.  Check your data for 
# missing bits like this. Perhaps your collaborator is forgetting to tell you that they didn't survey
# for new juveniles from 1996-2003.  It happens!

#####################################################################################  
####  3.  GENERATE A TRANSITION MATRIX
#####################################################################################
# this matrix links each individual to its fate in the next year/cycle/season
# this is why you need each individual to be "dead" for a year, but no longer.
# If your raw data is like mine and only taken on live plants, "dead" might be something you have to add.

# make columns for year2, fate, and seeds2 for the whole census
trans<-subset(merge(andre, andre, by = "PLANT_UNQ", sort = FALSE), YEAR.x == YEAR.y - 1)
head2(trans)

# rename rows and columns to improve clarity (I use the names used by popbio, which are similar to Morris and Doak)
# I think I worked out somewhat painfully that you need these particular column names.  
rownames(trans) <- 1:nrow(trans)
colnames(trans)[1:7] <- c("plant",  "year", "stage", "seeds",  "year2", "fate", "seeds2")
head2(trans)

#  add individual fertility estimates from the calculations above
seedlingtrans<-0.00305    # This is the rate at which a seed becomes a J individual (I estimated this elsewhere, see Appendix below)

# adding in the number of J individuals produced by each individual
trans$J<- trans$seeds * seedlingtrans   # note that J is not an integer, which is totally fine, its a rate of J production
head2(trans)

#####################################################################################  
####  3.  GENERATE ANNUAL MATRICES  --- THE SIMPLE WAY for 3 easy years
#####################################################################################

#################     NAME STAGES   ###########################
stages <- c("J", "A1", "A2", "A3", "A4")  
# you must have a vector named stages in this way for your own clases or stages used 
# in the "stage" vector of your data frame

#################   SET ITERATIONS  ############################
it<-100        # set the number of time steps for a deterministic model

# Make a demographic projection matrix for each year like so:

#     1994  
trans94 <-subset(trans, year == 1994, c(plant, stage, fate, J))
(proj94<-projection.matrix(trans94, stage, fate, J, sort = stages))  #this gives you a projection matrix for 1994

# you can do a simple deterministic projection of the matrix for just this year
(p94<-pop.projection(proj94, n, it))
(l94<-p94$lambda)   # wow! if we looked only at 1994 based on these estimates the population would be booming!

# stable.stage shows the proportion of the population in each stage class at the mythical equilibrium, 48% of plants are juveniles in 100 years
# note that the population has gone from around 200 individuals in 1994 to 2.8e+19 in 100 years, nice!

# Now make some matrices for other years

#     1995  
trans95 <-subset(trans, year == 1995, c(plant, stage, fate, J))
(proj95<-projection.matrix(trans95, stage, fate, J, sort = stages))

#     2011  
trans11 <-subset(trans, year == 2011, c(plant, stage, fate, J))
(proj11<-projection.matrix(trans11, stage, fate, J, sort = stages))

p95<-pop.projection(proj95, n, it)
(l95<-p95$lambda)  # lambda is much lower in 1995

(p11<-pop.projection(proj11, n, pi))
(l11<-p11$lambda)  # and based on 2011 alone extinction is eminent. The gist here is we need lots of years of data to make any decent estimation
                   # of what the population is really likely to do (ie more than the three here)

#####################################################################################  
####  4.  DETERMINISTIC ANALAYSIS OF THESE 3 WELL BEHAVED YEARS - basic non-stochastic PVA
#####################################################################################
thesearethemeanprojmats<-list(proj94, proj95, proj11)  # make a list of the three matrices
(meanxprojmat <-mean(thesearethemeanprojmats))   # make a mean of the three projection matrices for deterministic analysis

n         # recall what n is, our starting population vector, ie the # of individuals in each class at the start of the projection
(pprojme <- pop.projection(meanxprojmat, n))   # do the deterministic projection, lambda is the dominant left eigenvalue
(DetLamb<-pprojme$lambda)


#####################################################################################  
####  5.  RUN STOCHASTIC ANALYSES ON THESE 3 WELL BEHAVED YEARS
#####################################################################################

# for a stochastic analysis, include all of the annual matrices, then make a random draw with replacement
# for a series of time step, or until stable stage distribution is reached.

thesearethemeanprojmats # our list of projection matrices

stochme <- stoch.growth.rate(thesearethemeanprojmats, prob=NULL, maxt=50000, verbose=TRUE)

#  note that these stochastic approximations of lambda are in log form  (not immediately comparable to pop.project$lambda)
exp(stochme$approx) # is the analytic approximation of lambda via Tuljapakar's method
stochme$approx      # this is more accurate (perhaps) when there is a lot of covariation in matrix elements
exp(stochme$sim)    # gives stochastic growth rate by simulation, random draws of whole
stochme$sim

# Tuljapurkar’s approximation takes into account how stochastic variation in the matrix elements affects the 
# long-term stochastic growth rate (Caswell 2001). It can be more accurate in cases where there is 
# covariation between matrix elements within the same year but may not be as accurate when there is a high 
# level of temporal variation (Morris and Doak 2002, Stubben et al. 2012).

#### Fun with stochastic analyses!
# Its easy to give years "weights" in the stochastic model. For example, you can increase the 
# drought rate by weighting drought years (2011)

yearweight<-c(1,1,2)
moredrought <- stoch.projection(thesearethemeanprojmats, n, tmax=50, prob=yearweight,  nreps = 500)
# the output is population sizes, which are fun to graph when comparing models
yearweight<-c(1,1,0)
nodrought <- stoch.projection(thesearethemeanprojmats, n, tmax=50, prob=yearweight,  nreps = 500)

par(mfrow=c(2,1))
hist(log(apply(moredrought, 1, sum)),col="blue",density=50, ylim=c(0,150), 
     xlim=c(-1.3,25), xlab="", main='More drought')
abline(v=log10(200), lty=3)  # puts a line at the starting population size for reference
hist(log(apply(nodrought, 1, sum)),col="green3",density=50, ylim=c(0,150), 
     xlim=c(-1.3,25), xlab="", main='No drought')
abline(v=log10(200), lty=3)
# you can get fancy and put these on the same graph too to compare outcomes.
# y axis is frequency of final population size at tmax.  


#####################################################################################  
####  6.  QUASI EXTINCTION BASED ON THESE 3 WELL BEHAVED YEARS
#####################################################################################
# another useful way to think about populations. Since our ability to really estimate lambda
# is based on the assumption of equilibrium at stable stage, it might be more realistic
# to think about comparing extinction probabilities for different scenarios that comparing
# lambdas. 
# these are based on stochastic runs

obsd<-stoch.quasi.ext(thesearethemeanprojmats, n,prob=c(1,1,1), 
                   Nx=10, tmax = 50, maxruns = 10, nreps=500, sumweight=c(1,1,1))

drt<-stoch.quasi.ext(thesearethemeanprojmats, n,prob=c(1,1,2), 
                   Nx=10, tmax = 50, maxruns = 10, nreps=500, sumweight=c(1,1,1))

par(mfrow=c(2,1))
matplot(obsd, ylab="Quasi-extinction probability",  ylim=c(0,1.1),
        type="l", lty=1, col=rainbow(10), las=1, main='Observed climate',
        xlab="Years")
matplot(drt, ylab="Quasi-extinction probability",  ylim=c(0,1.1),
        type="l", lty=1, col=rainbow(10), las=1, main='Double drought',
        xlab="Years")

#####################################################################################  
####  7.   SENSITIVITY AND ELASTICITY
#####################################################################################
# SENSITIVITY is a measure of the amount of change is lambda give a small change in a matrix element.

# ELASTICITY is a measure of ``proportional'' effect, 
#i.e., the effect that a change in a given matrix element has as a proportional to the change in that element


meanxprojmat  # for an overall look at sensitivity and elasticity use the mean projection matrix
# you could do separate analyses by year or type of year too to examine how sensitivity and elasticity 
# vary among years

(eigout <- eigen.analysis(meanxprojmat))  # do the associated sensitivity analysis

colSums(eigout$sensitivities)  # this gives the cumulative sensitivity of each stage/class, fun to graph

# calculate fertility and survival sums
(projsums <- colSums(meanxprojmat))
(fert_row <- meanxprojmat[1,])
(surv_row <- projsums-fert_row)

# can do this for sens and elas too to make a bar chart 
# make x where columes are name, sens and elas
# for my data I'm summing elasticity and sensitivity for each class.  You could also pick the vital 
# rates that are most meaningful in your own analysis

(eigout <- eigen.analysis(meanxprojmat))
eigout$sensitivities

(fert_row_s <- eigout$sensitivities[1,])
(surv_row_s <-eigout$sensitivities[2,] + eigout$sensitivities[3,]
 +eigout$sensitivities[4,] + eigout$sensitivities[5,] )
sensme<-t(cbind(fert_row_s, surv_row_s))


(fert_row_e <- eigout$elasticities[1,])
(surv_row_e <-eigout$elasticities[2,] + eigout$elasticities[3,]
 +eigout$elasticities[4,] + eigout$elasticities[5,] )
(survtable<-t(rbind(surv_row_s, surv_row_e)))
(ferttable<-t(rbind(fert_row_s, fert_row_e)))

par(mfrow=c(2,2))
mynames <- c("Sensitivity", "Elasticity")
barplot(t(survtable[,1:2]), beside=TRUE, las=1, ylim=c(0, 3.5),xlab="Stage class",
        main="Growth and survival")
abline(h=0)
barplot(t(ferttable[,1:2]), beside=TRUE, las=1, ylim=c(0, .25), xlab="Stage class",
        main="Fertility")
abline(h=0)
legend("topright", mynames, fill=grey.colors(2))




#####################################################################################  
####  APPENDIX  A.   HOW I GENERATED A FECUNDITY ESTIMATE
#####################################################################################
# I'm including this for plant folks who might like to see how I made fecundity estimates from real data
# I also have developed a script to simulate juvenile numbers and transition rates for the years in my dataset 
# that are missing these data, and bootstraps of the whole model.  


#  Load up the raw-ish data

andre <-read.csv ("D__composite9_13_2012.csv")
str(andre)

#  BEST ESTIMATE OF SEED PRODUCTION:  average of seeds/fruit, weighed average of fruits/infl where 0.14 is 1/8 drought years
andre$SEEDS<-andre$MEDIAN_INFL * 14.35  # = 14.35 seeds/infl

# The big crux Part 2:  How many seeds makes a J plant the next year? 
# ie what is the transition rate or fecundity rate?

# observed juveniles in each year (sadly not all years have the same # of cohorts, so I adjust for that below) 
seedling_yr<-ddply (andre, c("YEAR"), function (df)
  return(c(sdlgs=sum(df$CLASS == "J"))))    # shows the OBSERVED # of J plants each year, 
# and some of these are not observation years and are omitted below
seedling_yr

#seedlings for 9 cohorts based on census from each year 
sdl94=seedling_yr$sdlgs[(seedling_yr$YEAR=="1994")]*9/2  # this year had only 2 cohorts
sdl95=seedling_yr$sdlgs[(seedling_yr$YEAR=="1995")]
sdl11=seedling_yr$sdlgs[(seedling_yr$YEAR=="2011")]
sdl12=0     # in 2012 observed seedlings were 0

seedling_pick=c(sdl94,sdl95,sdl11, sdl12)
seedling_pick

# estimate seedlings as the mean of the 4 OBSERVED years:
# in the other years seedlings where not surveyed for
(seedlings=mean(seedling_pick))  # = 125.375 seedlings/year on average 

# what's the annual total seed production rate?
# add up all the seeds estimated to be produced in each year
seed_yr<-ddply (andre, c("YEAR"), function (df)
  return(c(sumseeds=sum(df$SEEDS))))
seed_yr

# adjust so that 1994 has an estimate of all cohorts based on the observed 2 cohorts
seed_yr$sumseeds[(seedling_yr$YEAR=="1994")]<-(seed_yr$sumseeds[(seedling_yr$YEAR=="1994")]*9/2)

# get mean seeds/year
avg_seeds_p_yr<-mean(seed_yr$sumseeds)

# for each of the 4 years in which seedlings were observed, calculate an estimate of the transition rate from seed --> J
str(sdl.trans94<-sdl94/avg_seeds_p_yr)  
str(sdl.trans95<-sdl95/avg_seeds_p_yr)   
str(sdl.trans11<-sdl11/avg_seeds_p_yr)   
sdl.trans12<-0  # no seedlings observed this year

# looks at the options for transition rate
str(seedlingtrans_pick<-c(sdl.trans94,sdl.trans95,sdl.trans11,sdl.trans12))
# pick the mean for this analysis 
(seedlingtrans<-mean(seedlingtrans_pick))       # mean seedling transition rate  0.00305

### Now remove watered and caged plants from main dataset, these have different survival and transition rates.
# I included them thus far because we needed to get seed production estimates for them. Since they might have contributed
# to the observed juveniles 
andre<-andre[(andre$CAGED =="N"),]  
andre<-andre[(andre$WATERED =="N"),]
head2(andre)  

# reduce datafile to include only PLANT_UNQ, YEAR, CLASS and SEEDS
str(andre)
str(andre<-andre[,c(1:3,16)])

# Go back to step 2. 



