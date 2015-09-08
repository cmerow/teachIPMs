###########################################################################
###########################################################################
###########################################################################
# Introductory IPM exercises
###########################################################################
###########################################################################
###########################################################################

#Here, we introduce an extremely simple IPM for a the long- live alpine perennial plant Dracocephalum austriacum. The analyses continue from those in 'Beginner/Intro_to_IPMs/Exercises/Intro_to_IPMs_exercises.r'. After building a basic model, we illustrate ways to improve it, including, comparing alternative growth models, non-constant variance in growth, including a model for flowering probability.  We also illustrate a series of diagnostics, including correcting for eviction, determining asymptotic maximum size, and exploring the importance of transient dynamics. 

# OVERVIEW
# the document is organized as follows
# A. Rerun the core code from Intro_to_IPMs_exercises.r
# B. Non-constant variance in growth function
# C. Changing the mean growth
# D. Incorporating flowering probability

	# set up directory structure. we'll place a temp folder on your desktop to store some plots
if(!file.exists('~/Desktop/Temp_IPM_output')) dir.create('~/Desktop/Temp_IPM_output')
	# set this as the working directory
setwd('~/Desktop/Temp_IPM_output') 

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# A. rerun the core code from Intro_to_IPMs_exercises.r
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# just in case you had issues with Intro_to_IPMs_exercises.r, we begin by running just the core components of that code to give an IPM as a starting point.
  	# read in data. you'll have to set your own file path here. the data is a .csv file included in the workshop's dropbox folder
d=read.csv( '~/Documents/Work/teachIPMs/Beginner/Intro_to_IPMs/Exercises/Intro_to_IPMs_Exercises_Data.csv')

# 0. set up parameter list for regressions
	# this sets up a list of the model parameters. these parameters will be estimated and recorded below.
params=data.frame(
	surv.int=NA,
	surv.slope=NA,
	growth.int=NA,
	growth.slope=NA,
	growth.sd=NA,
	seed.int=NA,
	seed.slope=NA,
	recruit.size.mean=NA,
	recruit.size.sd=NA,
	establishment.prob=NA
)

# 1. survival regression
surv.reg=glm(surv~size,data=d,family=binomial())
params$surv.int=coefficients(surv.reg)[1]
params$surv.slope=coefficients(surv.reg)[2]

# 2. growth regression
growth.reg=lm(sizeNext~size,data=d)
params$growth.int=coefficients(growth.reg)[1]
params$growth.slope=coefficients(growth.reg)[2]
params$growth.sd=sd(resid(growth.reg))

# 3. seeds regression
# note that we are just pooling all individuals into this regression regardless of whether they flowered or not. a later exercise will be to explicitly model flowering probability. i.e. for the moment we're taking flowering_probability[z]=1.
seed.reg=glm(fec.seed~size,data=d,family=poisson())
params$seed.int=coefficients(seed.reg)[1]
params$seed.slope=coefficients(seed.reg)[2]

# 4. size distribution of recruits
	# in the dataframe, recruits are those individuals who have a value for sizeNext but not for size
params$recruit.size.mean=mean(d$sizeNext[is.na(d$size)])
params$recruit.size.sd=sd(d$sizeNext[is.na(d$size)])

# 5. establishment probability
	# these data represent a single year's worth of data, hence establishment probability can be estimated by dividing the number of observed recruits by the number of seeds. hence the growth/survival measurements were taken in year t which the recruit sizes were measured in year t+1.
params$establishment.prob=sum(is.na(d$size))/sum(d$fec.seed,na.rm=TRUE)

## vital rate functions

# 1. probability of surviving
s.x=function(x,params) {
	u=exp(params$surv.int+params$surv.slope*x)
	return(u/(1+u))
}

# 2. growth function
g.yx=function(xp,x,params) { 			
	dnorm(xp,mean=params$growth.int+params$growth.slope*x,sd=params$growth.sd)
}

# 3. reproduction function      
f.yx=function(xp,x,params) { 		
	params$establishment.prob*
	dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
	exp(params$seed.int+params$seed.slope*x)
}

# 1. boundary points b, mesh points y and step size h
	# integration limits - these limits span the range of sizes observed in the data set, and then some.
min.size=.9*min(c(d$size,d$sizeNext),na.rm=T)
max.size=1.1*max(c(d$size,d$sizeNext),na.rm=T)
	# number of cells in the discretized kernel
n=100 
	# boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y=0.5*(b[1:n]+b[2:(n+1)])
# width of the cells
h=y[2]-y[1]

# 2. make component kernels

G=h*outer(y,y,g.yx,params=params) 	# growth kernel
S=s.x(y,params=params) 							# survival 
P=G 																# placeholder;redefine P on the next line
for(i in 1:n) P[,i]=G[,i]*S[i]  		# growth/survival kernel  
F=h*outer(y,y,f.yx,params=params) 	# reproduction kernel 
K=P+F 															# full kernel

# 1. get lamda,v,w  
(lam=Re(eigen(K)$values[1])) # should be 1.013391


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# B.  Non-constant variance in growth function
# -------------------------------------------------------------------
# -------------------------------------------------------------------
  # Typically, the workflow for building IPMs is iterative; based on preliminary analyses, we determine which important factors are missing or superfluous in the models. In this section, we make three improvements to the model built thus far: (1) we make the variance in growth among individuals a function of size (i.e. growth is more variable in larger individuals); (2) add a quadratic term for size to the growth function to capture the plateauing pattern in the growth data (Fig. 2b); (3) modify the fecundity function to include flowering probability, thereby recognizing that not all individuals reproduce.
  
  # for all of the following exercises, code from above can be adapted with just a few minor changes.
  
  # make the variance of the growth regression a function of size and obtain the value of lambda. in the model above, we simply modeled the variance as constant, as seen in params$growth.sd. by looking at the growth data in figure 1 (section A) you might guess that the variance increases as a function of size. you can see this in the following plot, where we plot the absolute value of the residuals of the growth regression against size:
  plot(growth.reg$model$size,abs(resid(growth.reg)),xlab='size',ylab='residual')
  
  # incorporating size in the growth variance requires 4 steps:
  	# 1. build a regression that includes size as a predictor of the variance. this is most easily done using generalized least squares (gls). gls allows us to simultaneously fit a model for the expected size at t+1 and the variance. for this example, we'll model the variance in growth as an exponential function sigma.2(size)=(sigma^2)*exp[2*constant*size]. see ?varExp for more details. the gls model will estimate sigma^2 and the constant, in addition to the intercept and slope used in the mean of the growth regression. gls() is in the nlme package, which you may need to install and load using install.packages('nlme') and library(nlme).
			library(nlme)
			growth.reg=gls(sizeNext~size,weights=varExp(),na.action=na.omit, data=d)
			summary(growth.reg)
				# plot the model, just to be sure it's reasonable
			plot(d$size,d$sizeNext,main='Growth/Shrinkage/Stasis')	
			xx=seq(0,8,by=.01)
  		lines(xx,predict(growth.reg,data.frame(size=xx)),col='red',lwd=3)
 			
  	# 2. modify the params data frame (where we store all the coefficients used to build the ipm) to have coefficients called growth.sd.int and growth.sd.slope and set the values of those coefficients equal to the appropriate coefficients from part a.
			
			params$growth.int=coefficients(growth.reg)[1]
			params$growth.slope=coefficients(growth.reg)[2]
			params$growth.sigma2=summary(growth.reg)$sigma^2 
			params$growth.sigma2.exp=as.numeric(growth.reg$modelStruct$varStruct)

  	# 3. modify the growth function, g.xy, to allow the standard deviation argument, sd, to be a function of size. this will follow a similar pattern to the argument for mean. it's probably easiest not to rename the function g.yx, otherwise you'll have modify the code in section D.2. 
  	
			g.yx=function(xp,x,params) { 			
				dnorm(xp,mean=params$growth.int+params$growth.slope*x, sd=sqrt(params$growth.sigma2*exp(2*params$growth.sigma2.exp*x)))
			}

		# 4. rebuild the kernels and obtain lambda. this code is simple copied and pasted from section A above. the value should be 0.980247. plot the P (survival/growth) matrix and check that it reasonably predicts the growth observations. plot the results using code in section E and save the figure as 'IPM_output_v2.pdf' for later reference.

			# 2. make component kernels
			G=h*outer(y,y,g.yx,params=params) 	# growth kernel
			S=s.x(y,params=params) 							# survival 
			P=G 																# placeholder;redefine P on the next line
			for(i in 1:n) P[,i]=G[,i]*S[i]  		# growth/survival kernel
			F=h*outer(y,y,f.yx,params=params) 	# reproduction kernel	
			K=P+F 															# full kernel

			# 1. get lamda,v,w  
			(lam=Re(eigen(K)$values[1])) 
			w.eigen=Re(eigen(K)$vectors[,1])
			stable.dist=w.eigen/sum(w.eigen) 
			v.eigen=Re(eigen(t(K))$vectors[,1])
			repro.val=v.eigen/v.eigen[1] 
	
			# compute elasticity and sensitivity matrices
			v.dot.w=sum(stable.dist*repro.val)*h
			sens=outer(repro.val,stable.dist)/v.dot.w
			elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
		
			# plot results 
			# pdf('IPM_output_v2.pdf',h=8,w=12)
			par(mfrow=c(2,3)) 
			image(y,y,t(K), xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Kernel")
			contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
			plot(y,stable.dist,xlab="Size",type="l",main="Stable size distribution")
			plot(y,repro.val,xlab="Size",type="l",main="Reproductive values") 
			image(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity")
			image(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity")
  		image(y,y,t(P),main="Growth/Survival Kernel+Data")
  		points(d$size,d$sizeNext)	
  		# dev.off()
      # system('open IPM_output_v2.pdf')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# C.  Changing the mean growth
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# The next improvement to the IPM is to add a quadratic term to the growth regression to capture the slight downward curving pattern in growth (Fig. 2).  For simplicity, we keep the size dependence of the growth variance from the previous section in the model, but note that we rerun the variance regressions since the residuals will change when the quadratic term is added to the model.
 
 # add a quadratic term to growth regression build in question 1 and obtain the value of lambda (it should be 0.9777275). you should follow steps analogous to those in question 1. to tell the growth regression to use a quadratic term, use I(size^2) as a predictor variable. as in question 1d, plot the P matrix to check that the model is reasonable. i'll start you off with the basic structure of the code.
		growth.reg=gls(sizeNext~size+I(size^2),weights=varExp(), na.action=na.omit, data=d)
		summary(growth.reg)
		params$growth.int=coefficients(growth.reg)[1]
		params$growth.slope=coefficients(growth.reg)[2]
		params$growth.sqrd=coefficients(growth.reg)[3]
		params$growth.sigma2=summary(growth.reg)$sigma^2 
		params$growth.sigma2.exp=as.numeric(growth.reg$modelStruct$varStruct)  
	
		# now, following the example for g.yx above, modify it to include the quadratic term here
		g.yx=function(xp,x,params) { 			
			??
		}
		
		# insert code to rebuild kernels from above
				# make component kernels
		G=h*outer(y,y,g.yx,params=params) 	# growth kernel
		S=s.x(y,params=params) 							# survival 
		P=G 																# placeholder;redefine P on the next line
		for(i in 1:n) P[,i]=G[,i]*S[i]  		# growth/survival kernel
		F=h*outer(y,y,f.yx,params=params) 	# reproduction kernel	
		K=P+F 															# full kernel

		#  get lamda,v,w  
		(lam=Re(eigen(K)$values[1])) 
		w.eigen=Re(eigen(K)$vectors[,1])
		stable.dist=w.eigen/sum(w.eigen) 
		v.eigen=Re(eigen(t(K))$vectors[,1])
		repro.val=v.eigen/v.eigen[1] 

		# compute elasticity and sensitivity matrices
		v.dot.w=sum(stable.dist*repro.val)*h
		sens=outer(repro.val,stable.dist)/v.dot.w
		elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)

		# plot results 
		# pdf('IPM_output_v3.pdf',h=8,w=12)
		par(mfrow=c(2,3)) 
		image(y,y,t(K), xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Kernel")
		contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
		plot(y,stable.dist,xlab="Size",type="l",main="Stable size distribution")
		plot(y,repro.val,xlab="Size",type="l",main="Reproductive values") 
		image(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity")
		image(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity")
		image(y,y,t(P))
		points(d$size,d$sizeNext)	
		# dev.off()
		# system('open IPM_output_v3.pdf')
		
	# compare this figure to the first model (you may have saved it as 'IPM_output_v1') Compared to the model in Fig. IPM_output_v2, the ridge in the matrix corresponding to growth in Fig. IPM_output_v3 now has a parabolic shape. The new model prevents individuals from growing continuously and sets an (asymptotic) upper limit on size where the growth regression crosses the 1:1 line. Because individuals can not growth continuously in the model in Fig. IPM_output_v3, there is now more mass in the ’adult’ portion of the stable stage distribution (sizes > 4), because there are no longer very large individuals in the population, who previously were producing a very large number of offspring in the model in Fig. IPM_output_v2. Not surprisingly, there is slightly greater sensitivity/elasticity to survival of individuals reaching sizes near 6, corresponding to the largest individuals in current model.

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# D.  Incorporating flowering probability
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# The final improvement to the IPM is to incorporate a model for flowering probability. There are two issues in play here. One is that we know the biology and expect larger individuals to flower more. The other is that the data don’t support the form of the original model. More generally, it is not a good idea to combine multiple processes into a single function, as this makes it hard to model the data. Consider the case where the plant does not flower each year, perhaps dependent on stress deriving from some environmental factor. We expect that larger plants are more likely to flower because they can maintain a sufficient amount of stored resources to buffer against environmental variation. We might not know what triggers the decision to flower, but we can at least describe its size dependence by including the probability of flowering in the model. The flowering probability function will enter the reproduction kernel (called f.xy()), defined in section 4.3. The version of f.xy() used above simply assumes that all plants flower, so including the flowering probability function will reduce this number. (Alternatively, one could use a zero-inflated model, but since we have information on flowering, it is preferable to split up the distinct processes of flowering and reproductive output conditional on flowering.) Including a flowering probability model requires 5 steps:

# incorporate flowering probability and calculate lambda. your task is as follows:

	# 1. write the flowering probability function. flowering probability is most easily modeled using logistic regression, so this will be very similar to the survival function (s.xy) above in section C.1. call your function p.flower.x(). when writing this function, note that you'll store the slope and intercept of the flowering regression as params$flower.int and param$flower.slope, respectively, as discussed below.
 
 # probability of flowering
	 p.flower.x=function(x,params) {
			??
	 }
	 
# 2. modify the reproduction function (f.xy) to include the flowering probability function. this just amounts to multiplying the argument in f.xy by p.flower.x. 
	 # reproduction function      
	 f.yx=function(xp,x,params) { 		
		 p.flower.x(x,params)* 
		 params$establishment.prob*
		 dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
		 exp(params$seed.int+params$seed.slope*x)
	 }

# 3. build a logistic regression for flowering probability. see the data frame for binary flowering data under d$fec.flower. the code for estimating the survival regression can be easily modified to accomodate the flowering probability model since both should be logistic regressions on size. from this regression, you should obtain a slope and intercept, which should be stored in the param vector as as params$flower.int and param$flower.slope, respectively.  
	 flower.reg=glm(fec.flower~size,data=d,family=binomial())
	 summary(flower.reg)
	 params$flower.int=coefficients(flower.reg)[1]
	 params$flower.slope=coefficients(flower.reg)[2]

# 4. update the regression for seed number to include only the individuals that flowered. this is easily done by using the code for the seed regression in section B.3. and changing the data argument to data=d[d$fec.flower==1,]
	 seed.reg=glm(fec.seed~size,data=d[d$fec.flower==1,],family=poisson())
	 summary(seed.reg)
	 params$seed.int=coefficients(seed.reg)[1]
	 params$seed.slope=coefficients(seed.reg)[2]

# 5. build the kernel, and plot a picture of it using your new function for f.yx. this should use code from section D and E without modification. does the peak of the fecundity kernel increase or decrease or decrease compared to your kernel from exercise 2?

		# insert code from sections D and E here
			# 2. make component kernels
	G=h*outer(y,y,g.yx,params=params) 	# growth kernel
	S=s.x(y,params=params) 							# survival 
	P=G 																# placeholder;redefine P on the next line
	for(i in 1:n) P[,i]=G[,i]*S[i]  		# growth/survival kernel
	F=h*outer(y,y,f.yx,params=params) 	# reproduction kernel	
	K=P+F 															# full kernel

	# get lamda,v,w  
	(lam=Re(eigen(K)$values[1])) 
	w.eigen=Re(eigen(K)$vectors[,1])
	stable.dist=w.eigen/sum(w.eigen) 
	v.eigen=Re(eigen(t(K))$vectors[,1])
	repro.val=v.eigen/v.eigen[1] 

	# compute elasticity and sensitivity matrices
	v.dot.w=sum(stable.dist*repro.val)*h
	sens=outer(repro.val,stable.dist)/v.dot.w
	elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)

	# plot results 
	# pdf('IPM_output_v4.pdf',h=8,w=12)
	par(mfrow=c(2,3)) 
	image(y,y,t(K), xlab="Size (t)",ylab="Size (t+1)",col=topo.colors(100), main="Kernel")
	contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
	plot(y,stable.dist,xlab="Size",type="l",main="Stable size distribution")
	plot(y,repro.val,xlab="Size",type="l",main="Reproductive values") 
	image(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity")
	image(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity")
	image(y,y,t(P))
	points(d$size,d$sizeNext)	
	# dev.off()
	# system('open IPM_output_v4.pdf')
	
	# Note that λ decreases slightly compared to the model in Fig. IPM_output_v3, because fewer individuals reproduce. Compared to Fig. IPM_output_v3, the peak in the matrix and sensitivity is shifted toward individuals reaching a size near 6, where the flowering probability function asymptotes to 1 (Fig. IPM_output_v4). The stable size distribution shifts toward having more mass in established individuals (size>4) compared to Fig. IPM_output_v3. In the model in Fig. IPM_output_v3, small individuals reproduced 100% of the time; because fewer small individuals reproduce in the model in Fig. IPM_output_v4, a greater number of large individuals compensate for the lost reproductive output in the stable size distribution.