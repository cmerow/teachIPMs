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
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# Here's the piece you need !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			g.yx=function(xp,x,params) { 			
				dnorm(xp,
					mean=params$growth.int+params$growth.slope*x+params$growth.sqrd*x^2,
					sd=sqrt(params$growth.sigma2*exp(2*params$growth.sigma2.exp*x)))
			}
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

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
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# Here's the piece you need !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   # probability of flowering
	   p.flower.x=function(x,params) {
		   u=exp(params$flower.int+params$flower.slope*x)
		   return(u/(1+u))
	   }
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

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