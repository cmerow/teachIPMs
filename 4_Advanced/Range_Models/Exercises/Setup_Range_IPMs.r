# Note - The following functions aren't pretty and are not designed to work generally for other analyses. They are custom tailored for the garlic mustard models, often based on however impatient I was feeling at the time...

library(lme4)
library(MuMIn)
library(boot)
library(ez)
library(MCMCglmm)
library(IPMpack)
library(gplots)
library(fields)
library(coda)
library(Hmisc)
library(foreign)
library(sp)
library(raster)
library(maptools)
library(foreach)
library(doSNOW)
library(Matrix)
library(latticeExtra)
library(ROCR)
library(MCMCpack)
library(rasterVis)
library(rgdal)
library(reshape2)
library(mgcv)
library(dismo)
library(doParallel)

#== map of predictors
(load('ne.env.standardized.2.18.rdata'))

#== for plotting
cols1=function(x,bias=1) { colorRampPalette(c('steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x) 
}

#== group similar environmental conditions to speed up predictions. 
#== Predictions across large landscapes can be very slow, and there's no need to repeat demographic simulations for different cells with the same environment. This hack finds just the unique environments, and then redistributes the predictions back to the appropriate cells.
best.var=c('PAR','N','Ph.ave','mp.may','n.droughts.gs','mt.warm.month')
not.nas=complete.cases(values(ne.env)[,c("mp.may","n.droughts.gs","mat")]) 
a=values(ne.env)[not.nas,c('mt.warm.month','mp.may')]
# plot(a[,1],a[,2])
n.matrix1= 50; minSize1=min(a[,1])-.0001; maxSize1=max(a[,1])+.0001
b1=minSize1+c(0:n.matrix1)*(maxSize1-minSize1)/n.matrix1 
y1=0.5*(b1[1:n.matrix1]+b1[2:(n.matrix1+1)])

n.matrix2= 50; minSize2=min(a[,2])-.0001; maxSize2=max(a[,2])+.0001
b2=minSize2+c(0:n.matrix2)*(maxSize2-minSize2)/n.matrix2
y2=0.5*(b2[1:n.matrix2]+b2[2:(n.matrix2+1)])

# i didn't quite do this right; it gets assigned the largest value that its greater than and not the midpoint of an interval
ct=findInterval(a[,1],b1)
cp=findInterval(a[,2],b2)
cells=unique(cbind(ct,cp))
cells.index=which(!duplicated(cbind(ct,cp)))

# this is the thing to supply to predict functions
unique.env=data.frame(y1[cells[,1]],y2[cells[,2]])
names(unique.env)=c('mt.warm.month','mp.may')
tmp1=paste(ct,cp,sep='.')
tmp2=paste(cells[,1],cells[,2],sep='.')
# this is the thing to put the predicted values back on the map
env.cell.id=match(tmp1,tmp2)


#=======================================================================
#=======================================================================
#== PLOTTING FUNCTIONS
#=======================================================================
#=======================================================================

#== Growth =============================================================

#=======================================================================
gr.diagnostics=function(gr.reg,d1,...){
  par(mfrow=c(1,3),mar=c(5, 4, 4, .5))

	#== size(t) vs size(t+1) 
	plot(d1$size,d1$sizeNext, xlab='size(t)', ylab='size(t+1)',main='Growth Data',las=1,pch=19,cex=.7,xlim=c(min(d1$size,na.rm=T), max(d1$sizeNext,na.rm=T)),ylim=c(min(d1$size,na.rm=T), max(d1$sizeNext,na.rm=T)))
	abline(0,1,lty=2,col='red',lwd=3)

	#== predicted vs observed 
  gr.pred=predict(gr.reg,marginal=NULL,interval='confidence')
	plot(d1$sizeNext,gr.pred[,'fit'],xlab='observed', ylab='predicted',main='Growth Predictions',cex=.7, xlim=c(minSize, maxSize),ylim=c(minSize, maxSize),las=1) 
	errbar(d1$sizeNext,gr.pred[,'fit'],gr.pred[,'upr'],gr.pred[,'lwr'], add=T)
	abline(0,1,lty=2,col='red',lwd=3)
	
	#== mean growth curves at each site
  plot(0,0,col='white',xlab='size(t)',ylab='size(t+1)', main='Growth curves for each plot',xlim=c(minSize, maxSize),ylim=c(minSize, maxSize),las=1,cex=.7)
  for(i in 1:nrow(mean.envs)){
    gr.pred.new=gr.mean(gr.reg,data.frame(size= seq(minSize,maxSize,by=.1),mean.envs[i,]),...)
    lines(seq(minSize,maxSize,by=.1), gr.pred.new,col=cols1(nrow(mean.envs))[i],lwd=2)
  } 
  abline(0,1,lty=2,col='grey30',lwd=3)
  
  points(d1$size, d1$sizeNext, pch = 19)
  abline(0,1,lty=2,col='grey30',lwd=3)

}

#=========================================================================
gr.map=function(gr.reg,d,main='',df,zlims=range(values(pred.map),na.rm=T)){
 pred.map=ne.env[[1]] # dummy
 not.nas=complete.cases(values(ne.env)[,c("mp.may","n.droughts.gs","mat")]) 
 not.nas[is.na(not.nas)]=FALSE
 values(pred.map)[!not.nas]=NA
 preds=gr.mean(gr.reg,df)
 # make increment
 preds=preds-mean(df$size) 
 preds[preds<zlims[1]]=zlims[1]
 preds[preds>zlims[2]]=zlims[2]
 values(pred.map)[not.nas]=preds

 image.plot(pred.map, col=c(cols1(100,1)),xaxt='n', yaxt='n',xlab=' ',ylab=' ', bty='n',main=' ', zlim=zlims, legend.args=list(text="growth\nincrement",cex=.8))
  mtext(main,3,line=-2.8,at=-72.2)
}
  
#== Survival============================================================
#=========================================================================
#== survival curves 
sv.diagnostics=function(sv.reg,d1,marginal='none',form){
	plot(0,0, type='l',col='white',lwd=2,xlab='size',ylab='probability of survival',main='Survival curves for each plot',ylim=c(0,1),xlim=c(minSize,maxSize))
	for(i in 1:nrow(mean.envs)){
		sv.pred.new=sv(sv.reg,covariates=data.frame(size= seq(minSize,maxSize,by=0.1),mean.envs[i,],surv=1e4),form=form)
		points(seq(minSize,maxSize,by=0.1), sv.pred.new, type = "l", col=cols1(nrow(mean.envs))[i],lwd=2)
	} 
	rug(d1$size[d1$surv==0])
	rug(d1$size[d1$surv==1],side=3)
	# add points for whole data set
	ncuts=15
	os <- order(d1$size)
	os.surv <- (d1$surv)[os]
	os.size <- (d1$size)[os]
	psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
	ps <- tapply(os.surv, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
	points(as.numeric(psz), as.numeric(ps), pch = 19)
}

#=========================================================================
#== map survival
sv.map=function(reg,main='',new.df,random=FALSE,label="survival\nprobability",zlims=NA,form){
 pred.map=ne.env[[1]] # dummy
 not.nas=complete.cases(values(ne.env)[,c("mp.may","n.droughts.gs","mat")]) 
 not.nas[is.na(not.nas)]=FALSE
 values(pred.map)[!not.nas]=NA
 preds=sv(reg,new.df,random=FALSE,form)
 values(pred.map)[not.nas]=preds
 if(is.na(zlims[1])){zlims=range(values(pred.map),na.rm=T)}
 
 image.plot(pred.map, col=c(cols1(100,1)),xaxt='n', yaxt='n',xlab=' ',ylab=' ', bty='n',main=' ', zlim=zlims, legend.args=list(text=label,cex=.8))
  mtext(main,3,line=-2.8,at=-72.2)
}

#== Fecundity ============================================================
#=========================================================================
#== predicted vs observed plots
fec.diagnostics=function(shrub=FALSE){
  # predicted vs observed plots seed number
  par(mfrow=c(3,2),mar=c(5, 4, 4, 2) )
  fec1.pred=seed2(d.seed$size,apply(fec1.reg,2,mean),d.seed,seed.form)
  plot(d.seed$fec1,fec1.pred, xlab='observed', ylab='predicted', main='Seed Number (Log Scale)', ylim=c(10,max(c(fec1.pred,d$fec1),na.rm=T)), xlim=c(10,max(c(fec1.pred,d$fec1),na.rm=T)),log='xy')
  abline(0,1,lty=2,col='red',lwd=3)

  fec2.pred=flower(10000,apply(fec2.reg,2,mean),d[!is.na(d$fec3),])
  plot(d$fec2[!is.na(d$fec2)],fec2.pred, xlim=c(0,.7),ylim=c(0,.7), xlab='observed',ylab='predicted',main='Germination Probability')

  abline(0,1,lty=2,col='red',lwd=3)

  if(!shrub){
		fec3.pred=flower(10000,apply(fec3.reg,2,mean),d[!is.na(d$fec3),])
		plot(d$fec3[!is.na(d$fec3)],fec3.pred, xlim=c(0,1),ylim=c(0,1), xlab='observed',ylab='predicted',main='Germinant Survival')
		abline(0,1,lty=2,col='red',lwd=3)
	}

  fec1.pred=seed2(d.seed$size,apply(fec1.reg,2,mean),d.seed,seed.form)

  plot(d.seed$size,d.seed$fec1, xlab='size',ylab='log(seed #)',main='Seed number curves\nfor each plot',log='y',pch=19,cex=.7,las=1,cex.lab=1.5,cex.axis=1.5)
  for(i in 1:nrow(mean.envs)){
    fec1.pred.new=seed2(size=seq(minSize,maxSize,by=.1),seed.params= apply(fec1.reg,2,mean), covariates= mean.envs[i,],seed.form,random=FALSE)
    points(seq(minSize,maxSize,by=.1), fec1.pred.new, type = "l", col=cols1(nrow(mean.envs))[i],lwd=2)
  } 
  
    # size distribution of recruits
  hist(d$size[d$Age==1],breaks=20,freq=F,xlab='size',main='Recruit Size Distribution',las=1,ylim=c(0,.5),cex.lab=1.5,cex.axis=1.5)
  lines(seq(minSize,maxSize,by=.1),dnorm(seq(minSize,maxSize,by=.1), mean(offspr.reg$Sol),mean(sqrt(offspr.reg$VCV))),col='red',lwd=2)
  
}

#=========================================================================
fec23.preds=function(shrub=FALSE){
  # germination by plot
  par(mfrow=c(2,1),mar=c(4, 4, 2, 2) )
  fec2.pred.new=pred.mc(fec2.reg,newdata=mean.envs, interval='confidence', marginal=0)
  x.loc=1:nrow(mean.envs)+.1
  plot(x.loc,fec2.pred.new[,1],pch=20,col='red', xlab='plot',ylab='germination probability',main='Germination Probability',ylim=c(0,.6))
  abline(h=mean(fec2.pred.new[,1]),lty=2,lwd=3,col='grey60')
  arrows(x.loc, fec2.pred.new[,3], x.loc, fec2.pred.new[,2], angle = 90, code = 3, length = 0.08)
  indiv.means=aggregate(d[!is.na(d$fec2),c('fec2','wplot')], list(d[!is.na(d$fec2), 'wplot']),mean)
  points(x.loc-.2,indiv.means[,2],pch=20,col='blue')	
  text(4,.5,'random means',col='red')
  text(4,.45,'plot means',col='blue')
  text(4,.4,'group mean',col='grey60')
    
  # summer survival plot
  if(!shrub){
		fec3.pred.new=pred.mc(fec3.reg,newdata=mean.envs, interval='confidence', marginal=0)
		x.loc=1:nrow(mean.envs)+.1
		plot(x.loc,fec3.pred.new[,1],pch=20,col='red', xlab='plot',ylab='summer survival probability',main='Germinant survival through summer',ylim=c(0,1))
		arrows(x.loc, fec3.pred.new[,3], x.loc, fec3.pred.new[,2], angle = 90, code = 3, length = 0.08)
		indiv.means=aggregate(d[!is.na(d$fec3),c('fec3','wplot')], list(d[!is.na(d$fec3), 'wplot']),mean)
		abline(h=mean(indiv.means[,2]),lty=2,lwd=3,col='grey60')
		points(x.loc-.2,indiv.means[,2],pch=20,col='blue')	
		text(4,.9,'random means',col='red')
		text(4,.8,'plot means',col='blue')
		text(4,.7,'group mean',col='grey60')
	}  
}

#=========================================================================
#== map fec models
seed.map=function(reg,main='',df,label,zlims=NA,seed.form){
 pred.map=ne.env[[1]] # dummy
 not.nas=complete.cases(values(ne.env)[,c("mp.may","n.droughts.gs","mat")]) 
 not.nas[is.na(not.nas)]=FALSE
 values(pred.map)[!not.nas]=NA

 preds=seed2(size=df$size,seed.params=apply(reg,2,mean), covariates= df[,-which(names(df)=='size')],seed.form,random=FALSE)
 if(!is.na(zlims[1])){preds[preds>zlims[2]]=zlims[2]}
 values(pred.map)[not.nas]=preds
 if(is.na(zlims[1])){zlims=range(values(pred.map),na.rm=T)}
 
 image.plot(pred.map, col=c(cols1(100,1)),xaxt='n', yaxt='n',xlab=' ',ylab=' ', bty='n',main=' ', zlim=zlims, legend.args=list(text=label,cex=.8))
  mtext(main,3,line=-2.8,at=-72.2)
}
 
#=========================================================================

fl.diagnostics=function(...){ 
  par(mfrow=c(1,1),mar=c(5, 4, 4, 2) )
  plot(0,0, type='l',col='white',lwd=2,xlab='size',ylab='probability of flowering',main='Flowering curves for each plot',ylim=c(0,1),xlim=c(minSize,maxSize))
  for(i in 1:nrow(mean.envs)){
  	fl.pred.new=sv(fl.reg,covariates=data.frame(size= seq(minSize,maxSize,by=0.1),mean.envs[i,],flowering=1e4),...)
    points(seq(minSize,maxSize,by=.1), fl.pred.new, type = "l", col=cols1(nrow(mean.envs))[i],lwd=2)
    rug(d$size[d$flowering==0 & d$wplot==i])
    rug(d$size[d$flowering==1 & d$wplot==i],side=3)
    rug(d$size[d$flowering==1 & is.na(d$wplot==i)],side=3)
    rug(d$size[d$flowering==0 & is.na(d$wplot==i)])
  } 
  # add points for whole data set
  ncuts=15
  os.flowering <- (d$flowering)[order(d$size)]
  os.size <- (d$size)[order(d$size)]
  psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
  ps <- tapply(os.flowering, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
  points(as.numeric(psz), as.numeric(ps), pch = 19)	
} 
 
#=========================================================================
#== nice versions of response curves 
response.summary=function(mean.envs,grad='mt.warm.month',lab=letters[1:6],titl='',shrub=FALSE){
	cols=colorRampPalette(c('steelblue4','steelblue1','red1','red4'))(20)
	temp.mean.envs=mean.envs[order(mean.envs[,grad]),]
	col2=cols[round(rank(temp.mean.envs[,grad]))]

	par(mfrow=c(7,1),mar=c(4,5.5,1,1),oma=c(1,1,1,0))
# a - growth
  plot(0,0,col='white',xlab='size',ylab='size(t+1)', main='  ',xlim=c(minSize, maxSize),ylim=c(minSize, maxSize),las=1,cex=.7,bty='n',cex.lab=1.4)
  mtext(lab[1],3,at=minSize)
  for(i in 1:nrow(temp.mean.envs)){
    gr.pred.new=gr.mean(gr.reg,data.frame(size= seq(minSize,maxSize,by=.1),temp.mean.envs[i,]),random=FALSE)
    lines(seq(minSize,maxSize,by=.1), gr.pred.new,col=col2[i],lwd=2,lty=temp.mean.envs$habitat[i]+1)
  } 
 
# b - surv
	plot(0,0, type='l',col='white',lwd=2,xlab='size',ylab='survival probability',main=' ', cex.lab=1.4,ylim=c(0,1),xlim=c(minSize,maxSize),las=1,bty='n')
	mtext(lab[2],3,at=minSize)
	for(i in 1:nrow(temp.mean.envs)){
		sv.pred.new=sv(sv.reg,covariates=data.frame(size= seq(minSize,maxSize,by=0.1),temp.mean.envs[i,],surv=1e4),form=sv.form)
		points(seq(minSize,maxSize,by=0.1), sv.pred.new, type = "l", col=col2[i],lwd=2,lty=temp.mean.envs$habitat[i]+1)
	} 
	
  
#c flower
  plot(0,0, type='l',col='white',lwd=2,xlab='size',ylab='flowering\n probability',main=' ',cex.lab=1.4 ,ylim=c(0,1),xlim=c(minSize,maxSize),bty='n',las=1)
  mtext(lab[3],3,at=minSize)
  for(i in 1:nrow(mean.envs)){
    fl.pred.new=sv(fl.reg,covariates=data.frame(size= seq(minSize,maxSize,by=0.1),temp.mean.envs[i,],flowering=1e4),form=fl.form)
    points(seq(minSize,maxSize,by=.1), fl.pred.new, type = "l", col=col2[i],lwd=2,lty=temp.mean.envs$habitat[i]+1)
  } 
  
#d seed  
  plot(0,0, type='l',col='white',lwd=2,xlab='size',ylab='seed number',main=' ', cex.lab=1.4,ylim=c(10,600),xlim=c(minSize,maxSize),bty='n',las=1)
  mtext(lab[4],3,at=minSize)
  for(i in 1:nrow(mean.envs)){
    fec1.pred.new=seed2(size= seq(minSize,maxSize,by=0.1),apply(fec1.reg,2,mean),temp.mean.envs[i,],seed.form)
    points(seq(minSize,maxSize,by=0.1), fec1.pred.new, type = "l", col=col2[i],lwd=2,lty=temp.mean.envs$habitat[i]+1)
  } 

#e germ
	germ.pred.new=flower(1e4,apply(fec2.reg,2,mean),temp.mean.envs)
	barplot(germ.pred.new,width=.1,col=col2,beside=T, density=(temp.mean.envs$habitat+1)*50,main=' ',ylab='germination\nprobability',las=1,xlab='',cex.lab=1.4)
	mtext('plot',1,line=2)
	mtext(lab[5],3,at=.05)
	
if(!shrub){
#f germ survival
	sum.sv.pred.new=flower(10000,apply(fec3.reg,2,mean),temp.mean.envs)
	barplot(sum.sv.pred.new,width=.1,col=col2,beside=T, density=(temp.mean.envs$habitat+1)*50,main=' ',ylab='germinant\nsurvival probability',las=1,xlab='',cex.lab=1.4)#,names.arg=temp.mean.envs$Plot)
	mtext('plot',1,line=2)
	mtext(lab[6],3,at=.05)
}	
	# add color bar  
	plot(0,0, type='l',col='white',xlab=' ',ylab='',main=' ',ylim=c(10,600),xlim=c(minSize,maxSize),bty='n',las=1,xaxt='n', yaxt='n')
  tmp=values(ne.env)[,grad]
  tmp[tmp< min(temp.mean.envs[,grad])]=min(temp.mean.envs[,grad])
  tmp[tmp>max(temp.mean.envs[,grad])]=max(temp.mean.envs[,grad])
  tmp2=ne.env[[grad]]
  values(tmp2)=tmp
  image.plot(tmp2,add=TRUE,zlim=range(temp.mean.envs[,grad]), axis.args = list(cex.axis = 1.2),smallplot= c(.2,.9,.8,.95),legend.lab=titl,legend.line=3,col=col2,horizontal=TRUE)
}	


#========================================================================  
#========================================================================  
#== predict functions
#========================================================================  
#========================================================================  

#=========================================================================
# need to specify which type of random effects
gr.mean=function(reg,covariates,random=TRUE){
	temp.form=as.formula(paste('~',as.character(reg$Fixed$formula)[3]))
	design=model.matrix(temp.form, data=covariates)
	if(random) design= as.matrix(cbind(design,covariates[,grep('wplot.',names(covariates))]))
	design%*%apply(reg$Sol,2,mean)
} 

#=========================================================================
#== for making kernels 
gr2=function(size,size.next,mean.params,sd.param,covariates,gr.form, random=FALSE){
	covariates=data.frame(Intercept=1,size=size,covariates)
  design=model.matrix(gr.form, data=covariates)
  if(random) design= as.matrix(cbind(design,covariates[,grep('wplot.',names(covariates))]))
  dnorm(size.next,mean=design%*%mean.params,sd=sqrt(sd.param))
}

#=========================================================================
#== for mean prediction
sv=function(reg,covariates,random=FALSE,form){
	design=model.matrix(form, data=covariates)
	mu=design%*%apply(reg,2,mean)
	plogis(mu)

}  

#=========================================================================
#== for making kernels. 
#== note: also used for flowering
sv2=function(size,sv.params,covariates,sv.form, VCV,random=FALSE){
	covariates=data.frame(Intercept=1,size=size,covariates)
  design=model.matrix(sv.form, data=covariates)
  if(random){ # for different types of random effects
  	if( length(grep('wplot.',names(sv.params)))>0 )
  	design= as.matrix(cbind(design,covariates[,grep('wplot.',names(covariates))]))
  	if( length(grep('wregion.',names(sv.params)))>0  )
  	design= as.matrix(cbind(design,covariates[, grep('wregion.',names(covariates))]))
  }
	mu=design%*%sv.params
	plogis(mu)
}

#=========================================================================
#== for making kernels. 
#== note: also used for flowering
sv.shrub=function(size,sv.params,sv.a.params,covariates,sv.form, VCV,random=FALSE){
	covariates=data.frame(Intercept=1,size=size,covariates)
  design=model.matrix(sv.form, data=covariates)
  if(random){ # for different types of random effects
  	if( length(grep('wplot.',names(sv.params)))>0 )
  	design= as.matrix(cbind(design,covariates[,grep('wplot.',names(covariates))]))
  	if( length(grep('wregion.',names(sv.params)))>0  )
  	design= as.matrix(cbind(design,covariates[, grep('wregion.',names(covariates))]))
  }
	mu=design%*%sv.params
	mu=pmin(mu,sv.a.params)
	plogis(mu)
}

#=========================================================================
#== for making kernels 
offspr2=function(size,mean.params,sd.param,covariates,offspr.form){ 
  covariates=data.frame(Intercept=1,size=size,covariates)
	design=model.matrix(offspr.form, data=covariates)
  dnorm(size,mean=design%*%mean.params,sd=sqrt(sd.param))
}

#=========================================================================
#== seedhead function used to build IPM     
#== for making kernels
seed2=function(size,seed.params,covariates,seed.form,random=TRUE,VCV){ 
	covariates=data.frame(Intercept=1,size=size,covariates,fec1=1e4)
  design=model.matrix(seed.form, data=covariates)
  exp(design%*%seed.params)
} 

#=========================================================================
#== flowering function, used for building IPM below
flower=function(size,params,covariates,flower.form=NULL) { 
	if(is.null(flower.form)) flower.form=as.formula(paste('~',paste(names(params)[-1], collapse='+')))		
	design=model.matrix(flower.form, data=data.frame(size,covariates))
	plogis(design%*%params)
}

#=========================================================================
#== fast germ function used to build IPM     
germ2=function(germ.params,covariates,germ.form,random=TRUE){ 
	#covariates=data.frame(Intercept=1,size=size,covariates)
  design=model.matrix(germ.form, data=covariates)
  if(random){ #haven't checked this, as it's not used
  		if( length(grep('wplot.',names(germ.params)))>0 )
  		design= as.matrix(cbind(design,covariates[,grep('wplot.',names(covariates))]))
  		if( length(grep('wregion.',names(germ.params)))>0  )
  		design= as.matrix(cbind(design,covariates[, grep('wregion.',names(covariates))]))
  }
  u=exp(design%*%germ.params)
  u/(1+u)
} 

#=========================================================================
#== fast germinant survival function used to build IPM 
sum.sv2=function(sum.sv.params,covariates,sum.sv.form){ 
  design=model.matrix(sum.sv.form, data=covariates)
  u=exp(design%*%sum.sv.params)
	u/(1+u)
} 

#=========================================================================
#=========================================================================
# KERNEL FUNCTIONS ==========================================================
#=========================================================================
#=========================================================================

#=========================================================================
#=========================================================================
shrink.matrix=function(m,min=1e-5){
	m[m<min]=0
	m=as(m,'sparseMatrix')
}

#=========================================================================
#=========================================================================
chop.tasks=function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
#=========================================================================
#=========================================================================
#== make p kernel for mustard
p.mustard.kernel=function(task.split1){ post.temp=list(p=lapply(1:length(task.split1), function(i) 1))
	for(k in 1:length(task.split1)){
		ind=task.split1[k]
		tenv=data.frame(t(env.subset[ind,]))
		# growth kernel
		P=h*t(outer(y,y,gr2,gr.params.mean,gr.params.sd,tenv,gr.form, random=FALSE))
		# survival & not flowering
		S=sv2(y,sv.params,tenv,sv.form) 	
		# growth/survival kernel   
		Ps=matrix(rep(apply(P,2,sum),n.matrix),byrow=T,nrow=n.matrix)
		Ss=matrix(rep(S,n.matrix),byrow=T,nrow=n.matrix)
		P=(Ss*P)/Ps #correct eviction
		post.temp[[k]]=shrink.matrix(P)
	}
	post.temp
}
#=========================================================================
#=========================================================================
#== make f kernel for mustard
f.mustard.kernel=function(task.split,manage='',random=TRUE){
	if(manage=='') manage=1
	 post.temp=list(F=lapply(1:length(task.split), function(x) 1))
	 for(k in 1:length(task.split)){
		 ind=task.split[k]
		 tenv=data.frame(t(env.subset[ind,]),new.1.0=1e4)
		 offspr.size=offspr2(size=y,mean.params=offspr.params.mean, sd.param=offspr.param.sd, tenv,offspr.form)
		 tf=manage*matrix(rep(sv2(y,fl.params,tenv,fl.form,random=random), n.matrix), byrow=T,nrow=n.matrix)
		 tsh=matrix(rep(seed2(y,seed.params,tenv,seed.form,random=random), n.matrix), byrow=T,nrow=n.matrix)
		 germ=germ2(germ.params,tenv,germ.form,random=random)
		 sum.sv=sum.sv2(sum.sv.params,tenv,sum.sv.form)
		 F=h*offspr.size*germ*sum.sv*tf*tsh
		 post.temp[[k]]=shrink.matrix(F,1e-5)
	 }
	 post.temp
}

#=========================================================================
#=========================================================================
#== dynamics for mustard
dynamics.mustard=function(task.split,useful.calcs=FALSE,maybe.useful.calcs=FALSE){
	#==define things to save
	lam=r0=sens1=elas1=passage.time=stable=r.val=le= list(lapply(1:length(task.split), function(i) 1))
	for(i in 1:length(task.split)){
		ind=task.split[i]  
		#== for age 
		IPM=matrix(0,2*n.matrix,2*n.matrix)
		IPM[(n.matrix+1):(2*n.matrix), 1:n.matrix ]=as.matrix((p.post.mean[[ind]]))
		IPM[  1:n.matrix, (n.matrix+1):(2*n.matrix)] =as.matrix((f.post.mean[[ind]]))
		lam[[i]]=max(Re(eigen(IPM)$values))
		if(useful.calcs){
			r0[[i]]=R0Calc(as.matrix(p.post.mean[[ind]]),as.matrix(f.post.mean[[ind]]))
			sens1[[i]]=sens(IPM)
			elas1[[i]]=elas(IPM)
		}
		if(maybe.useful.calcs){
      stable[[i]]=abs(Re(eigen(IPM)$vectors[, 1]))
 			r.val[[i]]=abs(Re(eigen(t(IPM))$vectors[, 1]))
 			le[[i]]=meanLifeExpect(as.matrix(p.post.mean[[ind]]))
		}
	}
	list(lam=lam,r0=r0,sens1=sens1,elas1=elas1,passage.time=passage.time, stable=stable,r.val=r.val,le=le)
}

#=========================================================================
#=========================================================================
# plot kernels
plot.kernel=function(p.post.mean,flatten=1,main='',round.dig=2){
	par(mfrow=c(4,6),mar=c(1,1,3.5,1),oma=c(4,4,2,0))
	mat=lapply(p.post.mean,as.matrix)
	mat2=lapply(p.post.mean,as.matrix)

	for(i in 1:length(p.post.mean)) { 
		IPM.rotated.plot(mat[[i]],y,titles= paste('plot',i,' region=',mean.envs[i,1],'\nplot=',mean.envs[i,2],'  ', ifelse(mean.envs[i,'habitat']==1,'closed','open')),flatten=flatten,zlim=c(0,max(unlist(mat2))^flatten),minSize=minSize,maxSize=maxSize)
	}
	mtext(main,3,outer=TRUE)
	mtext(paste('flatten=',flatten),3,outer=TRUE,at=.75,cex=.7)
	mtext('Size (t+1)',2,outer=TRUE,line=2)
	mtext('Size (t)',1,outer=TRUE,line=2)
	par(oma=c(3.5,0,0,7))
	image.plot(legend.only=T,zlim=c(0,max(unlist(mat2))^flatten), col=cols1(100), legend.width=4,axis.args = list(cex.axis = 1,at=seq(0,max(unlist(mat2))^flatten,length=10),labels=round(seq(0,max(unlist(mat2)),length=10),round.dig))) 
}

#=========================================================================
#=========================================================================
pop.map=function(preds,main='',label='',cutoff=1,max.val=Inf,zlims=NA,pts=NA,cols=c(grey(seq(.7,.1,length=100)),cols1(100,1)),year=''){
 preds[preds>max.val]=max.val
 if(year=='') land=ne.env
 if(year==2041) land=ne.env.2041
 if(year==2090) land=ne.env.2090
 pred.map=land[[1]] # dummy
 not.nas=complete.cases(values(land)[,c("mp.may","mt.warm.month")]) 
 not.nas[is.na(not.nas)]=FALSE
 values(pred.map)[!not.nas]=NA
 values(pred.map)[not.nas]=preds
 if(is.na(zlims[1])) zlims=range(values(pred.map),na.rm=T)
 image.plot(pred.map, col=cols,xaxt='n', yaxt='n',xlab=' ',ylab=' ', bty='n',main=' ', zlim=zlims, legend.args=list(text=label,cex=.8),breaks=c(seq(0,cutoff-1e-4,len=length(cols)/2),seq(cutoff,zlims[2],len=length(cols)/2+1)))
 if(!is.na(pts)){
 	 points(pts,pch=21,col='grey95',lwd=1.3,cex=.5)
 	 valid=extract(pred.map,pts)
 }
  mtext(main,3,line=-3.8,at=-72)
}
#=========================================================================
#=========================================================================
# make binomial count data binary for use with logistic regression packages
bernoullize=function(dat,col.name1,col.name0){
	#takes a data frame and duplicates rows to turn binomial counts to binary data
	dat.tmp=dat[which(dat[,col.name1]>0),]
	out=dat.tmp[rep(row.names(dat.tmp), dat.tmp[,col.name1]), ]
	out[,'new.1.0']=1
	dat.tmp=dat[which(dat[,col.name0]>0),]
	out2=dat.tmp[rep(row.names(dat.tmp), dat.tmp[,col.name0]), ]
	out2[,'new.1.0']=0
	rr=rbind(out,out2)
}
	

#=========================================================================
#=========================================================================
#== geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
