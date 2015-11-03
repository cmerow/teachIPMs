setwd('~/Dropbox/Projects/ipms/teachIPMs/Advanced/Range_Models/Exercises/') # set this to your teachIPMs directory

# Note that a number of functions that simplify the following exercise are contained in this file. You don't really need to look at it unless you're going to run these analyses on your own data.
source("Setup_Range_IPMs.r")

############################################################################
#== OUTLINE
# 1. Set up data
# 2. Growth Model
# 3. Survival Model
# 4. Fecundity Models
# 5. Flowering Probability Model
# 6. Summary Plots
# 7. Make Mean Kernels
# 8. IPM diagnostics 

############################################################################
# 1. Set up data -----------------------------------------------------------

d=read.csv('garlic_mustard_data_frame_2_18.csv')
#== mean.envs are the average environmental conditions (over individuals) at each of the 21 common gardens where data were collected.
mean.envs=read.csv('garlic_mustard_plot_attributes.csv')
minSize=4 # of individuals in the IPM
maxSize=15

#== you may later find it useful to reload saved output from the models below if you need to close your R session. But for now, it's commented out.
# (load('Model_Output/AP_growth_v5.post'))
# (load('Model_Output/AP_surv_v5.post'))
# (load('Model_Output/AP_fec_v5.post'))
# (load('Model_Output/AP_flowering_v5.post'))
# gr.form= as.formula(paste('~',as.character(gr.reg$Fixed$formula)[3]))
# sv.form= as.formula(paste('~',paste0(dimnames(sv.reg)[[2]][-1],collapse='+')))
# seed.form= as.formula(paste('~',paste0(dimnames(fec1.reg)[[2]][-1],collapse='+')))
# germ.form= as.formula(paste('~',paste0(dimnames(fec2.reg)[[2]][-1],collapse='+')))
# sum.sv.form= as.formula(paste('~',paste0(dimnames(fec3.reg)[[2]][-1],collapse='+')))
# fl.form= as.formula(paste('~',paste0(dimnames(fl.reg)[[2]][-1],collapse='+')))
# offspr.form=as.formula(paste('~',as.character(offspr.reg$Fixed$formula)[3]))


############################################################################
#== 2. Growth Model --------------------------------------------------------
#== subset the data to just the cases that enter the growth model
d.temp=d[complete.cases(d[,c('sizeNext','size')]),]

gr.form= sizeNext~size+PAR+N+Ph.ave+mt.warm.month+mp.may
gr.reg=MCMCglmm(gr.form, data=d.temp, verbose=FALSE,burnin=3000, nitt=13000,thin=10)
summary(gr.reg)

#== save for later use
save(gr.reg,file='Model_Output/AP_growth_v5.post')
# (load('Model_Output/Posteriors/AP_growth_v2.post'))

#== explore output
#== if you're familiar with Bayesian models, check these out. If not the summary() call above tells you what you need to know about the model for this exercise.
posterior.mode(gr.reg$Sol)
coda::HPDinterval(gr.reg$Sol, 0.95) 
autocorr.plot(gr.reg$Sol)
plot((gr.reg$Sol)) 

#== diagnostics 
#== plot data and predictions 
pdf('Figures/AP_growth_diagnostics.pdf',h=3.5,w=10)
gr.diagnostics(gr.reg,d.temp,random=FALSE)
dev.off()
system("open Figures/AP_growth_diagnostics.pdf")

# map of growth predictions
pdf('Figures/AP_growth_map.pdf',h=5,w=10)
par(mfrow=c(1,2),oma=c(0,0,2,4))
which.subset=mean.envs$habitat==0 # choose closed canopy habitat, since the garlic mustard performs very differently in closed vs open habitat.
newdata=data.frame(size=quantile(d$size,.01,na.rm=T), PAR=mean(mean.envs$PAR[which.subset]),N=mean(mean.envs$N[which.subset]),Ph.ave=mean(mean.envs$Ph.ave[which.subset]),sm.t1=mean(mean.envs$sm.t1[which.subset]),values(ne.env)[not.nas,]) # data frame describing the landscape where you're predicting
gr.map(gr.reg,d,'AP\nclosed\ncanopy',newdata,zlims=c(0,2))
which.subset=mean.envs$habitat==1
newdata=data.frame(size=quantile(d$size,.01,na.rm=T), PAR=mean(mean.envs$PAR[which.subset]),N=mean(mean.envs$N[which.subset]),Ph.ave=mean(mean.envs$Ph.ave[which.subset]),sm.t1=mean(mean.envs$sm.t1[which.subset]),values(ne.env)[not.nas,])
gr.map(gr.reg,d,'AP\nopen\ncanopy',newdata,zlims=c(0,2))
dev.off()
system("open Figures/AP_growth_map.pdf")

############################################################################	
#== 3. Survival Model ------------------------------------------------------
d.sv=d[complete.cases(d[,c('surv')]),]

#== checking for correlated predictors, as this model has trouble converging.
round(cor(d.sv[,best.var]),2)

sv.form=surv~size+PAR+N+Ph.ave+mt.warm.month+mp.may
sv.reg=MCMClogit(sv.form, data=d.sv,b0=0, B0=.001,mcmc=5e4,thin=50)
summary(sv.reg)
save(sv.reg,file='Model_Output/AP_surv_v5.post')  
# (load('Model_Output/AP_surv_v3.post'))

#== explore output
# posterior.mode(sv.reg$Sol)
coda::HPDinterval(sv.reg, 0.95)
autocorr.plot(sv.reg)
plot((sv.reg)) 

#== diagnostics
#== plot data and predictions 
pdf('Figures/AP_sv_diagnostics.pdf',h=4,w=5)
sv.diagnostics(sv.reg,d.sv,form=sv.form)
dev.off()
system("open Figures/AP_sv_diagnostics.pdf")

#== map of surv predictions
pdf('Figures/AP_surv_map.pdf',h=5,w=10)
par(mfrow=c(1,2),oma=c(0,0,2,4))
which.subset=mean.envs$habitat==0
newdata=data.frame(size=quantile(d$size,.5,na.rm=T), PAR=mean(mean.envs$PAR[which.subset]),N=mean(mean.envs$N[which.subset]),Ph.ave=mean(mean.envs$Ph.ave[which.subset]),sm.t1=mean(mean.envs$sm.t1[which.subset]),values(ne.env)[not.nas,],surv=1e4)
 sv.map(sv.reg,'AP\nclosed\ncanopy',newdata,form=sv.form,zlims=c(0,1))
which.subset=mean.envs$habitat==1
newdata=data.frame(size=quantile(d$size,.5,na.rm=T), PAR=mean(mean.envs$PAR[which.subset]),N=mean(mean.envs$N[which.subset]),Ph.ave=mean(mean.envs$Ph.ave[which.subset]),sm.t1=mean(mean.envs$sm.t1[which.subset]),values(ne.env)[not.nas,],surv=1e4)
 sv.map(sv.reg,'AP\nopen\ncanopy',newdata,form=sv.form,zlims=c(0,1))
dev.off()
system("open Figures/AP_surv_map.pdf")
  
############################################################################
#== 4. Fecundity Model -----------------------------------------------------
#== The fecundity model has multiple components: (A) seeds per plant, (B) germination probability, (C) germinant survival from spring to summer, and (D) germinant size distribution.	
#== A. seed number
d.seed=d[!is.na(d$fec1) & !is.na(d$size),]

seed.form=fec1~size+PAR+Ph.ave+mt.warm.month+mp.may
fec1.reg=MCMCpoisson(seed.form, data=d.seed,b0=0, B0=.001,mcmc=5e4,thin=50)  
summary(fec1.reg)

#== B. % germination
d.germ=bernoullize(d[!is.na(d$fec2),],col.name1='n.germ.1', col.name0='n.germ.0')
germ.form=new.1.0~Ph.ave+light+mt.warm.month+mp.may
fec2.reg=MCMClogit(germ.form, data=d.germ,b0=0, B0=.001,mcmc=5e4,thin=50)  
summary(fec2.reg)

#== C. germinant survival
d.germ.s=bernoullize(d[!is.na(d$fec3),],col.name1='n.germ.surv.1', col.name0='n.germ.surv.0') # this turns counts of the number of 1s and 0s to a unique row for each 1 and 0, for use with MCMClogit
sum.sv.form=new.1.0~Ph.ave+light+mt.warm.month+mp.may
fec3.reg=MCMClogit(sum.sv.form, data=d.germ.s,b0=0, B0=.001,mcmc=5e4,thin=50)
summary(fec3.reg)

#== D. germinant size
d.off=subset(d, (d$Year.planted==d$Year.size) &  is.na(d$flowering))
offspr.reg=MCMCglmm(size~1, data=d.off, burnin=3000,nitt=13000,thin=10,verbose=T)
summary(offspr.reg)

save(fec1.reg,fec2.reg,fec3.reg,offspr.reg, file='Model_Output/AP_fec_v5.post')
# (load('Model_Output/AP fec 10-3-12.post'))

#== explore output
plot((fec1.reg))
plot((fec2.reg)) 
plot((fec3.reg)) 
  
#== diagnostics ------------------------------------------------------------
#== plot data and predictions 
pdf('Figures/AP_fec_diagnostics.pdf',h=9,w=6)
	fec.diagnostics()
dev.off()
system("open Figures/AP_fec_diagnostics.pdf")

# map of fec predictions
pdf('Figures/AP_fec_map.pdf',h=5,w=15)
par(mfrow=c(1,3),oma=c(0,0,2,4))
 newdata=data.frame(size=quantile(d$size,.75,na.rm=T), PAR=1.5,N=-1.5,Ph.ave=1.5,light=1.5,values(ne.env)[not.nas,],new.1.0=1e5)
 seed.map(fec1.reg,'AP',newdata,"seed\nnumber",seed.form=seed.form, zlims=c(0,2000))
 sv.map(fec2.reg,'AP',newdata,label="germination\nprobability", form=germ.form)
 sv.map(fec3.reg,'AP',newdata,label="germinant\nsurvival\nprobability", form=sum.sv.form)
dev.off()
system("open Figures/AP_fec_map.pdf")

############################################################################
#== 5. Flowering Model -----------------------------------------------------
#== flowering probability is a logistic regression, so i just borrow the machinery from the survival functions. 
d.fl=d[!is.na(d$flowering) & !is.na(d$size),]

fl.form=flowering~size+I(size^2)+I(size^3)+PAR+N+Ph.ave+mp.may+mt.warm.month
fl.reg=MCMClogit(fl.form, data=d.fl,b0=0, B0=.001,mcmc=5e4,thin=50)
summary(fl.reg)

save(fl.reg,file='Model_Output/AP_flowering_v5.post')
# (load('Model_Output/AP_flowering_v2.post'))

# explore output
autocorr.plot(fl.reg)
plot((fl.reg)) 

# diagnostics  ----------------------------------------------------------
#== plot data and predictions 
pdf('Figures/AP_fl_diagnostics.pdf',h=6,w=6)
	fl.diagnostics(random=FALSE,form=fl.form)
dev.off()
system("open Figures/AP_fl_diagnostics.pdf")

#== map of flowering predictions
pdf('Figures/AP_fl_map.pdf',h=5,w=5)
 newdata=data.frame(size=quantile(d$size,.75,na.rm=T), PAR=1.5,N=-1.5,Ph.ave=1.5,values(ne.env)[not.nas,],flowering=1e4)
 sv.map(fl.reg,'AP\nopen canopy',newdata,label='flowering\nprobability',form=fl.form)
dev.off()
system("open Figures/AP_fl_map.pdf")

#########################################################################
#== 6. Summary Plots-----------------------------------------------------

#== plot all response curves to temperature together
pdf('Figures/AP_all_response_curves_mt_warm_month.pdf',w=4,h=12)
	response.summary(mean.envs,grad='mt.warm.month',lab=letters[1:6],'Mean Temperature Warmest Month')
dev.off()
system("open Figures/AP_all_response_curves_mt_warm_month.pdf")

#########################################################################
#== 7.Make Mean Kernels--------------------------------------------------
#== these loops can be used to calculate the mean predicitons across the region

#== set up splitting for parallelizing
ntasks=7
#== open or closed habitat
habitat=1 #1=open; 0=closed
#== management
manage='' # for bakcward compatibility
#manage=.02 # needed for south open and closed
#== file label
label=paste0('NE_habitat',habitat,manage)
version='v6'
#== specify the environmental conditions
env.subset=as.matrix(data.frame(unique.env, PAR=ifelse(habitat==1,1.5,-1),N=1.5,Ph.ave=1.5,light=ifelse(habitat==1,1.5,-1)))
cell.trick=env.cell.id
#== set up indices for cells being used
use=1:nrow(env.subset)
(ncells=length(use)) # note that this is 
#== set up parallel
task.split=chop.tasks(use,ntasks)          
#== number of IPM cells
n.matrix = 50
b=minSize+c(0:n.matrix)*(maxSize-minSize)/n.matrix 
#== mesh points (midpoints of the cells)
y=0.5*(b[1:n.matrix]+b[2:(n.matrix+1)])
#== width of the cells
h=y[2]-y[1]
 
#== specify formulas for use in vital rate functions
gr.form=as.formula(paste('~',as.character(gr.reg$Fixed$formula)[3]))
sv.form=as.formula(paste('~',paste0(dimnames(sv.reg)[[2]][-1],collapse='+')))
seed.form=as.formula(paste('~',paste0(dimnames(fec1.reg)[[2]][-1],collapse='+')))
germ.form=as.formula(paste('~',paste0(dimnames(fec2.reg)[[2]][-1],collapse='+')))
sum.sv.form=as.formula(paste('~',paste0(dimnames(fec3.reg)[[2]][-1],collapse='+')))
fl.form=as.formula(paste('~',paste0(dimnames(fl.reg)[[2]][-1],collapse='+')))
offspr.form=as.formula(paste('~',as.character(offspr.reg$Fixed$formula)[3]))

#========================================================================
#== SET MEAN PARAMETERS
gr.params.mean=apply(gr.reg$Sol,2,mean) # growth
gr.params.sd=mean(gr.reg$VCV[,'units']) # growth variance
sv.params=apply(sv.reg,2,mean)     
seed.params=apply(fec1.reg,2,mean)
germ.params=apply(fec2.reg,2,mean)
sum.sv.params=apply(fec3.reg,2,mean)
offspr.params.mean=apply(offspr.reg$Sol,2,mean) # recruit size
offspr.param.sd=mean(offspr.reg$VCV[,'units']) # recruit size variance
fl.params=apply(fl.reg,2,mean)

#========================================================================
#== BUILD IPMS
registerDoParallel(ntasks)
#== GROWTH KERNEL  
p.post.mean=foreach(i = 1:ntasks,.packages="Matrix") %dopar% {	
	p.mustard.kernel(task.split[[i]])
}
p.post.mean=unlist(p.post.mean,recursive=F)

#== FECUNDITY KERNEL 
f.post.mean=foreach(i = 1:ntasks,.packages="Matrix") %dopar% {	
	f.mustard.kernel(task.split[[i]],manage=manage)
}
f.post.mean=unlist(f.post.mean,recursive=F)

#== SIMULATE DYNAMICS 
sim=foreach(i = 1:ntasks,.packages=c("Matrix","IPMpack")) %dopar% { dynamics.mustard(task.split[[i]])
}
n.lam=unlist(lapply(sim,function(x) x[['lam']]),recursive=FALSE)

#== save for later use
save(n.lam,file=paste0('Model_Output/AP_mean_',label,'_pop_stats_',version,'.pred'))

#==========================================================================
#==========================================================================
#== 8. IPM diagnostics  ---------------------------------------------------

#== map lambda
pdf(paste0('Figures/AP_lambda_map_',label,'_',version,'.pdf'),h=7,w=7)
 pop.map( unlist(n.lam)[cell.trick], paste0('AP\n',ifelse(habitat==1,'open','closed'), ' habitat'),'lambda',max.val=30)
dev.off()
system(paste0('open Figures/AP_lambda_map_',label,'_',version,'.pdf'))

#== plot lambda vs env
pdf(paste0('Figures/AP_lambda_response_curves_',label,'_' ,version,'.pdf'),h=3,w=6) 
par(mfrow=c(1,2),mar=c(5,4,3,1))

plot(values(ne.env)[not.nas,'mt.warm.month'],unlist(n.lam)[cell.trick], pch=19,col= rgb(100,100,100,40,maxColorValue=255), las=1,bty='n',xlab='Mean Temp. Warmest Month',ylab='Lambda')
d.tmp=data.frame(lam=c(unlist(n.lam)[cell.trick]),values(ne.env)[not.nas, c('mt.warm.month','mp.may','n.droughts.gs')])
m1=gam(lam~s(mt.warm.month),data=d.tmp)
xx=seq(min(values(ne.env)[,'mt.warm.month'],na.rm=T),max(values(ne.env)[, 'mt.warm.month'],na.rm=T),by=.05)
lines(xx,predict(m1,data.frame(mt.warm.month=xx)),col='red1',lwd=3)

plot(values(ne.env)[not.nas,'mp.may'],unlist(n.lam)[cell.trick],pch=19, col=rgb(100,100,100,40,maxColorValue=255), las=1,bty='n',xlab='Mean May Precip.',ylab='Lambda')
m1=gam(lam~s(mp.may),data=d.tmp)
xx=seq(min(values(ne.env)[,'mp.may'],na.rm=T),max(values(ne.env)[, 'mp.may'],na.rm=T),by=.05)
lines(xx,predict(m1,data.frame(mp.may=xx)),col='red1',lwd=3)

dev.off()
system(paste0('open Figures/AP_lambda_response_curves_',label,'_',version,'.pdf'))



