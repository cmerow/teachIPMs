# page 8 
rm(list = ls())
require(IPMpack)
data(dataIPMpackSuccisa)
Sp<-dataIPMpackSuccisa
head(Sp)

# page 10
survModelComp(dataf = Sp, makePlot = TRUE, legendPos = "bottomright", mainTitle = "Survival") 
so<-makeSurvObj(dataf=Sp,Formula="surv~size")

# page 11
growthModelComp(dataf = Sp, makePlot = TRUE, legendPos = "bottomright", mainTitle = "Growth") 
go<-makeGrowthObj(dataf=Sp,Formula=sizeNext~size+size2)

# page 12
Pmatrix<-makeIPMPmatrix(survObj=so, growObj=go, 	minSize=min(Sp$size,na.rm=T), maxSize=max(Sp$size,na.rm=T), 	correction="constant")

# page 13
diagnosticsPmatrix(Pmatrix, growObj=go, survObj=so, correction="constant")

# page 15
Pmatrix<-makeIPMPmatrix(survObj=so, growObj=go, minSize=min(Sp$size,na.rm=T)-1, maxSize=max(Sp$size,na.rm=T)+2, correction="constant")

# page 16
require(fields)
image.plot(Pmatrix@meshpoints,
           Pmatrix@meshpoints,
           t(Pmatrix),
           main = "Pmatrix: survival and growth",
           xlab = "Size at t",
           ylab = "Size at t+1")
abline(0,1,lty=2,lwd=3)

# page 18
diagnosticsPmatrix(Pmatrix, growObj=go, survObj=so, correction="constant")

# page 19
fo1<-makeFecObj(Sp, Formula=fec1Bolt~1, Family = "binomial") # Intercept only model 
fo2<-makeFecObj(Sp, Formula=fec1Bolt~size, Family = "binomial") 
fo3<-makeFecObj(Sp, Formula=fec1Bolt~size+size2, Family = "binomial")

fs <- order(Sp$size)
fs.fec <- (Sp$fec1Bolt)[fs]
fs.size <- (Sp$size)[fs]
pfz <- tapply(fs.size, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE)
ps <- tapply(fs.fec, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE)
plot(as.numeric(pfz), as.numeric(ps), pch = 19, cex=2, col="blue", ylim=c(0,1), xlab="size", ylab="proportion flowering", main="")
x<-seq(from=0,to=10,length=1001)
x0<-data.frame(size=x,size2=x*x)
y0<-predict(fo1@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="red")
y0<-predict(fo2@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="green")
y0<-predict(fo3@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="blue")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"),c(AIC(fo1@fitFec[[1]]),AIC(fo2@fitFec[[1]]),AIC(fo3@fitFec[[1]]))), col = c(2:4),lty = 1, xjust = 1, bg = "white")

# page 21
fo <- makeFecObj(Sp, Formula=fec1Bolt~size+size2, Family = "binomial")

# page 22
fo <- makeFecObj(Sp, Formula=c(fec1Bolt~size+size2, fec2Stem~1), Family = c("binomial","poisson"), Transform = c("none","-1"))

# page 23
fo <- makeFecObj(Sp, Formula=c(fec1Bolt~size+size2,fec2Stem~1, fec3Head~size), Family = c("binomial","poisson","poisson"), Transform = c("none","-1","none"))

# page 24
fo <- makeFecObj(Sp, 
                 Formula=c(fec1Bolt~size+size2,fec2Stem~1,fec3Head~size), 
                 Family = c("binomial","poisson","poisson"), 
                 Transform = c("none","-1","none"), 
                 fecConstants = data.frame(seedsPerHead=50, 
                                           seedlingEstablishmentRate= 0.02))

# page 27
Fmatrix <- makeIPMFmatrix(fecObj=fo, 
                          minSize=min(Sp$size,na.rm=T)-1,
                          maxSize=max(Sp$size,na.rm=T)+2,
                          correction="constant")

# page 29
co <- makeClonalObj(Sp, Formula = cloning~1, Family = c("binomial"), Transform=c("none"))

# page 30
co <- makeClonalObj(Sp, Formula = c(cloning~1, clonesNext~1), Family = c("binomial","poisson"), Transform=c("none","-1"))

# page 31
co <- makeClonalObj(Sp, offspringSizeExplanatoryVariables = "size", Formula = c(cloning~1, clonesNext~1), Family = c("binomial","poisson"), 	Transform=c("none","-1"))

# page 32
Cmatrix <- makeIPMCmatrix(clonalObj=co, minSize=min(Sp$size,na.rm=T)-1,
                          maxSize=max(Sp$size,na.rm=T)+2, 
                          correction="constant")

# page 33
IPM <- Pmatrix + Fmatrix + Cmatrix
Re(eigen(IPM)$value[1])