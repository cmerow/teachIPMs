##############################################################################
#####             Building an IPM for a complex life cycle             #####
#####                       Hypericum cumulicola                         #####
#####          										 #####
##############################################################################

####################   Code developed by IPMpack team   ######################

#The data for this exercise contains a single annual transition (1997-1998) for the endemic, fire dependent, herbaceous perennial plant /Hypericum cumulicola/ in a single population in the Florida scrubs. The data were collected by Pedro Quintana-Ascencio (UCF) and Eric Menges (Archbold Biological Station) Please contact them for permission on usage for publication purposes.

#CREATED: July 20th, 2012
#LAST TIME MODIFIED: Aug 9th, 2015

#This routine is structured in three modules:

#     Module 0. Basic description of the life-cycle of Hypericum cumulicola
#     Module 1. Builds an IPM on a complex lifecycle, that of Hypericum cumulicola, a perennial herb with a seedbank for a single annual transition. The steps are to do a bit of re-organizing in the data, determine some important constant vital rates, do some basic model selection on survival, growth and reproduction, and finally put the IPM together!
#     Module 2. Obtain some basic demographic output from the IPM


#############################################################################
#############################################################################
#####                 Module 0: Hypericum's life-cycle                  #####
#############################################################################
#############################################################################

# Hypericum cumulicola (Clusiaceae) is an endangered, fire-dependent, perennial herb endemic to the Florida scrub ecosystem. Individuals of this species are iteroparous and are short-lived. Their seeds need not germinate immediately, having the possibility of forming a permanent soil seedbank. More information about the biology and demography of the species can be found in Quintana-Ascencio P.F., Menges, E.S. & Weekley C.W. (2003) A fire-explicit population viability analysis of /Hypericum cumulicola/ in Florida rosemary scrub. Conservation Biology 17 (2), 433-449.

#The life-cycle of /H. cumulicola/ is modelled here on an annual time-step, as described in Metcalf C.J.E. et al (2012) IPMPack: an R package for Integral Projection Models. Methods in Ecology and Evolution, DOI: 10.111/2041-210X.12001. Briefly, photosynthetically active individuals where sampled every May at Archbold Biological Station, Florida (USA). For each individual, the size (stem maximum height) was measured every year, as well as the number of fruits produced per capita. The number of seeds contained in each fruit, as well as the number of seeds that germinate the same year of their production, that go into / survive in / emerge from the soil seedbank where estimated via side experiments described in Quintana-Ascencio et al. (2003).


#############################################################################
#############################################################################
#####              Module 1: Building IPMs with seedbank                #####
#############################################################################
#############################################################################

#In this module you will learn to:
#   (i)   Upload the data of Hypericum cumulicola
#   (ii)  Carry out model selection for the vital rates of survival, growth and fecundity and constructing P and F matrices
#   (iii) Building and saving the IPM and its output: kernel and the parameters of the vital rates functions


  #############################################################################
  #####              (i) Dataset uploading and formating                  #####
  #############################################################################

  #Clean memory
    rm(list=ls(all=TRUE))

  #Load IPMpack (Make sure you are working with version 2.1!)
    library(IPMpack)

  #Information on the subset of information made available in IPMpack 2.1 can be accessed through the help manual for "hyperDataCovSubset"
    data(dataIPMpackHypericum)
    
  #Take a look at the help file of this data in IPMpack to understand how it is structured
    help(dataIPMpackHypericum)
    d1 <- dataIPMpackHypericum

  #Due to the sampling design described in the help manual, here we consider only individuals for which we are certain about their recruitment origin 
    d1 <- subset(d1,is.na(d1$size)==FALSE | d1$ontogenyNext==1)

  #Side experiments carried out by Quintana-Ascencio and Menges estimated the following vital rates for the studied annual transition:
    #Number of seeds produced per fruit:
      fec2 <- 13.78
    #Probability of seedling establishment
      fec3 <- 0.001336
    #Probability of seedling survival half a year after germinating, corresponding to the next annual census
      fec4 <- 0.14
    #Probability of a seed going into the seed bank
      goSB <- 0.08234528
    #Probability of a seed staying in the seed bank
      staySB <- 0.671

    #The following part does a simple re-organization of the data, getting rid of non-critical information for this exercise
      d1 <- d1[,c("surv","size","sizeNext","fec0","fec1")]

    #The following lines of code state the continuous (max height of individual plant) part of the IPM. Note that the IPM to be constructed here contains a discrete stage: seedbank.
      d1$stageNext <- d1$stage <- "continuous"
      d1$stage[is.na(d1$size)] <- NA

    #If an individual did not survive, it is labelled as dead to t+1.
      d1$stageNext[d1$surv==0] <- "dead"

    #The following lines of command compile a dataframe that contains the population dynamics of the seedbank. Namely, seeds in the seedbank can remain there from one year to the next (seedbank -> seedbank), be recruited as a photosynthetically active individual (seedbank -> continuous) or be recruited into the seedbank from seed (continuous -> seedbank). Of course, there is an implicitly modelled probability of mortality in the seedbank such that the survival + mortality = 1 in that stage. 
      d1$number <- 1
      sb1 <- data.frame(stage=c("seedbank","seedbank","continuous"),
                      stageNext=c("seedbank","continuous","seedbank"),
                      surv=1,size=NA,sizeNext=NA,fec0=NA,fec1=NA,number=c(staySB,(1-staySB)*fec3*fec4,1))
      sb1

    #The following lines simply put together the parts of the dataset for continuous and discrete parts of the lifecycle
      d1 <- rbind(d1,sb1)
      d1$stage <- as.factor(d1$stage)
      d1$stageNext <- as.factor(d1$stageNext)
    #Note that the variables stage and stageNext must be factors to work with IPMpack (common mistake to forget to specify these data as such)


  #############################################################################
  #####                 (ii)  Model comparison and selection              #####
  #############################################################################

  ###                       Comparing survival models                       ###

  #The following plot will compare a model for survival that is independent of size, scales linearly with size, and scales with size linearly and squared
      soComparison <- survModelComp(d1,
                      expVars = c(surv~1,
                                  surv~size,
                                  surv~size + size2,
                                  surv~size + size2 + size3),
                                  testType = "AIC",
                                  makePlot = T)

  #Exercise 1. Choose from the numbered lines of code below the appropriate model for survival
  #Choice 1:
    so <- makeSurvObj(d1, Formula = surv~1)
  #Choice 2:
    so <- makeSurvObj(d1, Formula = surv~size)
  #Choice 3:
    so <- makeSurvObj(d1, Formula = surv~size +size2)
  #Choice 4:
    so <- makeSurvObj(d1, Formula = surv~size +size2+size3)



  #Answer: The polynomial of third degree has the lowest AIC score (choice 4). However, all models have very similar AIC scores and **biological knowledge of the species** should play an important role in model selection. The little bump in survival for intermediate sizes has been observed in the field (P. Quintana-Ascencio, pers. comm.)
      picSurv(d1, so)

  ###                       Comparing growth models                         ###

  #The following plot will compare growth models from t to t+1, just as you did before for survival:
      goComparison <- growthModelComp(d1,
                                      expVars = c(sizeNext~1,
                                                  sizeNext~size,
                                                  sizeNext~size + size2,
                                                  sizeNext~size + size2 + size3),
                                      testType = "AIC",
                                      makePlot = T,
                                      legendPos = "bottomright")

  #Exercise 2. Choose from the numbered lines of code below the appropriate model for growth
  #Choice 1:
    go <- makeGrowthObj(d1, Formula = sizeNext~1)
  #Choice 2:
    go <- makeGrowthObj(d1, Formula = sizeNext~size)
  #Choice 3:
    go <- makeGrowthObj(d1, Formula = sizeNext~size +size2)
  #Choice 4:
    go <- makeGrowthObj(d1, Formula = sizeNext~size +size2+size3)




  #Answer: The linear model has the lowest AIC score (choice 2). However, models in choice 3 and 4 have very similar AIC scores and biological knowledge of the species should play an important role in model selection.

  #The following will save in your writing directory the chosen model in a pdf graph:
      picGrow(d1, go)
      #Note that the black line depicts the 1:1 size to sizeNext ratio. Individuals on this line did not change in size between the studied years, 1997 and 1998

  #Based on the graph produced above, it should be straightforward to see that small individuals are more likely to grow, whereas larger individuals are more likely to shrink. Carefully examine the linear regression crossing over the 1:1 size:sizeNext relationship line for ca. 45 mm high individuals.


  ###                           Fecundity models                            ###

  #IPMpack 2.1 does not yet have the model selection options for fecundity as you saw above in survival and growth (stay tuned for version 2.0 coming out in just a few days!). In the mean time, let's compare the models from scratch with the following few lines of code.

  #First we look at the model that best describes the probability of reproduction (a.k.a. fec0). 

#Probability of reproduction (fec0)

  #The following are some of the options for fec0 that will be compared visually and quantitatively below
    fec01 <- makeFecObj(d1, Formula = fec0~1, Family = "binomial")
    fec02 <- makeFecObj(d1, Formula = fec0~size, Family = "binomial")
    fec03 <- makeFecObj(d1, Formula = fec0~size+size2, Family = "binomial")
    fec04 <- makeFecObj(d1, Formula = fec0~size+size2+size3, Family = "binomial")

  #The following are a few lines re-organize the data for plotting below
    #Reorganizes individuals by size in time t
      fs <- order(d1$size)
    #Values of reproduction organized by size of individuals
      fsFec0 <- (d1$fec0)[fs]
      fsSize <- (d1$size)[fs]
    #Means at cut-points to be plotted for size and fec0
      pfz <- tapply(fsSize, as.numeric(cut(fsSize, 21)), mean, na.rm = TRUE)
      ps0 <- tapply(fsFec0, as.numeric(cut(fsSize, 21)), mean, na.rm = TRUE)

  #Make dummy size axis
    x <- seq(from = 0, to = 100, length = 1001)
    x0 <- data.frame(size = x, size2 = x^2, size3 =x^3)

  #This is the actual plot produced to compare the models for fec0 (the probability of being reproductive at time t)
    plot(as.numeric(pfz), as.numeric(ps0), pch = 19, cex = 1, col = "black", ylim = c(0, 1),
      xlab = "size", ylab = "Proportion of reproductive individuals", main = "Probability flowering")
      y0 <- predict(fec01@fitFec[[1]], newdata = x0); y0 <- exp(y0)/(exp(y0)+1); lines(x, y0, col=2)
      y0 <- predict(fec02@fitFec[[1]], newdata = x0); y0 <- exp(y0)/(exp(y0)+1); lines(x, y0, col=3)
      y0 <- predict(fec03@fitFec[[1]], newdata = x0); y0 <- exp(y0)/(exp(y0)+1); lines(x, y0, col=4)
      y0 <- predict(fec04@fitFec[[1]], newdata = x0); y0 <- exp(y0)/(exp(y0)+1); lines(x, y0, col=5)
      legend("bottomright", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2","size+size2+size3"),
        c("AIC"), c(AIC(fec01@fitFec[[1]]), AIC(fec02@fitFec[[1]]), AIC(fec03@fitFec[[1]]), AIC(fec04@fitFec[[1]]))),
        col = c(2:5), lty = 1, xjust = 1, bg = "white")

  #Exercise 3: Choose (i.e. execute only that line) from the lines of code below the appropriate model for the probability of reproduction, fec0
  #Choice 1:
    fec0ChosenModel <- fec0~1
  #Choice 2:
    fec0ChosenModel <- fec0~size
  #Choice 3:
    fec0ChosenModel <- fec0~size+size2
  #Choice 4:
    fec0ChosenModel <- fec0~size+size2+size3




  #Answer: The cubic model has the lowest AIC score (choice 4) and it accurately describes the biological process depicted in the plot.


#Number of fruits produced per capita (fec1)
  #Values of fruit production organized by size of individuals
  fsFec1 <- (d1$fec1)[fs]
  #Means at cut-points to be plotted for size and fec0
  ps1 <- tapply(fsFec1, as.numeric(cut(fsSize, 21)), mean, na.rm = TRUE)

  #The following are some of the options for fec1 that will be compared visually and quantitatively below
    fec11 <- makeFecObj(d1, Formula = fec1~1, Family = "poisson")
    fec12 <- makeFecObj(d1, Formula = fec1~size, Family = "poisson")
    fec13 <- makeFecObj(d1, Formula = fec1~size+size2, Family = "poisson")
    fec14 <- makeFecObj(d1, Formula = fec1~size+size2+size3, Family = "poisson")

  #This is the actual plot produced to compare the models for fec1
    plot(as.numeric(pfz), as.numeric(ps1), pch = 19, cex = 1, col = "black", ylim = c(0, 1000),
      xlab = "size", ylab = "Per-capita fruit production", main = "Number of fruits per capita")
      y1 <- predict(fec11@fitFec[[1]], newdata = x0); y1 <- exp(y1)+1; lines(x, y1, col=2)
      y1 <- predict(fec12@fitFec[[1]], newdata = x0); y1 <- exp(y1)+1; lines(x, y1, col=3)
      y1 <- predict(fec13@fitFec[[1]], newdata = x0); y1 <- exp(y1)+1; lines(x, y1, col=4)
      y1 <- predict(fec14@fitFec[[1]], newdata = x0); y1 <- exp(y1)+1; lines(x, y1, col=5)
      legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2","size+size2+size3"),
        c("AIC"), c(AIC(fec01@fitFec[[1]]), AIC(fec02@fitFec[[1]]), AIC(fec03@fitFec[[1]]), AIC(fec04@fitFec[[1]]))),
        col = c(2:5), lty = 1, xjust = 1, bg = "white")

  #Exercise 4: Choose (i.e. execute only that line) from the lines of code below the appropriate model for fruit production, fec1
  #Choice 1:
    fec1ChosenModel <- fec1~1
  #Choice 2:
    fec1ChosenModel <- fec1~size
  #Choice 3:
    fec1ChosenModel <- fec1~size+size2
  #Choice 4:
    fec1ChosenModel <- fec1~size+size2+size3




  #Answer: Again, the cubic model has the lowest AIC score (choice 4).

  #Fecundity is almost always the most complicated part of an IPM, and things can get even more complicated when there is a discrete stage like seedbank involved in this case. Type in "?makeFecObj" to get a feeling for all the possible options with sexual reproduction.
    fo <- makeFecObj(d1, Formula = c(fec0ChosenModel, fec1ChosenModel),
                    Family = c("binomial", "poisson"),
                    Transform = c("none", "none"),
                    meanOffspringSize = mean(d1[is.na(d1$size)==TRUE & is.na(d1$sizeNext)==FALSE, "sizeNext"]),
                    sdOffspringSize = sd(d1[is.na(d1$size)==TRUE & is.na(d1$sizeNext)==FALSE, "sizeNext"]),
                    fecConstants = data.frame(fec2=fec2,fec3=fec3,fec4=fec4),
                    offspringSplitter = data.frame(seedbank = goSB,continuous=(1-goSB)),
                    vitalRatesPerOffspringType = data.frame(seedbank = c(1, 1, 1, 0, 0),
                                                  continuous = rep(1,5),
                                                  row.names = c("fec0", "fec1", "fec2", "fec3", "fec4")))

###                 Creating the slots to hold the seedbank              ###

  #The following part creates a matrix that will accommodate the probabilities for entering, staying and leaving the seedbank (the discrete stage), as well as the continuous stage. Type "?makeDiscreteTrans" for help on this command
    dto <- makeDiscreteTrans(d1)
    dummy <- as.matrix(fo@offspringRel$coefficients[1])
    dimnames(dummy) <- list(1, "seedbank")
    dto@meanToCont <- as.matrix(dummy, dimnames = c(1, "seedbank"))
    dummy <- as.matrix(fo@sdOffspringSize)
    dimnames(dummy) <- list(1,"seedbank")
    dto@sdToCont <- as.matrix(dummy, dimnames = c(1, "seedbank"))


  ###                        Creating the P matrix                         ###

  #Creating the P matrix with the space for the seedbank
    Pmatrix <- makeIPMPmatrix(growObj = go,
                              survObj = so,
                              discreteTrans = dto,
                              minSize = 0,
                              maxSize = 30,
                              nBigMatrix = 100,
                              correction = "constant")

  #Creating a P matrix reflecting only the continuous part of the model and check that binning, etc is adequate...
    PmatrixContinuousOnly <- makeIPMPmatrix(growObj = go,
                                            survObj = so,
                                            minSize = 0,
                                            maxSize = 30,
                                            nBigMatrix = 100,
                                            correction = "constant")

  #...by means of looking at the diagnostics of the P matrix
    diagnosticsPmatrix(PmatrixContinuousOnly,
                        growObj = go,
                        survObj = so,
                        dff = d1,
                        correction = "constant")

  #Exercise 5: Why is your diagnostics of the P matrix looking so funky?



  #Answer: Because the max size was fit to 30 and in doing so the model is losing a lot of individuals above that size limit. Go back to surv, growth and fec objects and fit it to the correct size range as it follows:

    Pmatrix <- makeIPMPmatrix(growObj = go,
                              survObj = so,
                              discreteTrans = dto,
                              minSize = 0,
                              maxSize = 70,
                              nBigMatrix = 100,
                              correction = "constant")

    PmatrixContinuousOnly <- makeIPMPmatrix(growObj = go,
                                            survObj = so,
                                            minSize = 0,
                                            maxSize = 70,
                                            nBigMatrix = 100,
                                            correction = "constant")

  #Now the diagnostics should look a lot better:
    diagnosticsPmatrix(PmatrixContinuousOnly,
                        growObj = go,
                        survObj = so,
                        dff = d1,
                        correction = "constant")



  ###                        Creating the F matrix                         ###
    Fmatrix <- makeIPMFmatrix(fecObj= fo,
                              minSize = 0,
                              maxSize = 70,
                              nBigMatrix = 100,
                              correction="constant")

  #############################################################################
  #####               (iii) Creating IPM and saving output                #####
  #############################################################################
  #Forming the IPM as a result of adding the P and F matrices
    IPM <- Pmatrix + Fmatrix

  #Visualizing the IPM kernel via the function IPMpack built-in function "contourPlot", which creates contour plots:
    image(c(0, Pmatrix@meshpoints),c(0, Pmatrix@meshpoints),t(IPM),
          xlab="Size (t)",ylab="Size (t+1)",col=heat.colors(30), main="Kernel",zlim=c(0,0.03))
    contour(c(0, Pmatrix@meshpoints),c(0, Pmatrix@meshpoints),t(IPM), add = TRUE, drawlabels = TRUE)
  
  #Exercise 6: Where is the seedbank represented in this IPM kernel? Try to identify the parts of the kernel that pertain to entering, staying and leaving the seedbank.
  #Hint: log the data for visual representation of survival-dependent and reproduction-dependent processes




  #Answer: The very bottom row of the IPM shows a very thin, white line that corresponds to seeds entering the seedbank. The bottom-left most corner cell represents stasis in the seedbank. The first column on your left-hand of the IPM kernel represents seeds leaving the seedbank.



#############################################################################
#############################################################################
#####       Module 2: Get some basic output out of your IPM             #####
#############################################################################
#############################################################################

  #############################################################################
  #####                         (i) Eigenstructure                        #####
  #############################################################################

  ###                   Deterministic population growth rate               ###

  #For the full IPM, that is, including the seedbank
    lambdaModelWithSeedbank <- Re(eigen(IPM)$value[1])
  #This population is at stationary equilibrium since its population growth rate = 1.03 (i.e. very close to ~1)

  #Exercise 7: How does the population growth rate calculated in this IPM differ from the actual ratio of changes in individuals between the studied years 1997 and 1998 (a.k.a. instantaneous population growth rate)?



  #Answer:
    #The instantaneaous change in individuals from one year to the next is the ration of surviving individuals to the next year plus new recruits in the next year over the individuals alive the previous year:
      surviving <- length(which(!is.na(d1$sizeNext) & !is.na(d1$size))==T)
      recruits <- length(which(!is.na(d1$sizeNext) & is.na(d1$size))==T)
      existing <- length(which(!is.na(d1$size))==T)

    lambdaInstantaneous <- (surviving + recruits)/existing
      
    barplot(c(lambdaModelWithSeedbank,lambdaInstantaneous),names=c("Model with seedbank","Instantaneous"),ylab="Population growth rate",ylim=c(0,2))
      abline(a=1, b=0, lty=2, col="black")

    #The model seems to have over-estimated the population growth rate. However, because these are only point-estimates, one should run bootstraps on the data and see if the distributions around both estimates overlap or not to test for significant differences beween modeled and actuarial changes in population growth rate. Furthermore, the modeled population growth rate is that achieved on the long-run under the assumption that the population is left under the current environmental conditions, which is highly unlikely.
   

  #Exercise 8: What effect does the seedbank have on the population growth rate?



  #Answer
    lambdaModelWithoutSeedbank <- Re(eigen(IPM[2:100, 2:100])$value[1])
    barplot(c(lambdaModelWithSeedbank, lambdaModelWithoutSeedbank),names=c("Model with seedbank", "Model without seedbank"), ylab = "Population growth rate", ylim = c(0, 2))
      abline(a=1, b=0, lty=2, col="black")
  #The effect seems negative because the lambda derived from the IPM without the seedbank is lower than that with the seedbank. This could be due to the fact that the survival probability involved with going into, staying in and potentially emerging from the seedbank is very low. However, again, in order to test whether these lambdas are statistically significant, boostraps should be ran on the data.


#Parameter level sensitivities and elasticities of population growth rate

res1 <- sensParams(growObj = go, survObj = so, fecObj = fo,
  nBigMatrix = 100, minSize = 0, maxSize = 70,  
  discreteTrans = dto, delta = 1e-04, 
  response="lambda")

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res1$sens, 
main = expression("Parameter sensitivity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 
barplot(res1$elas, 
main = expression("Parameter elasticity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 


#Parameter level sensitivities and elasticities of net reproductive output

res2 <- sensParams(growObj = go, survObj = so, fecObj = fo,
  nBigMatrix = 100, minSize = 0, maxSize = 70,  
  discreteTrans = dto, delta = 1e-04, 
  response="R0")

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res2$sens, 
main = expression("Parameter sensitivity of population growth rate Ro"), 
las = 2, cex.names = 0.5) 
barplot(res2$elas, 
main = expression("Parameter elasticity of population growth rate Ro"), 
las = 2, cex.names = 0.5) 




