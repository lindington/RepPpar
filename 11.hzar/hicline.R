#######################################################
#### Hybrid index geographic cline model selection ####
#######################################################
##adapted from Vitali
#I removed the neutral diffusion part, tell me if needed

library(tidyverse)
library(ggplot2)
library(data.table)
library(hzar)
library(dplyr)

###Probability of ancestry (HI) along distance graph, K=2

setwd("C:/Users/Jag/Documents/LMU/work/01.capture/07.hzar/hicline")

#file with individual;locality;Distance (csv)
#Distances <- read.csv(file = "../Distances_From_Gabas.csv", header = T, sep = ";", dec = ",")
Distances <- read.csv(file = "Distances_From_Gabas.csv", header = T, sep = ";", dec = ",")

#file with individuals HI (txt)
#HI<- read.table(file="sample_HI_p.txt", header=T)
HI<- read.table(file="hi_sample.txt", header=T)

geo_hi <- data.frame(Distances[,2], stringsAsFactors = T)
geo_hi[,2:3] <- data.frame(HI[,1:2], stringsAsFactors = T)
geo_hi[,4] <- Distances[,3]

#Observed cline width and center estimation

geo_hi[,5] <- 20.2 # average # of individuals across both transects x2 = 25.789 #effective number of samples observed at each locality, from vitali?
remove(Distances,HI)

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=
  list(A=c(597,527,127,977,547,97),
       B=c(521,121,971,541,91,591),
       C=c(122,972,542,92,592,522),
       D=c(978,598,528,128,98,548),
       E=c(99,979,599,529,129,549),
       G=c(544,974,594,524,124,94)
       )


## We are doing just the one allele at one locus, but it is
## good to stay organized.
hi <- list();
## Space to hold the observed data
hi$obs <- list();
## Space to hold the models to fit
hi$models <- list();
## Space to hold the compiled fit requests
hi$fitRs <- list();
## Space to hold the output data chains
hi$runs <- list();
## Space to hold the analysed data
hi$analysis <- list();


hi$obs <- hzar.doMolecularData1DPops(geo_hi$V4,
                                 geo_hi$HI,
                                 geo_hi$V5)

## Look at a graph of the observed data
hzar.plot.obsData(hi$obs);

## Make models

## Make a helper function
loadhiAmodel <- function(scaling,tails,
                         id=paste(scaling,tails,sep="."))
  hi$models[[id]] <<- hzar.makeCline1DFreq(hi$obs, scaling, tails)

loadhiAmodel("free","none","model1");
loadhiAmodel("none" ,"none","model2");
loadhiAmodel("none" ,"right","model3");
loadhiAmodel("none" ,"left","model4");
loadhiAmodel("none" ,"mirror","model5");
loadhiAmodel("none" ,"both","model6");

## Check the default settings
print(hi$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between 0 and 570 km.
hi$models <- sapply(hi$models,
                    hzar.model.addBoxReq,
                    -5, 30000,
                    simplify=FALSE)


## Check the updated settings
print(hi$models)

## Compile each of the models to prepare for fitting
hi$fitRs$init <- sapply(hi$models,
                        hzar.first.fitRequest.old.ML,
                        obsData=hi$obs,
                        verbose=FALSE,
                        simplify=FALSE)
## Update the settings for the fitter if desired.
hi$fitRs$init$model1$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model1$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model1$mcmcParam$seed[[1]] <-
  mainSeed$A


## Check fit request settings
print(hi$fitRs$init)


## Run just one of the models for an initial chain
hi$runs$init <- list()
hi$runs$init$model1 <-
  hzar.doFit(hi$fitRs$init$model1)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model1))



## Run another model for an initial chain
hi$fitRs$init$model2$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model2$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model2$mcmcParam$seed[[1]] <-
  mainSeed$B 


hi$runs$init$model2 <-
  hzar.doFit(hi$fitRs$init$model2)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model2))


## Run another model for an initial chain
hi$fitRs$init$model3$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model3$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model3$mcmcParam$seed[[1]] <-
  mainSeed$C 


hi$runs$init$model3 <-
  hzar.doFit(hi$fitRs$init$model3)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model3))


## Run another model for an initial chain
hi$fitRs$init$model4$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model4$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model4$mcmcParam$seed[[1]] <-
  mainSeed$D 

hi$runs$init$model4 <-
  hzar.doFit(hi$fitRs$init$model4)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model4))


## Run another model for an initial chain
hi$fitRs$init$model5$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model5$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model5$mcmcParam$seed[[1]] <-
  mainSeed$E 


hi$runs$init$model5 <-
  hzar.doFit(hi$fitRs$init$model5)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model5))


## Run another model for an initial chain
hi$fitRs$init$model6$mcmcParam$chainLength <-
  chainLength;                          #1e5
hi$fitRs$init$model6$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
hi$fitRs$init$model6$mcmcParam$seed[[1]] <-
  mainSeed$G 

hi$runs$init$model6 <-
  hzar.doFit(hi$fitRs$init$model6)

## Plot the trace
plot(hzar.mcmc.bindLL(hi$runs$init$model6))


## Compile a new set of fit requests using the initial chains 
hi$fitRs$chains <-
  lapply(hi$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
hi$fitRs$chains <-
  hzar.multiFitRequest(hi$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit

##center for model1 - model6
randcenter<-runif(18,-5,30000)
for (i in 1:18){
  hi$fitRs$chains[[i]]$modelParam$init["center"]=randcenter[[i]]
}
##width for model1 - model6
randwidth<-runif(18,0,30005)
for (i in 1:18){
  hi$fitRs$chains[[i]]$modelParam$init["width"]=randwidth[[i]]
}
##pMin for model1
randpmin<-runif(3,0,1)
for (i in 1:3){
  hi$fitRs$chains[[i]]$modelParam$init["pMin"]=randpmin[[i]]
}
##pMax for model1
randpmax<-runif(3,0,1)
for (i in 1:3){
  hi$fitRs$chains[[i]]$modelParam$init["pMax"]=randpmax[[i]]
}
##deltaR for model3 and 6
randdeltaR<-runif(18,0,30005)
for (i in c(8:10,16:18)){
  hi$fitRs$chains[[i]]$modelParam$init["deltaR"]=randdeltaR[[i]]
}
##tauR for model3 and 6
randtauR<-runif(18,0,1)
for (i in c(8:10,16:18)){
  hi$fitRs$chains[[i]]$modelParam$init["tauR"]=randtauR[[i]]
}
##deltaL for model4 and 6
randdeltaL<-runif(18,0,30005)
for (i in c(10:12,16:18)){
  hi$fitRs$chains[[i]]$modelParam$init["deltaL"]=randdeltaL[[i]]
}
##tauL for model4 and 6
randtauL<-runif(18,0,1)
for (i in c(10:12,16:18)){
  hi$fitRs$chains[[i]]$modelParam$init["tauL"]=randtauL[[i]]
}
##deltaM for model5
randdeltaM<-runif(18,0,30005)
for (i in c(13:15)){
  hi$fitRs$chains[[i]]$modelParam$init["deltaM"]=randdeltaM[[i]]
}
##tauM for model5
randtauM<-runif(18,0,1)
for (i in c(13:15)){
  hi$fitRs$chains[[i]]$modelParam$init["tauM"]=randtauM[[i]]
}

remove(randcenter,randwidth,randdeltaL,randdeltaM,randdeltaR,randpmax,randpmin,randtauL,randtauM,randtauR)

## Go ahead and run a chain of 3 runs for every fit request (16:29 - max 16:50 )
hi$runs$chains <-  hzar.doChain.multi(hi$fitRs$chains,
                                      doPar=FALSE,
                                      inOrder=TRUE,
                                      count=3)

saveRDS(hi, file="hi_new_hzar_obj")
hi<-readRDS("hi_new_hzar_obj")

## Did model1 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Yes it did.

## Did model2 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Yes it did.

## Did model3 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Yes it did.

## Did model4 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## warning message:
#In hzar.doFit(hzar.request) : Fitting failed.

#Error in attr(data, "mcpar") <- c(start, end, thin) : 
#  attempt to set an attribute on NULL

## Did model5 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Yes it did.

## Did model6 converge?
summary(do.call(mcmc.list,
                lapply(hi$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Yes it did.

## All three models have three convergent chains, so additional runs
## are not needed.


## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele
## frequency independent of distance along cline) to include in
## analysis.
hi$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(hi$obs))

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
hi$analysis$initDGs$model1 <-
  hzar.dataGroup.add(hi$runs$init$model1)
hi$analysis$initDGs$model2 <-
  hzar.dataGroup.add(hi$runs$init$model2)
hi$analysis$initDGs$model3 <-
  hzar.dataGroup.add(hi$runs$init$model3)
hi$analysis$initDGs$model4 <-
  hzar.dataGroup.add(hi$runs$init$model4)
hi$analysis$initDGs$model5 <-
  hzar.dataGroup.add(hi$runs$init$model5)
hi$analysis$initDGs$model6 <-
  hzar.dataGroup.add(hi$runs$init$model6)

## Create a hiDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (nullModel, model1,
## model2, model3).
hi$analysis$oDG <-
  hzar.make.obsDataGroup(hi$analysis$initDGs)
hi$analysis$oDG <-
  hzar.copyModelLabels(hi$analysis$initDGs,
                       hi$analysis$oDG)

## Convert all 27 runs to hzar.dataGroup objects, adding them to
## the hiDataGroup object. 17:00 -
hi$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(hi$runs$chains[-10:-12],
                                hzar.dataGroup.add),
                         hi$analysis$oDG);

saveRDS(hi, file="hi_new_hzar_obj")
hi<-readRDS("hi_hzar_obj")

saveRDS(hi, file="hi_hzar_obj")
hi<-readRDS("../hicline/hi_hzar_obj")


## Check to make sure that there are only four hzar.dataGroup
## objects named nullModel, model1, model2, and model3 in the
## hiDataGroup object.
print(summary(hi$analysis$oDG$data.groups))

## Compare the 3 cline models to the null model graphically
hzar.plot.cline(hi$analysis$oDG);

## Do model selection based on the AICc scores
print(hi$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(hi$analysis$oDG));

## Print out the model with the minimum AICc score
print(hi$analysis$model.name <-
        rownames(hi$analysis$AICcTable
        )[[which.min(hi$analysis$AICcTable$AICc)]])

## Extract the hzar.dataGroup object for the selected model
hi$analysis$model.selected <-
  hi$analysis$oDG$data.groups[[hi$analysis$model.name]]


## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(hi$analysis$model.selected,
                         names(hi$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(hi$analysis$model.selected))

## Get summary stats
hi_all_sum<-summary(do.call(mcmc.list,
                            lapply(hi$runs$chains[16:18],
                                   function(x) hzar.mcmc.bindLL(x[[3]]) )) )

hi_center<-hi_all_sum[[1]][1,1]
hi_center_CI<-1.96*hi_all_sum[[1]][1,2]
hi_width<-hi_all_sum[[1]][2,1]
hi_width_CI<-1.96*hi_all_sum[[1]][2,2]
#hi_center = 16581[15921-17241] +/- (1.96*336.7755)
#width = 5460[3400-7520] +/- (1.96*1050.9296)
w1<-hi_center-(0.5*hi_width)
w2<-hi_center+(0.5*hi_width)
w3<-hi_center-(0.5*(hi_width-hi_width_CI))
w4<-hi_center-(0.5*(hi_width+hi_width_CI))
w5<-hi_center+(0.5*(hi_width-hi_width_CI))
w6<-hi_center+(0.5*(hi_width+hi_width_CI))


## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(hi$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hi_plot<-hzar.plot.fzCline(hi$analysis$model.selected,fzCol = "gray");

abline(v=c(hi_center,w1,w2), col=c("black","black", "black"), lty=c(1,2,2), lwd=c(1,1, 1))
rect(xleft = w3, xright = w4, ybottom =  par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("lightgray", alpha = 0.3))
rect(xleft = w5, xright = w6, ybottom =  par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("lightgray", alpha = 0.3))
abline(v = c(hi_center-7133.989/2,hi_center + 7133.989/2), col = c("deeppink","deeppink"))
abline(v = c(hi_center-9209.94/2,hi_center + 9209.94/2), col = c("blue","blue"))

  #15 000 generations since contact
  #  geom_vline(xintercept = 12.68288 - 9.20994/2, linetype = 5, colour = "blue", size = 0.7) +
  #  geom_vline(xintercept = 12.68288 + 9.20994/2, linetype = 5, colour = "blue", size = 0.7) +

## plot the other models for good measure:
hzar.plot.cline(hi$analysis$oDG$data.groups$nullModel)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$nullModel);

hzar.plot.cline(hi$analysis$oDG$data.groups$model1)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model1);

hzar.plot.cline(hi$analysis$oDG$data.groups$model2)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model2);

hzar.plot.cline(hi$analysis$oDG$data.groups$model3)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model3);

hzar.plot.cline(hi$analysis$oDG$data.groups$model4)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model4);

hzar.plot.cline(hi$analysis$oDG$data.groups$model5)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model5);

hzar.plot.cline(hi$analysis$oDG$data.groups$model6)
hzar.plot.fzCline(hi$analysis$oDG$data.groups$model6);
## End Analysis



dev.off()

################################################
### Neutral Diffusion Model & Ancestry cline ###
################################################

###Neutral diffusion model

#From Gay et al. 2008:

#If a barrier to gene flow is absent or weak and involves few
#genes, allelic clines at neutral loci are expected to decay after the
#secondary contact. In this simple scenario, the width of a neutral
#cline (w) depends only on the dispersal rate (sig) and the number
#of generations since contact (t) (Endler 1977)

#w = sqrt(2*pi)*sig*sqrt(t) ||/sqrt(2*pi)*sig

#Then:
#sqrt(t) = w/(sqrt(2*pi)*sig) ||^2
t = (w/(sqrt(2*pi)*sig))^2

w = 4.15329 #why?
sig = 0.030 #(Hewitt 1993)
t =

sqrt(2*pi)*sig*sqrt(9000)

#9000 gen
#7.133989

#15000 gen
#9.20994

(w/(sqrt(2*pi)*sig))^2

#Optional vertical lines to input into plot script to show NDM expected width for
#9000 generations since contact
  geom_vline(xintercept = 12.68288 - 7.133989/2, linetype = 5, colour = "deeppink", size = 0.7) +
  geom_vline(xintercept = 12.68288 + 7.133989/2, linetype = 5, colour = "deeppink", size = 0.7) +
#15 000 generations since contact
#  geom_vline(xintercept = 12.68288 - 9.20994/2, linetype = 5, colour = "blue", size = 0.7) +
#  geom_vline(xintercept = 12.68288 + 9.20994/2, linetype = 5, colour = "blue", size = 0.7) +