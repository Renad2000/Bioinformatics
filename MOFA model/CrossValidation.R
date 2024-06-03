install.packages("MOFA2")
install.packages("caret")
install.packages("lattice")
install.packages("ggplot2")
install.packages('caTools')
install.packages("reticulate")
install.packages("basilisk")
install.packages("ggplot2movies")

install.packages("dplyr",dependencies=TRUE)
library(ggplot2movies)
library(MOFA2)
library(dplyr)
library(MOFAdata)
library(data.table)
library(tidyverse)
library(magrittr)
library(e1071)
library(ggplot2)
library(lattice)
library(caret)
library(patchwork)
library(caTools)
library(reticulate)
library(basilisk)
library(stats)

install.packages("xfun", type = "binary")
set.seed(123)
devtools::install_github("bioFAM/MOFA", subdir="mofapy2")
devtools::install_github("bioFAM/MOFA", subdir="mofapy2", dependencies = TRUE)

reticulate::use_python("C:\\Users\\Renad\\AppData\\Local\\Programs\\Python\\Python312")
if (!require("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

install.packages("remotes")

remotes::install_github("bioFAM/MOFAdata")
BiocManager::install("MOFAdata")

#manipulate dataset 
utils::data("CLL_data")
Drugs<-write.csv(CLL_data$Drugs,"Drugs.csv")
Methylation<-write.csv(CLL_data$Methylation,"Methylation.csv")
mRNA<-write.csv(CLL_data$mRNA,"mrna.csv")


Drugs <- read.csv("Drugs.csv")
Methylation <- read.csv("Methylation.csv")
mRNA <-read.csv("mrna.csv")


row.names(Drugs)<-Drugs$X
Drugs<-Drugs[ ,-1]
row.names(Methylation)<-Methylation$X
Methylation<-Methylation[ ,-1]
row.names(mRNA)<-mRNA$X
mRNA<-mRNA[ ,-1]

Drugs<-as.matrix(Drugs)
Methylation<-as.matrix(Methylation)
mRNA<-as.matrix(mRNA)

pcdata<-Drugs
pcdata[is.na(pcdata)]<- 0 
gsub(NA,0, pcdata)
pcdata<-data.matrix(pcdata)
pcdata<-as.numeric(as.factor(pcdata))


kFold=4
pc1=prcomp(pcdata)
trainIndex=groupKFold(pc1$x[,1],k=kFold)
pcdata<-t(pcdata)

foldMem=rep(0,ncol(pcdata))
for(i in 1:kFold){
  foldMem=foldMem+i*!(1:ncol(pcdata) %in% trainIndex[[i]])
}

foldMem

#run mofa
modata<-list(Drugs,Methylation,mRNA)
x1=t(modata[[1]])
x2=t(modata[[2]])
x3=t(modata[[3]])


object <- create_mofa(modata)
object
data_opts <- get_default_data_options(object)
data_opts
model_opts <- get_default_model_options(object)
model_opts
train_opts <- get_default_training_options(object)
train_opts


object <- prepare_mofa(object,
                       data_options = data_opts,
                       model_options = model_opts,
                       training_options = train_opts
)

objectcv <- run_mofa(object, outfile="/Users/ricard/Downloads/objectcv.hdf5",use_basilisk =TRUE)
saveRDS(objectcv,"MOFACV.rds")
objectcv



w1<-get_weights(objectcv , views = "view_1", factors = 1)$view_1
w2<-get_weights(objectcv, views = "view_2", factors = 1)$view_2
w3<-get_weights(objectcv, views = "view_3", factors = 1)$view_3

fullContribution_1_MOFA=x1%*%w1
fullContribution_2_MOFA=x2%*%w2
fullContribution_3_MOFA=x3%*%w3

fullContribution_MOFA=list("Drugs"=t(fullContribution_1_MOFA),
                           "Methylation"=t(fullContribution_2_MOFA),
                           "mRNA"=t(fullContribution_3_MOFA)
                          )


#part3
cvContribution_1_MOFA_Mat=matrix(rep(NA,kFold*nrow(x1)),ncol=kFold,nrow=nrow(x1))
cvContribution_2_MOFA_Mat=matrix(rep(NA,kFold*nrow(x2)),ncol=kFold,nrow=nrow(x2))
cvContribution_3_MOFA_Mat=matrix(rep(NA,kFold*nrow(x3)),ncol=kFold,nrow=nrow(x3))


for(i in 1:kFold){
  id1=which(foldMem!=i)
  valid1_id1 <- id1[id1 <= nrow(x1)]
  valid2_id1 <- id1[id1 <= nrow(x2)]
  valid3_id1 <- id1[id1 <= nrow(x3)]
  train1=x1[valid1_id1,]
  train2=x2[valid2_id1,]
  train3=x3[valid3_id1,]
  
  
  
  
  moDataTemp=list(t(train1),t(train2),t(train3))
  
  
  mofaTemp=create_mofa(moDataTemp)
  data_opts <- get_default_data_options(mofaTemp)
  model_opts <- get_default_model_options(mofaTemp)
  train_opts <- get_default_training_options(mofaTemp)
  
  
  mofaTemp <- prepare_mofa(mofaTemp,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
  )
  
  
  
  mofaTemp <- run_mofa(mofaTemp, outfile="/Users/ricard/Downloads/mofaTemp.hdf5",use_basilisk =TRUE)
  saveRDS(mofaTemp,"mofaTemp.rds")
  mofaTemp
  
  

  weight1<-get_weights(mofaTemp,views = "view_1", factors = 1)$view_1
  weight2<-get_weights(mofaTemp, views = "view_2", factors = 1)$view_2
  weight3<-get_weights(mofaTemp, views = "view_3", factors = 1)$view_3
  
  
  
  
  fit1=x1[-id1,]%*%weight1
  fit2=x2[-id1,]%*%weight2
  fit3=x3[-id1,]%*%weight3
  
  
  cvContribution_1_MOFA_Mat[-id1,i]=fit1
  cvContribution_2_MOFA_Mat[-id1,i]=fit2
  cvContribution_3_MOFA_Mat[-id1,i]=fit3
  
  
}

cvContribution_Drugs_MOFA=rowSums(cvContribution_1_MOFA_Mat,na.rm=T)
cvContribution_Methylation_MOFA=rowSums(cvContribution_2_MOFA_Mat,na.rm=T)
cvContribution_mRNA_MOFA=rowSums(cvContribution_3_MOFA_Mat,na.rm=T)

cvContribution_MOFA=list("Drugs"=t(cvContribution_Drugs_MOFA),
                         "Methylation"=t(cvContribution_Methylation_MOFA),
                         "mRNA"=t(cvContribution_mRNA_MOFA))




#run movie
for(i in seq_along(moDataTemp)) {
  n_cols <- ncol(moDataTemp[[i]])
  colnames(moDataTemp[[i]]) <- paste0("Sample", 1:n_cols)
}

library(devtools)
devtools::install_github("mccabes292/movie")
library(movie)



movieOb=makeMovieObject(fullContributions =fullContribution_MOFA,
                        cvContributions = cvContribution_MOFA,foldMem = foldMem,scaleType="SD")




plot(movieOb,plotType = "SideBySide",xAxisPlot=1,yAxisPlot=2,colorVar = as.factor(foldMem),grid=FALSE,colorVarLabel="Fold")
plot(movieOb,plotType = "Comparison",grid=TRUE)
plot(movieOb,plotType = "Comparison",xAxisPlot = 1,yAxisPlot = 2,grid=FALSE,colorVar=as.factor(foldMem),colorVarLabel="Fold")
plot(movieOb,plotType = "Full",grid=TRUE)
plot(movieOb,plotType = "Full",xAxisPlot = 1,yAxisPlot = 2,grid=FALSE,colorVar=as.factor(foldMem),colorVarLabel="Fold")
plot(movieOb,plotType = "CV",grid=TRUE)
plot(movieOb,plotType = "CV",xAxisPlot = 1,yAxisPlot = 2,grid=FALSE,colorVar=as.factor(foldMem),colorVarLabel="Fold")

plot(movieOb,plotType="Line")
plot(movieOb,plotType="Line",leg=FALSE)



