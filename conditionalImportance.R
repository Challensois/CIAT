library(caret)

source('~/RWorkspace/PB/profilePlot.R')

## Read the data, removes outliers and unwanted variables
data <- read.csv("original.csv",sep=";",na.strings="ND", dec=",",row.names=1)
data <- subset(data, DIBUJO_SIEMBRA_PLATANO!="EN_BARRERA")
data <- subset(data, d.interno!="LENTO A MUY LENTO")
data$ANALISIS_QUIMICO_PLATANO <- NULL
data$d.interno <- NULL
data$bio_3 <- NULL
data$bio_7 <- NULL

## Only analyses intercropping
data <- subset(data, CULT_ASOCIADO_PLATANO=="SI")
data$CULT_ASOCIADO_PLATANO <- NULL

data <- droplevels(data)

## Change the seed to get another split
set.seed(0)
## Splits the data in 70% training and 30% testing
inTrain <- createDataPartition(y=data$RDT_HA, p=0.7, list=F)
training <- data[inTrain,]
testing <- data[-inTrain,]

## Forces caret to test several mtries values
grid <- expand.grid(mtry=round(seq(from=2, to=ncol(training)-1, length=4)))

## Trains the first model with a first seed and 1000 trees
set.seed(1)
model1 <- train(training[,-ncol(training)], training[,ncol(training)],
                method="cforest", tuneGrid=grid,
                controls=cforest_unbiased(ntree=1000))

performance1 <- R2(predict(model1, testing), testing$RDT_HA) * 100

## Computes the conditional variable importance
currentVarImp1 <- varImp(model1, scale=F, conditional=T)

plot(currentVarImp1)

## Trains the second model with another seed and 1000 trees
set.seed(2)
model2 <- train(training[,-ncol(training)], training[,ncol(training)],
                method="cforest", tuneGrid=grid,
                controls=cforest_unbiased(ntree=1000))

performance2 <- R2(predict(model2, testing), testing$RDT_HA) * 100

currentVarImp2 <- varImp(model2, scale=F, conditional=T)

plot(currentVarImp2)

## Here, look at the most important variables and plot their profiles
