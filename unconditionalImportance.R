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
data <- droplevels(data)

## Only analyses intercropping
data <- subset(data, CULT_ASOCIADO_PLATANO=="SI")
data$CULT_ASOCIADO_PLATANO <- NULL

## Separates in different groups
# management <- data[,c(1:3,26)]
# climate <- data[,c(4:20,26)]
# soil <- data[,c(21:25,26)]
# mgtsoil <- data[,c(1:3,21:25,26)]
# clisoil <- data[,c(4:25,26)]
# mgtcli <- data[,c(1:20,26)]

## Will store the variable importance of each model
v <- integer()

## Will store the profiles of each variable in each model
profiles <- list()
length(profiles) <- length(names(data)[-ncol(data)])
names(profiles) <- names(data)[-ncol(data)]

## Will store the performance of each model
perf <- NULL

nb.it <- 100

## Change the seed to have a different result
set.seed(0)
for(i in 1:nb.it) {
  print(i)
  inTrain <- createDataPartition(y=data$RDT_HA, p=0.7, list=F)
  training <- data[inTrain,]
  testing <- data[-inTrain,]
  
  grid <- expand.grid(mtry=round(seq(from=2, to=ncol(training)-1, length=4)))
  model <- train(training[,-ncol(training)], training[,ncol(training)],
                 method="cforest", tuneGrid=grid,
                 controls=cforest_unbiased(ntree=2000))
  
  performance <- R2(predict(model, testing), testing$RDT_HA) * 100
  
  if(length(perf)==0) {
    perf <- performance
  } else {
    perf <- c(perf, performance)
  }
  
  currentVarImp <- t(varImp(model, scale=F)$importance)
  
  ## Scales the performance for comparison
  ## Here we make it so that the sum of v.i. is equal to the perf.
  scale <- performance / sum(currentVarImp)
  scaledVarImp <- scale * currentVarImp
  
  if(length(v) == 0) {
    v <- scaledVarImp
  } else {
    v <- rbind(v, scaledVarImp)
  }
  
  ordered <- sort(sapply(as.data.frame(scaledVarImp), median), decreasing=T)
  
  ## Stores the profiles of the 5 most important variables
  for(n in names(ordered)[1:5]) {
    profile <- profilePlot(model, n, data, F)
    if(length(profiles[[n]]) == 0) {
      profiles[[n]] <- data.frame(profile$y)
      row.names(profiles[[n]]) <- profile$x
    } else {
      profiles[[n]] <- cbind(profiles[[n]], profile$y)
    }
  }
}

v <- as.data.frame(v)

ordered <- sort(sapply(v, median), decreasing=T)

meanperf <- signif(median(perf, 5))

## Computes the 95% confidence interval of the mean
se <- apply(v, 2, function(x) 1.96*sd(x, na.rm=TRUE)/sqrt(nrow(v)))
mean <- as.data.frame(ordered)
mean <- cbind(mean, names(ordered))
names(mean) <- c("Mean", "Variable")
mean$Variable <- factor(names(ordered), levels= names(ordered))
limits <- aes(ymax = Mean + se, ymin=Mean - se)
m <- ggplot(mean, aes(x=Variable, y=Mean))
m + geom_bar(stat="identity", width=0.5, fill="blue") +
  ylab("Mean importance")+ geom_errorbar(limits, width=0.25) +
  theme_bw() + ggtitle(paste("Importance of variables (with a mean R2 of",
                             meanperf, "%)")) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))

## Here you can use the multiprofile function
