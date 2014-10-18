viewData <- function(data, dummify=T, predict=T, clustering=F, baseplot=T, old=F) {
  ## Libraries
  library(kohonen)
  library(caret)
  library(RColorBrewer)
  library(deldir)
  library(fields)
  
  computeUMatrix <- function(som, predict) {
    ## Handles the values from outside the matrix by returning NA
    getUMatrixValue <- function(umatrix, row, col) {
      if(row > 0 && row <= nrow(umatrix) && col > 0 && col <= nrow(umatrix)) {
        neighbor <- umatrix[row,col]
      } else {
        neighbor <- NA
      }
      neighbor
    }
    
    ## Only constructs the U-matrix on X-space
    if(predict) {
      dists <- som$codes$X
    } else {
      dists <- som$codes
    }
    
    ## Creates a Delaunay triangulation and Dirichlet tessellation
    d <- deldir(x=som$grid$pts[,1], y=som$grid$pts[,2])
    
    ## Removes the border effects
    corrected <- d$delsgs[order(d$delsgs$y1, d$delsgs$y2, d$delsgs$x1, d$delsgs$x2),5:6]
    trues <- abs(corrected$ind1 - corrected$ind2) <= som$grid$xdim + 1
    corrected <- corrected[trues,]
    rownames(corrected) <- 1:nrow(corrected)
    
    ## Creates the U-matrix and computes the "non-neurons" units
    u.matrix <- matrix(nrow=2*som$grid$ydim-1, ncol=2*som$grid$xdim-1)
    index <- 1
    for (row in 1:nrow(u.matrix)) {
      for (col in 1:ncol(u.matrix)) {
        if(row%%2 == 0 || col%%2 == 0) { ## This is a non-neuron unit
          u.matrix[col, row] <- dist(rbind(dists[corrected[index,1],],
                                           dists[corrected[index,2],]))
          index <- index+1
        } else { ## This is a neuron unit
          u.matrix[row, col] <- 0
        }
      }
    }
    ## The values in the neurons units are the mean of the
    ## surrounding (non-neurons) cells
    for (row in 1:nrow(u.matrix)) {
      for (col in 1:ncol(u.matrix)) {
        if(row%%2 == 0 || col%%2 == 0) { ## This is a non-neuron unit
          ## Do nothing
        } else { ## This is a neuron unit
          neighbors <- numeric(6)
          if(((row+1)/2)%%2 == 0) { ## Left-shifted row, neighbor cols are col and col-1
            neighbors.cols <- c(col-1, col, col-1, col+1, col-1, col)
          } else { ## Right-shifted row, neighbor cols are col and col+1
            neighbors.cols <- c(col, col+1, col-1, col+1, col, col+1)
          }
          ## Neighbors rows are row-1, row and row+1
          for(neighbor in 1:6) {
            neighbor.row <- row - (ceiling(neighbor/2) - 2)
            neighbors[neighbor] <- getUMatrixValue(u.matrix, neighbor.row, neighbors.cols[neighbor])
          }
          u.matrix[row, col] <- mean(neighbors, na.rm=T)
        }
      }
    }
    u.matrix
  }
  
  plotSOM <- function(som, title, matrix, type, palette) {
    if(baseplot) {
      plotSOMbase(som, title, matrix, type, palette)
      return(0)
    } else {
      return(plotSOMGG(som, title, matrix, type))
    }
  }
  
  plotSOMbase <- function(som, title, matrix, type, palette) {
    ## Creates empty plot with title and good dimensions
    dummyPlot(som, title)
    
    ## Sets the correct color options
    ## Either equal bins
    if(type == "Interval") {
      if(missing(palette)) {
        color.ramp <- designer.colors(n=50, c("grey75", "black"))
      } else {
        color.ramp <- rev(designer.colors(50, brewer.pal(9, palette)))
      }
      bins <- seq(min(matrix), max(matrix), length=length(color.ramp))
    }
    ## Or quantile bins
    if(type == "Quantile") {
      if(missing(palette)) {
        palette <- "Spectral"
      }
      color.ramp <- rev(designer.colors(50, brewer.pal(9, palette)))
      if(round(max(matrix)) == 1 & round(min(matrix)) == 0) {
        bins <- seq(min(matrix), max(matrix), length=length(color.ramp))
      } else {
        bins <- quantile(x=matrix, probs=cumsum(rep(1/length(color.ramp), length(color.ramp))), na.rm = T)
      }
    }
    color.codes <- rep("#FFFFFF", length(matrix))
    for (i in 1:length(matrix))
      if (!is.na(matrix[i])) color.codes[i] <- color.ramp[which.min(abs(bins-matrix[i]))] 
    
    ## Plots the hexagonal som and the legend
    plotHexagons(som, color.codes)
    image.plot(legend.only=TRUE, col=color.ramp, zlim=c(min(matrix), max(matrix)))
  }
  
  dummyPlot <- function(grid, title) {
    plot(0, 0, xlab="", ylab="", type="n", axes=F, main=title, asp=sqrt(3)/2,
         xlim=c(min(grid$pts[,1])-1/2,
                max(grid$pts[,1])+1/2),
         ylim=c(min(grid$pts[,2])-(sqrt(3)/3),
                max(grid$pts[,2])+(sqrt(3)/3))
    )
  }
  
  plotHexagons <- function(grid, color.codes) {
    Hexagon <- function (x, y, unitcell = 1, col = "grey", border=NA) {
      polygon(c(x-unitcell/2, x-unitcell/2, x, x+unitcell/2, x+unitcell/2, x),
              c(y-(sqrt(3)/6), y+(sqrt(3)/6), y+(sqrt(3)/3), y+(sqrt(3)/6),
                y-(sqrt(3)/6), y-(sqrt(3)/3)), 
              col = col, border=border)
    }
    
    for(i in 1:nrow(grid$pts)) {
      Hexagon(grid$pts[i,1], grid$pts[i,2], col = color.codes[i])
    }
  }
  
  plotSOMGG <- function(grid, title, matrix, type) {
    if(type == "Interval") {
      bins <- seq(min(matrix), max(matrix), length=50)
    } else {
      palette="Spectral"
      color.ramp <- rev(designer.colors(50, brewer.pal(9, palette)))
      if(length(unique(round(matrix)))==2) {
        bins <- seq(min(matrix), max(matrix), length=50)
      } else {
        bins <- quantile(x=matrix, probs=cumsum(rep(1/length(color.ramp), length(color.ramp))), na.rm = T)
      }
    }
    
    fill <- rep(0, length(matrix))
    for (i in 1:length(matrix))
      if (!is.na(matrix[i])) fill[i] <- which.min(abs(bins-matrix[i]))
    
    x <- as.vector(sapply(grid$pts[,1], function(x) c(x-1/2, x-1/2, x, x+1/2, x+1/2, x)))
    y <- as.vector(sapply(grid$pts[,2], function(y) c(y-(sqrt(3)/6), y+(sqrt(3)/6),
                                                      y+(sqrt(3)/3), y+(sqrt(3)/6),
                                                      y-(sqrt(3)/6), y-(sqrt(3)/3))))
    id <- rep(1:nrow(grid$pts), each=6)
    polygons <- data.frame(id=id, x=x, y=y, value=rep(matrix,each=6), fill=rep(fill,each=6))
    
    myplot <- ggplot(polygons, aes(x=x, y=y)) +
      geom_polygon(aes(fill=value, color=NULL, group=id)) +
      ggtitle(title) + coord_fixed(ratio=sqrt(3)/2) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),
            panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            plot.background=element_blank())
    if(type=="Quantile") {
      myplot <- myplot + scale_colour_manual(breaks = c(1:50),
                                             values = color.ramp)
    } else {
      myplot <- myplot + scale_colour_grey(start=0.2)
    }
    
    return(myplot)
  }
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  if(class(data)!="data.frame") {
    if(class(data)=="matrix") {
      data <- as.data.frame(data)
    } else {
      stop("Use viewData with a named matrix or data.frame")
    }
  }
  if(length(names(data)) != ncol(data)) {
    stop("Use viewData with a named matrix or data.frame")
  }
  
  ## Transforms factors in binaries
  if(dummify) {
    try(dummies <- dummyVars(~., data, levelsOnly = T), silent = T)
    if(!length(dummies) > 0) {
      dummies <- dummyVars(~., data)
    }
    data <- as.data.frame(predict(dummies, data))
  }
  
  ## Preprocess data accordingly
  data <- data[,sapply(data, is.numeric)]
  ## Kohonen package does NOT accept NA's
  data <- data[complete.cases(data),]
  ## Remove constant variables
  data <- data[,apply(data, 2, function(x) if(sd(x)!=0) T else F)]
  
  if(predict) {
    ## Don't separate into testing and training from now on
    if(old) {
      ## Separate in training and test sets
      inTrain <- createDataPartition(y=data[,ncol(data)], p=0.7, list=F)
      training <- data[inTrain,]
      testing <- data[-inTrain,]
      
      ## Scales the testing and the training set based only on the training set
      scaling <- preProcess(x=training)
      data.scale <- predict(scaling, training)
      testing <- predict(scaling, testing)
      
      #Build supervised SOM
      sqrt <- floor(sqrt(nrow(training)))
    }
    
    ## Scales the data
    scaling <- preProcess(x=data)
    data.scale <- predict(scaling, data)
    sqrt <- floor(sqrt(nrow(data)))
    
    aGrid <- somgrid(xdim = sqrt, ydim = sqrt, topo="hexagonal")
    som <- xyf(as.matrix(data.scale[,-c(ncol(data.scale))]),
               as.matrix(data.scale[,ncol(data.scale)]),
               grid=aGrid)
    
    if(old) {
      ## Computes the performance of the SOM
      predictions <- predict(som, as.matrix(testing))$prediction
      unscaled.output <- ((testing[,ncol(training)] * scaling$std[ncol(training)]) +
                            scaling$mean[ncol(training)])
      r2.te <- R2(predictions, testing[,ncol(training)])
      
      ## Plots the results of the SOM predictions on the test set
        X11()
        plot(unscaled.output, type="l", col="green", lwd=5,
             main=paste("Test predictions with SOM on", names(data)[ncol(data)],
                        "is", round(100*r2.te),
                        "%"), ylab=names(training)[ncol(training)])
        unscaled.predictions <- ((predictions * scaling$std[ncol(training)]) +
                                   scaling$mean[ncol(training)])
        points(unscaled.predictions, col="red", pch=19)
        legend("bottomright", c("Real","Predicted"), col=c("green","red"),
               lwd=5, lty=c(1,NA), pch=c(NA,19), cex=0.8)
    }
  } else {
    ## Scales all the data
    scaling <- preProcess(x=data)
    data.scale <- predict(scaling, data)
    
    #Build unsupervised SOM
    sqrt <- floor(sqrt(nrow(data)))
    aGrid <- somgrid(xdim = sqrt, ydim = sqrt, topo="hexagonal")
    som <- som(as.matrix(data.scale), grid=aGrid)
  }
  
  ## Creates a plotting region with correct dimension
  windows(width=16, height=9)
  sup.plots <- ifelse(clustering, 2, 1)
  nb.plots <- ncol(data)+sup.plots
  
  x = array(list(), nb.plots)
  
  nb.plots.x <- ceiling(sqrt(nb.plots))
  nb.plots.y <- ifelse(nb.plots.x^2-nb.plots.x>=nb.plots,nb.plots.x-1,nb.plots.x)
  par(mar=c(0,0,1,0), mfrow=c(nb.plots.y, nb.plots.x), pty = "m", oma=c(1,1,1,1))
  
  ## Computes U-Matrix before unscaling
  u.matrix <- computeUMatrix(som, predict)
  
  ## Creates a custom grid relevant for visualizing the U-Matrix
  u.grid <- somgrid(xdim=ncol(u.matrix), ydim=nrow(u.matrix), topo="hexagonal")
  for(i in 1:sqrt(length(u.matrix))) {
    for(j in 1:sqrt(length(u.matrix))) {
      if(((i+1)/2)%%2 == 0) {
        u.grid$pts[(i-1)*nrow(u.matrix)+j,1] <- u.grid$pts[(i-1)*nrow(u.matrix)+j,1]-1
      }
    }
  }
  x[[1]] <- plotSOM(u.grid, title="U-Matrix", matrix=u.matrix, type="Interval")
  
  ## Unscale SOM
  if(predict) {
    old.codes <- som$codes$X
    data.unscale <- cbind(som$codes$X, som$codes$Y)
    colnames(data.unscale) <- names(data)
  } else {
    old.codes <- som$codes
    data.unscale <- som$codes
  }
  for(i in 1:ncol(data.unscale)) {
    data.unscale[,i] <- ((data.unscale[,i] *
                            scaling$std[i]) +
                           scaling$mean[i])
  }
  som$codes <- data.unscale
  
  i = 2
  ## Plots the variables component planes
  for(p in colnames(som$codes)) {
    intensities <- som$codes[,p]
    x[[i]] <- plotSOM(som$grid, title=p, matrix=intensities, type="Quantile", palette="Spectral")
    i <- i+1
  }
  
  if(clustering) {
    ## Plots clusters according to pamk
    library(fpc)
    pamk <- pamk(dist(old.codes), krange = 2:10)
    color.ramp <- rainbow(pamk$nc)
    color.codes <- color.ramp[pamk$pamobject$clustering]
    dummyPlot(som$grid, "Clusters")
    plotHexagons(som$grid, color.codes)  
  }
  
  if(!baseplot) {
    multiplot(plotlist = x, cols = ceiling(sqrt(nb.plots)))
  }
  invisible(som)
}
