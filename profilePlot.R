profilePlot <- function(model, var, pred.data=model$trainingData, plot=T, ...) 
{
  xname <- ifelse(is.character(var), var,
                  ifelse(is.name(var), deparse(var), eval(var)))
  
  main = paste0("Partial Dependence on \"", xname, "\"")
  
  n.pt = min(length(unique(pred.data[, xname])), 50)
  xlab = xname
  ylab = names(pred.data)[ncol(pred.data)]
  
  xv <- pred.data[, xname]
  n <- nrow(pred.data)
  
  if (is.factor(xv) && !is.ordered(xv)) {
    x.pt <- levels(xv)
    y.pt <- numeric(length(x.pt))
    for (i in seq(along = x.pt)) {
      x.data <- pred.data
      x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
      y.pt[i] <- mean(predict(model, x.data), na.rm = TRUE)
    }
    
    if(plot) {
      library(fields)
      
      windows(width=16, height=9)
      m <- mean(y.pt)
      colors <- sapply(x.pt, function(x) sum(xv==x, na.rm=T))
      print(colors)
      palette <- two.colors(n=max(colors), start="lightgreen", end="darkgreen",
                            middle="green", alpha=1.0)
      colors <- palette[colors]
      barplot(y.pt - m, width = rep(1, length(y.pt)), col = colors,
              xlab = paste0("Mean difference on \"", ylab, "\""),
              ylab = xlab, main = main, names.arg = x.pt, horiz=TRUE, ...)
      image.plot(legend.only=TRUE, col=palette, zlim=c(0, length(xv)))
      abline(v = 0)
    }
  }
  else {
    if (is.ordered(xv)) 
      xv <- as.numeric(xv)
    x.pt <- seq(min(xv), max(xv), length = n.pt)
    if(is.integer(xv))
      x.pt <- as.integer(x.pt)
    y.pt <- numeric(length(x.pt))
    for (i in seq(along = x.pt)) {
      x.data <- pred.data
      x.data[, xname] <- rep(x.pt[i], n)
      predicted <- predict(model, x.data)
      y.pt[i] <- mean(predicted, na.rm = TRUE)
    }
    
    if(plot) {
      windows(width=16, height=9)
      par(mfrow=c(1,2), pty="s")

      plot(pred.data[, xname], pred.data[, ncol(pred.data)], pch=19,
           col=rgb(0,0,0,0.2), xlab=xlab, ylab=ylab, main=main)
      xl <- seq(min(x.pt),max(x.pt), (max(x.pt) - min(x.pt))/1000)
      lines(x.pt, y.pt, lwd=4, col="red")
      
      plot(x.pt, y.pt, type = "l", xlab = xlab, ylab = ylab,
           main = "Profile zoom", col="green", lwd=5)
      rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
    }
  }
  invisible(list(x = x.pt, y = y.pt))
}
