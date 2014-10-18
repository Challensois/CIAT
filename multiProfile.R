## Given a dataset and a list of profiles, plot profiles of a variable
multiProfile <- function(data, profiles, variable) {
  var <- data[,variable]
  if(is.factor(var)) {
    p <- as.data.frame(profiles[variable])
    ## Gets the 95% confidence interval of the mean
    se <- apply(p, 1, function(x) 1.96*sd(x, na.rm=TRUE)/nrow(p))
    mean <- rowMeans(p)
    ## Stores the mean effect of each category
    mean <- as.data.frame(mean - mean(mean))
    ## Stores the number of observations in each category
    mean <- cbind(mean, table(var))
    ## Stores their 95% c.i.
    mean <- cbind(mean, se)
    names(mean) <- c("Mean", "Class", "Freq", "SE")
    limits <- aes(ymax = Mean + SE, ymin=Mean - SE)
    m <- ggplot(mean, aes(x=Class, y=Mean, fill=Freq))
    m + geom_bar(stat="identity") + ylab("Mean effect on output")+
      coord_flip() + geom_errorbar(limits, width=0.25) +
      theme_bw() + ggtitle(paste("Individual influence of",
                                 variable,
                                 "(with", ncol(p), "profiles)")) +
      scale_fill_gradient2("Count", low = "red", high = "green",
                           midpoint=0)
  } else {
    p <- as.data.frame(profiles[variable])
    ## Stores all the points of interest for the variable in a column
    total <- data.frame(rep(as.numeric(row.names(p)), ncol(p)+1))
    ## Stores all the profiles and the mean profile in a column
    total <- cbind(total, c(unlist(p), rowMeans(p)))
    ## Stores the information of which run we are in
    total <- cbind(total, rep(1:(ncol(p)+1), each=nrow(p)))
    ## Stores if it is a "run" profile or the "mean" profile
    total <- cbind(total, c(rep("Runs", nrow(p)*ncol(p)),
                            rep("Mean",nrow(p))))
    names(total) <- c("Variable", "Outcome", "Run", "Group")
    
    ggplot() + geom_point(data=data,
                          aes_string(x=variable,
                                     y=names(data)[ncol(data)])) +
      geom_line(data=total,
                aes(x=Variable,
                    y=Outcome,
                    group=Run,
                    colour=Group,
                    size=Group,
                    alpha=Group),
      ) +
      theme_bw() + ggtitle(paste("Individual influence of",
                                 variable,
                                 "(with", ncol(p), "profiles)")) +
      ylab(names(data)[ncol(data)]) + xlab(variable) +
      scale_alpha_manual(values=c(1, 0.1)) +
      scale_size_manual(values = c(1, 1)) +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))
    }
}
