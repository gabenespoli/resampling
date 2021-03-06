# Resampling Package ----------------------------------------------------------
# This package contains functions for the resampling procedure from our paper.
# 
# Written by Gabriel A. Nespoli

# load required packages ------------------------------------------------------
library(Hmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(RColorBrewer)
library(psych)

# function RSSampleLoop -------------------------------------------------------
RSSampleLoop <- function(Rating, Trait,
                         Stimulus=factor(),
                         iterations=100,
                         repetitions=500,
                         replacement=TRUE,
                         confidence=0.95,
                         ci.method=2,
                         savedir="") {

  time.total <- proc.time() # start time of script


  if (length(Stimulus) != 0) {
      if (length(Stimulus) != length(Rating)) {
        stop("The length of Rating and Stimulus must be equal.")
      }
      if (length(Stimulus) != length(Trait)) {
        stop("The length of Trait and Stimulus must be equal.")
      }
  }
  if (length(Rating) != length(Trait)) {
    stop("The length of Rating and Trait must be equal.")
  }

  # create empty data frames to store data for output
  sampled <- data.frame()
  ci      <- data.frame()

  # loop through all Traits
  counter <- 0 # for displaying progress
  traits <- unique(Trait) # traits to loop through
  num.traits <- length(traits)
  for (trait in traits) {
    time.trait <- proc.time() # start time of script

    # display progress
    # e.g., "1/24 (competent): Sampling 100 value(s) 500 time(s)... "
    cat("\n")
    cat(paste(counter, "/", num.traits, " (", trait, "): ",
              "Sampling 1 to ", iterations, " value(s) ",
              repetitions, " time(s)... ",
              sep=""))

    cat(paste("\n", "Sampling 1 to ", iterations,
              " values from ", trait,
              sep=""))
    cat("\n")
    counter <- counter + 1

    # only consider current trait
    # restrict Rating and Stimulus to rows for the current trait
    Rating.trait   <- Rating  [Trait==trait]
    if (length(Stimulus) != 0) {
      Stimulus.trait <- Stimulus[Trait==trait]
    } else {
      Stimulus.trait <- factor()
    }

    # call RSSample to loop through Stimuli and do sampling
    temp <- RSSample(Rating.trait,
                        Stimulus.trait,
                        iterations=iterations,
                        repetitions=repetitions,
                        replacement=replacement,
                        confidence=confidence,
                        ci.method=ci.method)
    sampled.trait <- temp$data
    ci.trait <- temp$ci

    # prepend column with current trait
    sampled.trait <- cbind(Trait=rep(trait, nrow(sampled.trait)),
                           sampled.trait)
    ci.trait      <- cbind(Trait=rep(trait, nrow(ci.trait)),
                           ci.trait)

    # save this trait's sampling data to disk
    # also save a settings variable
    if (savedir != "") {
      settings <- list(iterations=iterations,
                       repetitions=repetitions,
                       replacement=replacement,
                       confidence=confidence,
                       ci.method=ci.method)
      fname <- file.path(savedir, paste("sampled_", trait, ".Rda", sep=""))
      cat("Saving '", fname, "'... ", sep="")
      save(sampled.trait, ci.trait, file=fname)
      cat("Done.\n")
    }

    # add to master data frame that includes all traits
    sampled <- rbind(sampled, sampled.trait)
    ci      <- rbind(ci,      ci.trait)

    # display trait processing time
    time.trait <- proc.time() - time.trait
    cat(paste("Total time for ", trait, ": ",
              round(time.trait["elapsed"] / 60, 1), " minutes.",
              sep=""))
    cat("\n")
  } # end looping traits

  # display total processing time
  time.total <- proc.time() - time.total
  cat(paste("Total time: ",
            round(time.total["elapsed"] / 60 / 60, 2), " hours.",
            sep=""))
  cat("\n")

  cat("Done sampling loop.")
  cat("\n")

  return(list(data=sampled, ci=ci, settings=settings))

}

# function RSSample -----------------------------------------------------------
RSSample <- function(Rating,
                     Stimulus=factor(),
                     iterations=100,
                     repetitions=500,
                     replacement=TRUE,
                     confidence=0.95,
                     ci.method=2) {
  # Sample data
  # Sample data. For each subset defined by unique combinations of the values
  #   of Stimulus, randomly sample an increasing number of values
  #   from Rating. This can be repeated many times at each iteration to make
  #   the sampling more robust.
  #
  # Args:
  #   Rating:       Vector of ratings to sample.
  #   Stimulus:     Factor with stimulus tokens to cycle through.
  #   iterations:   Maximum number of values to sample. Sampling will begin
  #                 with 1 value, and continue to increase by 1 until n.
  #                 value. Default 100.
  #   repetitions:  The number of times to sample at each iteration.
  #                 Default 500.
  #   replacement:  Whether to sample with replacement. Default TRUE.
  #   confidence:   Confidence of confidence interval. Default 0.95.
  #   ci.method:    1 = normal distribution, 2 = percentile
  #
  # Returns:
  #   A named list (sampled) with the following variables:
  #   data: [data frame] A data frame with columns Stimulus, and a column
  #         for each iteration (n1, n2, n3, etc.)
  #   ci:   [data frame] A data frame with columns, n, LL (lower limit of
  #         the confidence interval), and UL (upper limit).

  # check that inputs are the same length
  if ( (length(Stimulus) != 0) && (length(Stimulus) != length(Rating)) ) {
    stop("The length of Rating and Stimulus must be equal.")
  }

  # split Rating into a list of vectors by Stimulus so we can more
  #   do things separately for each stimulus
  if (length(Stimulus) != 0) {
    Rating.by.Stimulus <- split(Rating, Stimulus, drop=TRUE)
  } else {
    Rating.by.Stimulus <- list(data=Rating)
  }

  # center the ratings separately for each stimulus
  Rating.by.Stimulus <- lapply(Rating.by.Stimulus, RSMeanCenter)

  # loop iterations
  ci <- data.frame() # create container for ci data
  for (n in 1:iterations) {
    time.n <- proc.time() # start time for this iteration

    # display progress
    # e.g., "1/24 (competent): Sampling 100 value(s) 500 time(s)... "
    cat(paste("  Sampling ", n, " value(s) ", repetitions, " time(s)... ",
              sep=""))

    # do sampling for this n; get all repetitions
    # this returns a named list of stimuli, each item is `repetitions` long
    sampled.n.reps <- lapply(Rating.by.Stimulus, RSDoSampling,
                             n, replacement, repetitions)

    # calculate CI for this n and add to ci dataframe
    thisCI <- RSGetCI(unlist(sampled.n.reps), confidence, ci.method)
    ci[n, "n"]  <- n
    ci[n, "LL"] <- thisCI[1]
    ci[n, "UL"] <- thisCI[2]

    # convert to a data frame, add repetition column
    sampled.n <- stack(sampled.n.reps)
    sampled.n$rep <- 1:repetitions
    sampled.n <- sampled.n[, c(2, 3, 1)] # reorder cols
    names(sampled.n) <- c("Stimulus", "Repetition", as.character(n))

    # add this n to the other n's sampled data
    if (n > 1) {
      sampled <- merge(sampled, sampled.n)
    } else if (n == 1) {
      sampled <- sampled.n
    }

    # display processing time for this iteration
    time.n <- proc.time() - time.n
    cat(paste(round(time.n["elapsed"]), " seconds.", sep=""))
    cat("\n")
  } # end looping iterations

  cat("Done sampling.\n")
  return(list(data=sampled, ci=ci))
}

# function RSDoSamplng --------------------------------------------------------
RSDoSampling <- function(x, n, replacement=TRUE, repetitions=500) {
  # Repeatedly sample values from a vector.
  #
  # Args:
  #   x:            Vector to sample.
  #   n:            [numeric] Number of values to sample and average at each
  #                 repetition. Default 1.
  #   replacement:  [boolean] TRUE to sample with replacement; FALSE to sample
  #                 without replacement. Default TRUE.
  #   repetitions:  [numeric] The number of times to sample. Default 500.
  #
  # Returns:
  #   A vector the size of the number of repetitions.
  y <- rep(NA, repetitions) # container for repetitions
  for (i in 1:repetitions) {
    y[i] <- mean(sample(x, n, replace=replacement))
  }
  return(y)
}

# function RSMeanCenter -------------------------------------------------------
RSMeanCenter <- function(x) {
  # Center a vector on its mean by subtracting the mean from all values.
  x <- x - mean(x)
}

# function RSGetCI ------------------------------------------------------------
RSGetCI <- function(x, y=0.95, method=2) {
  # caller function for different confidence interval methods
  # Args:
  #   x:      Vector of values to get CI.
  #   y:      Numberic input to function being called. Default 0.95.
  #   method: Which function to call. See those functions for detail about
  #           their input. Use 1 for RSGetNormalRange or 2 for RSGetPercentile.
  #
  # Returns:
  #   CI <- c(LL, UL), where LL is the lower limit and UL is the upper limit.
  if (method == 1) {
    CI <- RSGetNormalRange(x, confidence=y)
  } else if (method == 2) {
    CI <- RSGetPercentile(x, percentile=y)
  } else {
    stop("Invalid method for getting CI.")
  }
  return(CI)
}

# function RSGetPercentile ----------------------------------------------------
RSGetPercentile <- function(x, percentile=0.95) {
  # Get the upper and lower percentile such that a 'percentile' amount of the
  #   values in x are contained within them.
  #
  # Args:
  #
  #
  # Returns:
  #   A 1-by-2 vector containing the lower and upper limit of the percentil CI.
  n     <- length(x)
  pct   <- (1 - percentile) / 2
  ind   <- floor(n * pct)
  x     <- sort(x)
  LL    <- x[ind]
  UL    <- x[n - ind]
  return(c(LL, UL))
}

# function RSGetNormalRange ---------------------------------------------------
RSGetNormalRange <- function(x, confidence=0.95) {
  # Get a confidence interval (CI) for a vector. Uses the following formula:
  #   mean +- ((z*sd)/sqrt(n)), where z is the critical value.
  #
  # Args:
  #   x:          Vector from which to calculate a CI.
  #   confidence: How much confidence the interval should have. Enter as a
  #               probability between 0 and 1. Default 0.95.
  #
  # Returns:
  #   A 1-by-2 vector containing the lower and upper limit of the CI.

  # split the confidence in half, one each for LL and UL
  # e.g., for 0.95 this will be 1.96; confidence = 0.975 for each of LL and UL
  z     <- qnorm(confidence + ((1 - confidence) / 2))
  n     <- length(x)
  meanx <- mean(x)
  sdx   <- sd(x)
  delta <- z * (sdx / sqrt(n))
  LL    <- meanx - delta # lower limit of CI
  UL    <- meanx + delta # upper limit of CI
  return(c(LL, UL))
}

# function RSGetPointOfStability ----------------------------------------------
RSGetPointOfStability <- function(ci, threshold, inarow=1) {
  # Get the point of stability.
  #
  # Args:
  #   ci:         Data frame of confidence intervals on which to base point of
  #               stability, as output from the RSSample function.
  #   threshold:  Threshold below which the CI must go before the data are
  #               considered stable. The absolute value of both the upper limit
  #               (UL) and the lower limit (LL) must be below the absolute value
  #               of the threshold.
  #   inarow:     Number of consecutive data points that must stay below the
  #               threshold to be considered stable.
  #
  # Returns:
  #   The n value at which the data are stable.
  threshold <- abs(threshold)
  ci.min <- apply(abs(cbind(ci$LL, ci$UL)), 1, min)
  stable <- ci.min < threshold
  stable.inarow <- logical()
  offset <- inarow - 1
  len <- length(stable) - offset
  for (i in 1:len) {
    j <- i + offset
    stable.inarow[i] <- all(stable[i:j])
  }
  pos <- min(which(stable.inarow==TRUE))
  return(pos)
}

# function RSFigure -----------------------------------------------------------
RSFigure <- function(sampled,
                     ci,
                     title="",
                     figPNGname="",
                     threshold=NA,
                     inarow=1) {
  # Create plots for resampling project.
  #
  # Args:
  #   sampled:    [data.frame] The sampled data from RSSample.
  #   ci:         [data.frame] The confidence intervals from RSSample.
  #   title:      [string] Title for the figure.
  #   figPNGname: [string] Filename to save the figure to disk.
  #   threshold:  [numeric] Plot vertical line at point of stability based on
  #               this threshold. See RSGetPointOfStability for more info.
  #   inarow:     [numeric] See RSGetPointOfStability for more info.
  #
  # Returns:
  #   Plot object.

  # reshape for plotting
  sampled <- reshape(sampled,
                     direction="long",
                     varying=list(3:ncol(sampled)))
  names(sampled) <- c("Stimulus", "Repetition", "n", "Rating", "id")

  # colours
  cols<-colorRampPalette(brewer.pal(9, "YlGnBu"))
  cols<-colorRampPalette(brewer.pal(11, "Spectral"))
  myPals <-cols(length(unique(sampled$Stimulus)))
  ci.colour <- "#000000"
  ci.linetype <- "solid"
  ci.size <- 1

  # plot each iteration and repetition
  p <- ggplot(sampled, aes(x=n, col=Stimulus)) +
    geom_line(aes(y=Rating), size=1, alpha=.2) +
    theme_minimal() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_line(colour="#e8e8e8"),
          panel.grid.minor = element_line(colour="#e8e8e8")) +
    scale_color_manual(values=myPals) +
    theme(legend.position="none") +
    xlab("n") +
    ylab("Rating") +
    ggtitle(title) +
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=24, face="bold"),
          plot.title=element_text(hjust=0.5, size=20)) +
    coord_cartesian(xlim=c(1, 100), ylim=c(-3.00, 3.00)) +
    scale_y_continuous(breaks=c(-3, -1, 0, 1, 3)) +
    scale_x_continuous(breaks=c(1, 20, 40, 60, 80, 100))

  # add lines for CI
  p <- p +
    geom_line(data=ci,
              aes(x=as.numeric(ci$n), y=ci$UL),
              colour=ci.colour,
              linetype=ci.linetype,
              size=ci.size) +
    geom_line(data=ci,
              aes(x=as.numeric(ci$n), y=ci$LL),
              colour=ci.colour,
              linetype=ci.linetype,
              size=ci.size)

  # add point of stability
  if (!is.na(threshold)) {
    pos <- RSGetPointOfStability(ci, threshold, inarow)
    p <- p +
      geom_vline(xintercept=pos,
                 colour=ci.colour,
                 linetype=ci.linetype,
                 size=ci.size)
  }

  # references for filesizes: Medium: h=22/3, w=30/3. Large: h=35/3, w=38/3
  if (figPNGname != ""){
    ggsave(plot = p, figPNGname, h = 27/3, w = 30/3, type = "cairo-png", dpi=72)
  }

  return(p)
}
