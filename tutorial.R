# Sampling Tutorial

# Setup -----------------------------------------------------------------------
# Place the below script from the OSF page into your R library. These will
# create pre-made functions. Although these functions have been imported, you
# will still have to run some code (which is below) to execute them.
source("ResamplingPackage.R")

# Load your own data or use our tutorial data, also available on the OSF page
df <- read.csv("TutorialData.csv")

# Filter the data to select the data to resample
trait <- "dom"
df.to.sample <- filter(df, Trait == trait)

# Do the resampling -----------------------------------------------------------

# The RSSample function expects at least a list of values to resample as the
# first argument. If the values come from different stimuli, a second argument
# can be used to specify the simulus each value belongs to. It is important
# that, for each separate stimulus, the values are centered around zero. This
# is done by default and can be controlled with the center.data argument.

# The optional arguments to control the resampling are:
#
#   iterations:  Maximum number of values to sample. Sampling will begin with 1
#                value, and continue to increase by 1 until n. value. Default
#                100.
#   repetitions: The number of times to sample at each iteration. Default 500.
#   replacement: Whether to sample with replacement. Default TRUE.
#   confidence:  Confidence of confidence interval. Default 0.95.
#   ci.method:   Method to use to calculate the confidence interval. Default is
#                2 to use percentile. Enter 1 to use normal distribution.

# The function will return a list with the resampled data ($data) and the
# confidence intervals ($ci) for each iteration. Use these to plot the results
# with the RSFigure function.

trait.sampled <- RSSample(
    df.to.sample$Rating,
    df.to.sample$Stimulus,
    iterations = 100,
    repetitions = 500,
    replacement = TRUE,
    confidence = 0.95,
    ci.method = 2,
    center.data = TRUE
)

# Plot the result -------------------------------------------------------------

# Pass the result of the RSSample to RSFigure to plot the results. The first
# argument is the resampled data ($data), the second is the confidence
# intervals ($ci).

# The threshold argument can be used to plot a vertical line at the point of
# stability based on this threshold. The CI must go below this threshold before
# the data are considered stable. The absolute value of both the upper limit
# (UL) and the lower limit (LL) must be below the absolute value of the
# threshold. Use the inarow argument to specify how many iterations the CI must
# be below the threshold before the data are considered stable.

# You can use the title argument to specify the title of the plot, and the
# figPNGname argument to specify a filename save the plot as a PNG file.

RSFigure(
    trait.sampled$data,
    trait.sampled$ci,
    threshold = NA,
    title = "",
    figPNGname = "",
    inarow = 1
)

