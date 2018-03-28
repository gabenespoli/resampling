# Resampling Tool

Across diverse areas of research, it is common to average a series of observations, and to use these averages in subsequent analyses. Research using this approach faces the challenge of knowing when these averages are stable. Meaning, to what extent do these averages change when additional observations are included? Using averages that are not stable introduces a great deal of error into any analysis. The current research develops a tool, implemented in R, to assess when averages are stable. Using a sequential sampling approach, it determines how many observations are needed before additional observations would no longer meaningfully change an average.

## How it works

Given a list of values, an increasing number of items are randomly sampled for the given number of iterations. For example, on the first iteration 1 value is sampled, on the second iteration 2 values are sampled, and so on. The number of iterations can be controlled with the `iterations` input to the `RSSample` and `RSSampleLoop` functions.

At each iteration, the sampling is repeated and the different samplings are averaged together. The number of repetitions can be controlled with the `repetitions` input to the `RSSample` and `RSSampleLoop` functions.

Also at each iteration, a confidence interval (CI) is created around the repeated samplings. Two different CI methods are available.

Once the sampling has been completed, the CIs can be used to obtain the point of stability of the average using the `RSGetPointOfStability` and a specified threshold. The results can also be plotted using `RSFigure`.

See the comments within each function in `ResamplingTool.R` for a complete list of parameters and their default values.

## Tutorial

First, source the tool:

```R
source("ResamplingTool.R")
```

To sample from a vector, use the `RSSample` function:

```R
sampled <- RSSample(rating)
```

The output is a list of two data frames: 

`data` (`sampled$data` in the above example) has a column for each iteration containing the average of all of the samplings of values.

`ci` (`sampled$ci` in the above example) has three columns for the iteration, the lower limit of the CI, and the upper limit of the CI.

If you would like to do the sampling separately for different groups of values, you can pass the groupings as a factor. In our case, the data were ratings of different stimuli, and we wanted to do the sampling separately for each stimulus.

```R
sampled <- RSSample(rating, stimulus)
```

The tool also provides a loop wrapper for `RSSample`, creatively called `RSSampleLoop`. In our case, we had many different ratings for each stimulus, and wanted to sample separately for each rating-type and stimulus. Simply pass a factor as the second argument to `RSSampleLoop` to have the tool loop through the unique values in that factor. You can also provide a directory where the sampled data for each will be saved as an .Rda file.

```R
sampled <- RSSampleLoop(rating, trait, stimulus)
```

To find out how many iterations are needed before the average is stable, the `RSGetPointOfStability` function checks when the range of the CI falls below a specified threshold. Note that the units of the CI are the same as your data. In our case, these were ratings on a 7-point Likert scale (i.e., integers from 1 to 7).

```R
pos <- RSGetPointOfStability(sampled$ci, 1)
```

Use the `RSFigure` function to create a plot with iterations on the x-axis and the average value on the y-axis. Each stimulus as well as the CIs are plotted. YOu can specify a threshold (like in `RSGetPointOfStability`) so that a vertical line will be plotted at the point of stability. Save the figure as a .png by specifying a filename.

```R
RSFigure(sampled$data, sampled$ci, threshold=1, figPNGname="myFigure.png")
```

## Acknowledgements

This code was developed in collaboration with [Eric Hehman](http://erichehman.com), Sally Xie, and Eugene Ofosu.

This project can also be found on the [Open Science Framework](https://osf.io/82dsj/) and on the [Hehman Lab's Toolbox page](http://hehmanlab.org/toolbox).
