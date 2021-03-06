% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/object.R
\name{plot.nimue_simulation}
\alias{plot.nimue_simulation}
\title{nimue simulation plot}
\usage{
\method{plot}{nimue_simulation}(
  x,
  var_select = NULL,
  replicates = FALSE,
  summarise = TRUE,
  ci = TRUE,
  q = c(0.025, 0.975),
  summary_f = mean,
  x_var = "t",
  date_0 = NULL,
  particle_fit = FALSE,
  ...
)
}
\arguments{
\item{x}{An squire_simulation object}

\item{var_select}{Vector of variable names to plot (default is all)}

\item{replicates}{Plot replicates}

\item{summarise}{Logical, add summary line}

\item{ci}{logical add confidence interval ribbon}

\item{q}{Quantiles for upper and lower of interval ribbon}

\item{summary_f}{Function to summarise each compartment
passed to the \code{fun} argument of \code{\link[ggplot2]{stat_summary}}}

\item{x_var}{X variable to use for plotting (default is \code{"t"},
but can be set to, \code{"date"}, if \code{date_0} provided), which will
cause the date to be plotted rather than time.}

\item{date_0}{Date of time 0 (e.g. "2020-03-01"), if specified a date column
will be added}

\item{particle_fit}{If the squire_simulation provided is the result of
running the particle filter, do we want to just plot the fit. Default =
FALSE}

\item{...}{additional arguments affecting the plot produced.}
}
\description{
nimue simulation plot
}
