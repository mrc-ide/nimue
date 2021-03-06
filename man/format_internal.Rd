% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{format_internal}
\alias{format_internal}
\title{Internals of Format vaccine model output as data.frame}
\usage{
format_internal(x, compartments, summaries, reduce_age, index, time, replicate)
}
\arguments{
\item{x}{squire_simulation object}

\item{compartments}{Vector of compartment names, e.g. \code{c("S", "R")}, or sub-compartment names, e.g. \code{c("S", "E1", "E2")}}

\item{summaries}{Vector of summary names, which may be:
\itemize{
      \item{"deaths"}{ Deaths per day }
      \item{"infections"}{ Infections per day. New infections (note this is currently a slightly different definitionto the main Squire mode)}
      \item{"hospitilisations"}{ Hospitalisations per day (Note this takes into account hospital capacity)}
      \item{"hospital_occupancy"}{ Occupied Hospital Beds }
      \item{"ICU_occupancy"}{ Occupied ICU Beds }
      \item{"hospital_demand}{ Required Hospital Beds }
      \item{"ICU_demand}{ Required ICU Beds }
      \item{"vaccinated"}{ Vaccines administered per day}
      }}

\item{reduce_age}{Collapse age-dimension, calculating the total in the
compartment.}

\item{index}{odin ouput index}

\item{time}{time vector}

\item{replicate}{outpu replicate number}
}
\description{
Internals of Format vaccine model output as data.frame
}
