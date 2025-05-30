% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiCure.R
\name{test.k.cure}
\alias{test.k.cure}
\title{Wald test across K groups for cure probabilities}
\usage{
test.k.cure(df, time, status, group)
}
\arguments{
\item{df}{A data frame in long format.}

\item{time}{Unquoted name of the time column.}

\item{status}{Unquoted name of the event indicator column (0/1).}

\item{group}{Unquoted name of the grouping factor (with >= 2 levels).}
}
\value{
A list with:
\describe{
\item{test_stat}{Wald statistic.}
\item{p_value}{p-value of the \\(\\chi^2_{K-1}\\) test.}
\item{df}{Degrees of freedom (\\(K-1\\)).}
\item{invertible}{Logical, whether contrast covariance was invertible.}
}
}
\description{
Tests equality of cure probabilities across \\(K\\) groups using a multivariate Wald statistic.
}
