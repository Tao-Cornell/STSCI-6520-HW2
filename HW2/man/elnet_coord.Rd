% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elnet_coord.R
\name{elnet_coord}
\alias{elnet_coord}
\title{elnet_coord Function}
\usage{
elnet_coord()
}
\value{
The function will return a matrix of regression coefficients for each \eqn{\lambda} collected as column vectors, the simulated design matrix \eqn{X}, the simulated response matrix \eqn{Y}.
 and together with a solution path plot shown in the plots panel.
}
\description{
This function fits elastic net to data using coordinate descent algorithm. The simulation setting is specified in homework 1 problem 3
part 2.The function will automatically generate the simulation setting specified in homework 1 problem 3 part 2 with \eqn{\lambda} running over a equally sapced grid of size \eqn{100} on \eqn{[0,3]}, so it does not take any arguments.
}
