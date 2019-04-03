\name{strength}
\alias{strength}
\title{Generate paper strength values}
\description{Generates paper strength values for an experiment with different
             temperatures.}
\usage{strength(nodays, noruns, temperature, ident)}
\arguments{
 \item{nodays}{The number of days over which the experiment is to be run.}
 \item{noruns}{The number of runs to be performed on each day of the experiment.}
 \item{temperature}{A \code{\link{factor}} that encapsulates the layout by giving the
                    temperature to be investigated for each run on each day.
                    These must be ordered so that the temperatures for the
                    first day are given in the order in which they are to be
                    investigated on that day. These must be followed by the
                    noruns temperatures for the second day and so on.
                    Consequently, the factor temperature will have nodays*noruns
                    values.}
 \item{ident}{The digits of your student identity number. That is, leave out
              any letters.}
}
\value{A \code{data.frame} object containing the \code{factors} day, run and 
       temperature and a \code{vector} of the generated strengths.}
\author{Chris Brien}
\examples{
## Here temperature is a factor with 4*3 = 12 values whose
## first 3 values specify the temperatures to be applied in
## the 3 runs on the first day, values 4 to 6 specify the
## temperatures for the 3 runs on day 2, and so on.
temperature <- factor(rep(c(80,85,90), 4))
exp.strength <- strength(nodays = 4, noruns = 3,
                         temperature = temperature, ident = 0123456)

## In this second example, a completely randomized design is generated 
## for the same 3 temperatures replicated 4 times. The layout is stored 
## in the data.frame called Design.
Design <- designRandomize(allocated = temperature, 
                          recipient = list(runs = 12), 
                          seed = 5847123)
## eradicate the unrandomized version of temperature
remove("temperature")

## The 12 temperatures in Design are to be regarded as being assigned to 
## days and runs in the same manner as for the first example.
exp.strength <- strength(nodays = 4, noruns = 3,
                         temperature = Design$temperature, ident = 0123456)
}
\keyword{datagen}
\keyword{design}
