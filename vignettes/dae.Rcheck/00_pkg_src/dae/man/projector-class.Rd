\name{projector-class}
\Rdversion{1.1}
\docType{class}
\alias{projector-class}
\alias{coerce,projector,matrix-method}
\alias{coerce<-,projector,matrix-method}

\title{Class projector}
\description{The class "\code{\link[=projector-class]{projector}}" is the 
             subclass of matrices that are square, symmetric and idempotent.
                          
             \code{\link{is.projector}} is the membership function for this class.
             
             \code{\link{degfree}} is the extractor function for the degrees of freedom and 
             \code{\link[=degfree]{degfree<-}} is the replacement function.
             
             \code{\link{correct.degfree}} checks whether the stored degrees of freedom are correct.}
\section{Objects from the Class}{
An object of class "\code{\link[=projector-class]{projector}}" consists of a 
square, symmetric, idempotent matrix along with its degrees of freedom (rank).

Objects can be created by calls of the form \code{new("projector", data, nrow, ncol, byrow, dimnames, ...)}.
However, this does not add the degrees of freedom to the object. These can be 
added using the replacement function \code{\link[=degfree]{degfree<-}}.  
Alternatively, the function \code{\link{projector}} creates the new object 
from a \code{\link{matrix}}, adding its degrees of freedom at the same time.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}}
    \item{\code{degfree}:}{Object of class \code{"integer"}}
  }
}
\section{Extends}{
Class "\code{\link[=matrix-class]{matrix}}", from data part.
Class "\code{\link[=array-class]{array}}", by class "matrix", distance 2.
Class "\code{\link[=structure-class]{structure}}", by class "matrix", distance 3.
Class "\code{\link[=vector-class]{vector}}", by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "projector", to = "matrix")} }
    \item{print}{\code{signature(x = "projector")}}
    \item{show}{\code{signature(object = "projector")}}
	 }
}
\author{Chris Brien}

\seealso{\code{\link{projector}}, \code{\link{degfree}}, \code{\link{correct.degfree}} 
in package \pkg{dae}.
}
\examples{
showClass("projector")

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## check that it is a valid projector
is.projector(proj.m)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom and print the projector
degfree(proj.m) <- proj.m
}
\keyword{classes}
\keyword{array}
\keyword{projector}
