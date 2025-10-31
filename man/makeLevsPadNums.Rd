\name{makeLevsPadNums}
\alias{makeLevsPadNums}
\title{make a set of levels that consist of a prefix combined with padded numeric values}
\description{Makes a set of levels by combining a single \code{string} with \code{string}s 
             formed by either left- or right-padding of \code{numeric} values so that the 
             width of the padded \code{string} is at least \code{ndigits}. The \code{prefix} 
             and the padded string are combined using the function \code{paste0}. 
             If \code{pad.char} is \code{"0"},  \code{pad.direction} is \code{"right"}, 
             and there is a mixture of integers and reals, full stops will be added to the 
             integers as part of the padding process. Factor levels with a single-character 
             \code{prefix} and a sequence of padded integers from 1 to the number of levels are 
             often used as generic, abbreviated levels.
             
             A warning is generated if there is more than one \code{prefix} \code{string} or 
             non-numeric values are supplied to \code{nums}.}
\usage{makeLevsPadNums(prefix, nums, ndigits = 2, pad.char = "0", pad.direction = "left")}
\arguments{
 \item{prefix}{A \code{string} to use as the prefix to a set of padded numeric values.}
 \item{nums}{A \code{vector} of \code{numeric} values that are to be left- or right-padded 
               using \code{str_pad_left} or \code{str_pad_right} from the package 
               \code{stringi}. The width of the padded \code{string} will be at least 
               \code{ndigits}.}
 \item{ndigits}{An \code{integer} specifying the minimum width of the padded \code{string} 
                that is to be appended to the \code{prefix}.}
 \item{pad.char}{A \code{string} to use in padding the numeric values.}
 \item{pad.direction}{A \code{string} to specify whether to pad on the \code{left} or 
       the \code{right}.}
}
\value{A \code{\link{vector}} of the same length as \code{nums}, provided only a 
       single \code{prefix} is provided.}
\author{Chris Brien}
\seealso{\code{stringi::stri_pad_left}, \code{stringi::stri_pad_right}.}
\examples{
treat.levs <- makeLevsPadNums("T", 1:12)

treat.levs <- makeLevsPadNums("T", seq(0.0, 2.5, 0.5), ndigits = 3, pad.dir = "r")

levs <- makeLevsPadNums("", 1:12)
}
\keyword{datagen}
