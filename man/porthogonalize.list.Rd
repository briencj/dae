\name{porthogonalize.list}
\alias{porthogonalize.list}
\alias{porthogonalize}
\title{Takes a list of \code{\link{projector}s} and constructs a \code{\link{pstructure.object}} that includes projectors, each of which has been orthogonalized to all projectors preceding it in the list.}
\description{Constructs a \code{\link{pstructure.object}} that includes a 
             set of mutually orthogonal projectors, one for each of the 
             \code{\link{projector}s} in the \code{\link{list}}. 
             These specify a structure, or an orthogonal decomposition of the 
             data space. This function externalizes the process previously performed 
             within \code{\link{pstructure.formula}} to orthogonalize 
             \code{\link{projector}}s. There are three methods 
             available for carrying out orthogonalization: \code{differencing}, 
             \code{eigenmethods} or the default \code{hybrid} method. 

             It is possible to use this function to find out what sources 
             are associated with the terms in a model and to determine the 
             marginality between terms in the model. The marginality matrix 
             can be saved.
}
\usage{
\method{porthogonalize}{list}(projectors, formula = NULL, keep.order = TRUE, 
               grandMean = FALSE, orthogonalize = "hybrid", labels = "sources", 
               marginality = NULL, check.marginality = TRUE, 
               omit.projectors = FALSE, 
               which.criteria = c("aefficiency","eefficiency","order"), 
               aliasing.print = TRUE, ...)}
\arguments{
 \item{projectors}{A \code{\link{list}} each of whose components is a \code{\link{projector}}.}
 \item{formula}{An object of class \code{\link{formula}} from which the 
                \code{\link{projector}s} have been obtained. If \code{NULL}, then 
                the \code{differencing} option of \code{orthogonalize} is not available.}
 \item{keep.order}{A \code{\link{logical}} indicating whether the terms should keep 
                   their position in the expanded \code{formula} projector, or 
                   reordered so that main effects precede two-factor 
                   interactions, which precede three-factor interactions and 
                   so on.}
 \item{grandMean}{A \code{\link{logical}} indicating whether the projector for the 
                       grand mean is to be included in the set produced.}
 \item{orthogonalize}{A \code{\link{character}} vector indicating the method for orthogonalizing a 
                      projector to those for terms that occurred previously 
                      in the formula. Three options are available:  
                      \code{hybrid}; \code{differencing}; \code{eigenmethods}, 
                      unless \code{formula} is \code{NULL} in which case 
                      \code{differencing} is not available.
                      The \code{hybrid} option is the most general and uses the 
                      relationships between the projection operators for the 
                      terms in the \code{formula} to decide which 
                      \code{\link{projector}}s to substract and which to 
                      orthogonalize using eigenmethods. The \code{differencing} 
                      option subtracts, from the current \code{\link{projector}}, 
                      those previously orthogonalized \code{\link{projector}}s for 
                      terms whose factors are a subset of the current 
                      \code{\link{projector}}'s factors. The \code{eigenmethods} 
                      option recursively orthogonalizes the \code{\link{projector}}s 
                      using an eigenanalysis of each \code{\link{projector}} 
                      with previously orthogonalized \code{\link{projector}}s.}
 \item{labels}{A \code{\link{character}} nominating the type of labels to be used in 
                      labelling the projectors, and which will be used also 
                      in the output tables, such the tables of the aliasing 
                      in the structure. The two alternatives are \code{terms} and 
                      \code{sources}. Terms have all factors/variables in it 
                      separated by colons (\code{:}). Sources have factors/variables 
                      in them that represent interactions separated by hashes 
                      (\code{#}); if some factors are nested within others, the 
                      nesting factors are surrounded by square brackets 
                      (\code{[} and \code{]}) and separated by colons (\code{:}). 
                      If some generalized, or combined, factors have no marginal 
                      terms, the constituent factors are separated by colons 
                      (\code{:}) and if they interact with other factors in 
                      the source they will be parenthesized.}
 \item{marginality}{A square \code{\link{matrix}} that can be used to supply the 
                       marginality \code{\link{matrix}} when it is desired to overwrite 
                       the calculated marginality \code{\link{matrix}} or when it is not 
                       being calculated. It should consist of zeroes and ones that 
                       gives the marginalites of the terms in the formula. It must have 
                       the row and column names set to the terms from the expanded 
                       \code{formula}, including being in the same order as these terms. 
                       
                       The entry in the ith row and jth column will be one if the 
                       ith term is marginal to the jth term i.e. the column space of the 
                       ith term is a subspace of that for the jth term and so the source 
                       for the jth term is to be made orthogonal to that for the ith term. 
                       Otherwise, the entries are zero. A row and column should not be 
                       included for the grand mean even if \code{grandMean} is \code{TRUE}.}
 \item{check.marginality}{A \code{\link{logical}} indicating whether the marginality matrix, 
                       when it is supplied, is to be checked against that computed 
                       by \code{\link{porthogonalize.list}}. It is ignored when 
                       \code{orthogonalize} is set to \code{eigenmethods}.}
\item{omit.projectors}{A \code{\link{logical}}, which, if \code{TRUE}, 
           results in the \code{\link{projector}}s in the \code{Q} of the 
           \code{pstructure.object} being replaced by their degrees of freedom. 
           These will be the degrees of freedom of the sources. 
           This option is included a device for saving storage when the 
           \code{\link{projector}}s are not required for further analysis.}
 \item{which.criteria}{A character \code{vector} nominating the efficiency 
                       criteria to be included in the summary of aliasing 
                       between terms. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency}, 
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}. 
                       If \code{none}, no summary is printed.}
 \item{aliasing.print}{A \code{\link{logical}} indicating whether the aliasing between sources 
                       within the structure is to be printed.}
 \item{...}{further arguments passed to \code{terms}.}
}
\value{A \code{\link{pstructure.object}}.}
\details{It is envisaged that the \code{\link{projector}s} in the \code{\link{list}} 
        supplied to the \code{projectors} argument correspond to the terms in a 
        linear model. One way to generate them is to obtain the design matrix \bold{X} 
        for a term and then calculate its projector as 
        \eqn{\mathbf{X(X'X)^-X'}}{X(X'X)^-X'}, There are three methods available for 
        orhtogonalizing the supplied projectors: \code{differencing}, 
        \code{eigenmethods} or the default \code{hybrid} method. 

         \code{Differencing} relies on 
         comparing the factors involved in two terms, one previous to the 
         other, to identify whether to subtract the orthogonalized projector 
         for the  previous term from the primary projector of the other. It 
         does so if factors/variables for the previous term are a subset of 
         the factors/variablesfor for the other term. This relies on ensuring that all 
         projectors whose factors/variables are a subset of the current 
         projector occur before it in the expanded formula. It is checked that 
         the set of matrices are mutually orthogonal. If they are not then 
         a warning is given. It may happen that differencing does not produce 
         a projector, in which case \code{eigenmethods} must be used.

         \code{Eigenmethods} forces each projector to be orthogonal 
         to all terms previous to it in the expanded formula. It uses 
         equation 4.10 of James and Wilkinson (1971), which involves 
         calculating the canonical efficiency factors for pairs of primary 
         projectors. It produces a 
         table of efficiency criteria for partially aliased terms. Again,
         the order of terms is crucial. This method has the disadvantage that 
         the marginality of terms is not determined and so sources names are set 
         to be the same as the term names, unless a \code{marginality} matrix 
         is supplied.
         
         The \code{hybrid} method is the most general and uses the relationships 
         between the projection operators for the terms in the \code{formula} 
         to decide which projectors to subtract and which to orthogonalize using 
         eigenmethods. If \eqn{\mathbf{Q}_i}{Q_i} and \eqn{\mathbf{Q}_j}{Q_j} are 
         two projectors for two different terms, with \eqn{i < j}{i < j}, then
         \enumerate{
         \item if \eqn{\mathbf{Q}_j\mathbf{Q}_i \neq \mathbf{0}}{Q_j Q_i != 0} then 
               have to orthogonalize \eqn{\mathbf{Q}_j}{Q_j} to \eqn{\mathbf{Q}_i}{Q_i}.
         \item if \eqn{\mathbf{Q}_j\mathbf{Q}_i = \mathbf{Q}_j}{Q_j Q_i = Q_j}  
               then,  if \eqn{\mathbf{Q}_i = \mathbf{Q}_j}{Q_i = Q_j}, they are equal 
               and \eqn{\mathbf{Q}_j}{Q_j} will be removed from the list of terms; 
               otherwise they are marginal and \eqn{\mathbf{Q}_i}{Q_i} is subtracted 
               from \eqn{\mathbf{Q}_j}{Q_j}.
         \item if have to orthogonalize and 
               \eqn{\mathbf{Q}_j\mathbf{Q}_i = \mathbf{Q}_i}{Q_j Q_i = Q_i} then 
               \eqn{\mathbf{Q}_j}{Q_j} is aliased with previous terms and will be 
               removed from the list of terms; otherwise \eqn{\mathbf{Q}_i}{Q_i} is 
               partially aliased with \eqn{\mathbf{Q}_j}{Q_j} and 
               \eqn{\mathbf{Q}_j}{Q_j} is orthogonalized to \eqn{\mathbf{Q}_i}{Q_i} 
               using eigenmethods.
         }
         The order of projections matrices in the \code{\link{list}} is crucial in this process.
         
        Of the three methods, \code{eigenmethods} is least likely to fail, but it 
        does not establish the marginality between the terms. It is often needed 
        when there is nonorthogonality between terms, such as when there are several 
        linear covariates. It can also be more efficeint in these circumstances.
        
        The process can be computationally expensive, particularly for a large data set 
        (500 or more observations) and/or when many terms are to be orthogonalized. 
             
        If the error \code{Matrix is not idempotent} should occur then, especially 
        if there are  many terms, one might try using \code{\link{set.daeTolerance}} 
        to reduce the tolerance used in determining if values are either the same 
        or are zero; it may be necessary to lower the tolerance to as low as 0.001. 
        Also, setting \code{orthogonalize} to \code{eigenmethods} is worth a try.
}
\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{pstructure.formula}}, \code{\link{proj2.efficiency}}, 
         \code{\link{proj2.combine}}, \code{\link{proj2.eigen}}, \cr
         \code{\link{projs.2canon}} in package \pkg{dae}, 
         \code{\link{eigen}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)
## manually obtain projectors for units
Q.G <- projector(matrix(1, nrow=24, ncol=24)/24)                         
Q.B <- projector(fac.meanop(PBIBD2.lay$Block))
Q.BU <- projector(diag(1, nrow=24))

## manually obtain projector for trt
Q.T <- projector(fac.meanop(PBIBD2.lay$trt) - Q.G)

##Orthogonalize the projectors using porthogonalize.list
Qs <- list(Mean = Q.G, Block = Q.B, "Block:Unit" = Q.BU)
struct <- porthogonalize(Qs, grandMean = TRUE)
Qs <- struct$Q
(lapply(Qs, degfree))

#Add a linear covariate
PBIBD2.lay <- within(PBIBD2.lay,
                     {
                       cBlock <- as.numfac(Block)
                       cBlock <- cBlock - mean(unique(cBlock))
                     })
X <- model.matrix(~ cBlock, data = PBIBD2.lay)
Q.cB <- projector(X \%*\% mat.ginv(t(X) \%*\% X) \%*\% t(X))
Qs <- list(cBlock = Q.cB, Block = Q.B, "Block:Unit" = Q.BU)
struct <- porthogonalize(Qs, grandMean = FALSE)
Qs <- struct$Q
(lapply(Qs, degfree))
}
\keyword{array}
\keyword{design}
\keyword{projector}
