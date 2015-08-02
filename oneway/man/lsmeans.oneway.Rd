\name{lsmeans.oneway}
\alias{lsmeans.oneway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lsmeans.oneway(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    object <- summary(object)
    if (object$P > 0.05) {
        warning("F-test is not significant at alpha=0.05.")
    }
    compare <- function(i, j) {
        d <- object$ls.means[i] - object$ls.means[j]
        s.e <- sqrt(object$MS[2] * (1/object$n.vec[i] + 1/object$n.vec[j]))
        t.val <- d/s.e
        round(2 * pt(abs(t.val), object$N - object$a, lower.tail = FALSE), 
            4)
    }
    p.vals <- pairwise.table(compare.levels = compare, level.names = object$groups, 
        p.adjust.method = "none")
    result <- list(p.value = p.vals, call = match.call())
    class(result) <- "lsmeans"
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
