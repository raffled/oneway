\name{oneway.default}
\alias{oneway.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oneway.default(z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
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
function (z, ...) 
{
    if (!is.list(z) | is.data.frame(z)) {
        stop("data structure is not one of:\nlist, factor & vector, or formula")
    }
    a <- length(z)
    N <- length(unlist(z))
    TT <- sum(unlist(z)^2)
    A <- sum(sapply(z, function(x) sum(x)^2/length(x)))
    CF <- sum(unlist(z))^2/N
    ss.a <- A - CF
    ss.e <- TT - A
    df.a <- a - 1
    df.e <- N - a
    ybar.vec <- unlist(lapply(z, mean))
    n.vec <- unlist(lapply(z, length))
    if (is.null(names(z))) {
        groups <- as.character(1:length(z))
    }
    else {
        groups <- names(z)
    }
    names(ybar.vec) <- groups
    res <- list(df = c(df.a, df.e), SS = c(ss.a, ss.e), groups = groups, 
        call = match.call(), data = z, ls.means = ybar.vec, n.vec = n.vec, 
        N = N, a = a)
    class(res) <- "oneway"
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
