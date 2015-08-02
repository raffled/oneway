\name{summary.oneway}
\alias{summary.oneway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.oneway(object, ...)
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
    attach(object)
    ss.t <- SS[1] + SS[2]
    df.t <- df[1] + df[2]
    ms.a <- SS[1]/df[1]
    ms.e <- SS[2]/df[2]
    F <- ms.a/ms.e
    p <- pf(F, df[1], df[2], lower.tail = FALSE)
    tab <- with(object, cbind(DF = c(df, df.t), SS = c(SS, ss.t), 
        MS = c(ms.a, ms.e, NA), F = c(F, NA, NA), `Pr(>F)` = c(p, 
            NA, NA)))
    rownames(tab) <- c("Among Groups", "Within Groups (error)", 
        "Total")
    res <- list(call = call, tab = tab, groups = groups, ls.means = ls.means, 
        P = p, MS = c(ms.a, ms.e), n.vec = n.vec, N = N, a = a)
    class(res) <- "summary.oneway"
    detach(object)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
