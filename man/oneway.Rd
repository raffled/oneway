\name{oneway}
\alias{oneway}
\alias{oneway.default}
\alias{oneeway.formula}
\alias{oneway.factor}
\alias{print.oneway}
\alias{summary.oneway}
\alias{print.summary.oneway}
\alias{lsmeans}
\alias{lsmeans.default}
\alias{print.lsmeans}
\alias{plot.oneway}
\title{
  Oneway Analysis of Variance
}
\description{
  \code{oneway} performs a oneway analysis of variance (ANOVA) for data
  with one factor.  Fisher's LSD is implemented with the \code{lsmeans}
  function.
}
\usage{
  oneway(z, ...)
  \method{oneway}{default}(z, ...)
  \method{oneway}{factor}(z, y, ...)
  \method{oneway}{formula}(formula, data=list(), ...)
  \method{print}{oneway}(x, ...)
  \method{summary}{oneway}(object, ...)
  \method{print}{summary.oneway}(x, ...)
  \method{lsmeans}{default}(object, ...)
  \method{lsmeans}{oneway}(object, ...)
  \method{plot}{oneway}(x, names=x$groups, xlab="Group", ylab="Response",
  main=capture.output(x$call), ...)
}

\arguments{
  \item{z}{
    a (possibly named) list where each element is a vector of
    observations for a level of the factor or a factor, depending on how
    the function is called.
  }
  \item{y}{ a vector of observations corresponding to the factor z }
  \item{formula}{ the model to fit of the form \code{response~factor} }
  \item{data}{ an optional dataframe containing the variables in
    the model represented by \code{formula} }
  \item{x}{ a oneway object }
  \item{object}{ a oneway object }
  \item{names}{ default value to use the names from the list as
    labels in the boxplot }
  \item{xlab}{ default value to used \"Group\" as the xlabel of the
    plot }
  \item{ylab}{ default option to use \"Response\" as the ylabel
    of the plot }
  \item{main}{ default option to use the function
    all that created the oneway object as the title of the plot }
  \item{\dots}{ only used by \code{plot()} to pass optional graphics
    parameters. }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \code{oneway()} returns an object of class \code{oneway}, which contains:
  \item{df}{a vector containing the degrees of freedom (among, within)}
  \item{SS}{a vector containing the sums of squares (among, within)}
  \item{groups}{a vector of the group (treatment) level names}
  \item{call}{the function call that created the object}
  \item{data}{a list of the original data (see \code{arguments})}
  \item{ls.means}{a vector of the group (treatment) means}
  \item{n.vec}{a vector of the group (treatment) sample sizes}
  \item{N}{the total number of observations}
  \item{a}{the number of groups (treatments)}


  \code{print()} gives a succinct summary of the data (degrees of freedom
  and sums of squares).

  \code{summary()} creates/displays the \code{ANOVA} table and group means.

  \code{plot()} displays a side-by-side boxplot of the groups (treatments).

  \code{lsmeans()} performs and displays Fisher's LSD procedure on a
  \code{oneway} object.  This function will throw an error if called on
  any object of any other class.
}
\author{
Doug Raffle
}
\examples{
data(coagulation)
attach(coagulation)

##From a list:
coag.list <- split(coag, diet)
(coag.oneway <- oneway(coag.list))
summary(coag.oneway)
lsmeans(coag.oneway)
plot(coag.oneway)

## From a factor & vector:
oneway(diet, coag)

##From a formula:
oneway(coag~diet)

detach(coagulation)
}
