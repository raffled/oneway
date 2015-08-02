######################################################################
#####################           oneway.R          ####################
######################################################################
## contains source code for oneway S3 class.                        ##
######################################################################

######################################################################
####################          Methods           ######################
######################################################################
oneway <- function(z, ...) UseMethod("oneway")

####################      oneway.default      ########################
## Default method, takes list & returns oneway object.
oneway.default <- function(z, ...) {
    ## Make sure we're the right data structure
    if(!is.list(z) | is.data.frame(z)){
        stop("data structure is not one of:\nlist, factor & vector, or formula")
    }
    ## computing formulas
    a <- length(z)
    N <- length(unlist(z))
    TT <- sum(unlist(z)^2)
    A <- sum(sapply(z, function(x) sum(x)^2/length(x)))
    CF <- sum(unlist(z))^2/N
    ## computes SS's and DF's
    ss.a <- A - CF ## within
    ss.e <- TT - A ## between
    df.a <- a - 1
    df.e <- N - a
    ## get LS Means & sample sizes
    ybar.vec <- unlist(lapply(z, mean))
    n.vec <- unlist(lapply(z, length))
    ## find group names
    if(is.null(names(z))){
        groups <- as.character(1:length(z))
    }
    else{
        groups <- names(z)
    }
    names(ybar.vec) <- groups
    ## create object, add class, & return
    res <- list(df = c(df.a, df.e), SS = c(ss.a, ss.e),
                groups = groups, call = match.call(), data = z,
                ls.means = ybar.vec, n.vec = n.vec, N = N, a = a)
    class(res) <- "oneway"
    return(res)
}

####################      oneway.factor      #########################
## takes a vector & factor, calls oneway.default
oneway.factor <- function(z, y, ...) {
    foo <- oneway.default(split(y, z))
    foo$call <- match.call()
    foo
}

####################      oneway.formula      ########################
## takes formula & extracts vector & factor, calls oneway.factor
oneway.formula <- function(formula, data=list(), ...) {
    mf <- model.frame(formula, data)
    foo <- oneway.factor(mf[,2], mf[,1])
    foo$call <- match.call()
    foo
}

######################################################################
####################          Utility           ######################
######################################################################

####################       print.oneway       ########################
## prints basic summary
print.oneway <- function(x, ...) {
   print(x$call)
   cat("\nWithin SS:", x$SS[1], "on", x$df[1],
       "degrees of freedom.\n")
   cat("Between SS:", x$SS[2], "on", x$df[2],
       "degrees of freedom.\n")
}

###################       summary.oneway       #######################
## creates summary object
summary.oneway <- function(object, ...) {
    attach(object)
    ## Get total SS & df
    ss.t <- SS[1] + SS[2]
    df.t <- df[1] + df[2]
    ## Calculate mean squares
    ms.a <- SS[1]/df[1]
    ms.e <- SS[2]/df[2]
    ## get F stat & p-val
    F <- ms.a/ms.e
    p <- pf(F, df[1], df[2], lower.tail = FALSE)

    ## construct AOV table
    tab <- with(object, cbind(DF = c(df, df.t),
                         SS = c(SS, ss.t),
                         MS = c(ms.a, ms.e, NA),
                         F = c(F, NA, NA),
                         "Pr(>F)" = c(p, NA, NA)))
    rownames(tab) <- c("Among Groups", "Within Groups (error)",
                       "Total")
    res <- list(call=call, tab=tab, groups=groups, ls.means=ls.means,
                P = p, MS = c(ms.a, ms.e), n.vec = n.vec, N = N, a = a)
    class(res) <- "summary.oneway"
    detach(object)
    return(res)
}

################       print.summary.oneway       ####################
## prints the summary object
print.summary.oneway <- function(x, ...) {
    ## function call
    cat("Call:\n\t")
    print(x$call)
    ## LS means
    cat("\nMeans:\n")
    print(x$ls.means)
    cat("\n")
    # AOV Table
    cat("\nAnalysis of Variance Table:\n")
    printCoefmat(x$tab, P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars=TRUE, na.print="")
}

####################       print.oneway       ########################
## prints side-by-sidee boxplot of data
plot.oneway <- function(x, names=x$groups, xlab="Group", ylab="Response", main=capture.output(x$call), ...){
    boxplot(x=x$data, names=names, xlab=xlab, ylab=ylab, main=main,
            ...)
}

######################################################################
####################          lsmeans           ######################
######################################################################
## perform's Fisher's LSD test for a oneway object

####################         methods          ########################
lsmeans <- function(object, ...) UseMethod("lsmeans")
lsmeans.default <- function(object, ...){
    if(!(class(object)=="oneway")){
        stop("lsmeans only accepts class \"oneway\"")
    }
}
lsmeans.oneway <- function(object, ...) {
    object <- summary(object)
    if(object$P > 0.05){
        warning("F-test is not significant at alpha=0.05.")
    }
    compare <- function(i, j){
        d <- object$ls.means[i] - object$ls.means[j]
        s.e <- sqrt(object$MS[2]*(1/object$n.vec[i] + 1/object$n.vec[j]))
        t.val <- d/s.e
        round(2*pt(abs(t.val), object$N-object$a, lower.tail=FALSE),4)
    }
    p.vals <- pairwise.table(compare.levels=compare,
                             level.names=object$groups,
                             p.adjust.method="none")
    result <- list(p.value=p.vals, call=match.call())
    class(result) <- "lsmeans"
    result
}
#####################         print          #########################
print.lsmeans <- function(x, ...){
    cat("Call:\n\t")
    print(x$call)
    cat("\nFisher's LSD Table\n")
    cat("\nP-Values:\n")
    print.table(x$p.value, na.print="-")
}
