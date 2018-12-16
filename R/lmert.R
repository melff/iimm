#' @title Inference based on a t-distribution for 'lmer' objects.
#' 
#' @description Function \code{lmer_t} provides confidence intervals and p-values
#' based on the Kenward-Roger method for models fitted with \code{lmer}.
#'
#' @return The function returns a list with two elements \code{"lmer"} - the original
#' object - and \code{"coefTab"} which is a matrix with estimates, standard errors,
#' t-statistics, degrees of freedom, lower and upper confidence intervals, and p-values.
#'
#' @param object an object representing a linear mixed effects model, as returned by \code{lmer}
#' @param method a character string that selects the method used to determine the degrees of freedom.
#' @param level a number between 0 and 1, the level of the confidence intervals
#' @export
#' @import lme4 pbkrtest
#' @import Matrix 
lmer_t <- function(object,
                   method = c("Heuristic","Kenward-Roger"),
                   level = .95){
    alpha <- (1-level)/2
    method <- match.arg(method)

    coefs <- fixef(object)

    if(method=="Heuristic"){
        se <- sqrt(Matrix::diag(vcov(object)))
        tval <- coefs/se
        X <- getME(object,"X")
        grps <- getME(object,"flist")
        ddf <- getDF(X,grps)
        pval <- 2*pt(abs(tval),df=ddf,lower.tail=FALSE)
        wdth2 <- qt(1-alpha,df=ddf)*se
        ci.coefs <- coefs + cbind(-wdth2,wdth2)
    }
    else if(method=="Kenward-Roger") {
        vcov.raw <- vcov(object)
        vcov.adj <- vcovAdj(object)
        se.raw <- sqrt(Matrix::diag(vcov.raw))
        se <- sqrt(Matrix::diag(vcov.adj))
        tval <- coefs/se
        ddf <- numeric(length=length(coefs))
        for(i in 1:length(coefs)){
            l <- numeric(length(coefs))
            l[i] <- 1
            
            ddf[i] <- ddf_Lb(vcov.adj,l,vcov.raw)
        }
        pval <- 2*pt(abs(tval),df=ddf,lower.tail=FALSE)
        wdth2 <- qt(1-alpha,df=ddf)*se
        ci.coefs <- coefs + cbind(-wdth2,wdth2)
    }
    else{
        stop("method='",method,"' not yet supported.")       
    }

    coefTab <- cbind(
        Estimate=coefs,
        "Std.Err"=se,
        "t value"=tval,
        "  Df"=ddf,
        "Lower"=ci.coefs[,1],
        "Upper"=ci.coefs[,2],
        "Pr(>|t|)"=pval
    )
  
    ans <- list(
        lmer=object,
        coefTab=coefTab)
    class(ans) <- "lmer_t"
    ans
}

#' @export
print.lmer_t <- function(x,...) print(x$lmer,...)

#' @export
summary.lmer_t <- function(object, ...){
    ans <- list(
        summary.lmer = summary(object$lmer),
        coefTab = object$coefTab
    )
    class(ans) <- "summary.lmer_t"
    ans
}

#' @export
print.summary.lmer_t <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                #symbolic.cor = x$symbolic.cor, 
                                signif.stars = getOption("show.signif.stars"),
                                ranef.comp = c("Variance", "Std.Dev."),
                                show.resids = TRUE,
                                ...){
    coefTab <- x$coefTab

    s <- x$summary.lmer

    lme4:::.prt.methTit(s$methTitle, s$objClass)
    lme4:::.prt.family(s)
    lme4:::.prt.call(s$call); cat("\n")
    lme4:::.prt.aictab(s$AICtab); cat("\n")
    if (show.resids)
        ## need residuals.merMod() rather than residuals():
        ##  summary.merMod has no residuals method
        lme4:::.prt.resids(s$residuals, digits = digits)
    lme4:::.prt.VC(s$varcor, digits = digits, useScale = s$useScale,
            comp = ranef.comp, ...)
    lme4:::.prt.grps(s$ngrps, nobs = s$devcomp$dims[["n"]])
    
    
    nc <- ncol(coefTab)
    zap.i <- 1L:(nc-1)
    #coefTab[,4] <- round(coefTab[,4],1)
    pcoefTab <- cbind(
        format(round(coefTab[,1],digits),digits=digits,scientific=FALSE),
        format(round(coefTab[,2],digits),digits=digits,scientific=FALSE),
        format(round(coefTab[,3],digits),digits=digits,scientific=FALSE,justify="centre",width=8),
        format(round(coefTab[,4],1),scientific=FALSE),
        format(round(coefTab[,5:6],digits),digits=digits,scientific=FALSE),
        format.pval(coefTab[,7],digits=digits)
    )
    dimnames(pcoefTab) <- dimnames(coefTab)

    if(signif.stars){
        Signif <- symnum(coefTab[,7], corr = FALSE, na = FALSE,
                         cutpoints = c(0,  .001,.01,.05, .1, 1),
                         symbols   =  c("***","**","*","."," "))
        pcoefTab <- cbind(pcoefTab, format(Signif))
    }

    cat("\nCoefficients:\n")
    
    print.default(pcoefTab,quote=FALSE,
                  na.print="NA")
    
    if(signif.stars){
        sig.legend <- attr(Signif,"legend")
        w <- getOption("width")
        if(nchar(sig.legend) > w)
            sig.legend <- strwrap(sig.legend, width = w - 2,
                                  prefix = "  ")
        cat("---\nSignif. codes:  ", sig.legend, sep = "",
            fill = w + 4 + max(nchar(sig.legend,"bytes") - nchar(sig.legend)))
    }

    if(length(s$fitMsgs) && any(nchar(s$fitMsgs) > 0)) {
        cat("fit warnings:\n"); writeLines(s$fitMsgs)
    }
    lme4:::.prt.warn(s$optinfo,summary=FALSE)
    
    invisible(x)    
}

# Heuristically determine the degrees of freedom for all variables that are collumns in
# matrix X.
# Returns a matrix with one column for each random effects level and one row for
# each column in X.
getDF <- function(X,grps) pmin(sapply(grps,getDF1,X))

getDF1 <- function(g,X){
    res <- numeric()
    m <- nlevels(g)
    n <- nrow(X)
    p <- ncol(X)
    l <- 0
    has.interc <- FALSE
    machine.eps <- .Machine$double.eps
    has.within <- rep(FALSE,p)
    namesX <- colnames(X)
    names(has.within) <- namesX
    interpos <- match("(Intercept)",namesX,nomatch=0L)
    for(i in 1:ncol(X)){
        x <- X[,i]
        if(var(x) <= machine.eps){
            has.interc <- TRUE
        } else {
            x.bar <- ave(x,g)
            var.within <- var(x - x.bar)
            var.between <- var(x.bar)
            if(var.within > var.between)
                has.within[i] <- TRUE
            else
                l <- l + 1
        }
    }
    within.df  <- n - (p - (m - l))
    between.df <- m - l
    res <- ifelse(has.within,within.df,between.df)
    if(interpos > 0){
        res <- res - 1
        res[interpos] <- m - 1
    }
    res
}

#' @import memisc 
#' @export
getSummary.lmer_t <- function (obj, alpha = 0.05, ...) {

    smry <- summary(obj$lmer)
    ctab <- obj$coefTab

    est  <- ctab[,1]
    se   <- ctab[,2]
    tval <- ctab[,3]
    ddf  <- ctab[,4]
    pval <- ctab[,7]
    
    lower <- est + se*qt(p = alpha/2,df=ddf)
    upper <- est + se*qt(p = 1 - alpha/2,df=ddf)

    coef <- cbind(est,se,tval,pval,lower,upper)

    dn.cf <- list(
        rownames(coef),
        c("est","se","stat","p","lwr","upr"),
        names(obj$lmer@frame)[1]
    )
    dim(coef) <- c(dim(coef)[1],dim(coef)[2],1)
    dimnames(coef) <- dn.cf
    
    varcor <- smry$varcor

    VarPar <- list()
    VarPar.names <- c()
    
    for(i in seq_along(varcor)){
        vc.i <- varcor[[i]]
        lv.i <- names(varcor)[i]
        vr.i <- diag(vc.i)
        cv.i <- vc.i[lower.tri(vc.i)]
        nms.i <- rownames(vc.i)
        nms.i <- gsub("(Intercept)","1",nms.i,fixed=TRUE)
        vrnames.i <- paste0("Var(~",nms.i,"|",lv.i,")")
        cvnames.i <- t(outer(nms.i,nms.i,FUN=paste,sep=":"))
        cvnames.i <- cvnames.i[lower.tri(cvnames.i)]
        if(length(cvnames.i))
            cvnames.i <- paste0("Cov(~",cvnames.i,"|",lv.i,")")
        vp.i <- matrix(NA,nrow=length(vr.i)+length(cv.i),ncol=6)
        vp.i[,1] <- c(vr.i,cv.i)
        dim(vp.i) <- c(dim(vp.i),1)
        dimnames(vp.i) <- list(c(vrnames.i,cvnames.i),
                               c("est","se","stat","p","lwr","upr"),
                               names(obj$lmer@frame)[1])
        VarPar <- c(VarPar,list(vp.i))
        VarPar.names <- c(paste0("Var(",lv.i,")"),
                          VarPar.names)
    }
    if(smry$sigma>1){
        vp.i <- matrix(NA,nrow=1,ncol=6)
        vp.i[1] <- smry$sigma^2
        dim(vp.i) <- c(dim(vp.i),1)
        dimnames(vp.i) <- list("Var(residual)",
                               c("est","se","stat","p","lwr","upr"),
                               names(obj$lmer@frame)[1])
        VarPar <- c(list(vp.i),VarPar)
        VarPar.names <- c("Var(residual)",VarPar.names)
    }
    names(VarPar) <- VarPar.names
    
    ## Factor levels.
    xlevels <- list()
    Contr <- names(attr(model.matrix(obj$lmer), "contrasts"))
    for (c in Contr) xlevels[[c]] <- levels(obj$lmer@frame[,c])

    ## Model fit statistics.
    ll <- logLik(obj$lmer)[1]
    isREML <- inherits(obj$lmer@resp,"lmerResp") && obj$lmer@resp$REML > 0
    if(!isREML)
        deviance <- deviance(obj$lmer)
    else 
        deviance <- lme4::REMLcrit(obj$lmer)
    AIC <- AIC(obj$lmer)
    BIC <- BIC(obj$lmer)
    
    N <- c(Total=nobs(obj$lmer))
    G <-as.integer(smry$ngrps)
    names(G) <- names(smry$ngrps)
    G <- c(N,G)
    
    sumstat <- c(logLik = ll,
                 deviance = deviance,
                 AIC = AIC,
                 BIC = BIC)
    ## Return model summary.
    
    ans <- list(coef= coef)

    ans <- c(ans,VarPar)
    
    ans <- c(ans,       
             list(Groups = G,
                  sumstat = sumstat,
                  contrasts = Contr, ## Reuse 'Contr'
                  xlevels = xlevels,
                  call = obj$lmer@call))
    return(ans)
}
