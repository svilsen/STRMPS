.ibeta <- function(k, lambda, theta) {
    pbeta(lambda/(theta + lambda), k, theta)*beta(k, theta)
}

.dLambdaIbeta <- function(k, lambda, theta) {
    numerator <- lambda^(k-1)*theta^theta
    denominator <- (lambda + theta)^(k + theta)
    numerator/denominator
}

.dThetaIbeta <- function(k, lambda, theta) {
    numerator <- lambda^k * theta^(theta - 1)
    denominator <- (lambda + theta)^(k + theta)
    integrand <- function(x) x^(k - 1) * (1 - x)^(theta - 1) * log(1 - x)
    integral.ibeta <- integrate(integrand, lower = 0, upper = lambda/(lambda + theta), subdivisions = 1000L, stop.on.error = F)$value
    return(-numerator/denominator + integral.ibeta)
}

kiNegBin <- function(parms, X, K, Y, Y0, Y1, weights, offsetx, offsetz, kx, kz, kInflation) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(exp(K %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(kInflation, size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    loglik
}

gradientkiNegBin <- function(parms, X, K, Y, Y0, Y1, weights, offsetx, offsetz, kx, kz, kInflation) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(K %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- exp(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    clogdens0 <- dnbinom(kInflation, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

    wres_count <- ifelse(Y1,
                         Y/mu - (Y + theta)/(mu + theta),
                         exp(-log(dens0) + log(1 - muz) + clogdens0)*((theta + kInflation)/(theta + mu) + kInflation/mu))

    wres_zero <- ifelse(Y1,
                        -1/(1 - muz) * exp(etaz),
                        (exp(etaz) - exp(clogdens0) * exp(etaz))/dens0)

    wres_theta <- theta * ifelse(Y1,
                                 digamma(Y + theta) - digamma(theta) + log(theta) -
                                     log(mu + theta) + 1 - (Y + theta)/(mu + theta),
                                 exp(-log(dens0) + log(1 - muz) + clogdens0) *
                                     (log(theta) - log(mu + theta) + 1 - (theta + Y)/(mu + theta) +
                                          digamma(Y + theta) - digamma(theta))
    )

    colSums(cbind(wres_count * weights * X, wres_zero * weights * K, wres_theta))
}

kiNegBinTruncated <- function(parms, X, K, Y, Y0, Y1, weights, offsetx, offsetz, kx, kz, kInflation) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(exp(K %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(kInflation, size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE)) - suppressWarnings(pbeta((mu/(mu + theta)), kInflation, theta, log.p = TRUE))

    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    loglik
}

gradientkiNegBinTruncated <- function(parms, X, K, Y, Y0, Y1, weights, offsetx, offsetz, kx, kz, kInflation) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(K %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- exp(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    clogdens0 <- dnbinom(kInflation, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

    wres_count <- ifelse(Y1,
                         Y/mu - (Y + theta)/(mu + theta) - 1/.ibeta(kInflation + 1, mu, theta)*.dLambdaIbeta(kInflation + 1, mu, theta),
                         exp(-log(dens0) + log(1 - muz) + clogdens0)*((theta + kInflation)/(theta + mu) + kInflation/mu)
    )

    wres_zero <- ifelse(Y1,
                        -1/(1 - muz) * exp(etaz),
                        (exp(etaz) - exp(clogdens0) * exp(etaz))/dens0
    )

    wres_theta <- theta * ifelse(Y1,
                                 digamma(Y + theta) - digamma(theta + kInflation + 1) + log(theta) -
                                     log(mu + theta) + 1 - (Y + theta)/(mu + theta) -
                                     1/.ibeta(kInflation + 1, mu, theta)*.dThetaIbeta(kInflation + 1, mu, theta),
                                 exp(-log(dens0) + log(1 - muz) + clogdens0) *
                                     (log(theta) - log(mu + theta) + 1 - (theta + Y)/(mu + theta) +
                                          digamma(Y + theta) - digamma(theta))
    )

    colSums(cbind(wres_count * weights * X, wres_zero * weights * K, wres_theta))
}

.modelOffset2 <- function (x, terms = NULL, offset = TRUE) {
    if (is.null(terms))
        terms <- attr(x, "terms")
    offsets <- attr(terms, "offset")
    if (length(offsets) > 0) {
        ans <- if (offset)
            x$"(offset)"
        else NULL
        if (is.null(ans))
            ans <- 0
        for (i in offsets) ans <- ans + x[[deparse(attr(terms,
                                                        "variables")[[i + 1]])]]
        ans
    }
    else {
        ans <- if (offset)
            x$"(offset)"
        else NULL
    }
    if (!is.null(ans) && !is.numeric(ans))
        stop("'offset' must be numeric")
    ans
}

#' K-Inflated Negative Binomial Control Function
#'
#' @return Returns a list of default objects for the \link{kinb}-function.
#' @export
kinb.control <- function (method = "BFGS", maxit = 10000, trace = FALSE, EM = FALSE,
                          start = NULL, return.pars.only = FALSE, ...) {
    rval <- list(method = method, maxit = maxit, trace = trace,
                 EM = EM, start = start, return.pars.only = return.pars.only)
    rval <- c(rval, list(...))
    if (!is.null(rval$fnscale))
        warning("fnscale must not be modified")
    rval$fnscale <- -1
    if (!is.null(rval$hessian))
        warning("hessian must not be modified")
    rval$hessian <- TRUE
    if (method == "L-BFGS-B") {
        if(is.null(rval$factr))
            rval$factr <- 1e-8
        if(is.null(rval$pgtol))
            rval$pgtol <- 0
    }
    else {
        if (is.null(rval$reltol))
            rval$reltol <- .Machine$double.eps^(1/1.6)
    }
    if (is.null(rval$reltol))
        rval$reltol <- .Machine$double.eps^(1/1.6)
    rval
}

#' K-Inflated Negative Binomial Models
#'
#' @description A function fitting k-inflated negative binomial model, (k-1)-truncated and non-truncated.
#'
#' @param truncated a TRUE/FALSE statement indicating whether the model is (k-1)-truncated.
#' @param kInflation the point of inflation. The variable defaults to 1.
#' @param control a \link{kinb.control} object.
#'
#' @return A list of parameter estimates (can return a shorter list if needed).
#'
#' @export
kinb <- function(formula, data, subset, na.action, weights, offset, link = "log",
                 control = kinb.control(), model = TRUE, y = TRUE, x = FALSE,
                 truncated = TRUE, kInflation = 1L, ...) {
    if (truncated) {
        loglikfun <- kiNegBinTruncated
        gradfun <- gradientkiNegBinTruncated
    }
    else {
        loglikfun <- kiNegBin
        gradfun <- gradientkiNegBin
    }

    linkstr <- match.arg(link)
    linkobj <- make.link(linkstr)
    linkinv <- linkobj$linkinv
    if (control$trace)
        cat("K-inflated Count Model\n",
            paste("count model:", dist, "with log link\n"),
            paste("k-inflation model: negative binomial with", linkstr, "link\n"), sep = "")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights",
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
        ff <- formula
        formula[[3]][1] <- call("+")
        mf$formula <- formula
        ffc <- . ~ .
        ffz <- ~.
        ffc[[2]] <- ff[[2]]
        ffc[[3]] <- ff[[3]][[2]]
        ffz[[3]] <- ff[[3]][[3]]
        ffz[[2]] <- NULL
    } else {
        ffz <- ffc <- ff <- formula
        ffz[[2]] <- NULL
    }
    if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
        ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])),
                                         deparse(ffz))))
    }
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    mtX <- terms(ffc, data = data)
    X <- model.matrix(mtX, mf)
    mtK <- terms(ffz, data = data)
    mtK <- terms(update(mtK, ~.), data = data)
    K <- model.matrix(mtK, mf)
    Y <- model.response(mf, "numeric")
    if (length(Y) < 1)
        stop("empty model")
    if (all(Y > kInflation))
        stop("invalid dependent variable, minimum count is not k")
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
        stop("invalid dependent variable, non-integer values")
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0))
        stop("invalid dependent variable, negative counts")
    if (control$trace) {
        cat("dependent variable:\n")
        tab <- table(factor(Y, levels = 0:max(Y)), exclude = NULL)
        names(dimnames(tab)) <- NULL
        print(tab)
    }
    n <- length(Y)
    kx <- NCOL(X)
    kz <- NCOL(K)
    Y0 <- Y <= kInflation
    Y1 <- Y > kInflation
    weights <- model.weights(mf)
    if (is.null(weights))
        weights <- 1
    if (length(weights) == 1)
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)
    offsetx <- .modelOffset2(mf, terms = mtX, offset = TRUE)
    if (is.null(offsetx))
        offsetx <- 0
    if (length(offsetx) == 1)
        offsetx <- rep.int(offsetx, n)
    offsetx <- as.vector(offsetx)
    offsetz <- .modelOffset2(mf, terms = mtK, offset = FALSE)
    if (is.null(offsetz))
        offsetz <- 0
    if (length(offsetz) == 1)
        offsetz <- rep.int(offsetz, n)
    offsetz <- as.vector(offsetz)
    start <- control$start
    if (!is.null(start)) {
        valid <- TRUE
        if (!("count" %in% names(start))) {
            valid <- FALSE
            warning("invalid starting values, count model coefficients not specified")
            start$count <- rep.int(0, kx)
        }
        if (!("k" %in% names(start))) {
            valid <- FALSE
            warning("invalid starting values, k-inflation model coefficients not specified")
            start$k <- rep.int(0, kz)
        }
        if (length(start$count) != kx) {
            valid <- FALSE
            warning("invalid starting values, wrong number of count model coefficients")
        }
        if (length(start$k) != kz) {
            valid <- FALSE
            warning("invalid starting values, wrong number of k-inflation model coefficients")
        }
        if (!("theta" %in% names(start)))
            start$theta <- 1
        start <- list(count = start$count, k = start$k,
                      theta = as.vector(start$theta[1]))
        if (!valid)
            start <- NULL
    }
    if (is.null(start)) {
        if (control$trace)
            cat("generating starting values...")
        model_count <- glm.fit(X, Y, family = poisson(), weights = weights,
                               offset = offsetx)
        model_k <- glm.fit(K, as.integer(Y0), weights = weights,
                           family = binomial(link = linkstr), offset = offsetz)
        start <- list(count = model_count$coefficients, k = model_k$coefficients, theta = 1)
        if (control$EM) {
            mui <- model_count$fitted
            probi <- model_k$fitted
            probi <- probi/(probi + (1 - probi) * dnbinom(kInflation, size = start$theta, mu = mui))
            probi[Y1] <- 0
            ll_new <- loglikfun(c(start$count, start$k, log(start$theta)), X, K, Y, Y0, Y1, weights, offsetx, offsetz, kInflation)
            ll_old <- 2 * ll_new
            while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
                ll_old <- ll_new
                model_count <- suppressWarnings(glm.nb(Y ~ 0 + X + offset(offsetx),
                                                       weights = weights * (1 - probi),
                                                       start = start$count, init.theta = start$theta))
                model_k <- suppressWarnings(glm.fit(K, probi, weights = weights, offset = offsetz,
                                                    family = binomial(link = "log"),
                                                    start = start$k))
                start <- list(count = model_count$coefficients,
                              k = model_k$coefficients, theta = model_count$theta)
                mui <- model_count$fitted
                probi <- model_k$fitted
                probi <- probi/(probi + (1 - probi) * dnbinom(kInflation, size = start$theta, mu = mui))
                probi[Y1] <- 0
                ll_new <- loglikfun(c(start$count, start$k, log(start$theta)), X, K, Y, Y0, Y1, weights, offsetx, offsetz, kInflation)
            }
        }
        if (control$trace)
            cat("done\n")
    }
    if (control$trace)
        cat("calling optim() for ML estimation:\n")
    method <- control$method
    hessian <- control$hessian
    ocontrol <- control
    control$method <- control$hessian <- control$EM <- control$start <- NULL

    fit <- suppressWarnings(optim(fn = loglikfun, gr = gradfun,
                                  par = c(start$count, start$k, log(start$theta)),
                                  X = X, K = K, Y = Y, Y0 = Y0, Y1 = Y1, weights = weights,
                                  offsetx = offsetx, offsetz = offsetz, kx = kx, kz = kz, 
                                  kInflation = kInflation,
                                  method = method, hessian = hessian, control = control))
    
    if (fit$convergence > 0)
        warning("optimization failed to converge")
    coefc <- fit$par[1:kx]
    names(coefc) <- names(start$count) <- colnames(X)
    coefz <- fit$par[(kx + 1):(kx + kz)]
    names(coefz) <- names(start$k) <- colnames(K)
    vc <- -solve(as.matrix(fit$hessian))
    np <- kx + kz + 1
    theta <- as.vector(exp(fit$par[np]))
    SE.logtheta <- as.vector(sqrt(diag(vc)[np]))
    vc <- vc[-np, -np, drop = FALSE]
    colnames(vc) <- rownames(vc) <- c(paste("count", colnames(X),
                                            sep = "_"), paste(kInflation, colnames(K), sep = "_"))
    mu <- exp(X %*% coefc + offsetx)[, 1]
    phi <- exp(K %*% coefz + offsetz)[, 1]
    Yhat <- (1 - phi) * mu
    res <- sqrt(weights) * (Y - Yhat)
    nobs <- sum(weights > 0)

    if(control$return.pars.only) {
        rval <- list(size = theta, mu = coefc, pi = mean(phi))
    }
    else {
        rval <- list(coefficients = list(count = coefc, k = coefz),
                     residuals = res, fitted.values = Yhat, optim = fit, method = method,
                     control = ocontrol, start = start, weights = if (identical(as.vector(weights), rep.int(1L, n))) NULL else weights,
                     offset = list(count = if (identical(offsetx, rep.int(0, n))) NULL else offsetx,
                                   k = if (identical(offsetz, rep.int(0, n))) NULL else offsetz),
                     n = nobs, df.null = nobs - 2, df.residual = nobs - (kx + kz + 1),
                     terms = list(count = mtX, k = mtK, full = mt), theta = theta, size = theta, phi = phi,
                     SE.logtheta = SE.logtheta, loglik = fit$value, vcov = vc,
                     dist = dist, link = "log", linkinv = exp, converged = fit$convergence < 1,
                     call = cl, formula = ff, levels = .getXlevels(mt, mf),
                     contrasts = list(count = attr(X, "contrasts"), k = attr(K, "contrasts")))
        if (model)
            rval$model <- mf
        if (y)
            rval$y <- Y
        if (x)
            rval$x <- list(count = X, k = K)
        if(!truncated) class(rval) <- "kinb" else class(rval) <- "kinb.truncated"
    }
    return(rval)
}

setClass("kinb")
setClass("kinb.truncated")

.createNoiseThreshold <- function(stringCoverageListObject, predefinedProportions, dist = c("kinb", "geometric"), numberOfStandardDeviations = 3, numberOfThreads) {
    dist <- match.arg(dist)
    thresholdList <- structure(mclapply(seq_along(stringCoverageListObject), function(s) {
        df <- data.frame(Coverage = stringCoverageListObject[[s]]$Coverage)
        if (dist == "kinb") {
            kinb_fit <- kinb(formula = Coverage ~ 1, data = df, kInflation = 1, control = kinb.control(EM = TRUE))
            m <- kinb_fit$coefficients$count
            theta <- kinb_fit$theta
            sd <- sqrt(m + m^2/theta)

            thresholds <- qnbinom(predefinedProportions[s], size = theta, mu = m) + numberOfStandardDeviations*sd
        }
        else if (dist == "geometric") {
            geometric_fit <- fitdistr(df$Coverage - 1, densfun = "geometric")
            prob <- as.numeric(geometric_fit$estimate)
            sd <- sqrt((1 - prob)/prob^2)

            thresholds <- (qgeom(predefinedProportions[s], prob) + 1L) + numberOfStandardDeviations*sd
        }
        else  {
            thresholds <- NULL
        }
        return(thresholds)
    }, mc.cores = numberOfThreads), .Names = names(stringCoverageListObject))

    return(do.call(c, thresholdList))
}

#' Sets Noise Threshold
#'
#' @description The function creates noise thresholds based on the provided stringCoverageList-object, and a set of pre-defined proportions, using either a \link{kinb} or a geometric distribution.
#'
#' @param stringCoverageListObject
#' @param predefinedProportions
#' @param dist
#' @param numberOfStandardDeviations
#' @param numberOfThreads
#'
#' @return A set of thresholds for noise identification.
#'
#' @export
setGeneric("noiseThresholds", signature = "stringCoverageListObject",
           function(stringCoverageListObject, predefinedProportions, numberOfStandardDeviations, numberOfThreads)
               standardGeneric("noiseThresholds")
)

setMethod("noiseThresholds", "stringCoverageList",
          function(stringCoverageListObject, predefinedProportions, numberOfStandardDeviations = 3, numberOfThreads)
              .createNoiseThreshold(stringCoverageListObject, predefinedProportions, numberOfStandardDeviations, numberOfThreads)
)
