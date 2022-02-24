#' This implements the stochastic search based on Simulated Anneling strategy.
#'
#' Highest posterior model is widely accepted as a good model among available
#' models. In terms of variable selection highest posterior model is often the true
#' model. Our stochastic search process SAHPM based on simulated annealing maximization method
#' tries to find the highest posterior model by maximizing the model space with
#' respect to the posterior probabilities of the models. This function currently
#' contains the SAHPM method only for linear models. The codes for GLM will be
#' added in future.
#'
#' The model is:
#' \deqn{y= \alpha + X\beta+\epsilon, \epsilon \sim N(0,\sigma^2)}
#' The Zellner's \eqn{g} prior is used with default \eqn{g = n}.
#'
#' @references Maity, A., K., and Basu, S. Highest Posterior Model Computation and
#' Variable Selection via the Simulated Annealing
#'
#' @export
#'
#' @param formula an object of class \code{\link{formula}} (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in the
#' model. If not found in data, the variables are taken from environment(formula),
#' typically the environment from which \code{\link{lm}} is called.
#' @param na.action a function which indicates what should happen when the data contain
#' \code{NA}s.  The default is set by the \code{na.action} setting of
#' \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.
#' The \dQuote{factory-fresh} default is \code{\link{na.omit}}.  Another possible
#' value is \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.
#' @param g value of \eqn{g} for \eqn{g} prior. Default is sample size \eqn{n}.
#' @param nstep maximum number of steps for simulated annealing search.
#' @param abstol desired level of difference of marginal likelihoods between two steps.
#' @param replace logical. If \code{TRUE} the replce step is considered in the search. Default is FALSE.
#' @param burnin logical. If \code{TRUE} the burnin is added. Default is FALSE. Number of burnin is specified by the next input.
#' @param nburnin Number of burnin (required if burnin = TRUE). Default is 50.
#'
#'
#' @return \item{final.model}{A column vector which corresponds to the original
#' variable indices.}
#' \item{history}{A history of the search process. By columns: Step number,
#' temperature, current objective function value, current minimal objective
#' function value, current model, posterior probability of current model.}
#'
#' @examples
#' require(mvtnorm)     # for multivariate normal distribution
#' n <- 100             # sample size
#' k <- 40              # number of variables
#' z <- as.vector(rmvnorm(1, mean = rep(0, n), sigma = diag(n)))
#' x <- matrix(NA, nrow = n, ncol = k)
#' for(i in 1:k)
#' {
#' x[, i] <- as.vector(rmvnorm(1, mean = rep(0, n), sigma = diag(n))) + z
#' }                    # this induce 0.5 correlation among the variables
#' beta <- c(rep(0, 10), rep(2, 10), rep(0, 10), rep(2, 10))
#'                      # vector of coefficients
#' sigma <- 1
#' sigma.square <- sigma^2
#' linear.pred <- x %*% beta
#' y <- as.numeric(t(rmvnorm(1, mean = linear.pred, sigma = diag(sigma.square, n))))
#'                      # response
#' answer <- sahpmlm(formula = y ~ x)
#' answer$final.model
#' answer$history
#'
#'
#' \dontrun{
#' # With small effect size
#' beta <- c(rep(0, 10), rep(1, 10), rep(0, 10), rep(1, 10))
#'                      # vector of coefficients
#'
#' linear.pred <- x %*% beta
#' y <- as.numeric(t(rmvnorm(1, mean = linear.pred, sigma = diag(sigma.square, n))))
#'                      # response
#' answer <- sahpmlm(formula = y ~ x)
#' answer$final.model  # Might miss some of the true predictors
#' answer$history
#'
#' # Able to recover all the predictors with 50 burnin
#' answer <- sahpmlm(formula = y ~ x, burnin = TRUE, nburnin = 50)
#' answer$final.model  # Misses some of the true predictors
#' answer$history
#' }


sahpmlm <- function(formula, data, na.action, g = n, nstep = 200, abstol = 0.0000001,
                    replace = FALSE, burnin = FALSE, nburnin = 50)
{
  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- stats::model.response(mf, "numeric")
  if (length(dim(y)) == 1) {
    nm <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nm))
      names(y) <- nm
  }

  if (stats::is.empty.model(mt)) {
    X <- NULL
  }
  else {
    X <- stats::model.matrix(mt, mf, stats::contrasts)
  }

  x     <- X[, -1]
  k     <- ncol(x)   # number of variables
  n     <- nrow(x)   # sample size
  error <- 1


  # Fixing the parameters
  T0 <- 1
  Tf <- 0.1
  C  <- 0.9  # Cooling Constant


  current.model <- NULL
  current.T     <- T0  # fixes the temperature
  current.gamma <- rep(0, k)  # gamma (indicator) vector
  # current.gamma <- c(rep(1, 0), rep(0, 10), rep(1, 10), rep(0, 10))
  current.step  <- 1  # Number of steps
  model.matrix  <- NULL
  trace.matrix <- matrix(NA, nstep, 4)
  # A matrix which contains the history of the algorithm. (By columns: Step number,
  # temperature, current objective function value, current minimal objective
  # function value, current model, current highest posterior probability model)


  model <- function(current.gamma)
  {
    return(paste0(which(current.gamma == 1), collapse = "-"))
  }



if(burnin == TRUE)
  {

  marginal.likelihood <- function(current.gamma)
  { # Returns the log of marginal likelihood given in the article

    current.model <- model(current.gamma)

    if(sum(current.model == model.matrix[, 1]) >= 1)
    {

      result <- as.numeric(model.matrix[which(current.model == model.matrix[, 1]), 2])[1]

    } else {

      subset <- which(current.gamma == 1)
      x.sub <- x[, subset]

      if(is.vector(x.sub)) {
        x.sub <- matrix(x.sub, ncol=1)
        k <- 1} else {
          k <- ncol(x.sub)
        }

      if(length(subset) == 0)
      {
        regression.fit <- stats::lm(y ~ 1)  # Intercept Regression Model
      } else {

        regression.fit <- stats::lm(y ~ x.sub)  # Regression Model
      }

      SSE <- regression.fit$residuals %*% regression.fit$residuals


      regression.fit0 <- stats::lm(y ~ 1)  # Intercept Regression Model
      SSE0            <- stats::anova(regression.fit0)[["Sum Sq"]]
      # Sum of squared errors for Intercept model


      logB0 <- (-(n - 1)/2) * log(1 + g * (SSE/SSE0)) +
        ((n - k - 1)/2) * log(1 + g)

      result <- exp(logB0)
    }

    return(as.numeric(result))  # Highest Posterior model

  }


  marg.like <- marginal.likelihood(current.gamma)


  current.marg.like <- marg.like


  model.matrix <- cbind(99, current.marg.like)
  # This will contain the model and the corresponding marginal likelihood


  trace.matrix[1, 1] <- current.step
  trace.matrix[1, 2] <- current.T
  trace.matrix[1, 3] <- current.marg.like
  trace.matrix[1, 4] <- paste0(current.model, collapse = "-")

  for(i in 1:nburnin)  # burn-in nburnin
  {

    old.marg.like <- current.marg.like
    set.new.gamma <- matrix(current.gamma, nrow = k, ncol = k)
    # This will contain full set of next gamma's


    for(j in 1:k)
    {
      set.new.gamma[j, j] <- 1 - current.gamma[j]
    }

    # Following will create the set \gamma_0 described in the Hans et al article
    # The other 2 sets are created in the previous loop.


    p <- length(current.model)

    if(p != 0)
    {
      model.minus <- utils::combn(current.model, (p - 1))
      full.model <- 1:k
      complement.current.model <- full.model[is.na(pmatch(full.model, current.model))]
      temp.model <- t(utils::combn(complement.current.model, 1))
      model0 <- NULL
      for(j in 1:ncol(model.minus))
      {
        for(i in 1:nrow(temp.model))
        {
          model0 <- cbind(model0, c(temp.model[i, ], model.minus[, j]))
        }
      }


      for(j in 1:ncol(model0))
      {
        temp.gamma <- rep(0, k)
        temp.gamma[model0[, j]] <- 1
        set.new.gamma <- cbind(set.new.gamma, temp.gamma)
      }

    }

    set.new.gamma <- cbind(set.new.gamma, current.gamma)
    # We add the current model as well

    proposed.model <- apply(set.new.gamma,
                            2, model)
    proposed.marginal.likelihood <- as.numeric(apply(set.new.gamma,
                                                     2, marginal.likelihood))



    model.matrix <- rbind(model.matrix,
                          cbind(proposed.model, proposed.marginal.likelihood))


    proposed.probability <- proposed.marginal.likelihood
    # We shall jump to the next model with these proposed probabilities
    new.gamma.index <- sample(1:ncol(set.new.gamma), size = 1,
                              prob = proposed.probability)

    new.gamma <- set.new.gamma[, new.gamma.index]
    new.model <- which(new.gamma == 1)
    new.model <- paste0(new.model, collapse = "-")


    marg.like <- marginal.likelihood(new.gamma)
    new.marg.like <- marg.like

    model.matrix <- rbind(model.matrix,
                          cbind(new.model, new.marg.like))


    if(new.marg.like >= current.marg.like)
    {
      current.gamma     <- new.gamma
      current.marg.like <- new.marg.like
    } else {
      if(stats::runif(1, 0, 1) <=
         exp((log(new.marg.like) - log(current.marg.like))/current.T))
      {
        current.gamma     <- new.gamma
        current.marg.like <- new.marg.like
      }

    }



    error             <- abs(current.marg.like - old.marg.like)
    current.model     <- which(current.gamma == 1)
    current.step      <- current.step + 1  # Number of steps
    # current.T        <- ((Tf - T0)/(nrep - 1)) * (current.step) + T0  # fixes the temperature
    current.T         <- current.T * C  # fixes the temperature
    # current.T         <- (k * (log(old.marg.like) -
    #                            log(current.marg.like)))/log(current.step)
    current.marg.like <- current.marg.like


    trace.matrix[current.step, 1] <- current.step
    trace.matrix[current.step, 2] <- current.T
    trace.matrix[current.step, 3] <- current.marg.like
    trace.matrix[current.step, 4] <- paste0(current.model, collapse = "-")


  }


  while((current.step < nstep) && (error > abstol))
  {

    old.marg.like <- current.marg.like
    set.new.gamma <- matrix(current.gamma, nrow = k, ncol = k)
    # This will contain full set of next gamma's


    for(j in 1:k)
    {
      set.new.gamma[j, j] <- 1 - current.gamma[j]
    }

    # Following will create the set \gamma_0 described in the Hans et al article
    # The other 2 sets are created in the previous loop.


    p <- length(current.model)

    if(p != 0)
    {
      model.minus <- utils::combn(current.model, (p - 1))
      full.model <- 1:k
      complement.current.model <- full.model[is.na(pmatch(full.model, current.model))]
      temp.model <- t(utils::combn(complement.current.model, 1))
      model0 <- NULL
      for(j in 1:ncol(model.minus))
      {
        for(i in 1:nrow(temp.model))
        {
          model0 <- cbind(model0, c(temp.model[i, ], model.minus[, j]))
        }
      }


      for(j in 1:ncol(model0))
      {
        temp.gamma <- rep(0, k)
        temp.gamma[model0[, j]] <- 1
        set.new.gamma <- cbind(set.new.gamma, temp.gamma)
      }

    }

    set.new.gamma <- cbind(set.new.gamma, current.gamma)
    # We add the current model as well

    proposed.model <- apply(set.new.gamma,
                            2, model)
    proposed.marginal.likelihood <- as.numeric(apply(set.new.gamma,
                                                     2, marginal.likelihood))



    model.matrix <- rbind(model.matrix,
                          cbind(proposed.model, proposed.marginal.likelihood))


    proposed.probability <- proposed.marginal.likelihood
    # We shall jump to the next model with these proposed probabilities
    new.gamma.index <- sample(1:ncol(set.new.gamma), size = 1,
                              prob = proposed.probability)

    new.gamma <- set.new.gamma[, new.gamma.index]
    new.model <- which(new.gamma == 1)
    new.model <- paste0(new.model, collapse = "-")


    marg.like <- marginal.likelihood(new.gamma)
    new.marg.like <- marg.like

    model.matrix <- rbind(model.matrix,
                          cbind(new.model, new.marg.like))


    if(new.marg.like >= current.marg.like)
    {
      current.gamma     <- new.gamma
      current.marg.like <- new.marg.like
    } else {
      if(stats::runif(1, 0, 1) <=
         exp((log(new.marg.like) - log(current.marg.like))/current.T))
      {
        current.gamma     <- new.gamma
        current.marg.like <- new.marg.like
      }

    }



    error             <- abs(current.marg.like - old.marg.like)
    current.model     <- which(current.gamma == 1)
    current.step      <- current.step + 1  # Number of steps
    # current.T        <- ((Tf - T0)/(nrep - 1)) * (current.step) + T0  # fixes the temperature
    current.T         <- current.T * C  # fixes the temperature
    # current.T         <- (k * (log(old.marg.like) -
    #                            log(current.marg.like)))/log(current.step)
    current.marg.like <- current.marg.like


    trace.matrix[current.step, 1] <- current.step
    trace.matrix[current.step, 2] <- current.T
    trace.matrix[current.step, 3] <- current.marg.like
    trace.matrix[current.step, 4] <- paste0(current.model, collapse = "-")

}  # end of while loop
  } else {


    marginal.likelihood <- function(current.gamma)
    { # Returns the log of marginal likelihood given in the article

      current.model <- model(current.gamma)

      if(sum(current.model == model.matrix[, 1]) >= 1)
      {

        result <- as.numeric(model.matrix[which(current.model == model.matrix[, 1]), 2])[1]

      } else {

        subset <- which(current.gamma == 1)
        x.sub <- x[, subset]

        if(is.vector(x.sub)) {
          x.sub <- matrix(x.sub, ncol=1)
          k <- 1} else {
            k <- ncol(x.sub)
          }

        if(length(subset) == 0)
        {
          regression.fit <- stats::lm(y ~ 1)  # Intercept Regression Model
        } else {

          regression.fit <- stats::lm(y ~ x.sub)  # Regression Model
        }

        SSE <- regression.fit$residuals %*% regression.fit$residuals


        regression.fit0 <- stats::lm(y ~ 1)  # Intercept Regression Model
        SSE0            <- stats::anova(regression.fit0)[["Sum Sq"]]
        # Sum of squared errors for Intercept model


        logB0 <- (-(n - 1)/2) * log(1 + g * (SSE/SSE0)) +
          ((n - k - 1)/2) * log(1 + g)

        result <- -exp(logB0)
      }

      return(as.numeric(result))  # Highest Posterior model

    }


    marg.like <- marginal.likelihood(current.gamma)


    current.marg.like <- marg.like


    model.matrix <- cbind(99, current.marg.like)
    # This will contain the model and the corresponding marginal likelihood


    trace.matrix[1, 1] <- current.step
    trace.matrix[1, 2] <- current.T
    trace.matrix[1, 3] <- current.marg.like
    trace.matrix[1, 4] <- paste0(current.model, collapse = "-")


  while((current.step < nstep) && (error > abstol))
  {

    old.marg.like <- current.marg.like
    set.new.gamma <- matrix(current.gamma, nrow = k, ncol = k)
    # This will contain full set of next gamma's


    for(j in 1:k)
    {
      set.new.gamma[j, j] <- 1 - current.gamma[j]
    }

    # Following will create the set \gamma_0 described in the Hans et al article
    # The other 2 sets are created in the previous loop.
    if(replace == TRUE)
    {
      if(current.step %% 5 == 0)
      {
        p <- length(current.model)

        if((p > 0) & (p < k))
        {
          model.minus <- utils::combn(current.model, (p - 1))
          full.model <- 1:k
          complement.current.model <- full.model[is.na(pmatch(full.model, current.model))]
          temp.model <- t(utils::combn(complement.current.model, 1))
          model0 <- NULL
          for(j in 1:ncol(model.minus))
          {
            for(i in 1:nrow(temp.model))
            {
              model0 <- cbind(model0, c(temp.model[i, ], model.minus[, j]))
            }
          }


          for(j in 1:ncol(model0))
          {
            temp.gamma <- rep(0, k)
            temp.gamma[model0[, j]] <- 1
            set.new.gamma <- cbind(set.new.gamma, temp.gamma)
          }

        }

      }

    }


    set.new.gamma <- cbind(set.new.gamma, current.gamma)
    # We add the current model as well

    proposed.model <- apply(set.new.gamma, 2, model)
    proposed.marginal.likelihood <- as.numeric(apply(set.new.gamma,
                                                     2, marginal.likelihood))



    model.matrix <- rbind(model.matrix,
                          cbind(proposed.model, proposed.marginal.likelihood))


    proposed.probability <- exp(-proposed.marginal.likelihood -
                                  max(-proposed.marginal.likelihood))
    # We shall jump to the next model with these proposed probabilities
    new.gamma.index <- sample(1:ncol(set.new.gamma), size = 1,
                              prob = proposed.probability)

    new.gamma <- set.new.gamma[, new.gamma.index]
    new.model <- which(new.gamma == 1)
    new.model <- paste0(new.model, collapse = "-")


    marg.like <- marginal.likelihood(new.gamma)
    new.marg.like <- marg.like
    #current.T         <- (k * (new.marg.like - old.marg.like))/log(current.step + 1)

    model.matrix <- rbind(model.matrix,
                          cbind(new.model, new.marg.like))


    if(new.marg.like <= current.marg.like)
    {
      current.gamma     <- new.gamma
      current.marg.like <- new.marg.like
    } else {
      temp.prob <- exp(-(new.marg.like - old.marg.like)/current.T)
      if(stats::rbinom(1, size = 1, temp.prob) == 1)
      {
        current.gamma     <- new.gamma
        current.marg.like <- new.marg.like
      }
    }




    error             <- abs(current.marg.like - old.marg.like)
    current.model     <- which(current.gamma == 1)
    current.step      <- current.step + 1  # Number of steps
    current.T         <- current.T * C  # fixes the temperature
    # current.T <- current.T/current.step
    current.marg.like <- current.marg.like


    trace.matrix[current.step, 1] <- current.step
    trace.matrix[current.step, 2] <- current.T
    trace.matrix[current.step, 3] <- current.marg.like
    trace.matrix[current.step, 4] <- paste0(current.model, collapse = "-")

  }

}
  return(list(final.model = current.model, history = trace.matrix[1:current.step, ]))

}






