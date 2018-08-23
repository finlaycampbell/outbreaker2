## #' Shape MCMC samples into outputs for outbreaker
## #'
## #' This function shapes MCMC samples for outbreaker into a data.frame with the class \code{outbreaker_chains}.
## #'
## #' @author Thibaut Jombart (\email{thibautjombart@@gmail_com})
## #'
## #' @param param a list of output items as returned by \code{create_param}
## #'
## #' @param data a list of data items as returned by \code{outbreaker_data}
## #'
## #' @export
## #'
outbreaker_mcmc_shape <- function(param, data) {
  ## unfold ancestries ##
  if (!all(vapply(param$alpha, length, integer(1))==data$N)) {
    stop("some ancestries are missing in the param")
  }
  param$alpha <- matrix(unlist(param$alpha), ncol = data$N, byrow = TRUE)
  colnames(param$alpha) <- paste("alpha", seq_len(data$N), sep=".")

  ## unfold infection dates ##
  if (!all(vapply(param$t_inf, length, integer(1))==data$N)) {
    stop("some infection dates are missing in the param")
  }
  param$t_inf <- matrix(unlist(param$t_inf), ncol = data$N, byrow = TRUE)
  colnames(param$t_inf) <- paste("t_inf", seq_len(data$N), sep=".")

  ## unfold t_onw dates ##
  if(data$move_t_onw) {
    if (!all(vapply(param$t_onw, length, integer(1))==data$N)) {
      stop("some onward infection dates are missing in the param")
    }
    t_onw <- matrix(unlist(param$t_onw), ncol = data$N, byrow = TRUE)
    t_onw[t_onw == -1000] <- NA
    colnames(t_onw) <- paste("t_onw", seq_len(data$N), sep=".")
  }

  if(data$move_t_onw) {
    pi2 <- param$pi2
  }
  
  ## unfold number of generations ##
  if (!all(vapply(param$kappa, length, integer(1))==data$N)) {
    stop("some ancestries are missing in the param")
  }
  param$kappa <- matrix(unlist(param$kappa), ncol = data$N, byrow = TRUE)
  colnames(param$kappa) <- paste("kappa", seq_len(data$N), sep=".")

  ## shape data.frame and convert ##
  param <- data.frame(step = param$step,
                      post = param$post, like = param$like, prior = param$prior,
                      mu = param$mu, pi = param$pi, eps = param$eps,
                      lambda = param$lambda, param$alpha, param$t_inf,
                      param$kappa)
  if(data$move_t_onw) param$pi2 <- pi2
  if(data$move_t_onw) param <- cbind(param, t_onw)
  names(param) <- gsub("[.]", "_", names(param))

  ## output is a data.frame containing all parameters and augmented data, with a dedicated
  ## class (for summary, plotting, etc.)
  class(param) <- c("outbreaker_chains","data.frame")
  return(param)
}
