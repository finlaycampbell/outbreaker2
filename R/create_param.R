#' Initializes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param data A list of data items as returned by \code{outbreaker_data}, or
#' arguments passed to this function.
#'
#' @param config A list of settings as returned by \code{create_config}, or
#' arguments passed to this function.
#'
#' @export
#'
#' @aliases outbreaker_store
#'
#' @aliases outbreaker_param
#'
#' @return
#'
#' A list containing two components \code{$store} and
#' \code{$current}. \code{store} is a list with the class
#' \code{outbreaker_store}, used for storing 'saved' states of the
#' MCMC. \code{current} is a list with the class \code{outbreaker_param}, used
#' for storing 'current' states of the MCMC. \cr \cr
#'
#' \code{outbreaker_store} class content:
#' \itemize{
#'
#'  \item \code{size}: The length of the list, corresponding to the number of
#' samples saved from the MCMC.
#'
#'  \item \code{step}: A vector of integers of length \code{size}, storing the
#' steps of the MCMC corresponding to the saved samples.
#'
#'  \item \code{post}: A numeric vector of length \code{size}, storing
#' log-posterior values.
#'
#'  \item \code{like}: A numeric vector of length \code{size}, storing
#' log-likelihood values.
#'
#'  \item \code{prior}: A numeric vector of length \code{size},
#' storing log-prior values.
#'
#'  \item \code{alpha}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing indices (from 1 to N) of
#' infectors for each case.
#'
#'  \item \code{t_inf}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing dates of infections for
#' each case.
#'
#'  \item \code{mu}: A numeric vector of length \code{size}, storing values of
#' the mutation rate.
#'
#'  \item \code{kappa}: A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing the number of generations
#' before the last sampled ancestor for each case.
#'
#'  \item \code{pi}: A numeric vector of length \code{size}, storing values of
#' the reporting probability.
#'
#'  \item \code{eps}: A numeric vector of length \code{size}, storing values of
#' the contact reporting coverage.
#'
#'  \item \code{eta}: A numeric vector of length \code{size}, storing values of
#' the contact sensitivity.
#'
#'  \item \code{lambda}: A numeric vector of length \code{size}, storing values of
#' the non-infectious contact rate.
#'
#'  \item \code{counter}: A counter used to keep track of the current iteration
#' of the MCMC (used internally).
#'
#' }
#'
#'
#' \code{outbreaker_param} class content:
#' \itemize{
#'
#'  \item \code{alpha}: An integer vector of length \code{data$N}, storing
#' indices (from 1 to N) of infectors for each case.
#'
#'  \item \code{t_inf}: An integer vector of length \code{data$N}, storing dates
#' of infections for each case.
#'
#'  \item \code{mu}: The value of the mutation rate.
#'
#'  \item \code{kappa}: An integer vector of length \code{data$N}, storing the
#' number of generations before the last sampled ancestor for each case.
#'
#'  \item \code{pi}: The value of the reporting probability.
#'
#'  \item \code{eps}: The value of the contact reporting coverage.
#'
#'  \item \code{eta}: The value of the contact sensitivity.
#' 
#'  \item \code{lambda}: The value of the non-infectious contact rate.
#'
#' }
#'
#' @examples
#'
#' ## load data
#' x <- fake_outbreak
#' data <- outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#'
#' ## modify config settings
#' config <- create_config(move_alpha = FALSE, n_iter = 2e5, sample_every = 1000)
#'
#' ## create param object
#' param <- create_param(data = data, config = config)
#'
create_param <- function(data = outbreaker_data(),
                         config = create_config()) {
  ## CREATE EMPTY OUTPUT VECTORS ##
  size <- round(config$n_iter/config$sample_every)
  step <- integer(size)
  post <- prior <- like <- mu <- pi <- tau <- double(size)
  eps <- as.list(integer(size))
  eta <- as.list(integer(size))
  lambda <- as.list(integer(size))
  alpha <- as.list(integer(size))
  t_inf <- as.list(integer(size))
  t_onw <- as.list(integer(size))
  ward <- as.list(integer(size))
  kappa <- as.list(integer(size))
  
  ## SET CURRENT VALUES AND COUNTER ##
  step[1] <- 1L
  current_mu <- mu[1] <- config$init_mu
  current_alpha <- alpha[[1]] <- config$init_alpha
  current_kappa <- kappa[[1]] <- config$init_kappa
  current_pi <- pi[1] <- config$init_pi
  current_tau <- tau[1] <- config$init_tau
  current_eps <- eps[[1]] <- config$init_eps
  current_eta <- eta[[1]] <- config$init_eta
  current_lambda <- lambda[[1]] <- config$init_lambda
  if (is.null(config$init_t_inf)) {
    current_t_inf <- t_inf[[1]] <- data$dates - which.max(data$f_dens) + 1L
  } else {
    current_t_inf <- t_inf[[1]] <- config$init_t_inf
  }
  if (is.null(config$init_t_onw)) {
    tmp_t_onw <- round(config$init_t_inf - sum(data$w_dens*seq_along(data$w_dens))/2)
    tmp_t_onw[config$init_kappa == 1] <- -1000
    current_t_onw <- t_onw[[1]] <- as.integer(tmp_t_onw)
  } else {
    tmp_t_onw <- config$init_t_onw
    tmp_t_onw[config$init_kappa == 1] <- -1000
    current_t_onw <- t_onw[[1]] <- as.integer(tmp_t_onw)
  }

  ## Make sure kappa = 1 matches with ward = 0
  if (is.null(config$init_ward)) {
    tmp_ward <- ifelse(config$init_kappa > 1, 1, 0)
    tmp_ward[is.na(config$init_alpha)] <- 0
    current_ward <- ward[[1]] <- as.integer(tmp_ward)
  } else {
    tmp_ward <- config$init_ward
    tmp_ward[tmp_ward == 0 & config$init_kappa > 1] <- 1
    tmp_ward[tmp_ward != 0 & config$init_kappa == 1] <- 0
    current_ward <- ward[[1]] <- as.integer(tmp_ward)
  }

  ## Calculate initial ward transition probabilities
  if (!is.null(data$wards)) {
    ward_mat <- get_ward_p(data$p_ward, config$init_eps, config$init_tau, config$max_kappa)
    ward_mat_1 <- get_ward_p(data$p_ward, config$init_eps, 1, 1)
  } else {
    ward_mat <- NULL
    ward_mat_1 <- NULL
  }
  
  counter <- 1L

  store <- list(
    size = size, step = step,
    post = post, like = like, prior = prior,
    alpha = alpha, t_inf = t_inf, t_onw = t_onw,
    mu = mu, kappa = kappa, pi = pi, tau = tau,
    eps = eps, eta = eta, lambda = lambda, ward = ward,
    counter = counter
  )
  class(store) <- c("outbreaker_store", "list")

  current  <- list(alpha = current_alpha, t_inf = current_t_inf,
                   t_onw = current_t_onw, mu = current_mu,
                   kappa = current_kappa, pi = current_pi, tau = current_tau,
                   eps = current_eps, eta = current_eta,
                   lambda = current_lambda, ward = current_ward,
                   ward_mat = ward_mat, ward_mat_1 = ward_mat_1)
  class(current) <- c("outbreaker_param", "list")

  tmp <- cpp_swap_cases(current, 1, FALSE)
  
  ## Count the number of contacts
  if(!is.null(data$ctd_timed)) {
    current$n_contacts <- local_n_contacts(data, current, seq.int(1, data$N))
  } else {
    current$n_contacts <- matrix(, ncol = 0, nrow = 0)
  }
  ## SHAPE CHAIN ##
  out <- list(store = store,
              current = current)
  return(out)
}

