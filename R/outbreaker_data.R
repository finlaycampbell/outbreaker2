#' Process input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.  It
#' takes a list of named items as input, performs various checks, set defaults
#' where arguments are missing, and return a correct list of data input. If no
#' input is given, it returns the default settings.
#'
#' Acceptables arguments for ... are:
#' \describe{
#' \item{dates}{dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date.}
#'
#' \item{dna}{the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.}
#'
#' \item{ctd}{the contact tracing data provided as a matrix or dataframe of two
#' columns, indicating a reported contact between the two individuals whose ids
#' are provided in a given row of the data.}
#'
#' \item{w_dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t = 1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f_dens}{similar to \code{w_dens}, except that this is the distribution
#' of the colonization time, i_e. time interval during which the pathogen can
#' be sampled from the patient.}
#'
#'}
#'
#' @param ... a list of data items to be processed (see description)
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker_data}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @export
#'
#' @examples
#'
#' x <- fake_outbreak
#' outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#'
outbreaker_data <- function(..., data = list(...)) {

  ## SET DEFAULTS ##
  defaults <- list(dates = NULL, w_dens = NULL, f_dens = NULL, dna = NULL,
                   ctd = NULL, wards = NULL, ids = NULL, N = 0L, L = 0L,
                   D = NULL, max_range = NA, can_be_ances = NULL,
                   log_w_dens = NULL, log_f_dens = NULL, contacts = NULL,
                   C_combn = NULL, C_nrow = NULL, contacts_timed = NULL,
                   C_ind = NULL, has_dna = logical(0), id_in_dna = integer(0))

  ## MODIFY DATA WITH ARGUMENTS ##
  data <- modify_defaults(defaults, data, FALSE)


  ## CHECK DATA ##
  ## CHECK DATES
  if (!is.null(data$dates)) {
    min_date <- min(data$dates)
    if (inherits(data$dates, "Date")) {
      data$dates <- data$dates - min_date
    }
    if (inherits(data$dates, "POSIXct")) {
      data$dates <- difftime(data$dates, min_date, units="days")
    }
    data$dates <- as.integer(round(data$dates))
    data$N <- length(data$dates)
    data$max_range <- diff(range(data$dates))
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    ## Calculate the serial interval from w_dens and f_dens
    .get_SI <- function(w_dens, f_dens) {
      wf <- stats::convolve(w_dens, rev(f_dens), type = 'open')
      conv <- stats::convolve(rev(f_dens), rev(wf), type = 'open')
      lf <- length(f_dens)
      lw <- length(w_dens)
      return(data.frame(x = (-lf + 2):(lw + lf - 1), d = conv))
    }
    ## Check if difference in sampling dates falls within serial interval
    ## This allows for i to infect j even if it sampled after (SI < 0)
    .can_be_ances <- function(date1, date2, SI) {
      tdiff <- date2 - date1
      out <- sapply(tdiff, function(i) if(i %in% SI$x) return(TRUE) else return(FALSE))
      return(out)
    }
    SI <- .get_SI(data$w_dens, data$f_dens)
    data$can_be_ances <- outer(data$dates,
                               data$dates,
                               FUN=.can_be_ances,
                               SI = SI) # strict < is needed as we impose w(0)=0
    diag(data$can_be_ances) <- FALSE
  }

  ## CHECK ID
  if(!is.null(data$ids)) {

  } else {
    data$ids <- seq_len(data$N)
  }

  ## CHECK W_DENS
  if (!is.null(data$w_dens)) {
    if (any(data$w_dens<0)) {
      stop("w_dens has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$w_dens))) {
      stop("non-finite values detected in w_dens")
    }


    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(data$w_dens[length(data$w_dens)] < 1e-15) {
      final_index <- max(which(data$w_dens > 1e-15))
      data$w_dens <- data$w_dens[1:final_index]
      warning("Removed trailing zeroes found in w_dens")
    }

    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if (length(data$w_dens) < data$max_range) {
      length_to_add <- (data$max_range-length(data$w_dens)) + 10 # +10 to be on the safe side
      val_to_add <- stats::dexp(seq_len(length_to_add), 1)
      val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
      data$w_dens <- c(data$w_dens, val_to_add)
    }

    ## standardize the mass function
    data$w_dens <- data$w_dens / sum(data$w_dens)
    data$log_w_dens <- matrix(log(data$w_dens), nrow = 1)
  }

  ## CHECK F_DENS
  if (!is.null(data$w_dens) && is.null(data$f_dens)) {
    data$f_dens <- data$w_dens
  }
  if (!is.null(data$f_dens)) {
    if (any(data$f_dens<0)) {
      stop("f_dens has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$f_dens))) {
      stop("non-finite values detected in f_dens")
    }

    data$f_dens <- data$f_dens / sum(data$f_dens)
    data$log_f_dens <- log(data$f_dens)
  }


  ## CHECK DNA

  if (!is.null(data$dna)) {
    if (!inherits(data$dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if (!is.matrix(data$dna)) data$dna <- as.matrix(data$dna)

    ## get matrix of distances

    data$L <- ncol(data$dna) #  (genome length)
    data$D <- as.matrix(ape::dist.dna(data$dna, model="N")) # distance matrix
    storage.mode(data$D) <- "integer" # essential for C/C++ interface

    ## get matching between sequences and cases

    if (is.null(rownames(data$dna))) {
      if (nrow(data$dna) != data$N) {
        msg <- sprintf(paste("numbers of sequences and cases differ (%d vs %d):",
                             "please label sequences"),
                       nrow(data$dna), data$N)
        stop(msg)
      }

      rownames(data$dna) <- rownames(data$D) <- colnames(data$D) <- seq_len(data$N)
    }

    data$id_in_dna <- match(as.character(seq_len(data$N)), rownames(data$dna))
    if(all(is.na(data$id_in_dna))) {
      stop("DNA sequence labels don't match case ids")
    }

  } else {
    data$L <- 0L
    data$D <- matrix(integer(0), ncol = 0, nrow = 0)
    data$id_in_dna <- rep(NA_integer_, data$N)
  }
  data$has_dna <- !is.na(data$id_in_dna)


  ## CHECK CTD
  if (!is.null(data$ctd)) {
    
    if (!inherits(data$ctd, c("matrix", "data.frame"))) {
      stop("ctd is not a matrix or data.frame")
    }
    
    if(ncol(data$ctd) == 2) {
      if (!is.matrix(data$ctd)) data$ctd <- as.matrix(data$ctd)
      not.found <- data$ctd[any(!data$ctd %in% 1:data$N)]
      if (length(not.found) != 0) {
        not.found <- sort(unique(not.found))
        stop(paste("Individual(s)", paste(not.found, collapse = ", "),
                   "are unknown cases (idx < 1 or > N")
             )
      }
      data$contacts <- matrix(0, data$N, data$N)
      data$contacts[data$ctd] <- data$contacts[data$ctd[,c(2, 1)]] <- 1
      data$C_combn <- data$N*(data$N - 1)/2
      data$C_nrow <- nrow(data$ctd)
    } else if(ncol(data$ctd) == 3) {
      
      if(class(min_date) != class(data$ctd[,3])) {
        stop("Sampling dates and contacts dates are not in the same format")
      }
      
      ## Set min(data$date) as day = 1 and adjust contact dates accordingly
      if (inherits(data$ctd[,3], "Date")) {
        data$ctd[,3] <- data$ctd[,3] - min_date
      }
      if (inherits(data$ctd[,3], "POSIXct")) {
        data$ctd[,3] <- difftime(data$ctd[,3], min_date, units="days")
      }
      data$ctd[,3] <- as.integer(round(data$ctd[,3]))

      ## Set up an array to cover all transmission pairs and times
      data$contacts_timed <- array(0, c(data$N, data$N, diff(range(data$ctd[,3])) + 1))
      data$C_ind <- -min(data$ctd$date)

      ## Assign a unique identifier to each combination of ID | ID | Date,
      ## defined as the coordinate of that point in a 3D matrix
      get.coor <- function(x) data$N*data$N*(x[3] + data$C_ind) + data$N*(x[1] - 1) + x[2] - 1
      
      ## Fill in reported contacts as 1 in the contact array
      ##tmp <- as.matrix(cbind(data$ctd[, 1:2], data$ctd[,3] + data$C_ind + 1))
      ##data$contacts_timed[tmp] <- data$contacts_timed[tmp[,c(2, 1, 3)]] <- 1

      data$contacts_timed <- c(apply(data$ctd, 1, get.coor),
                               apply(data$ctd[,c(2, 1, 3)], 1, get.coor))
      
      
      ## The total number of contacts is N*(N - 1)/2 (ie pairwise contacts)
      ## times the total number of timesteps in our analysis
      data$C_combn <- (data$N*(data$N - 1)/2)*(diff(range(data$ctd$date)) + 1)
      data$C_nrow <- nrow(data$ctd)

      ## UPDATE CAN_BE_ANCES
      
    } else {
      stop("ctd must have two or three columns")
    }
  } else {
    data$contacts <- data$contacts_timed <- matrix(integer(0), ncol = 0, nrow = 0)
  }

  if(!is.null(data$wards)) {
    
    if(!is.null(data$ctd)) {
      stop("Cannot use ward data and contact data; remove one.")
    }
    if(class(min_date) != class(data$wards[,3]) |
       class(min_date) != class(data$wards[,4])) {
      stop("Sampling dates and ward dates are not in the same format.")
    }
    if(!inherits(data$wards[,1], c("numeric", "character", "integer"))) {
      stop("IDs in ward data must be numbers of characters")
    }
    
    ## Set min(data$date) as day = 1 and adjust contact dates accordingly
    if (inherits(data$wards[,3], "Date")) {
      data$wards[,3] <- data$wards[,3] - min_date
      data$wards[,4] <- data$wards[,4] - min_date
    }
    if (inherits(data$wards, "POSIXct")) {
      data$wards[,3] <- difftime(data$wards[,3], min_date, units="days")
      data$wards[,4] <- difftime(data$wards[,4], min_date, units="days")
    }
    data$wards[,3] <- as.integer(round(data$wards[,3]))
    data$wards[,4] <- as.integer(round(data$wards[,4]))

    if(any(data$wards[,4] - data$wards[,3] < 0) {
      stop("Ward discharge dates must be after admission dates.")
    }

    if(any(!data$wards[,1] %in% data$ids)) {
      stop("IDs in ward data not found in case IDs.")
    }
    
    data$wards[,1] <- match(data$wards[,1], data$ids)
    data$wards[,2] <- as.character(data$wards[,2])
    
    ward_dat <- matrix(0, nrow = data$N,
                       ncol = max(data$wards[,4]) - min(data$wards[,3]) + 1)

    data$C_ind <- -min(data$wards$adm)
    ward_order <- unique(data$wards[,2])

    for(i in 1:nrow(data$wards)) {
      ind <- (data$wards[i,3] + data$C_ind + 1):(data$wards[i,4] + data$C_ind + 1)
      ward_dat[data$wards[i, 1], ind] <- match(data$wards[i, 2], ward_order)
    }
    
  } else {
    data$ward_dat <- matrix(integer(0), ncol = 0, nrow = 0)
  }
  
  ## output is a list of checked data
  return(data)

}
