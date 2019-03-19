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
#' @importFrom magrittr %>%
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
  defaults <- list(dates = NULL,
                   w_dens = NULL,
                   f_dens = NULL,
                   dna = NULL,
                   ctd = NULL,
                   ctd_timed = NULL,
                   wards = NULL,
                   ids = NULL,
                   ctd_directed = FALSE,
                   w_unobs = NULL,
                   N = 0L,
                   L = 0L,
                   D = NULL,
                   max_range = NA,
                   can_be_ances = NULL,
                   log_w_dens = NULL,
                   log_w_unobs = NULL,
                   log_f_dens = NULL,
                   contacts = NULL,
                   C_combn = NULL,
                   C_nrow = NULL,
                   C_ind = NULL,
                   has_dna = logical(0),
                   move_t_onw = FALSE,
                   between_wards = FALSE,
                   N_ward = NULL,
                   p_ward = 1,
                   id_in_dna = integer(0))

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
    data$dates <- data$dates - min(data$dates)
    data$N <- length(data$dates)
    data$max_range <- diff(range(data$dates))
  }

  ## CHECK ID
  if(!is.null(data$ids)) {
    data$ids <- as.character(data$ids)
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
    if(data$w_dens[length(data$w_dens)] == 0) {
      final_index <- max(which(data$w_dens > 0))
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

  ## CHECK W_UNOBS
  if (!is.null(data$w_unobs)) {
    if (any(data$w_unobs<0)) {
      stop("w_unobs has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$w_unobs))) {
      stop("non-finite values detected in w_unobs")
    }

    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(data$w_unobs[length(data$w_unobs)] == 0) {
      final_index <- max(which(data$w_unobs > 0))
      data$w_unobs <- data$w_unobs[1:final_index]
      warning("Removed trailing zeroes found in w_unobs")
    }

    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if (length(data$w_unobs) < data$max_range) {
      length_to_add <- (data$max_range-length(data$w_unobs)) + 10 # +10 to be on the safe side
      val_to_add <- stats::dexp(seq_len(length_to_add), 1)
      val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
      data$w_unobs <- c(data$w_unobs, val_to_add)
    }

    ## standardize the mass function
    data$w_unobs <- data$w_unobs / sum(data$w_unobs)
    data$log_w_unobs <- matrix(log(data$w_unobs), nrow = 1)
  } else {
    data$w_unobs <- data$w_dens
    data$log_w_unobs <- data$log_w_dens
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

  ## Add temporal ordering constraints using Serial Interval
  if(!is.null(data$dates)) {
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

    data$id_in_dna <- match(as.character(data$ids), rownames(data$dna))
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
    
    if(! ncol(data$ctd) %in% c(2, 3)) {
      stop("ctd must have two columns")
    }

    ## Convert to character to prevent factors from interfering
    data$ctd[,1] <- as.character(data$ctd[,1])
    data$ctd[,2] <- as.character(data$ctd[,2])

    ## Ensure all cases found in linelist
    unq <- unique(unlist(data$ctd[,1:2]))
    not_found <- unq[!unq %in% data$ids]
    if (length(not_found) != 0) {
      stop(paste("Individual(s)", paste(not_found, collapse = ", "),
                 "are unknown cases (idx < 1 or > N")
           )
    }

    ## Determine the number of contact types
    if(ncol(data$ctd) == 3) {
      if(!inherits(data$ctd[,3], "factor")) {
      data$ctd[,3] <- factor(data$ctd[,3],
                             levels = unique(data$ctd[,3]))
      }
    } else {
      data$ctd$type <- factor('foo')
    }

    data$ctd_matrix <- list()

    ## Create a contact matrix for each contact type
    ## If contact data is not directed, fill the matrix in both directions
    for(i in seq_along(levels(data$ctd[,3]))) {
      to_keep <- data$ctd[,3] == levels(data$ctd[,3])[i]
      tmp <- subset(data$ctd, to_keep)
      mat <- matrix(0, nrow = data$N, ncol = data$N)
      for(j in 1:nrow(tmp)) {
        mtch_1 <- match(tmp[,1], data$ids)
        mtch_2 <- match(tmp[,2], data$ids)
        mat[cbind(mtch_2, mtch_1)] <- 1
        if(!data$ctd_directed) mat[cbind(mtch_1, mtch_2)] <- 1
      }
      data$ctd_matrix[[i]] <- mat
    }
    
    if(data$ctd_directed) {
      data$C_combn <- data$N*(data$N - 1)
    } else {
      data$C_combn <- data$N*(data$N - 1)/2
    }
    ## Count the number of each type of contact - these are ordered by their
    ## factor order, which we will use throughout the analysis
    data$C_nrow <- as.vector(table(data$ctd[,3]))

  } else {
    data$ctd_matrix <- list()
  }

  if (!is.null(data$ctd_timed)) {

    if (!inherits(data$ctd_timed, c("matrix", "data.frame"))) {
      stop("ctd_timed is not a matrix or data.frame")
    }
    if(!ncol(data$ctd_timed) %in% c(4, 5)) {
      stop(paste0("Timed contact data must have four columns (ID | Place",
                  " | Start date | End date), with an optional fifth column",
                  " specifying the type of contact"))
    }
    if(class(min_date) != class(data$ctd_timed[,3]) |
       class(min_date) != class(data$ctd_timed[,4])) {
      stop("Sampling dates and contact dates are not in the same format.")
    }
    if(inherits(data$ctd_timed[,1], "factor")) {
      data$ctd_timed[,1] <- as.character(data$ctd_timed[,1])
    }
    if(!inherits(data$ctd_timed[,1], c("numeric", "character", "integer"))) {
      stop("IDs in contact data must be numbers or characters")
    }

    if(any(!data$ctd_timed[,1] %in% data$ids)) {
      stop("IDs in contact data not found in case IDs.")
    }

    data$N_place <- length(unique(data$ctd_timed[,2]))

    ## Replace IDs with their numeric indices
    data$ctd_timed[,1] <- as.character(data$ctd_timed[,1])

    ## Set place names as characters to prevent accidental indicing with numbers
    data$ctd_timed[,2] <- as.character(data$ctd_timed[,2])
    
    ## Set min(data$date) as day = 0 and adjust contact dates accordingly
    if (inherits(data$ctd_timed[,3], "Date")) {
      data$ctd_timed[,3] <- data$ctd_timed[,3] - min_date
      data$ctd_timed[,4] <- data$ctd_timed[,4] - min_date
    }
    if (inherits(data$ctd_timed, "POSIXct")) {
      data$ctd_timed[,3] <- difftime(data$ctd_timed[,3], min_date, units="days")
      data$ctd_timed[,4] <- difftime(data$ctd_timed[,4], min_date, units="days")
    }
    data$ctd_timed[,3] <- as.integer(round(data$ctd_timed[,3]))
    data$ctd_timed[,4] <- as.integer(round(data$ctd_timed[,4]))

    if(any(data$ctd_timed[,4] - data$ctd_timed[,3] < 0)) {
      stop("End dates must be after start dates.")
    }
    if(ncol(data$ctd_timed) == 5) {
      data$ctd_timed[,5] <- factor(data$ctd_timed[,5],
                                   levels = unique(data$ctd_timed[,5]))
    } else {
      data$ctd_timed$type <- factor('foo')
    }

    data$ctd_timed_matrix <- list()
    data$C_ind <- -min(data$ctd_timed[,3])
    
    ## Pass places as numbers rather than strings to speed up comparison
    place_order <- unique(data$ctd_timed[,2])

    for(i in seq_along(levels(data$ctd_timed[,5]))) {
      tmp <- subset(data$ctd_timed,
                    data$ctd_timed[,5] == levels(data$ctd_timed[,5])[i])
      mat <- matrix(0, nrow = data$N,
                    ncol = max(data$ctd_timed[,4]) - min(data$ctd_timed[,3]) + 1)
      for(j in 1:nrow(tmp)) {
        ind <- (tmp[j,3] + data$C_ind + 1):(tmp[j,4] + data$C_ind + 1)
        mat[match(tmp[j, 1], data$ids), ind] <- match(tmp[j, 2], place_order)
      }
      data$ctd_timed_matrix[[i]] <- mat
    }
    
    ## Define the time of contacts as occuring between the first admission date
    ## and last discharge date
    data$C_combn <- data$N*(data$N - 1)/2
    data$C_nrow <- c(data$C_nrow, rep(NA, length(data$ctd_timed_matrix)))
    
    ## UPDATE CAN_BE_ANCES
    
  } else {
    data$ctd_timed_matrix <- list()
  }

  if(!is.null(data$wards)) {

    if(ncol(data$wards) != 4) {
      stop(paste0("Ward data must have four columns (ID | Ward",
                  " | Admission date | Discharge date)"))
    }
    if(!is.null(data$ctd)) {
      stop("Cannot use ward data and contact data; remove one.")
    }
    if(class(min_date) != class(data$wards[,3]) |
       class(min_date) != class(data$wards[,4])) {
      stop("Sampling dates and ward dates are not in the same format.")
    }
    if(inherits(data$wards[,1], "factor")) {
      data$wards[,1] <- as.character(data$wards[,1])
    }
    if(!inherits(data$wards[,1], c("numeric", "character", "integer"))) {
      stop("IDs in ward data must be numbers or characters")
    }
    if(any(!data$wards[,1] %in% data$ids)) {
      stop("IDs in ward data not found in case IDs.")
    }
    if(data$p_ward != 1 && sum(data$p_ward) != 1) {
      stop("p_ward must sum to 1")
    }

    data$N_ward <- length(unique(data$wards$ward))

    if(is.null(data$p_ward)) {
      data$p_ward <- rep(1/data$N_ward, data$N_ward)
    }
    
    ## Replace IDs with their numeric indices
    data$wards[,1] <- as.character(data$wards[,1])

    ## Set ward names as characters to prevent accidental indicing with numbers
    data$wards[,2] <- as.character(data$wards[,2])
    
    ## Set min(data$date) as day = 0 and adjust contact dates accordingly
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

    if(any(data$wards[,4] - data$wards[,3] < 0)) {
      stop("Ward discharge dates must be after admission dates.")
    }

    ward_matrix <- matrix(0, nrow = data$N,
                          ncol = max(data$wards[,4]) - min(data$wards[,3]) + 1)
    data$C_ind <- -min(data$wards$adm)
    ward_order <- unique(data$wards[,2])

    for(i in 1:nrow(data$wards)) {
      ind <- (data$wards[i,3] + data$C_ind + 1):(data$wards[i,4] + data$C_ind + 1)
      ward_matrix[match(data$wards[i, 1], data$ids), ind] <- match(data$wards[i, 2], ward_order)
    }

    ## Calculate the number of contacts from ward data
    n_contacts <- ward_matrix %>%
      apply(2, function(i) table(i[i != 0])) %>%
      lapply(function(i) i*(i - 1)/2) %>%
      lapply(sum) %>%
      unlist
    
    ## Define the time of contacts as occuring between the first admission date
    ## and last discharge date
    data$C_nrow <- sum(n_contacts)
    data$C_combn <- (data$N*(data$N - 1)/2)*ncol(ward_matrix)
    data$ward_matrix <- ward_matrix
    data$ward_ncol <- ncol(ward_matrix)
    
  } else {
    data$ward_matrix <- matrix(0, nrow = 0, ncol = 0)
  }
  
  ## output is a list of checked data
  return(data)

}
