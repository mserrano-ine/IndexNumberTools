#' Get pyp indices
#'
#' Computes the pyp index series from a chain-linked series.
#'
#' @param x A chain-linked series (with annual overlap) of class "ts".
#' @returns The pyp series of class "ts".
#' @export
pyp_from_chain <- function(x) {
  if (!("ts" %in% class(x))) {
    stop("<x> is not 'ts'.")
  }
  if ((frequency(x) != 1) & (start(x)[2] != 1)) {
    stop("<x> must start at the beginning of the year.")
  }
  x_a <- aggregate.ts(x, FUN = mean)
  x_a_lagged <- lag(x_a, k = -1)
  x_a_aux <- x_a_lagged[rep(1:length(x_a),
                                    times = rep(frequency(x), length(x_a)))] |>
    ts(start = start(x_a_lagged), frequency = frequency(x))
  x_pyp <- x/x_a_aux*100
  return(x_pyp)
}

#' Checks the reference period
#'
#' Takes the mean of the values in the reference period and checks if
#' it is equal to 100.
#'
#' @param x A chain-linked series of indices of class "ts".
#' @param ref_period Reference period to be tested.
#' @param tol Tolerance
#' @returns TRUE if the reference period is c
check_reference <- function(x, ref_period, tol = 0.01) {
  if (!is.null(ref_period)) {
    # check <old_ref> is correct.
    if (!methods::is(ref_period, "list")) {
      if (length(ref_period) == 1) {
        check <- mean(window(x, start = c(ref_period,1), end = c(ref_period, frequency(x))))
      } else {
        check <- window(x, start = ref_period, end = ref_period)
      }
    } else {
      if (length(ref_period) == 2) {
        check <- mean(window(x, start = ref_period[[1]], end = ref_period[[2]]))
      } else {
        check <- window(x, start = ref_period[[1]], end = ref_period[[1]])
      }
    }
    if (abs(check - 100) < tol) {
      warning("The average in the reference period doesn't equal 100.")
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}


#' Change reference year
#'
#' Changes the reference year of a chain-linked series (with annual
#' overlap).
#' @param x A chain-linked series (with annual overlap) of class "ts".
#' @param new_ref New reference year. Must be such that \code{start(x) <=
#' new_ref <= end(x)}.
#' @returns The re-referenced index series of class "ts".
#' @importFrom methods is
#' @importFrom stats frequency
#' @importFrom stats window
#' @examples
#' set.seed(23)
#' x <- sample(95:105, 12, TRUE) |>
#'   ts(start = 2005, frequency = 4) |>
#'   chain_from_pyp(2006)
#' aggregate(x, FUN = mean)
#' y <- change_ref_year(x, 2007)
#' aggregate(y, FUN = mean)
#' plot(x)
#' lines(y, col = "red")
#' @export
change_ref_year <- function(x, new_ref) {
  denom <- window(x, start = c(new_ref,1),
                  end = c(new_ref, frequency(x))) |>
    mean()
  y <- x / denom * 100
  return(y)
}

#' Concatenate two reference years
#'
#' Takes two overlapping sections of an index series, each with a
#' different reference year, and outputs a joined series with reference
#' that of the first section.
#' @param x1 First section of the series.
#' @param x2 Second section of the series.
#' @param tol Tolerance for checking the overlapping period.
#' @importFrom TimeSeriesTools concat
#' @export
#' @examples
#' x <- ts(sample(98:103, 20, TRUE), start = 2001, frequency = 4) |>
#'   chain_from_pyp(2001)
#' x1 <- window(x, end = 2004) |> change_ref_year(2002)
#' x2 <- window(x, start = 2003) |> change_ref_year(2004)
#' y <- concat_references(x1, x2) |> change_ref_year(2001)
#' plot(x)
#' lines(x1, col = "green")
#' lines(x2, col = "red")
#' points(y, col = "blue")
#'
concat_references <- function(x1, x2, tol = 1e-2) {
  if (!(end(x1)[1] > start(x2)[1])) {
    stop("Series don't overlap (enough). It must be a full year overlap.")
  }
  link <- window(x1, start = start(x2)[1],
                 end = c(start(x2)[1], frequency(x2))) |> mean()
  x2 <- change_ref_year(x2, start(x2)[1]) * link / 100
  if (any(abs(x1-x2) > tol)) {
    warning("The overlapping periods for the two sections
            don't exactly match.")
  }
  y <- concat(window(x1, end = c(start(x2)[1]-1, frequency(x1))), x2)
  return(y)
}

#' Get chain-linked indices
#'
#' Computes chain-linked index series from a pyp series.
#'
#' @param x A pyp series of class "ts".
#' @param ref_year Reference year ("num") for the chain-linked series.
#' @param normalize Make reference year = 100.
#' @returns The chain-linked series of class "ts".
#' @export
chain_from_pyp <- function(x, ref_year, normalize = T) {
  if (frequency(x) == 4) {
    result <- chain_from_pyp_q(x, ref_year, normalize)
  } else if (frequency(x) == 1) {
    result <- chain_from_pyp_annual(x, ref_year, normalize)
  }
  return(result)
}


#' Get chain-linked indices
#'
#' Computes the quarterly chain-linked index series from a pyp series.
#'
#' @param x A pyp series of class "ts".
#' @param ref_year Reference year ("num") for the chain-linked series.
#' @param normalize Make reference year = 100.
#' @returns The chain-linked series of class "ts".
#' @export
chain_from_pyp_q <- function(x, ref_year, normalize = T) {
  x_pyp_a <- aggregate.ts(x, FUN = mean)
  x_chain_a <- chain_from_pyp_annual(x_pyp_a, ref_year)
  x_chain_a_aux <- x_chain_a |> rep(times = rep(4,length(x_chain_a))) |>
    ts(start = start(x_chain_a)[1] + 1, frequency = 4)
  x_chain_q <- x * x_chain_a_aux/ 100
  if (normalize) {
    value_ref_year <- window(x_chain_a, start = ref_year, end = ref_year) |> c()
    x_chain_q <- x_chain_q / value_ref_year * 100
  }
  return(x_chain_q)
}

#' Get chain-linked indices
#'
#' Computes the quarterly chain-linked index series from a pyp series.
#'
#' @param x A pyp series of class "ts".
#' @param ref_year Reference year ("num") for the chain-linked series.
#' @param normalize Make reference year = 100.
#' @returns The chain-linked series of class "ts".
#' @export
chain_from_pyp_annual <- function(x, ref_year, normalize = T) {
  chain <- cumprod(c(100,x/100)) |> ts(start = start(x)[1] - 1)
  if (normalize) {
    value_ref_year <- window(chain, start = ref_year, end = ref_year) |> c()
    chain <- chain / value_ref_year * 100
  }
  return(chain)
}

#' Get index from current and constant prices
#'
#'
#' @param current Current prices series of class "ts".
#' @param constant Constant prices (pyp) series of class "ts".
#' @returns List of time series ("ts"): the quarterly IPs and IQs (pyp).
#' @export
pyp_from_money <- function(current, constant) {
  if (frequency(current) == 4) {
    curr_a <- aggregate.ts(current, FUN = sum)
    const_a <- aggregate.ts(constant, FUN = sum)
    curr_a_aux <- (curr_a/4) |> rep(times = rep(4,length(curr_a))) |>
      ts(start = start(curr_a)[1] + 1, frequency = 4)
    const_a_aux <- (const_a/4) |> rep(times = rep(4,length(const_a))) |>
      ts(start = start(const_a)[1] + 1, frequency = 4)
    iq_pyp <- constant / curr_a_aux * 100
    ip_pyp <- current / const_a_aux * 100
  } else if (frequency(current) == 1) {
    iq_pyp <- constant / lag(current, -1) * 100
    ip_pyp <- current / constant * 100
  }
  return(list(iq_pyp = iq_pyp, ip_pyp = ip_pyp))
}


