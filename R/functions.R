#' Aggregate, repeat and lag
#'
#' Helper function to repeat the aggregate annual value of a series on
#' each period, and possibly lag it.
#' @param x (ts) Any time series
#' @param fun (function) Aggregation function, mean by default
#' @param k (int) Units to lag.
#' @details
#' Applies \code{aggregate.ts} to the series to get the annual values
#' and then repeats those values for every subyear period.
#'
#' The \code{k} parameter is passed to \code{stats::lag}.
#' @returns description
#' @examples
#' aggr_rep_lag(gdp_volume) |> plot()
#' @export
aggr_rep_lag <- function(x, fun = mean, k = 0) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, aggr_rep_lag, fun, k)
  } else {
    f <- frequency(x)
    ratio <- floor(length(x)/f)
    aux <- aggregate.ts(x, FUN = fun) |> rep(times = rep(f, ratio))
    y <- ts(aux, start = start(x), frequency = f) |> stats::lag(k)
  }
  return(y)
}

#' Apply method to multivariate
#'
#' This function applies a function for univariate series ("ts") to a
#' multivariate series ("mts").
#' @param x (mts) A multivariate time series.
#' @param f (function) A function that takes an univariate series as input.
#' @param ... Arguments for \code{f}.
apply_to_columns <- function(x, f, ...){
  n <- ncol(x)
  y <- lapply(1:n, \(i) {
    z <- f(x[,i], ...)
    w <- ts(z, start = start(z), frequency = frequency(z))
    return(w)
  })
  result <- do.call(cbind, y)
  colnames(result) <- colnames(x)
  return(result)
}

#' Get pyp indices
#'
#' Computes the pyp index series from a chain-linked series.
#'
#' @param x (ts) Chain-linked series with annual overlap.
#' @param x_a (ts) Annual chain-linked series. If not given, it's computed by
#' taking the average of each year.
#' @returns The pyp series.
#' @details
#' The time series should start at (y,1) where y is any year.
#'
#' @examples
#' get_pyp(gdp_volume)
#'
#' @export
get_pyp <- function(x, x_a = NULL) {
  f <- stats::frequency(x)
  ratio <- floor(length(x)/f)
  if (methods::is(x, "mts")) {
    x_pyp <- apply_to_columns(x, get_pyp)
  } else {
    if (is.null(x_a)) {
      aux <- aggr_rep_lag(x) |> stats::lag(-f)
    } else {
      aux <- rep(x_a, times = rep(f, ratio)) |>
        stats::ts(start = stats::start(x), frequency = f) |>
        stats::lag(-f)
    }
    x_pyp <- x/aux*100
  }
  return(x_pyp)
}

#' Change reference year
#'
#' Changes the reference year of a chain-linked series (with annual
#' overlap).
#' @param x (ts) A chain-linked series with annual overlap.
#' @param new_ref_year (num) New reference year.
#' @returns The re-referenced index series.
#' @examples
#' change_ref_year(gdp_volume, 2015)
#' plot(gdp_volume)
#' lines(change_ref_year(gdp_volume, 2015))
#'
#' @export
change_ref_year <- function(x, new_ref_year) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, change_ref_year, new_ref_year)
  } else {
    denom <- stats::window(x, start = c(new_ref_year,1),
                    end = c(new_ref_year, stats::frequency(x))) |> mean()
    y <- x / denom * 100
  }
  return(y)
}

#' Get chain-linked indices
#'
#' Computes chain-linked index series from a pyp series.
#'
#' @param x (ts) A pyp series.
#' @param ref_year (num) Reference year for the chain-linked series.
#' @param x_a (ts) Annual pyp series. If not given, it's computed by taking the
#' average of each year.
#' @returns The chain-linked series.
#' @details
#' The chain-linked series x_chain is computed with the annual overlap method.
#' Suppose the x series runs from (y0, p0 = 0) to (y1, p1), where pi is a subyear
#' period. Then the chain-linked series at (y2, p2) is given by the cumulative
#' product of the annual series from y0 to y2-1 times x at (y2, p2).
#' @examples
#' gdp_pyp <- get_pyp(gdp_volume)
#' get_chain_linked(gdp_pyp, 2015)
#'
#' @export
get_chain_linked <- function(x, ref_year, x_a = NULL) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, get_chain_linked, ref_year)
  } else {
    if (frequency(x) == 1) {
      x_chain <- cumprod(c(100,x/100)) |> ts(start = start(x)[1] - 1)
      value_ref_year <- c(window(x_chain, start = ref_year, end = ref_year))
      y <- x_chain / value_ref_year * 100
    } else {
      if (is.null(x_a)) {
        x_a <- aggregate.ts(x, FUN = mean)
      }
      x_chain_a <- get_chain_linked(x_a, ref_year)
      s <- frequency(x)
      x_chain_a_aux <- x_chain_a |> rep(times = rep(s,length(x_chain_a))) |>
        ts(start = start(x_chain_a)[1] + 1, frequency = s)
      x_chain <- x * x_chain_a_aux/ 100
      value_ref_year <- c(window(x_chain_a, start = ref_year, end = ref_year))
      y <- x_chain / value_ref_year * 100
    }
  }
  return(y)
}

#' Get quantity index
#'
#' Returns the series of quantity indices in previous year prices from
#' a current prices and
#' @param current (ts) Values at current prices.
#' @param constant (ts) Values at previous year prices.
#' @returns Series of quantity indices for previous year prices.
#' @examples
#' gdp_pyp <- get_pyp(gdp_volume)
#' gdp_constant <- gdp_current / gdp_pyp * 100
#' get_q_index(gdp_current, gdp_constant)
#' @export
get_q_index <- function(current, constant) {
  if (methods::is(current, "mts") & methods::is(constant, "mts")) {
    n <- ncol(current)
    y <- lapply(1:n, \(i) {
      z <- get_q_index(current[, i], constant[, i])
      w <- ts(z, start = start(z), frequency = frequency(z))
      return(w)
    })
    result <- do.call(cbind, y)
    colnames(result) <- colnames(current)
  } else if (!methods::is(current, "mts") & !methods::is(constant, "mts")) {
    aux <- aggr_rep_lag(current) |> stats::lag(-frequency(current))
    result <- constant / aux * 100
  }
  return(result)
}

#' Get value index
#'
#' Returns the (not chain-linked) series of value indices from a series
#' of current prices.
#' @param current (ts) Series of current prices series.
#' @returns (ts) Series of value indices.
#' @details
#' The value of the resulting series x at (y,s), where y is the year and s
#' is the subyear period, is current(y,s)/current(y)
#'
#' @examples
#' get_v_index(gdp_current)
#'
#' @export
get_v_index <- function(current) {
  if (methods::is(current, "mts")) {
    y <- apply_to_columns(current, get_v_index)
  } else {
    aux <- aggr_rep_lag(current) |> stats::lag(-frequency(current))
    y <- current / aux * 100
  }
  return(y)
}
