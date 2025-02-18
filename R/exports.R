#' Aggregate and repeat
#'
#' Helper function to repeat the aggregate annual value of a series on
#' each period.
#' @param x (ts) Any time series
#' @param fun (function) Aggregation function, mean by default
#' @details
#' Applies \code{aggregate.ts} to the series to get the annual values
#' and then repeats those values for every subyear period.
#' @returns description
#' @examples
#' rep_year_value(gdp_es_index) |> plot()
#' @export
aggr_and_rep <- function(x, fun = mean) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, aggr_and_rep)
  } else {
    f <- frequency(x)
    ratio <- floor(length(x)/f)
    aux <- aggregate.ts(x, FUN = fun) |> rep(times = rep(f, ratio))
    y <- ts(aux, start = start(x), frequency = f)
  }
  return(y)
}

#' Apply method to multivariate
#'
#' This function applies a function for univariate series ("ts") to a
#' multivariate series ("mts").
#' @param x A multivariate time series of class "mts".
#' @param f A function that takes an univariate series as input.
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
#' @param x_a (ts) Annual chain-linked series.
#' @returns The pyp series of class "ts".
#' @examples
#' get_pyp(gdp_es_volume)
#'
#' @export
get_pyp <- function(x, x_a = NULL) {
  if (methods::is(x, "mts")) {
    x_pyp <- apply_to_columns(x, get_pyp)
  } else {
    if (!("ts" %in% class(x))) {
      stop("<x> is not 'ts'.")
    }
    if ((frequency(x) != 1) & (start(x)[2] != 1)) {
      stop("<x> must start at the beginning of the year.")
    }
    if (is.null(x_a)) {
      aux <- aux <- aggr_and_rep(x)
    } else {
      f <- frequency(x)
      ratio <- floor(length(x)/f)
      aux <- rep(x_a, times = rep(f, ratio)) |>
        ts(aux, start = start(x), frequency = f)
    }
    x_pyp <- x/aux*100
  }
  return(x_pyp)
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
#' change_ref_year(gdp_es_volume, 2015)
#' plot(gdp_es_volume)
#' lines(change_ref_year(gdp_es_volume, 2015))
#'
#' @export
change_ref_year <- function(x, new_ref) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, change_ref_year_uni, new_ref)
  } else {
    denom <- window(x, start = c(new_ref,1),
                    end = c(new_ref, frequency(x))) |> mean()
    y <- x / denom * 100
  }
  return(y)
}

#' Get chain-linked indices
#'
#' Computes chain-linked index series from a pyp series.
#'
#' @param x A pyp series of class "ts".
#' @param ref_year Reference year ("num") for the chain-linked series.
#' @param x_a Annual pyp series. If not given, it is assumed that it's computed
#' like \code{x_a = aggregate.ts(x, FUN = mean)}.
#' @param normalize Make reference year = 100.
#' @returns The chain-linked series of class "ts".
#' @details
#' The chain-linked series x_chain is computed with the annual overlap method.
#' Suppose the x series runs from (y0, p0 = 0) to (y1, p1), where pi is a subyear
#' period. Then the chain-linked series at (y2, p2) is given by the cumulative
#' product of the annual series from y0 to y2-1 times x at (y2, p2).
#' @examples
#' gdp_es_pyp <- get_pyp(gdp_es_volume)
#' get_chain_linked(gdp_es_pyp)
#'
#' @export
get_chain_linked <- function(x, ref_year, x_a = NULL) {
  if (methods::is(x, "mts")) {
    y <- apply_to_columns(x, chain_from_pyp_uni, ref_year)
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
#' @param current Current prices series of class "ts".
#' @param constant Constant prices (pyp) series of class "ts".
#' @returns List of time series ("ts"): the quarterly IPs and IQs (pyp).
#' @examples
#' gdp_es_pyp <- get_pyp(gdp_es_volume)
#' gdp_es_constant <- gdp_es_current / gdp_es_constant * 100
#' get_q_index(gdp_es_current, gdp_es_constant)
#'
#' @export
get_q_index <- function(current, constant) {
  if (methods::is(current, "mts") & methods::is(constant, "mts")) {
    if(ncol(current) != ncol(constant)) {
      stop("current and constant don't have the same number of columns!")
    }
    n <- ncol(current)
    y <- lapply(1:n, \(i) {
      z <- get_q_index(current[, i], constant[, i])
      w <- ts(z, start = start(z), frequency = frequency(z))
      return(w)
    })
    result <- do.call(cbind, y)
    colnames(result) <- colnames(current)
  } else if (!methods::is(current, "mts") & !methods::is(constant, "mts")) {
    aux <- aggr_and_rep(current) |> stats::lag(-frequency(current))
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
#' get_v_index(gdp_es_current)
#'
#' @export
get_v_index <- function(current) {
  if (methods::is(current, "mts")) {
    y <- apply_multivariate(current, get_v_index)
  } else {
    aux <- aggr_and_rep(current) |> stats::lag(-frequency(current))
    y <- current / aux * 100
  }
  return(y)
}
