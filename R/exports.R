#' Apply method to multivariate
#'
#' This function applies a function for univariate series ("ts") to a multivariate
#' ("mts").
#' @param x A multivariate time series of class "mts".
#' @param f A function that takes an univariate series as input.
#' @param ... Arguments for \code{f}.
apply_multivariate <- function(x, f, ...){
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
#' @param x A chain-linked series (with annual overlap) of class "ts".
#' @returns The pyp series of class "ts".
#' @examples
#' set.seed(23)
#' x <- sample(95:105, 24, TRUE) |>
#'   ts(start = 2005, frequency = 4)
#' y <- x |> chain_from_pyp(2006)
#' (x - pyp_from_chain(y)) |> round(5)
#' set.seed(23)
#' x <- sample(95:105, 24, TRUE) |>
#'   matrix(ncol = 2) |>
#'   ts(start = 2005, frequency = 4)
#' y <- x |> chain_from_pyp(2006)
#' (x - pyp_from_chain(y)) |> round(5)
#' @export
pyp_from_chain <- function(x) {
  if (methods::is(x, "mts")) {
    x_pyp <- apply_multivariate(x, pyp_from_chain)
  } else {
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
#' set.seed(23)
#' x <- sample(95:105, 12, TRUE) |>
#'   ts(start = 2005, frequency = 4) |>
#'   chain_from_pyp(2006)
#' aggregate(x, FUN = mean)
#' y <- change_ref_year(x, 2007)
#' aggregate(y, FUN = mean)
#' set.seed(23)
#' x <- sample(95:105, 24, TRUE) |>
#'   matrix(ncol = 2) |>
#'   ts(start = 2005, frequency = 4) |>
#'   chain_from_pyp(2006)
#' aggregate(x, FUN = mean)
#' y <- change_ref_year(x, 2007)
#' aggregate(y, FUN = mean)
#' @export
change_ref_year <- function(x, new_ref) {
  if (methods::is(x, "mts")) {
    y <- apply_multivariate(x, change_ref_year_uni, new_ref)
  } else {
    denom <- window(x, start = c(new_ref,1),
                    end = c(new_ref, frequency(x))) |>
      mean()
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
#'
#' @export
chain_from_pyp <- function(x, ref_year, x_a = NULL) {
  if (methods::is(x, "mts")) {
    y <- apply_multivariate(x, chain_from_pyp_uni, ref_year)
  } else {
    if (frequency(x) == 1) {
      x_chain <- cumprod(c(100,x/100)) |> ts(start = start(x)[1] - 1)
      value_ref_year <- c(window(x_chain, start = ref_year, end = ref_year))
      y <- x_chain / value_ref_year * 100
    } else {
      if (is.null(x_a)) {
        x_a <- aggregate.ts(x, FUN = mean)
      }
      x_chain_a <- chain_from_pyp(x_a, ref_year)
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

#' Get volume index from current and pyp prices
#'
#' Returns the series of pyp volume indices given current prices
#' and pyp prices.
#' @param current Current prices series of class "ts".
#' @param constant Constant prices (pyp) series of class "ts".
#' @returns List of time series ("ts"): the quarterly IPs and IQs (pyp).
#' @export
iq_pyp_from_money <- function(current, constant) {
  if (methods::is(current, "mts")) {
    if(ncol(current) != ncol(constant)) {
      stop("current and constant don't have the same number of columns!")
    }
    n <- ncol(current)
    y <- lapply(1:n, \(i) {
      z <- iq_pyp_from_money(current[, i], constant[, i])
      w <- ts(z, start = start(z), frequency = frequency(z))
      return(w)
    })
    result <- do.call(cbind, y)
    colnames(result) <- colnames(current)
  } else {
    current_aux <- current |>
      aggregate.ts(FUN = mean) |>
      rep(times = rep(frequency(current),
                      times = floor(length(current)/frequency(current)))) |>
      ts(start = start(current)[1] + 1, frequency = frequency(current))
    result <- constant / current_aux * 100
  }
  return(result)
}

#' Get chain-linked volume measure
#'
#' Returns chain-linked volume given chain-linked volume indices and current prices.
#' @param iq_chain Chain-linked volume indices.
#' @param current Current prices.
#' @param ref_year Reference year for the chain-linked indices.
#' @returns The chain-linked volume measures.
#' @export
chain_vol_from_iq <- function(iq_chain, current, ref_year){
  s <- frequency(iq_chain)
  if (methods::is(current, "mts")) {
    vol_chain <- sapply(1:ncol(iq_chain), \(i) {
      c(sum(window(current[,i], start = ref_year, end = c(ref_year, s))))*iq_chain[,i]/400
    }) |> ts(start = start(iq_chain), frequency = s)
    colnames(vol_chain) <- colnames(current)
  } else {
    vol_chain <- c(sum(window(current, start = ref_year, end = c(ref_year, s))))*iq_chain/400
  }
  return(vol_chain)
}

#' Get value index from current prices
#'
#' Returns the series of pyp value indices given current prices.
#' @param current Current prices series of class "ts".
#' @returns the series of IVs (pyp).
#' @export
iv_pyp_from_money <- function(current) {
  if (methods::is(current, "mts")) {
    n <- ncol(current)
    y <- lapply(1:n, \(i) {
      z <- iv_pyp_from_money(current[, i])
      w <- ts(z, start = start(z), frequency = frequency(z))
      return(w)
    })
    result <- do.call(cbind, y)
    colnames(result) <- colnames(current)
  } else {
    aux <- current |>
      aggregate.ts(FUN = mean) |>
      rep(times = rep(frequency(current),
                      times = floor(length(current)/frequency(current)))) |>
      ts(start = start(current)[1] + 1, frequency = frequency(current))
    result <- current / aux * 100
  }
  return(result)
}
