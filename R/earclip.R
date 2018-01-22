#' Ear clipping triangulation
#'
#' Ear clipping is a relatively cheap method for constrained triangulation
#' of polygons.
#'
#' Originally this code was published in the `rgl` package, and written
#' by Duncan Murdoch. The code from rgl 0.99.9 was copied on 2018-01-22.
#' @author Duncan Murdoch, in the rgl package
#' @examples
#' data("minimal_mesh", package = "silicate")
#' rbind_na <- function(x) {
#' ## the inner head-1 is to remove the final closing/start point from simple features
#' head(do.call(rbind, lapply(x, function(a) rbind(head(a, -1), NA))), -1)
#' }
#' strip_sf <- function(x) {
#'   lapply(unlist(lapply(unclass(x), unclass), recursive = FALSE), rbind_na)
#' }
#' x <- strip_sf(minimal_mesh$geom)
#' lapply(x, earclip_rgl)
#' #ncsf <- sf::read_sf(system.file("shape/nc.shp", package="sf"))
#' #x <- strip_sf(ncsf$geometry)
#' #lapply(x, earclip_rgl)
#' ## ## 20s in rgl, 14s with col/rowSums
#' #rbenchmark::benchmark(lapply(x[1:12], earclip_rgl), replications = 40)
#' @export
earclip_rgl <- function (x, y = NULL, z = NULL, random = TRUE, plot = FALSE,
          partial = NA)
{
  xyz <- xyz.coords(x, y, z)
  if (xyz$xlab == "Index" && is.null(z) && (is.null(ncol(x)) ||
                                            ncol(x) == 2L)) {
    x <- xyz$y
    y <- xyz$z
  }
  else {
    x <- xyz$x
    y <- xyz$y
    if (!diff(range(x)))
      x <- xyz$z
    else if (!diff(range(y)))
      y <- xyz$z
  }
  nesting <- #nestPolys_rgl(x, y)
  verts <- nesting$verts
  nextvert <- rep(NA, length(x))
  processInside <- function(v) {
    for (i in nesting$nesting[[v]]) processOutside(i)
  }
  processOutside <- function(fwd) {
    fwd1 <- verts[[fwd]]
    nextvert[fwd1] <<- c(fwd1[-1], fwd1[1])
    reversed <- nesting$nesting[[fwd]]
    for (rev in reversed) {
      rev1 <- rev(verts[[rev]])
      nextvert[rev1] <<- c(rev1[-1], rev1[1])
      processInside(rev)
      done <- FALSE
      pairs <- expand.grid(seq_along(fwd1), seq_along(rev1))
      if (random)
        pairs <- pairs[sample(nrow(pairs)), ]
      for (p in seq_len(nrow(pairs))) {
        i <- fwd1[pairs[p, 1]]
        j <- rev1[pairs[p, 2]]
        seg <- cbind(c(x[i], y[i]), c(x[j], y[j]))
        clear <- TRUE
        for (q in seq_along(verts)) {
          i1 <- verts[[q]]
          if (!length(i1))
            next
          i2 <- c(i1[-1], i1[1])
          for (v in seq_along(i1)) if (length(intersect(c(i1[v],
                                                          i2[v]), c(i, j))) == 0 && intersectSegSeg_rgl(seg,
                                                                                                    cbind(c(x[i1[v]], y[i1[v]]), c(x[i2[v]],
                                                                                                                                   y[i2[v]])))) {
            clear <- FALSE
            break
          }
          if (!clear)
            break
        }
        if (clear) {
          i <- pairs[p, 1]
          j <- pairs[p, 2]
          ind <- c(fwd1[seq_len(i)], rev1[j:length(rev1)],
                   rev1[seq_len(j)])
          if (i < length(fwd1))
            ind <- c(ind, fwd1[i:length(fwd1)])
          verts[[fwd]] <<- fwd1 <- ind
          verts[[rev]] <<- integer(0)
          done <- TRUE
          break
        }
      }
      if (!done)
        stop("Cannot simplify polygon")
    }
    ind <- verts[[fwd]]
    tri <- triangulateSimple_rgl(x[ind], y[ind], random, plot,
                             partial = FALSE)
    if (is.null(tri))
      stop("Cannot triangulate polygon")
    dim <- dim(tri)
    tri <- ind[tri]
    dim(tri) <- dim
    subtri[[fwd]] <<- tri
  }
  subtri <- list()
  for (i in nesting$toplevel) processOutside(i)
  res <- matrix(nrow = 3, ncol = 0)
  for (i in seq_along(subtri)) res <- cbind(res, subtri[[i]])
  attr(res, "nextvert") <- nextvert
  res
}


nestPolys_rgl <-
function (x, y = NULL)
{
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  n <- length(x)
  nas <- c(which(is.na(x) | is.na(y)), n + 1L)
  prev <- 0L
  verts <- list()
  for (i in seq_along(nas)) {
    verts[[i]] <- ind <- (prev + 1L):(nas[i] - 1L)
    tri <- triangulateSimple_rgl(x[ind], y[ind], random = TRUE,
                             plot = FALSE, partial = FALSE)
    if (is.null(tri))
      verts[[i]] <- rev(verts[[i]])
    prev <- nas[i]
  }
  nesting <- rep(list(integer()), length(verts) + 1)
  place <- function(new, toplevel) {
    placed <- FALSE
    contains <- integer()
    if (length(nesting[[toplevel]])) {
      newverts <- rbind(x[verts[[new]]], y[verts[[new]]])
      for (j in nesting[[toplevel]]) {
        prev <- rbind(x[verts[[j]]], y[verts[[j]]])
        if (pointInPoly_rgl(prev, newverts[, 1])) {
          place(new, j)
          placed <- TRUE
          break
        }
        if (pointInPoly_rgl(newverts, prev[, 1]))
          contains <- c(contains, j)
      }
    }
    if (!placed) {
      nesting[[toplevel]] <<- c(setdiff(nesting[[toplevel]],
                                        contains), new)
      nesting[[new]] <<- contains
    }
  }
  for (i in seq_along(verts)) {
    place(i, length(verts) + 1)
  }
  list(verts = verts, nesting = nesting[-length(nesting)],
       toplevel = nesting[length(nesting)])
}


triangulateSimple_rgl <- function (x, y, random = TRUE, plot = FALSE, partial = NA)
{
  n <- length(x)
  stopifnot(n == length(y))
  stopifnot(n > 2)
  it <- matrix(NA_integer_, nrow = 3, ncol = n - 2)
  verts <- seq_len(n)
  while ((m <- length(verts)) > 3) {
    i1 <- 1:m
    i2 <- i1%%m + 1
    i3 <- i2%%m + 1
    theta3 <- atan2(y[verts[i3]] - y[verts[i1]], x[verts[i3]] -
                      x[verts[i1]])
    theta2 <- atan2(y[verts[i2]] - y[verts[i1]], x[verts[i2]] -
                      x[verts[i1]])
    diff <- ((theta3 - theta2)/pi + 4)%%2
    convex <- which(diff < 1)
    if (random && length(convex) > 1)
      convex <- sample(convex)
    good <- FALSE
    for (k in convex) {
      i <- c(i1[k], i2[k], i3[k])
      tri <- rbind(x[verts[i]], y[verts[i]])
      good <- TRUE
      for (j in 2:(m - 1)) {
        i4 <- (i1[k] + j - 1)%%m + 1
        i5 <- (i1[k] + j)%%m + 1
        j <- c(i4, i5)
        if (intersectTriSeg_rgl(tri, rbind(x[verts[j]], y[verts[j]]))) {
          good <- FALSE
          break
        }
      }
      if (good) {
        if (plot)
          polygon(x[verts[i]], y[verts[i]], col = m)
        it[, m - 2] <- verts[i]
        verts <- verts[-i2[k]]
        break
      }
    }
    if (!good)
      break
  }
  if (!good) {
    if (is.na(partial)) {
      warning("Triangulation is incomplete")
      partial <- TRUE
    }
    if (partial)
      it <- it[, seq_len(n - m) + m - 2, drop = FALSE]
    else it <- NULL
  }
  else {
    if (plot)
      polygon(x[verts], y[verts], col = 3)
    it[, 1] <- verts
  }
  it
}

intersectTriSeg_rgl <- function (tri, seg)
{

  coeffs <- try(solve(rbind(tri, 1), rbind(seg, 1)), silent = TRUE)
  if (inherits(coeffs, "try-error"))
    return(TRUE)
  coeffs <- zapsmall(coeffs)
 # browser()
  dm <- dim(coeffs)
  m <- dm[1L]
  n <- dm[2L]
  #if (any(apply(coeffs <= 0, 1, all)))
  if (any(.rowSums(coeffs <= 0, m, n) == dim(coeffs)[2L]))
    return(FALSE)
  #if (any(apply(coeffs > 0, 2, all)))
  if (any(.colSums(coeffs > 0, m, n) == dim(coeffs)[1L]))
    return(TRUE)
  up <- coeffs[, 1] < 0
  dn <- coeffs[, 2] < 0
  lb <- max(-coeffs[up, 1]/(coeffs[up, 2] - coeffs[up, 1]))
  ub <- 1 - max(-coeffs[dn, 2]/(coeffs[dn, 1] - coeffs[dn,
                                                       2]))
  lb <= ub
}

pointInPoly_rgl <- function (poly, pt)
{
  n <- ncol(poly)
  i1 <- seq_len(n)
  i2 <- i1%%n + 1
  x <- poly[1, i1] + (poly[1, i2] - poly[1, i1]) * (pt[2] -
                                                      poly[2, i1])/(poly[2, i2] - poly[2, i1])
  crossings <- ((poly[2, i1] < pt[2]) & (pt[2] <= poly[2, i2]) |
                  (poly[2, i2] < pt[2]) & (pt[2] <= poly[2, i1])) & pt[1] <
    x
  sum(crossings)%%2 == 1
}

intersectSegSeg_rgl <- function (seg1, seg2)
{
  coeffs <- try(solve(cbind(seg1[, 2] - seg1[, 1], seg2[, 1] -
                              seg2[, 2]), seg2[, 1] - seg1[, 1]), silent = TRUE)
  if (inherits(coeffs, "try-error"))
    return(FALSE)
  all(zapsmall(coeffs) >= 0) && all(zapsmall(1 - coeffs) >=
                                      0)
}
