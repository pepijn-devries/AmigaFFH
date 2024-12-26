#' @rdname IFFChunk
#' @method IFFChunk IFF.ILBM
#' @export
IFFChunk.IFF.ILBM <- function(x, ...) {
  if (inherits(x, "matrix")) stop(paste("This IFF chunk interpretation is probably based on",
                                         "an ANIM DLTA chunk. It can't be converted back into an IFFChunk object."))
  rasterToIFF(x, ...)@chunk.data[[1]]
}

#' @rdname IFFChunk
#' @method IFFChunk IFF.CMAP
#' @export
IFFChunk.IFF.CMAP <- function(x, ...) {
  result <- colourToAmigaRaw(x, n.bytes = "3", ...)
  return(new("IFFChunk", chunk.type = "CMAP", chunk.data = list(result)))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.BMHD <- function(x, ...) {
  compr <- which(c("cmpNone", "cmpByteRun1") == x$Compression) - 1
  ## if the compression type is unknown, set to 0
  if (length(compr) == 0) compr <- 0
  mask <- which(c("mskNone", "mskHasMask", "mskHasTransparentColour", "mskLasso") == x$Masking) - 1
  ## if the masking type is unknown, set to 0
  if (length(mask) == 0) mask <- 0
  result <- c(.amigaIntToRaw(c(x$w, x$h), 16, FALSE),
              .amigaIntToRaw(c(x$x, x$y), 16, TRUE),
              .amigaIntToRaw(c(x$nPlanes, mask, compr), 8, FALSE),
              as.raw(x$pad)[[1]],
              .amigaIntToRaw(x$transparentColour, 16, FALSE),
              .amigaIntToRaw(c(x$xAspect, x$yAspect), 8, FALSE),
              .amigaIntToRaw(c(x$pageWidth, x$pageHeight), 16, FALSE))
  return(new("IFFChunk", chunk.type = "BMHD", chunk.data = list(result)))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.CAMG <- function(x, ...) {
  return(.inverseViewPort(x$display.mode, x$monitor))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.CRNG <- function(x, ...) {
  flag <- which(c("RNG_OFF", "RNG_ACTIVE", "RNG_REVERSE") == x$flags) - 1
  ## if the flag is unknown, set to 0
  if (length(flag) == 0) flag <- 0
  dat <- c(
    x$padding,
    .amigaIntToRaw(round(x$rate*(2^14)/60), 16, FALSE),
    .amigaIntToRaw(flag, 16, FALSE),
    .amigaIntToRaw(c(x$low, x$high), 8, FALSE)
  )
  result <- new("IFFChunk", chunk.type = "CRNG", chunk.data = list(dat))
  return(result)
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.ANIM <- function(x, ...) {
  rasterToIFF(x, ...)@chunk.data[[1]]
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.ANHD <- function(x, ...) {
  oper <- which(x$operation == c("standard", "XOR", "LongDeltaMode", "ShortDeltaMode", "GeneralDeltamode", "ByteVerticalCompression", "StereoOp5", "ShortLongVerticalDeltaMode")) - 1
  ## if the flag is unknown, set to 0
  if (length(oper) == 0) oper <- 0
  result <- c(.amigaIntToRaw(oper, 8, FALSE),
              .bitmapToRaw(x$mask, FALSE, FALSE),
              .amigaIntToRaw(c(x$w, x$h), 16, FALSE),
              .amigaIntToRaw(c(x$x, x$y), 16, TRUE),
              .amigaIntToRaw(c(x$abstime, x$reltime), 32, TRUE),
              .amigaIntToRaw(x$interleave, 8, FALSE),
              x$pad0,
              .bitmapToRaw(x$flags, TRUE, TRUE),
              x$pad1
  )
  result <- new("IFFChunk", chunk.type = "ANHD", chunk.data = list(result))
  return(result)
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.DLTA <- function(x, ...) {
  return(new("IFFChunk", chunk.type = "DLTA", chunk.data = list(x)))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.DPAN <- function(x, ...) {
  result <- c(.amigaIntToRaw(c(x$version, x$nframes), 16, FALSE),
              .bitmapToRaw(x$flags, TRUE, TRUE)
  )
  result <- new("IFFChunk", chunk.type = "DPAN", chunk.data = list(result))
  return(result)
}

#' Convert AmigaFFH objects into grDevices raster images
#'
#' Convert AmigaFFH objects that contain bitmap images into grDevices raster
#' images.
#'
#' Images on the Amiga were stored as bitmap images with indexed colour
#' palettes. This was mainly due to hardware and memory limitations.
#' Bitmap images could also be embedded in several file types. This method
#' can be used to convert AmigaFFH objects read from such files into
#' grDevices raster images ([grDevices::as.raster()]).
#'
#' @rdname as.raster
#' @name as.raster
#' @param x Object that needs to be converted into a `grDevices` raster. It
#' can be an [IFFChunk()] containing an interleaved bitmap image
#' (ILBM) or animation (ANIM), a [hardwareSprite()], an
#' [AmigaBitmapFont()] object or an [AmigaBitmapFontSet()] object.
#' @param background Use the argument `background` to
#' specify a background colour in case `x` is a [hardwareSprite()].
#' @param selected When `x` is an object of class [AmigaIcon()], `selected` can be
#' used to select a specific state. When set to `TRUE`, the raster of the [AmigaIcon()]
#' will be based on the `selected' state of the icon. Otherwise it will be based on the
#' deselected state (default).
#' 
#' When `x` is an [AmigaBasicShape()] class object, `selected` can be used to select a
#' specific layer of the shape to plot, which can be one of `"bitmap"` (default), `"shadow"` or `"collision"`.
#' @param ... Currently ignored.
#' @returns Returns a `grDevices` raster image ([grDevices::as.raster()])
#' based on `x`. If `x` is an animation ([IFFChunk()] of type ANIM),
#' a `list` of raster objects is returned.
#' @examples
#' ## load an IFF file
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## The file contains an interleaved bitmap image that can be
#' ## converted into a raster:
#' example.raster <- as.raster(example.iff)
#' 
#' ## the raster can be plotted:
#' plot(example.raster)
#' 
#' ## note that the IFFChunk can also be plotted directly:
#' plot(example.iff)
#' 
#' ## Hardware sprites can also be converted into raster images.
#' ## Let's generate a 16x16 sprite with a random bitmap:
#' spr <- new("hardwareSprite",
#'            VStop = 16,
#'            bitmap = as.raw(sample.int(255, 64, replace = TRUE)))
#' 
#' ## now convert it into a raster image.
#' ## as the background colour is not specified for hardware
#' ## sprite, we can optionally provide it here.
#' spr.raster <- as.raster(spr, background = "green")
#'
#' ## AmigaBasicShape objects can also be converted into rasters:
#' ball <- read.AmigaBasicShape(system.file("ball.shp", package = "AmigaFFH"))
#' ball.rst <- as.raster(ball)
#' @family iff.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
as.raster.IFFChunk <- function(x, ...) {
  if (x@chunk.type == "FORM") {
    result <- lapply(x@chunk.data, function(y){
      if (y@chunk.type == "ILBM") return(as.raster(y))
      if (y@chunk.type == "ANIM") return(interpretIFFChunk(y))
    })
    if (length(result) == 1) result <- result[[1]]
    return (result)
  } else if (x@chunk.type == "ILBM") {
    sub.chunks <- unlist(lapply(x@chunk.data, function(y) y@chunk.type))
    if (!("BMHD" %in% sub.chunks)) stop("No bitmap header present. Can't interpret bitmap.")
    if (!("BODY" %in% sub.chunks)) stop("No BODY chunk present. Can't convert the bitmap into a raster.")
    bm.header  <- interpretIFFChunk(getIFFChunk(x, "BMHD"))
    if (!("CAMG" %in% sub.chunks)) {
      bm.vp.mode <- NULL
      warning("No Amiga viewport available, interpretation possibly incorrect.")
    } else {
      bm.vp.mode <- interpretIFFChunk(getIFFChunk(x, "CAMG"))
      bm.vp.mode <- .display.properties(bm.vp.mode$display.mode, bm.vp.mode$monitor)
    }
    ## default palette in case 'CMAP' chunk is missing:
    bm.palette <- grDevices::gray(round(seq(0, 15, length.out = 2^bm.header$nPlanes))/15)
    ## Colours are interpreted as 12 bit when al low bits are 0, otherwise as 24 bit
    try(bm.palette <- interpretIFFChunk(getIFFChunk(x, "CMAP")), silent = TRUE)
    ## extend the palette when the extra halfbright mode is used:
    if (!is.null(bm.vp.mode) && bm.vp.mode$is.halfbright) {
      bm.palette <- c(bm.palette,
                      substr(grDevices::adjustcolor(bm.palette, 1, 0.5, 0.5, 0.5), 1, 7))
    }
    if (bm.header$Compression == "cmpByteRun1") {
      bm <- unPackBitmap(interpretIFFChunk(getIFFChunk(x, "BODY")))
    } else if (bm.header$Compression == "cmpNone") {
      bm <- interpretIFFChunk(getIFFChunk(x, "BODY"))
    } else {
      stop("Bitmap data is compressed with unsupported algorithm.")
    }
    np <- bm.header$nPlanes
    if (length(bm.palette) < (2^np)) {
      if (bm.vp.mode$is.HAM) np <- ifelse(np == 8, 6, 5)
      bm.palette <- c(bm.palette, rep("#000000", (2^np) - length(bm.palette)))
    }
    rm(np)
    if (bm.header$Masking == "mskHasTransparentColour") {
      transparent <- bm.header$transparentColour + 1
      bm.palette[transparent] <- grDevices::adjustcolor(bm.palette[transparent],
                                                        alpha.f = 0)
    }
    attr.palette <- NULL
    if ("palette" %in% names(list(...)) && is.null(list(...)$palette)){
      attr.palette <- bm.palette
      bm.palette <- NULL
    }
    if (bm.vp.mode$is.HAM) {
      result <- bitmapToRaster(bm, bm.header$w, bm.header$h, bm.header$nPlanes, NULL)
      ## It is assumed that the image is HAM8 in case there are 8 bitplanes,
      ## HAM6 in all other cases...
      ## if bm.palette is null, we are dealing with an animation
      ## and we need to return the indices rather than the colours
      if (!is.null(bm.palette)) {
        result <- .indexToHAMraster(result, bm.header$nPlanes, bm.palette, bm.header$transparentColour)
      }
    } else {
      result <- bitmapToRaster(bm, bm.header$w, bm.header$h, bm.header$nPlanes, bm.palette)
    }
    if (bm.vp.mode$is.HAM) {
      attributes(result)[["mode"]] <- ifelse(bm.header$nPlanes == 8, "HAM8", "HAM6")
    }
    attributes(result)[["asp"]] <- bm.vp.mode$aspect.y/bm.vp.mode$aspect.x
    if (!is.null(attr.palette)) attributes(result)[["palette"]] <- attr.palette
    return(result)
  } else if (x@chunk.type == "ANIM") {
    return(interpretIFFChunk(x))
  } else {
    stop(sprintf("IFF chunk of type %s cannot be converted into a raster.", x@chunk.type))
  }
}

#' @rdname plot
#' @name plot
#' @export
plot.IFF.ILBM <- function(x, y, ...) {
  ## For ANIM frames, a matrix of palette indices is returned
  ## turn that into a raster object with a grayscale palette
  if (inherits(x, "matrix")) {
    pal <- grDevices::gray(seq(0, 1, length.out = max(round(abs(x)))))
    asp <- attributes(x)$asp
    x <- as.raster(apply(x, 2, function(z) pal[round(abs(z)) + 1]))
    attributes(x)$asp <- asp
  }
  class(x) <- "raster"
  if ("asp" %in% names(list(...)))
    graphics::plot(x, y, ...) else
      graphics::plot(x, y, asp = attributes(x)$asp, ...)
}

#' @rdname plot
#' @name plot
#' @export
plot.IFF.ANIM <- function(x, y, ...) {
  invisible(lapply(x, plot.IFF.ILBM, ...))
}

#' Convert a grDevices raster image into an IFF formated bitmap image
#'
#' Convert grDevices raster images ([grDevices::as.raster()])
#' into a formal [IFFChunk()] object, as an interleaved bitmap (ILBM)
#' image.
#'
#' Convert any modern image into a interleaved bitmap (image) conform
#' Interchange File Format (IFF) specifications. If your original image
#' is in true colour (i.e., a 24 bit colour depth) it will be converted
#' into a bitmap image with an indexed palette.
#'
#' @rdname rasterToIFF
#' @name rasterToIFF
#' @param x A raster object created with [grDevices::as.raster()] which
#' needs to be converted into an IFF formated bitmap image. It is also possible to let `x` be
#' a matrix of `character`s, representing colours.
#' @param display.mode Specify the Amiga display mode that should be used.
#' See [amiga_display_modes()] for all possible options.
#' "`LORES_KEY`" is used by default, this is the lowest resolution
#' possible on the Amiga.
#' @param monitor The Amiga monitor on which the needs to be displayed.
#' See [amiga_monitors()] for more details and posible options.
#' By default "`DEFAULT_MONITOR_ID`" is used.
#' @param anim.options Currently ignored. This argument will potentitally be implemented
#' in future versions of this package. Currently, animations are always encoded
#' with the "ByteVerticalCompression" in this package (when `x` is a list of
#' `raster` objects).
#' @param ... Arguments passed on to [rasterToBitmap()].
#' @returns Returns an [IFFChunk()] object holding an Interleaved
#' Bitmap (ILBM) image based on `x`.
#' @examples
#' ## first: Let's make a raster out of the 'volcano' data, which we can use in the example:
#' volcano.raster <- as.raster(t(matrix(terrain.colors(1 + diff(range(volcano)))[volcano -
#'   min(volcano) + 1], nrow(volcano))))
#' 
#' ## Turning the raster into an IFFChunk object is easy:
#' volcano.iff <- rasterToIFF(volcano.raster)
#' 
#' ## This object can be saved as an IFF file using write.iff
#' 
#' ## in special modes HAM6 and HAM 8 higher quality images
#' ## can be obtained. See 'rasterToBitmap' for more info on the
#' ## special HAM modes.
#' volcano.ham <- rasterToIFF(volcano.raster, "HAM_KEY", depth = "HAM8")
#' 
#' ## The result can be further improved by applying dithering
#' volcano.ham.dither <- rasterToIFF(volcano.raster, "HAM_KEY", depth = "HAM8",
#'   indexing = function(x, length.out) {
#'     index.colours(x, length.out, dither = "JJN", iter.max = 20)
#'   })
#' @family iff.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
rasterToIFF <- function(x,
                        display.mode = as.character(AmigaFFH::amiga_display_modes$DISPLAY_MODE),
                        monitor = as.character(AmigaFFH::amiga_monitors$MONITOR_ID),
                        anim.options,
                        ...) {
  display.mode <- match.arg(display.mode)
  monitor <- match.arg(monitor)
  pars <- list(...)
  if (is.null(pars$depth)) pars$depth <- 3
  if (is.null(pars$colour.depth)) pars$colour.depth <- "12 bit"
  if (grepl("EHB|EXTRAHALFBRITE", display.mode)) stop("Sorry, 'extra halfbrite' modes is currently not implemented")
  special.mode <- "none"
  if (pars$depth %in% c("HAM6", "HAM8")) {
    if (!grepl("HAM", display.mode)) warning("Display mode should be a HAM mode, when 'depth' is set to 'HAM6' or 'HAM8'. Display mode is corrected to 'HAM_KEY'.")
    display.mode <- "HAM_KEY"
  }
  if (grepl("HAM", display.mode)) {
    if (pars$depth %in% c("HAM6", "HAM8")) {
      special.mode <- pars$depth
      pars$colour.depth <- ifelse(special.mode == "HAM6", "12 bit", "24 bit")
      pars$depth <- ifelse(special.mode == "HAM6", 6, 8)
    } else {
      special.mode <- ifelse(pars$colour.depth == "24 bit", "HAM8", "HAM6")
    }
  }
  if (is.list(x)) {
    if (length(x) < 2) stop("When x is a list of rasters, it will be converted to an anim. x should have a length of at least 2.")
    if (any(unlist(lapply(x, function(y) !inherits(y, c("raster", "matrix")) || !all(.is.colour(y))))))
      stop("All elements of x should be a grDevices raster or a matrix of colours")
    if ("indexing" %in% names(list(...))) {
      x <- list(...)$indexing(x = x, length.out = ifelse(special.mode %in% c("HAM6", "HAM8"),
                                                         special.mode,
                                                         2^pars$depth))
    } else {
      x <- index.colours(x, length.out = ifelse(special.mode %in% c("HAM6", "HAM8"),
                                                special.mode,
                                                2^pars$depth))
    }
    pal <- attributes(x)$palette
    trans <- attributes(x)$transparent

    anhd <- lapply(1:(length(x) + 2), function(z) {
      anhdz <- list(
        operation       = "ByteVerticalCompression",
        mask            = rep(FALSE, 8),
        w               = dim(x[[1]])[[2]],
        h               = dim(x[[1]])[[1]],
        x               = 0,
        y               = 0,
        abstime         = z*2,
        reltime         = 2,
        interleave      = 0,
        pad0            = raw(1),
        flags           = rep(FALSE, 32),
        pad1            = raw(16)
      )
      class(anhdz) <- "IFF.ANHD"
      anhdz <- IFFChunk(anhdz)
      anhdz
    })
    
    frame1 <- .indexToBitmap(x[[1]], pars$depth, TRUE)
    frame1 <- .bitmapToILBM(frame1, dim(x[[1]]), display.mode, monitor, pars$depth, pars$colour.depth, pal, trans)
    ## DPaint needs this class to set the correct number of frames:
    dpan <- list(version = 4, nframes = length(x), flags = rep(FALSE, 32))
    class(dpan) <- "IFF.DPAN"
    frame1@chunk.data <- c(frame1@chunk.data[1],
                           anhd[[1]],
                           frame1@chunk.data[2],
                           IFFChunk(dpan),
                           frame1@chunk.data[3:4])
    frame1 <- new("IFFChunk", chunk.type = "FORM", chunk.data = list(frame1))

    ## .byteVerticalCompression returns a list of DLTA chunks, with the first frame missing (NULL)
    x <- .byteVerticalCompression(x, pars$depth)
    x[-1] <- lapply(2:length(x), function(y) {
      res <- new("IFFChunk", chunk.type = "ILBM", chunk.data = list(anhd[[y]], x[[y]]))
      res <- new("IFFChunk", chunk.type = "FORM", chunk.data = list(res))
      return(res)
    })
    x[[1]] <- frame1
    x <- new("IFFChunk", chunk.type = "ANIM", chunk.data = x)
    x <- new("IFFChunk", chunk.type = "FORM", chunk.data = list(x))
    return(x)
  }
  if (!inherits(x, c("raster", "matrix")) || !all(.is.colour(x))) stop("x should be a raster object or a matrix of colours.")
  if ("depth" %in% names(list(...))) {
    bm <- rasterToBitmap(x, ...)
  } else {
    bm <- rasterToBitmap(x, depth = ifelse(special.mode %in% c("HAM6", "HAM8"),
                                           special.mode,
                                           pars$depth), ...)
  }
  pal <- attributes(bm)$palette
  transparent <- attributes(bm)$transparent

  ## Create an ILBM chunk based on the bitmap data
  ilbm <- .bitmapToILBM(bm, dim(x), display.mode, monitor, pars$depth, pars$colour.depth, pal, transparent)
  ## Create a FORM chunk, encapsulating the bitmap info
  form <- new("IFFChunk", chunk.type = "FORM", chunk.data = list(ilbm))
  return(form)
}

## decompression for anim dlta frames
## dlta = raw data from dlta chunk
.byteVerticalDecompression <- function(dlta, w, h, interleave, use.xor, previous = NULL) {
  interleave[interleave == 0] <- 2
  pointers <- .rawToAmigaInt(dlta[1:(16*4)], 32, FALSE)[1:8]
  if (all(pointers == 0)) {
    if (is.null(previous) || length(previous) == 0) {
      result <- matrix(0, h, w)
      class(result) <- c("IFF.ILBM", "IFF.ANY", class(result))
      return(result)
    } else {
      prev <- length(previous) - interleave + 1
      prev[prev < 1] <- 1
      return(previous[[prev]])
    }
  }
  bitmap.layers <- which(pointers > 1)
  bitmap.layers <- 1:bitmap.layers[bitmap.layers == max(bitmap.layers)]
  bitmap.layers <- pointers[bitmap.layers]
  ## loop the bitmap depth dimensions:
  result <- lapply(1:length(bitmap.layers), function(y) {
    offs <- bitmap.layers[[y]]
    prev <- NULL
    if (!is.null(previous)) {
      prev <- length(previous) - (interleave - 1)
      if (prev < 1) prev <- 1
      prev <- previous[[prev]]
      prev <- apply(prev, 2, function(z) {
        as.logical(floor(z/(2^(y - 1))) %% 2)
      })
      prev <- cbind(prev, matrix(FALSE, h, w - ncol(prev)))
    }
    ## if the pointer is zero for a bitmap layer, there is no change. Just return previous frame
    if (offs == 0) {
      if (is.null(prev)) {
        return (matrix(0, h, w))
      } else return(prev)
    }
    ## loop columns:
    layer <- NULL
    for (i in 1:(w/8)) {
      op.count <- .rawToAmigaInt(dlta[1 + offs], 8, FALSE)
      offs <- offs + 1
      row.result <- matrix(FALSE, 0, 8)
      if (op.count > 0) {
        for (j in 1:op.count) {
          op <- .rawToAmigaInt(dlta[1 + offs], 8, FALSE)
          offs <- offs + 1
          if (op == 0) { ## RUN operation; followed by number of repetitions and the repeating byte
            rep.count <- .rawToAmigaInt(dlta[1 + offs], 8, FALSE)
            rep.dat <- matrix(rep(
              as.logical(.rawToBitmap(dlta[2 + offs], TRUE, FALSE)),
              rep.count
            ),
            ncol = 8, byrow = TRUE)
            if (use.xor)
              rep.dat <- xor(prev[nrow(row.result) + 1:nrow(rep.dat), i*8 + (-7:0)], rep.dat)
            row.result <- rbind(row.result, rep.dat)
            offs <- offs + 2
          } else if (op < 0x80) { ## SKIP operation; move the cursor op bytes forward
            if (is.null(previous) || (is.list(previous) && length(previous) == 0)) {
              row.result <- rbind(row.result, matrix(FALSE, nrow = op, ncol = 8))
            } else {
              nroff <- nrow(row.result) + 1
              if (length(nroff) == 0) nroff <- 1
              row.result <- rbind(row.result, prev[nroff + 0:(op - 1), i*8 + (-7:0)])
            }
          } else { ## DUMP operation; use op bytes literally
            op.cor <- op - 0x80
            temp <- as.logical(.rawToBitmap(dlta[(1 + offs):(offs + op.cor)], TRUE, FALSE))
            temp <- matrix(temp, ncol = 8, byrow = TRUE)
            if (use.xor)
              temp <- xor(prev[nrow(row.result) + 1:nrow(temp), i*8 + (-7:0)], temp)
            row.result <- rbind(row.result, temp)
            offs <- offs + op.cor
          }
        }
      }
      if (nrow(row.result) == 0 || op.count == 0) {
        if (is.null(previous)) {
          row.result <- matrix(FALSE, ncol = 8, nrow = h)
        } else {
          row.result <- prev[1:h, i*8 + (-7:0)]
        }
      }
      if (is.null(dim(row.result)) || !all(dim(row.result) == c(h, 8))) stop("Could not decode the bitmap correctly. If the bitmap was encoded with the AmigaFFH package, please contact the package author to get this fixed.")
      layer <- cbind(layer, row.result)
    }
    layer <- apply(layer, 2, function(z) (2^(y - 1))*as.numeric(z))
    return(layer)
  })
  result <- Reduce("+", result)
  class(result) <- c("IFF.ILBM", "IFF.ANY", class(result))
  return(result)
}

## vertical delta compression
## x should be a list of matrices with palette indices
## length of x should be > 1 (not checked here)
## the function returns a list of DLTA chunks, with the first element empty (NULL).
.byteVerticalCompression <- function(x, depth) {
  x <- c(x, x[1:2])
  h   <- dim(x[[1]])[[1]]
  w   <- dim(x[[1]])[[2]]
  x   <- lapply(x, function(y) cbind(y, matrix(0, nrow = h, ncol = (-w %% 16))))
  wbm <- w + (-w %%16)
  ## first frame should be ILBM with BMHD and BODY (constructed later and skipped here)
  result <- vector("list", length(x))
  for (i in 2:length(x)) {
    prev.id <- i - 2
    if (prev.id < 1) prev.id <- 1
    ## loop bitmap layers
    pointers <- NULL
    no.change <- rep(FALSE, 16)
    for (j in 1:depth) {
      ## loop columns
      dep.data <- NULL
      for (k in 1:(wbm/8)) {
        previous.column <- floor((x[[prev.id]][,k*8 + (-7:0)] - 1)/(2^(j - 1))) %% 2
        curr.column     <- floor((x[[i]][,k*8 + (-7:0)] - 1)/(2^(j - 1))) %% 2
        row.similarity <- apply(previous.column == curr.column, 1, all)
        ## when there are only short fragments the same, classify them as different
        ## this is to avoid the large overhead of the operation bytes
        run.lengths <- rle(row.similarity)$lengths
        run.lengths <- rep(run.lengths, run.lengths)
        row.similarity[row.similarity & run.lengths < 3] <- FALSE
        ## If the entire row is similar to the previous frame, add 0 as op.count to the result
        cursor <- 1
        op.count <- 0
        op.data <- NULL
        if (!all(row.similarity)) {
          while (cursor <= h) {
            if (row.similarity[[cursor]]) { ## when elements at the cursor are the same compared to the previous frame, the SKIP operation is in place
              skip <- which(diff(c(row.similarity[cursor:h], FALSE)) != 0)[[1]]
              repskip <- floor(skip/127)
              ## note that the op.count needs to preceed each column of data
              ## it will be added after all ops of the column have been processed
              op.data <- c(op.data,
                           .amigaIntToRaw(c(
                             rep(127, repskip),
                             skip %% 127), 8, FALSE))
              op.count <- op.count + repskip + 1
              cursor <- cursor + skip
            } else {
              dup <- c(TRUE, duplicated(curr.column[cursor:h,, drop = FALSE], fromLast = FALSE)[-1]) & !row.similarity[cursor:h]
              dup.run.length <- rle(dup)$lengths
              dup.run1 <- dup.run.length[[1]]
              dup.run1[dup.run1 > 255] <- 255
              dup.run.length <- c(1, dup.run.length[1] - 1, dup.run.length[-1])
              dup.run.length <- rep(dup.run.length, dup.run.length)
              if (dup.run1 > 3 && length(dup) > 3 && dup[[2]]){ ## When elements at the cursor are not the same as the previous frame, and they are repetative, the RUN operation is in place
                op.data <- c(op.data,
                             .amigaIntToRaw(c(0, dup.run1), 8, FALSE),
                             .bitmapToRaw(as.logical(t(curr.column[cursor,, drop = FALSE])),
                                                      TRUE, FALSE))
                op.count <- op.count + 1
                cursor <- cursor + dup.run1
              } else { ## When elements at the cursor are not the same as the previous frame, and they are not repetative, the DUMP operation is in place
                skip <- which(diff(c(row.similarity[cursor:h] |
                                       c(dup.run.length[-1] > 2 & dup[-1], FALSE), TRUE)) != 0)[[1]]
                skip[skip > 127] <- 127
                dat <- .bitmapToRaw(as.logical(t(curr.column[cursor:(cursor + skip - 1),, drop = FALSE])),
                                                TRUE, FALSE)
                op.data <- c(op.data,
                             .amigaIntToRaw(skip + 0x80, 8, FALSE),
                             dat)
                cursor <- cursor + skip
                op.count <- op.count + 1
              }
            }
          }
        }
        ## After each column, bind the result to the previous column. Each column should start with the count of
        ## the total number of operations in the column, followed by the encoded data.
        dep.data <- c(dep.data,
                      .amigaIntToRaw(op.count, 8, FALSE),
                      op.data)
      }
      if (all(dep.data == raw(1))) {
        no.change[[j]] <- TRUE
        dep.data <- NULL
      }
      pointers <- c(pointers, length(dep.data))
      if (!is.null(dep.data)) result[[i]] <- c(result[[i]], dep.data) 
    }
    pointers <- cumsum(pointers)
    ptrs <- rep(0, 16)
    ptrs[1:length(pointers)] <- 64 + c(0 , pointers[0:(length(pointers) - 1)])
    ptrs[no.change] <- 0
    ## combine the op-data from all bitmap layers, starting with pointers to the start of the data for each layer
    ## if there was no change to a bitmap layer, there is no data. The pointer is set to 0 in that case.
    result[[i]] <- c(
      .amigaIntToRaw(ptrs, 32, FALSE),
      result[[i]]
    )
    ## DLTA should be preceded by an ANHD chunk, this is not added by this function
    result[[i]] <- new("IFFChunk", chunk.type = "DLTA", chunk.data = list(result[[i]]))
  }
  return(result)
}

.bitmapToILBM <- function(bm, dim.x, display.mode, monitor, depth, colour.depth, pal, transparent) {
  bm <- .bitmapToRaw(bm, TRUE, FALSE)
  bm <- matrix(bm, nrow = 2*ceiling(dim.x[[2]]/16), byrow = FALSE)
  bm <- c(unlist(apply(bm, 2, packBitmap)))
  
  ## Create a BODY chunk based on the bitmap data
  body <- new("IFFChunk", chunk.type = "BODY", chunk.data = list(bm))

  ## Create a CAMG chunk, specifying the displaymode
  camg <- .inverseViewPort(display.mode, monitor)
  disp <- .amigaViewPortModes(camg@chunk.data[[1]])
  disp.prop <- .display.properties(disp$display.mode, disp$monitor)
  ## Create BMHD chunk, with bitmap header information
  BMHD <- list(
    w                 = dim.x[2],
    h                 = dim.x[1],
    x                 = 0,
    y                 = 0,
    nPlanes           = depth,
    Masking           = ifelse(!is.null(transparent) && !is.na(transparent), "mskHasTransparentColour", "mskNone"),
    Compression       = "cmpByteRun1",
    pad               = raw(1),
    transparentColour = ifelse(!is.null(transparent) && !is.na(transparent), transparent - 1, 0),
    xAspect           = disp.prop$aspect.x,
    yAspect           = disp.prop$aspect.y,
    pageWidth         = disp.prop$screenwidth,
    pageHeight        = disp.prop$screenheight
  )
  class(BMHD) <- "IFF.BMHD"
  hdr <- IFFChunk(BMHD)
  ## Create a CMAP chunk with palette information
  class(pal) <- "IFF.CMAP"
  cmap <- IFFChunk(pal, colour.depth = colour.depth)
  # cmap <- new("IFFChunk", chunk.type = "CMAP", chunk.data = list(colourToAmigaRaw(pal, pars$colour.depth, "3")))
  ## Create a ILBM chunk, which hold all required bitmap information
  ilbm <- new("IFFChunk", chunk.type = "ILBM", chunk.data = list(
    hdr, cmap, camg, body
  ))
  return(ilbm)
}

.indexToHAMraster <- function(x, depth, palette, transparentColour) {
  control.mask  <- bitwShiftL(3, depth - 2)
  max_color     <- ifelse(depth == 8, 255, 15)
  color_divisor <- ifelse(depth == 8, 1, 17)
  color_multi   <- ifelse(depth == 8, 255/63, 1)
  x <- apply(x, 1, function(y) {
    control.flags <- 3*bitwAnd(y, control.mask)/control.mask
    y.shift <- (y - control.mask*control.flags/3)
    z             <- rep(NA, length(y))
    z[control.flags == 0] <- palette[y.shift[control.flags == 0] + 1]
    for (i in 1:length(y)) {
      if (is.na(z[i])) {
        z0   <- z[i - 1]
        if (length(z0) == 0) z0 <- palette[transparentColour + 1]
        cl   <- grDevices::col2rgb(z0)
        z[i] <- grDevices::rgb(
          ifelse(control.flags[i] == 2, color_multi*y.shift[i], cl["red",]/color_divisor),
          ifelse(control.flags[i] == 3, color_multi*y.shift[i], cl["green",]/color_divisor),
          ifelse(control.flags[i] == 1, color_multi*y.shift[i], cl["blue",]/color_divisor),
          maxColorValue = max_color
        )
      }
    }
    z
  })
  as.raster(t(x))
}
