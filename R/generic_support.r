#' @rdname dither
#' @name dither
#' @export
dither <- function(x, method, ...) {
  ## I made this an S3 generic such that I could implement a dither function in the future
  ## for audio waves for downsampling audio...
  UseMethod("dither", x)
}

#' Convert colours to Amiga compatible raw data or vice versa
#'
#' Convert colours to Amiga compatible raw data or vice versa, such that
#' it can be used in graphical objects from the Commodore Amiga.
#'
#' On the original Commodore Amiga chipset, graphics used indexed palettes of
#' 12 bit colours. Colours are specified by their RGB (Red, Green and Blue)
#' values, each component requiring 4 bits (with corresponding values ranging
#' from 0 up to 15). Data structures on the Amiga were WORD (2 bytes) aligned.
#' Colours are therefore typically stored in either 2 bytes (skipping the
#' first four bits) or 3 bytes (one byte for each value).
#' 
#' These functions can be used to convert R colours into the closest matching
#' Amiga colour in a `raw` format, or vice versa. Note that later Amiga
#' models with the advanced (graphics) architecture (known as AA or AGA) allowed
#' for 24 bit colours.
#'
#' @param x In the case `amigaRawToColour` is called, `x` should
#' be a `vector` of `raw` data. The length of this vector should
#' be a multiple of 2 (when `n.bytes = "2"`) or 3 (when
#' `n.bytes = "3"`). When `colourToAmigaRaw` is called, `x`
#' should be a `character` strings representing a colour.
#' @param colour.depth A `character` string: `"12 bit"` (default) or
#' `"24 bit"`. The first should be used in most cases, as old Amigas
#' have a 12 bit colour depth.
#' @param n.bytes A `character` string: `"2"` or `"3"`. The
#' number of bytes that is used or should be used to store each colour.
#' @returns In the case `amigaRawToColour` is called, a (vector of)
#' colour `character` string(s) is returned. When `colourToAmigaRaw`
#' is called, `raw` representing the colour(s) specified in `x` is
#' returned.
#' 
#' @rdname colourToAmigaRaw
#' @name colourToAmigaRaw
#' @examples
#' ## Let's create some Amiga palettes:
#' colourToAmigaRaw(c("red", "navy blue", "brown", "#34AC5A"))
#' 
#' ## let's do the reverse.
#' ## this is white:
#' amigaRawToColour(as.raw(c(0x0f, 0xff)))
#' 
#' ## this is white specified in 3 bytes:
#' amigaRawToColour(as.raw(c(0xf0, 0xf0, 0xf0)), n.bytes = "3")
#' 
#' ## lower nybbles are ignored, you will get a warning when it is not zero:
#' # amigaRawToColour(as.raw(c(0xf0, 0xf0, 0x0f)), n.bytes = "3")
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
colourToAmigaRaw <- function(x, colour.depth = c("12 bit", "24 bit"), n.bytes = c("2", "3")) {
  colour.depth <- match.arg(colour.depth, c("12 bit", "24 bit"))
  n.bytes <- match.arg(n.bytes, c("2", "3"))
  if (colour.depth == "24 bit" && n.bytes == "2") stop("3 bytes are required to store 24 bit colours!")
  col <- grDevices::col2rgb(x)
  if (colour.depth == "12 bit") {
    col <- floor((col + 4)/16.5)
  }
  if (colour.depth == "24 bit") col <- col/16
  if (n.bytes == "3") {
    as.vector(apply(col, 2, function(y) .amigaIntToRaw(16*y, 8, FALSE)))
  } else {
    as.vector(apply(col, 2, function(y) as.raw(c(y[1], y[2]*16 + y[3]))))
  }
}

#' @rdname colourToAmigaRaw
#' @name amigaRawToColour
#' @export
amigaRawToColour <- function(x, colour.depth = c("12 bit", "24 bit"), n.bytes = c("2", "3")) {
  ## x = raw
  
  colour.depth <- match.arg(colour.depth, c("12 bit", "24 bit"))
  n.bytes <- match.arg(n.bytes, c("2", "3"))
  if (n.bytes == "2" && (length(x) %% 2) != 0) stop("x should be a vector of even length.")
  if (n.bytes == "3" && (length(x) %% 3) != 0) stop("x should be a vector with a multiple length of 3.")
  if (colour.depth == "24 bit" && n.bytes == "2") stop("3 bytes are required to store 24 bit colours!")
  hi <- .hiNybble(x)
  lo <- .loNybble(x)
  if (colour.depth == "24 bit" && n.bytes == "3") {
    sq <- seq(1, to = length(x), by = 3)
    x <- .rawToAmigaInt(x, 8, FALSE)
    return(grDevices::rgb(x[sq]/255, x[sq + 1]/255, x[sq + 2]/255))
  } else if (colour.depth == "12 bit" && n.bytes == "3") {
    sq <- seq(1, to = length(x), by = 3)
    hi <- .hiNybble(x)
    if (any(lo != 0)) warning("The low nybble is not zero for all colours.")
    return(grDevices::rgb(hi[sq]/15, hi[sq + 1]/15, hi[sq + 2]/15))
  } else {
    x <- as.vector(rbind(hi, lo))
    sq <- seq(1, to = length(x), by = 4)
    if (any(x[sq] != 0)) warning("The low nybble is not zero for all colours.")
    x <- x[-sq]
    sq <- seq(1, to = length(x), by = 3)
    return(grDevices::rgb(x[sq]/15, x[sq + 1]/15, x[sq + 2]/15))
  }
}

#' A routine to (un)pack bitmap data
#'
#' A very simplistic lossless routine to (un)pack repetitive bitmap data. Often
#' used in InterLeaved BitMap (ILBM) images in IFF containers ([IFFChunk()]).
#'
#' InterLeaved BitMap (ILBM) images on the Amiga often use a packing algorithm
#' referred to as `ByteRun1'. This routine was introduced first on
#' the Macintosh where it was called PackBits. It is a form of run-length encoding
#' and is very simple:
#' when a specific byte is repeated in a bitmap, it is replaced by
#' a (signed negative) byte telling how many times the following byte
#' should be repeated. When a series of bytes are not repetitive, it
#' is preceded by a (signed positive) byte telling how long the non
#' repetitive part is.
#' 
#' Not very complicated, but for most images some bytes can be shaved
#' off the file. This was very useful when everything had to be stored
#' on 880 kilobyte floppy disks with little CPU time to spare. Note
#' that the file size can also increase for (noisy) images.
#' 
#' This packing routine will pack the entire bitmap (`x`)
#' at once. The IFF file format requires packing of bitmap data per
#' scanline. This is done automatically by the [rasterToIFF()]
#' function, which calls this packing routine per scanline.
#'
#' @param x `raw` data, usually representing a (packed) bitmap.
#' @returns Returns packed or unpacked `raw` data, depending on
#' whether `packBitmap` or `unPackBitmap` was called.
#' 
#' @rdname packBitmap
#' @name packBitmap
#' @examples
#' ## generate some random raw data:
#' dat.rnd <- as.raw(sample.int(10, 100, TRUE))
#' 
#' ## try to pack it:
#' pack.rnd <- packBitmap(dat.rnd)
#' 
#' ## due to the random nature of the source data
#' ## the data could not be packed efficiently.
#' ## The length of the packed data is close to
#' ## the length of the original data:
#' length(pack.rnd) - length(dat.rnd)
#' 
#' ## Now generate similar data but sort it
#' ## to generate more repetitive data:
#' dat.srt  <- as.raw(sort(sample.int(10, 100, TRUE)))
#' pack.srt <- packBitmap(dat.srt)
#' 
#' ## This time the packing routing is more successful:
#' length(pack.srt) - length(dat.srt)
#' 
#' ## The original data can always be obtained
#' ## from the packed data:
#' all(dat.rnd == unPackBitmap(pack.rnd))
#' all(dat.srt == unPackBitmap(pack.srt))
#' @references <http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/node01C0.html>
#' @references <https://en.wikipedia.org/wiki/PackBits>
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
packBitmap <- function(x) {
  if (typeof(x) != "raw") stop("Argument 'x' should be raw data")
  n <- length(x)
  y <- x[-1L] != x[-length(x)]
  i <- c(which(y | is.na(y)), n)
  l <- diff(c(0L, i))
  while (any(l > 128)) {
    i <- sort(c(i, i[l > 128] - l[l > 128] + 128))
    l <- diff(c(0L, i))
  }
  ## Skip double repeats, as there is a large overhead from the packing byte:
  sel <- l > 1 & l < 4
  i <- c(i, rep(i[sel], l[sel] - 1))
  l <- c(l, rep(1, length(i) - length(l)))
  l[l > 1 & l < 4] <- 1
  l <- l[order(i)]
  i <- i[order(i)]
  while (any(duplicated(i))) {
    i[duplicated(i, fromLast = TRUE)] <- i[duplicated(i, fromLast = TRUE)] - 1
  }
  ## End skipping double repeats
  one.series.start <- which(diff(c(FALSE, l == 1, FALSE)) == 1)
  one.series.end <- which(diff(c(FALSE, l == 1, FALSE)) == -1) - 1
  if (length(one.series.start) != length(one.series.end)) stop("Unexpected error in packing the bitmap. Please report this error to the package author.")
  one.series <- mapply(function(start, end) {
    list(x[i[start[[1]]]:i[end[[1]]]])
  }, start = one.series.start,
  end   = one.series.end)
  one.series <- lapply(one.series, function(y) {
    yl <- length(y)
    result <- NULL
    while (yl > 128) {
      result <- c(result, .amigaIntToRaw(127, 8, TRUE), y[1:128])
      yl <- yl - 128
      y <- y[-1:-128]
    }
    return(c(result, .amigaIntToRaw(yl - 1, 8, TRUE), y))
  })
  result <- rep(list(raw(0)), length(l))
  result[one.series.start] <- one.series
  more.series <- mapply(function(y, dat, rep) {
    list(c(.amigaIntToRaw(-rep + 1, 8, TRUE), dat))
  }, dat = x[i[l > 1]], rep = l[l > 1])
  result[l > 1] <- more.series
  result <- unlist(result)
  return (result)
}

#' @rdname packBitmap
#' @name unPackBitmap
#' @export
unPackBitmap <- function(x) {
  if (typeof(x) != "raw") stop("Argument 'x' should be raw data")
  ## Very simple packing routine for bitmap images
  ## TODO this routine is very slow due to the while loop. See if this routine can be implemented more efficiently
  result <- raw(0)
  offset <- 0
  while (offset < length(x)) {
    n <- .rawToAmigaInt(x[offset + 1], 8, TRUE)
    if (n == -128) {
      offset <- offset + 1
    } else if (n < 0) {
      result <- c(result, rep(x[offset + 2], -n + 1))
      offset <- offset + 2
    } else {
      result <- c(result, x[offset + 2:(n + 2)])
      offset <- offset + 2 + n
    }
  }
  return(result)
}

#' Convert an Amiga bitmap image into a raster
#'
#' Amiga images are usually stored as bitmap images with indexed colours. This
#' function converts raw Amiga bitmap data into raster data
#' ([grDevices::as.raster()]).
#'
#' Bitmap images stored as raw data, representing palette index colours, can
#' be converted into raster data ([grDevices::as.raster()]). The latter
#' data can easily be plotted in R. It is usually not necessary to call this function
#' directly, as there are several more convenient wrappers for this function. Those
#' wrappers can convert specific file formats (such as IFF ILBM and Hardware Sprites,
#' see [AmigaFFH::as.raster()]) into raster objects. This function is
#' provided for completeness sake (or for when you want to search for images in an
#' amiga memory dump).
#'
#' @param x a `vector` of `raw` values, representing bitmap data.
#' @param w Width in pixels of the bitmap image. Can be any positive value. However,
#' bitmap data is `word' aligned on the amiga. This means that the width of the stored
#' bitmap data is a multiple of 16 pixels. The image is cropped to the width specified here.
#' @param h Height in pixels of the bitmap image.
#' @param depth The colour depth of the bitmap image (i.e., the number of bit planes).
#' The image will be composed of `2^depth` indexed colours.
#' @param palette A `vector` of `2^depth` colours, to be used for the indexed
#' colours of the bitmap image. By default, a grayscale palette is used.
#' When explicitly set to `NULL`, this function returns a matrix with palette index
#' values.
#' @param interleaved A `logical` value, indicating whether the bitmap is interleaved.
#' An interleaved bitmap image stores each consecutive bitmap layer per horizontal scanline.
#' @returns Returns a raster object ([as.raster()]) as specified in
#' the [grDevices()] package. Unless, `palette` is set to `NULL`,
#' in which case a `matrix` with `numeric` palette index values is returned.
#' 
#' @rdname bitmapToRaster
#' @name bitmapToRaster
#' @examples
#' ## first load an example image:
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## get the raw bitmap data, which is nested in the InterLeaved BitMap (ILBM)
#' ## IFF chunk as the BODY:
#' bitmap.data <- interpretIFFChunk(getIFFChunk(example.iff, c("ILBM", "BODY")))
#' 
#' ## In order to translate the bitmap data into a raster object we need
#' ## to know the image dimensions (width, height and colour depth). This
#' ## information can be obtained from the bitmap header (BMHD):
#' 
#' bitmap.header <- interpretIFFChunk(getIFFChunk(example.iff, c("ILBM", "BMHD")))
#' 
#' ## First the bitmap data needs to be unpacked as it was stored in a compresssed
#' ## form in the IFF file (see bitmap.header$Compression):
#' 
#' bitmap.data <- unPackBitmap(bitmap.data)
#' 
#' ## It would also be nice to use the correct colour palette. This can be obtained
#' ## from the CMAP chunk in the IFF file:
#' 
#' bitmap.palette <- interpretIFFChunk(getIFFChunk(example.iff, c("ILBM", "CMAP")))
#' 
#' example.raster <- bitmapToRaster(bitmap.data,
#'                                  bitmap.header$w,
#'                                  bitmap.header$h,
#'                                  bitmap.header$nPlanes,
#'                                  bitmap.palette)
#' 
#' ## We now have a raster object that can be plotted:
#' 
#' plot(example.raster, interpolate = FALSE)
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
bitmapToRaster <- function(x, w, h, depth, palette = grDevices::gray(seq(0, 1, length.out = 2^depth)), interleaved = TRUE) {
  if (!is.raw(x)) stop("x should be a vector of raw values.")
  w <- round(w)
  h <- round(h)
  if (w < 1 || h < 1) stop("Width and heigth should both be at least 1 pixel.")
  if (depth != round(depth) || depth < 1) stop("Depth should be a whole positive number.")
  if (!is.null(palette) && any(!.is.colour(palette))) stop("Palette should be composed of colours only.")
  if (!is.null(palette) && length(palette) != (2^depth)) stop("Palette should have a length of 2^depth.")
  if (length(interleaved) > 1) warning("More than 1 interleave value is given, only the first element of the vector is used.")
  interleaved <- as.logical(interleaved[[1]])
  ## invert bytes and longs is opposite to the defaults in adfExplorer.
  ## Does the user need to be able to change these values for bitmap images?
  x <- .rawToBitmap(x, invert.bytes = TRUE, invert.longs = FALSE)
  if (interleaved) {
    x <- array(x, c(16*ceiling(w/16), depth, h))
    x <- apply(x, c(1, 3), function(y) {
      sum(2^(0:(length(y) - 1)) * as.numeric(y))
    })
  } else {
    x <- array(x, c(16*ceiling(w/16), h, depth))
    x <- apply(x, c(1, 2), function(y) {
      sum(2^(0:(length(y) - 1)) * as.numeric(y))
    })
  }
  if (is.null(palette)) {
    x <- matrix(x, ncol = h, byrow = FALSE)
    x <- t(x)[, 1:w, drop = FALSE]
    return(x)
  } else {
    x <- matrix(palette[x + 1], ncol = h, byrow = FALSE)
    x <- t(x)[, 1:w, drop = FALSE]
    return(grDevices::as.raster(x))
  }
}

#' Convert a grDevices `raster` object into binary bitmap data
#'
#' Converts an image represented by a grDevices `raster` object into binary
#' (Amiga) bitmap data.
#'
#' Images represented by grDevices `raster` objects are virtually true colour (24 bit
#' colour depth) and an alpha layer (transparency). On the early Amiga's the chipset
#' (in combination with memory restrictions) only allowed images with indexed
#' palettes. The colour depth was 12 bit with the original chipset and the number
#' of colours allowed in a palette also depended on the chipset. This function
#' will allow you to convert a `raster` object into binary bitmap data with
#' an indexed palette. This means that the image is converted in a lossy way
#' (information will be lost). So don't expect the result to have the same quality as
#' the original image.
#'
#' With the `depth` argument, the raster can also be converted
#' to special mode bitmap images. One of these modes is the
#' \sQuote{hold and modify} (HAM). In this mode two of the bitplanes
#' are reserved as modifier switches. If the this switch equals
#' zero, the remainder of the bitplanes are used as an index for
#' colours in a fixed palette. If the switch equals 1, 2 or 3, the
#' red, green or blue component of the previous is modified, using the
#' number in the remainder of the bitplanes. So it holds the previous
#' colour but modifies one of the colour components (hence the term
#' \sQuote{hold and modify}.) Here only the HAM6 and
#' the HAM8 mode are implemented. HAM6 uses 6 bitplanes and a 12 bit
#' colour depth, HAM8 uses 8 bitplanes and a 24 bit colour depth.
#' 
#' The HAM mode was a special video modes supported by Amiga hardware.
#' Normal mode bitmap images with a 6 bit depth would allow for a
#' palette of 64 (2^6) colours, HAM6 can display 4096 colours with
#' the same bit depth.
#' 
#' In addition to HAM6 and HAM8, sliced HAM (or SHAM) was another
#' HAM variant. Using the coprocessor on the Amiga, it was possible
#' to change the palette at specific scanlines, increasing the number
#' of available colours even further. The SHAM mode is currently not
#' supported by this package.
#' @param x A raster object created with [grDevices::as.raster()] which
#' needs to be converted into bitmap data. It is also possible to let `x` be
#' a matrix of `character`s, representing colours.
#' @param depth The colour depth of the bitmap image. The image will be composed
#' of `2^depth` indexed colours.
#' 
#' `depth` can also be a `character` string "HAM6" or "HAM8"
#' representing special Amiga display modes (see details).
#' @param interleaved A `logical` value, indicating whether the bitmap needs to be
#' interleaved. An interleaved bitmap image stores each consecutive bitmap layer per
#' horizontal scanline.
#' @param indexing A function that accepts two arguments: `x` (a grDevices
#' `raster` object); `length.out`, a numeric value indicating the
#' desired size of the palette (i.e., the number of colours). It should return
#' a matrix with numeric palette indices (ranging from 1 up to the number of
#' colours in the palette). The result should have an attribute named `palette' that
#' contains the colours that correspond with the index numbers. The result should
#' also carry an attribute with the name `transparent', with a single numeric value
#' representing which colour in the palette should be treated as transparent (or
#' `NA` when no transparency is required). By default the
#' function [index.colours()] is used. You are free to provide
#' a customised version of this function (see examples).
#' @returns The bitmap is returned as a `vector` of `logical` values.
#' The `logical` values reflect the bits for each bitplane. The palette used
#' for the bitmap is returned as attribute to the `vector`. There will also be
#' an attribute called `transparent'. This will hold a numeric index corresponding
#' with the colour in the palette that will be treated as transparent. It will be
#' `NA` when transparency is not used.
#' 
#' @rdname rasterToBitmap
#' @name rasterToBitmap
#' @examples
#' ## first: Let's make a raster out of the 'volcano' data, which we can use in the example:
#' volcano.raster <- as.raster(t(matrix(terrain.colors(1 + diff(range(volcano)))[volcano -
#'   min(volcano) + 1], nrow(volcano))))
#' 
#' ## convert the raster into binary (logical) bitmap data:
#' volcano.bm <- rasterToBitmap(volcano.raster)
#' 
#' ## The palette for the indexed colours of the generated bitmap is returned as
#' ## attribute. There is no transparency is the image:
#' attributes(volcano.bm)
#' 
#' ## We can also include a custom function for colour quantisation. Let's include
#' ## some dithering:
#' volcano.dither <- rasterToBitmap(volcano.raster,
#'                                  indexing = function(x, length.out) {
#'                                    index.colours(x, length.out,
#'                                                  dither = "floyd-steinberg")
#'                                  })
#'
#' ## You can also use a custom indexing function to force a specified palette,
#' ## in this case black and white:
#' volcano.bw <- rasterToBitmap(volcano.raster,
#'                              indexing = function(x, length.out) {
#'                                index.colours(x, length.out,
#'                                              palette = c("black", "white"),
#'                                              dither = "floyd-steinberg")
#'                              })
#' 
#' ## Make a bitmap using a special display mode (HAM6):
#' volcano.HAM <- rasterToBitmap(volcano.raster, "HAM6")
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
rasterToBitmap <- function(x, depth = 3, interleaved = TRUE, indexing = index.colours) {
  special.mode <- "none"
  if (depth %in% c("HAM6", "HAM8")) {
    special.mode <- depth
    depth <- ifelse(depth == "HAM6", 6, 8)
  }
  depth <- round(depth[[1]])
  if (depth < 1) stop("Bitmap depth should be at least 1.")
  interleaved <- interleaved[[1]]
  if (!is.logical(interleaved)) stop("Interleaved should be a logical value.")
  if (!inherits(indexing, "function")) stop("'indexing' should be a function")
  if (!all(c("x", "length.out") %in% names(formals(indexing)))) stop("Function 'indexing' should require arguments 'x' and 'length.out'.")
  x <- as.matrix(x)
  x <- indexing(x = x, length.out = ifelse(special.mode %in% c("HAM6", "HAM8"),
                                           special.mode,
                                           2^depth))
  palette <- attributes(x)$palette
  transparent <- attributes(x)$transparent
  x <- .indexToBitmap(x, depth, interleaved)
  attributes(x) <- list(palette = palette, transparent = transparent)
  return (x)
}

#' Quantisation of colours and indexing a grDevices raster image
#'
#' Converts an image represented by a grDevices `raster` object into a
#' matrix containing numeric indices of a quantised colour palette.
#'
#' Determines the optimal limited palette by clustering colours in an image
#' with [stats::kmeans()]. The result of the optimisation routine
#' will depend on the randomly chosen cluster centres by this algorithm. This
#' means that the result may slightly differ for each call to this function. If
#' you want reproducible results, you may want to reset the random seed
#' ([set.seed()]) before each call to this function.
#' 
#' @param x A raster object ([grDevices::as.raster()]), or a `matrix`
#' containing `character` strings representing colours. `x` can also
#' be a `list` of such matrices or rasters. All elements of this list should
#' have identical dimensions. An overall palette will be generated for elements in the
#' list.
#' @param length.out A `numeric` value indicating the number of desired
#' colours in the indexed palette.
#' 
#' It can also be a `character` string indicating which special
#' Amiga display mode should be used when indexing colours.
#' \sQuote{`HAM6`} and \sQuote{`HAM8`} are supported.
#' See [rasterToBitmap()] for more details on these
#' special modes.
#' @param palette A vector of no more than `length.out` colours, to be used
#' for the bitmap image. When missing or set to `NULL`, a palette will be
#' generated based on the provided colours in raster `x`. In that case,
#' [stats::kmeans()] is used on the hue, saturation, brightness and
#' alpha values of the colours in `x` for clustering the colours. The cluster
#' centres will be used as palette colours.
#' @param background On the Amiga, indexed images could not be semi-transparent.
#' Only a single colour could be designated as being fully transparent. The
#' ``background`' argument should contain a background colour with which
#' semi-transparent colours should be mixed, before colour quantisation. It is
#' white by default.
#' @param dither Dither the output image using the algorithm specified here.
#' See the usage section for possible options. By default no dithering ("`none`")
#' is applied. See [dither()] for more details.
#' @param colour.depth A `character` string indicating the colour depth to be used.
#' Can be either "`12 bit`" (default, standard on an Amiga with original chipset),
#' or "`24 bit`".
#' 
#' This argument is overruled when `length.out` is set to \dQuote{`HAM6`}
#' or \dQuote{`HAM8`}. In that case the colour depth linked to that special mode
#' is used (12 bit for HAM6, 24 bit for HAM8).
#' @param ... Arguments that are passed onto [stats::kmeans()] (see
#' `palette` argument).
#' @returns Returns a `matrix` with the same dimensions as `x` containing
#' `numeric` index values. The corresponding palette is returned as attribute,
#' as well as the index value for the fully transparent colour in the palette.
#' When `x` is a `list` a `list` of matrices is returned.
#' 
#' @rdname index.colours
#' @name index.colours
#' @examples
#' ## first: Let's make a raster out of the 'volcano' data, which we can use in the example:
#' volcano.raster <- as.raster(t(matrix(terrain.colors(1 + diff(range(volcano)))[volcano -
#'   min(volcano) + 1], nrow(volcano))))
#'
#' ## This will create an image of the original raster using an indexed palette:
#' volcano.index <- index.colours(volcano.raster)
#' 
#' ## The index values can be converted back into colours, using the palette:
#' volcano.index <- as.raster(apply(volcano.index, 2,
#'                                  function(x) attributes(volcano.index)$palette[x]))
#' 
#' ## Create an indexed image using dithering
#' volcano.dith <- index.colours(volcano.raster, dither = "floyd-steinberg")
#' volcano.dith <- as.raster(apply(volcano.dith, 2,
#'                                 function(x) attributes(volcano.dith)$palette[x]))
#' 
#' ## plot the images side by side for comparison
#' par(mfcol = c(1, 3))
#' plot(volcano.raster, interpolate = FALSE)
#' plot(volcano.index, interpolate = FALSE)
#' plot(volcano.dith, interpolate = FALSE)
#' @family colour.quantisation.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
index.colours <- function(x, length.out = 8, palette = NULL, background = "#FFFFFF",
                          dither = c("none", "floyd-steinberg", "JJN", "stucki", "atkinson", "burkse", "sierra", "two-row-sierra", "sierra-lite"),
                          colour.depth = c("12 bit", "24 bit"), ...) {
  special.mode <- "none"
  x.is.list <- is.list(x)
  list.length <- 1
  if (x.is.list) list.length <- length(x)
  if (x.is.list) x <- lapply(x, as.matrix) else x <- as.matrix(x)
  if (!all(.is.colour(c(unlist(x))))) stop("x should be a matrix of colours or a grDevices raster object.")
  if (length.out %in% c("HAM6", "HAM8")) {
    special.mode <- length.out
    length.out <- ifelse(length.out == "HAM6", 16, 64)
    ## overrule the colour.depth argument when HAM6 or HAM8
    colour.depth <- ifelse(special.mode == "HAM6", "12 bit", "24 bit")
  } else {
    length.out <- round(length.out[[1]])
    if (length.out < 2) stop("length.out should be 2 or more.")
  }
  if (!is.null(palette) && !all(.is.colour(palette))) stop("palette should consist of colours.")
  if (!is.null(palette) && length(palette) < 2) stop("palette should consist of at least 2 colours.")
  background <- background[[1]]
  if (!.is.colour(background)) stop("background is not a valid colour.")
  
  colour.depth <- match.arg(colour.depth)
  if (colour.depth != "12 bit" && special.mode == "HAM6") stop("HAM6 required 12 bit colour depth")
  if (colour.depth != "24 bit" && special.mode == "HAM8") stop("HAM8 required 24 bit colour depth")
  dither <- match.arg(dither)

  background <- grDevices::col2rgb(background)

  if (x.is.list) {
    c.dim <- do.call(rbind, lapply(x, dim))
    if (any(!apply(c.dim, 2, function(y) all(y == y[[1]]))))
      stop("The dimensions of all elements in x should be equal")
    c.dim <- c.dim[1,]
    x <- unlist(x)
  } else {
    c.dim <- dim(x)
  }

  col.vals <- grDevices::col2rgb(x, TRUE)
  if (special.mode %in% c("HAM6", "HAM8")) col.vals.rgb <- col.vals
  alpha <- col.vals[4,]
  col.vals <- col.vals[-4,]
  col.vals <- (col.vals*rbind(alpha, alpha, alpha) +
                 rep(background, ncol(col.vals))*(255 - rbind(alpha, alpha, alpha)))/255
  col.vals <- grDevices::rgb2hsv(col.vals)
  col.vals[,alpha == 0] <- grDevices::rgb2hsv(background)
  alpha[alpha > 0] <- 255
  x <- apply(rbind(col.vals, alpha/255), 2,
             function(y) grDevices::hsv(y[1], y[2], y[3], y[4]))

  x <- array(x, c(c.dim, list.length))
  x <- lapply(1:list.length, function(y) as.raster(x[,,y]))
  col.vals <- rbind(col.vals, 1 - as.numeric(alpha == 0))
  current.unique.length <- length(unique(c(unlist(x))))
  current.total.length <- length(unlist(x))
  result <- NULL
  transparent <- NA
  if (is.null(palette)) {
    if (current.total.length <= length.out || current.unique.length < length.out) {
      palette <- rep("#000000", length.out)
      palette[1:current.unique.length] <- unique(c(unlist(x)))
      transparent <- which(substr(palette, 8, 9) == "00")[1]
      result <- lapply(x, function(y) apply(y, 2, match, table = palette))
    } else {
      if (special.mode %in% c("HAM6", "HAM8")) {
        col.diff <- array(col.vals.rgb[-4,], c(3, c.dim, list.length))
        col.diff <- c(apply(col.diff, 4, function(z) {
          z <- (z[,,-1] - z[,,-dim(z)[3]])^2
          z <- apply(z, c(2, 3), function(z2) {
            z2[which(z2 == max(z2))[[1]]] <- 0
            prod(1 + z2)/(256*256)
          })
          z <- cbind(rep(0, nrow(z)), z)
          z
        }))
        ## include information on where the image changes a lot in R, G and B value
        col.vals <- rbind(col.vals, col.diff)
        palette <- stats::kmeans(as.matrix(t(col.vals)), length.out, ...)
      } else {
        palette <- stats::kmeans(as.matrix(t(col.vals)), length.out, ...)
        result <- palette$cluster
        result <- array(palette$cluster, c(c.dim, list.length))
        result <- lapply(1:list.length, function(y) result[,,y])
      }
      transparent <- which(palette$centers[,4] == 0)[1]
      palette <- apply(palette$centers, 1, function(x) grDevices::hsv(x[1], x[2], x[3], x[4]))
    }
    # sort colours such that the most frequently occuring colours are listed first
    freqs   <- table(factor(unlist(result), as.character(1:length.out)))
    ord     <- order(-freqs)
    rnk     <- rank(-freqs, ties.method = "first")
    palette <- as.vector(palette[ord])
    transparent <- as.vector(rnk[transparent])
    if (!is.null(result)) {
      result <- lapply(result, function(y) as.vector(rnk)[y])
      result <- lapply(result, matrix, nrow = c.dim)
    }
  } else {
    palette <- grDevices::col2rgb(palette, TRUE)
    transparent <- which(palette[4,] == 0)[1]
    palette[4,palette[4,] > 0] <- 255
    palette <- grDevices::rgb(palette[1,], palette[2,], palette[3,], palette[4,], maxColorValue = 255)
    result <- lapply(x, function(y) apply(y, 2, match, table = palette))
  }

  if (dither != "none" || special.mode %in% c("HAM6", "HAM8")) { ## dithering should also be called in case of HAM modes
    if (x.is.list) {
      result <- lapply(x, function(y) dither(y, method = dither, palette = palette, mode = special.mode))
    } else {
      result <- dither(x[[1]], method = dither, palette = palette, mode = special.mode)
    }
  } else if (!x.is.list) {
    result <- result[[1]]
  }
  palette <- suppressWarnings(amigaRawToColour(colourToAmigaRaw(palette, "24 bit", "3"), colour.depth, "3"))
  attributes(result)[["palette"]] <- as.vector(palette)
  attributes(result)[["transparent"]] <- transparent
  return(result)
}

#' Image dithering
#'
#' Dither is an intentional form of noise applied to an image to avoid colour
#' banding when reducing the amount of colours in that image. This function
#' applies dithering to a grDevices `raster` image.
#'
#' The approaches implemented here all use error diffusion to achieve dithering.
#' Each pixel is scanned (from top to bottom, from left to right), where the actual
#' colour is sampled and compared with the closest matching colour in the palette.
#' The error (the differences between the actual and used colour) is distributed over
#' the surrounding pixels. The only difference between the methods implemented here
#' is the way the error is distributed. The algorithm itself is identical. For more
#' details consult the listed references.
#'
#' Which method results in the best quality image will depend on the original image
#' and the palette colours used for dithering, but is also a matter of taste. Note
#' that the dithering algorithm is relatively slow and is provided in this package
#' for your convenience. As it is not in the main scope of this package you should
#' use dedicated software for faster/better results.
#' @param x Original image data that needs to be dithered. Should be a raster object
#' ([grDevices::as.raster()]), or a matrix of `character` string
#' representing colours.
#' @param method A `character` string indicating which dithering method should
#' be applied. See usage section for all possible options (Note that the "JJN" is
#' the Jarvis, Judice, and Ninke algorithm). Default is "`none`", meaning that
#' no dithering is applied.
#' @param palette A palette to which the image should be dithered. It should be a
#' `vector` of `character` strings representing colours.
#' @param mode A `character` string indicating whether a special
#' Amiga display mode should be used when dithering. By default
#' \sQuote{`none`} is used (no special mode). In addition,
#' \sQuote{`HAM6`} and \sQuote{`HAM8`} are supported.
#' See [rasterToBitmap()] for more details.
#' @param ... Currently ignored.
#' @returns Returns a `matrix` with the same dimensions as `x` containing
#' `numeric` index values. The corresponding palette is returned as attribute,
#' as well as the index value for the fully transparent colour in the palette.
#' 
#' @rdname dither
#' @name dither
#' @aliases dither.raster
#' @examples
#' ## first: Let's make a raster out of the 'volcano' data, which we can use in the example:
#' volcano.raster <- as.raster(t(matrix(terrain.colors(1 + diff(range(volcano)))[volcano -
#'   min(volcano) + 1], nrow(volcano))))
#'
#' ## let's dither the image, using a predefined two colour palette:
#' volcano.dither <- dither(volcano.raster,
#'                          method = "floyd-steinberg",
#'                          palette = c("yellow", "green"))
#' 
#' ## Convert the indices back into a raster object, such that we can plot it:
#' volcano.dither <- as.raster(apply(volcano.dither, 2, function(x) c("yellow", "green")[x]))
#' par(mfcol = c(1, 2))
#' plot(volcano.raster, interpolate = FALSE)
#' plot(volcano.dither, interpolate = FALSE)
#' 
#' ## results will get better when a better matching colour palette is used.
#' ## for that purpose use the function 'index.colours'.
#' @references R.W. Floyd, L. Steinberg, *An adaptive algorithm for spatial grey scale*. Proceedings of the Society of Information Display 17, 75-77 (1976).
#' @references J. F. Jarvis, C. N. Judice, and W. H. Ninke, *A survey of techniques for the display of continuous tone pictures on bilevel displays*. Computer Graphics and Image Processing, 5:1:13-40 (1976).
#' @references <https://en.wikipedia.org/wiki/Floyd-Steinberg_dithering>
#' @references <https://tannerhelland.com/4660/dithering-eleven-algorithms-source-code/>
#' @family colour.quantisation.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
dither.raster <- function(x, method = c("none", "floyd-steinberg", "JJN", "stucki", "atkinson", "burkse", "sierra", "two-row-sierra", "sierra-lite"), palette, mode = c("none", "HAM6", "HAM8"), ...) {
  mode <- match.arg(mode)
  if (!all(.is.colour(c(x)))) stop("x should be a matrix of colours or a grDevices raster object.")
  if (!is.null(palette) && !all(.is.colour(palette))) stop("palette should consist of colours.")
  if (!is.null(palette) && length(palette) < 2) stop("palette should consist of at least 2 colours.")

  x <- matrix(x, nrow = dim(x))
  method <- match.arg(method)
  c.dim <- dim(x)
  
  ## create an array with width, height, r, g, b and alpha as separate dimensions
  x <- grDevices::col2rgb(x, TRUE)
  x <- lapply(split(x, row(x)), matrix, nrow = c.dim)
  x <- array(c(x[[1]], x[[2]], x[[3]], x[[4]]), dim = c(rev(c.dim), 4))
  
  pal.rgb <- col2rgb(palette, TRUE)

  result <- matrix(rep(NA, prod(c.dim)), nrow = c.dim)
  if (method == "floyd-steinberg") {
    e2 <- matrix(c(0, 3, -16, 5, 7, 1), nrow = c(2, 3))/16
    ir2 <- 0:1
    jr2 <- -1:1
  } else if (method == "JJN") {
    e2 <- matrix(c(0, 3, 1, 0, 5, 3, -48, 7, 5, 7, 5, 3, 5, 3, 1), nrow = c(3, 5))/48
    ir2 <- 0:2
    jr2 <- -2:2
  } else if (method == "stucki") {
    e2 <- matrix(c(0, 2, 1, 0, 4, 2, -42, 8, 4, 8, 4, 2, 4, 2, 1), nrow = c(3, 5))/42
    ir2 <- 0:2
    jr2 <- -2:2
  } else if (method == "atkinson") {
    e2 <- matrix(c(0, 1, 0, -8, 1, 1, 1, 1, 0, 1, 0, 0), nrow = c(3, 4))/8
    ir2 <- 0:2
    jr2 <- -1:2
  } else if (method == "burkse") {
    e2 <- matrix(c(0, 2, 0, 4, -32, 8, 8, 4, 4, 2), nrow = c(2, 5))/32
    ir2 <- 0:1
    jr2 <- -2:2
  } else if (method == "sierra") {
    e2 <- matrix(c(0, 2, 0, 0, 4, 2, -32, 5, 3, 5, 4, 2, 3, 2, 0), nrow = c(3, 5))/32
    ir2 <- 0:2
    jr2 <- -2:2
  } else if (method == "two-row-sierra") {
    e2 <- matrix(c(0, 1, 0, 2, -16, 3, 4, 2, 3, 1), nrow = c(2, 5))/16
    ir2 <- 0:1
    jr2 <- -2:2
  } else if (method == "sierra-lite") {
    e2 <- matrix(c(0, 1, -4, 1, 2, 0), nrow = c(2, 3))/4
    ir2 <- 0:1
    jr2 <- -1:1
  }
  if (method == "none" & !(mode %in% c("HAM6", "HAM8"))) {
    result <- apply(x, 2, function(a) {
      res <- apply(a, 1, function(b) {
        dst <- sqrt(colSums((pal.rgb - b)^2))
        which(dst == min(dst))[[1]]
      })
      res
    })
    result <- t(result)
  } else {
    color_multi   <- ifelse(mode == "HAM8", 255/63, 17)
    for(j in 1:dim(x)[2]) {
      if (mode %in% c("HAM6", "HAM8")) prev <- c(grDevices::col2rgb(palette[1]))
      for(i in 1:dim(x)[1]) {
        ## find the closest matching colour in the palette compared to the
        ## current pixel. This is the colour where the Euclidean distance
        ## in RGBA space is smallest compared to the actual colour:
        if (mode %in% c("HAM6", "HAM8")) {
          dst <- apply(pal.rgb, 2, function(z) {
            dst <- abs(x[i, j, ] - z)
            sqrt(sum(dst^2))
          })
          dst.diff <- abs(x[i, j, 1:3] - prev)
          control.flag = which(dst.diff == max(dst.diff))[[1]]
          dst.diff[control.flag] <- 0
          dst.diff <- sqrt(sum(dst.diff^2))
          if (all(dst.diff < dst)) {
            idx <- round(x[i,j,][control.flag]/ifelse(mode == "HAM6", 17, (255/63)))
            prev[control.flag] <- color_multi*idx
            control.flag <- c(2, 3, 1)[control.flag]
          } else {
            control.flag <- 0
            ## Possible improvement for future versions:
            ## When multiple colours in the palette match best with the current
            ## pixel, now the first matching colour is selected.
            ## it is better to also look ahead to see if the pixel to the
            ## right matches best with this colour in the palette.
            idx <- which(dst == min(dst))[[1]] - 1
            prev <- c(grDevices::col2rgb(palette[[idx + 1]]))
          }
          result[j, i] <- idx + bitwShiftL(control.flag, ifelse(mode == "HAM6", 4, 6))
        } else {
          dst <- sqrt(colSums((pal.rgb - x[i, j, ])^2))
          result[j, i] <- which(dst == min(dst))[[1]]
        }
        if (method != "none" && !(j == dim(x)[[2]] && i == dim(x)[[1]])) {
          if (mode %in% c("HAM6", "HAM8")) {
            ## seems to create a slight horizontal stripes artifact in HAM modes.
            ## See if this can be avoided
            P            <- c(prev, 255)
          } else {
            P            <- pal.rgb[,result[j, i]]
          }
          
          ## calculate the error (difference) between the actual colour and the colour
          ## from the palette:
          e            <- x[i, j, ] - P
          
          ## get the proper row and column indices for the error distribution matrix.
          ## This is necessary when we are close to the edge of the image:
          sel.i        <- i + ir2
          ir           <- ir2[sel.i %in% 1:dim(x)[1]]
          sel.j        <- j + jr2
          jr           <- jr2[sel.j %in% 1:dim(x)[2]]
          
          ## Distribute the error (e) over the surrounding pixels using the error
          ## distribution matrix (e2) for the selected method:
          repl <- x[i + ir, j + jr, ]
          repl <- repl + (e2 %o% e)[ir - min(ir2) + 1, jr -min(jr2) + 1,]
          ## Put some constrains on the error:
          repl[repl < 0] <- 0
          repl[repl > 255] <- 255
          x[i + ir, j + jr, ] <- repl
        }
      }
    }
    if (mode %in% c("HAM6", "HAM8")) result <- result + 1
  }
  return(result)
}

#' @rdname dither
#' @name dither
#' @aliases dither.matrix
#' @export
dither.matrix <- function(x, method = c("none", "floyd-steinberg", "JJN", "stucki", "atkinson", "burkse", "sierra", "two-row-sierra", "sierra-lite"), palette, mode = c("none", "HAM6", "HAM8"), ...) {
  dither.raster(grDevices::as.raster(x), method, palette, mode, ...)
}

#' (De)compress 8-bit continuous signals.
#'
#' Use a lossy delta-Fibonacci (de)compression to continuous 8-bit signals.
#' This algorithm was used to compress 8-bit audio wave data on the Amiga.
#'
#' This form of compression is lossy, meaning that information and quality will get lost.
#' 8-bit audio is normally stored as an 8-bit signed value representing the amplitude
#' at specific time intervals. The delta-Fibonacci compression instead stores the
#' difference between two time intervals (delta) as a 4-bit index. This index in turn
#' represents a value from the Fibonacci series (hence the algorithm name). The compression
#' stores small delta values accurately, but large delta values less accurately.
#' As each sample is stored as a 4-bit value instead of an 8-bit value, the amount of
#' data is reduced with almost 50\% (the exact compression ratio is (4 + n)/(2n)).
#' 
#' The algorithm was first described by Steve Hayes and was used in 8SVX audio stored in
#' the Interchange File Format (IFF). The quality loss is considerable (especially
#' when the audio contained many large deltas) and was even in
#' the time it was developed (1985) not used much. The function is provided here for
#' the sake of completeness. The implementation here only compresses 8-bit data, as
#' for 16-bit data the quality loss will be more considerable.
#' @param x A `vector` of `raw` data that needs to be (de)compressed.
#' @param ... Currently ignored.
#' @returns Returns a `vector` of the resulting (de)compressed `raw` data.
#' @rdname deltaFibonacciCompress
#' @name deltaFibonacciCompress
#' @examples
#' ## Let's get an audio wave from the ProTrackR package, which we
#' ## can use in this example:
#' buzz     <- ProTrackR::PTSample(ProTrackR::mod.intro, 1)
#' 
#' ## Let's convert it into raw data, such that we can compress it:
#' buzz.raw <- as.integer(ProTrackR::waveform(buzz) - 128) |>
#'   bitwAnd(0xFF) |>
#'   as.raw()
#' 
#' ## Let's compress it:
#' buzz.compress <- deltaFibonacciCompress(buzz.raw)
#' 
#' ## Look the new data uses less memory:
#' length(buzz.compress)/length(buzz.raw)
#' 
#' ## The compression was lossy, which we can examine by decompressing the
#' ## sample again:
#' buzz.decompress <- deltaFibonacciDecompress(buzz.compress)
#' 
#' ## And turn the raw data into numeric data:
#' buzz.decompress <-
#'   ifelse(buzz.decompress > 0x7f, as.integer(buzz.decompress) - 256L,
#'          as.integer(buzz.decompress))
#' 
#' ## Plot the original wave in black, the decompressed wave in blue
#' ## and the error in red (difference between the original and decompressed
#' ## wave). The error is actually very small here.
#' plot(ProTrackR::waveform(buzz) - 128, type = "l")
#' lines(buzz.decompress, col = "blue")
#' buzz.error <- ProTrackR::waveform(buzz) - 128 - buzz.decompress
#' lines(buzz.error, col = "red")
#' 
#' ## this can also be visualised by plotting the orignal wave data against
#' ## the decompressed data (and observe a very good correlation):
#' plot(ProTrackR::waveform(buzz) - 128, buzz.decompress)
#' 
#' ## Let's do the same with a sample of a snare drum, which has larger
#' ## delta values:
#' snare.drum <- ProTrackR::PTSample(ProTrackR::mod.intro, 2)
#' 
#' ## Let's convert it into raw data, such that we can compress it:
#' snare.raw <- as.integer(ProTrackR::waveform(snare.drum) - 128L) |>
#'   bitwAnd(0xFF) |>
#'   as.raw()
#' 
#' ## Let's compress it:
#' snare.compress <- deltaFibonacciCompress(snare.raw)
#' 
#' ## Decompress the sample:
#' snare.decompress <- deltaFibonacciDecompress(snare.compress)
#' 
#' ## And turn the raw data into numeric data:
#' snare.decompress <-
#'   ifelse(snare.decompress > 0x7f, as.integer(snare.decompress) - 256L,
#'          as.integer(snare.decompress))
#' 
#' ## Now if we make the same comparison as before, we note that the
#' ## error in the decompressed wave is much larger than in the previous
#' ## case (red line):
#' plot(ProTrackR::waveform(snare.drum) - 128, type = "l")
#' lines(snare.decompress, col = "blue")
#' snare.error <- ProTrackR::waveform(snare.drum) - 128 - snare.decompress
#' lines(snare.error, col = "red")
#' 
#' ## this can also be visualised by plotting the orignal wave data against
#' ## the decompressed data (and observe a nice but not perfect correlation):
#' plot(ProTrackR::waveform(snare.drum) - 128, snare.decompress)
#' @references <https://en.wikipedia.org/wiki/Delta_encoding>
#' @references <http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/node02D6.html>
#' @author Pepijn de Vries
#' @export
deltaFibonacciCompress <- function(x, ...) {
  ## Steve Hayes' Fibonacci Delta sound compression technique
  ## algorithm results in slightly different compression than
  ## achieved with Audiomaster IV. But the total error is smaller
  ## in this implementation
  result <- c(raw(1), x[1])
  x <- .rawToAmigaInt(x, 8, TRUE)
  fibonacci <- rev(c(-34,-21,-13,-8,-5,-3,-2,-1,0,1,2,3,5,8,13,21))
  fib.deltas <- rep(NA, length(x))
  new.wave <- fib.deltas
  new.wave[1] <- x[1]
  for (i in 1:length(x)) {
    target.value <- x[i] + 128
    achieved.value <- (c(x[1], new.wave)[i] + fibonacci + 128) %% 256
    value.dif <- target.value - achieved.value
    fib.deltas[i] <- 16 - (which(abs(value.dif) == min(abs(value.dif)))[[1]])
    
    new.wave[i] <- c(x[1], new.wave)[i] + fibonacci[16 - fib.deltas[i]]
  }
  fib.even <- as.raw(fib.deltas)
  fib.odd <- fib.even[seq(1, length(fib.even), by = 2)]
  if (length(fib.even) == 1)
    fib.even <- as.raw(8)
  else
    fib.even <- fib.even[seq(2, length(fib.even), by = 2)]
  if (length(fib.odd) < length(fib.even)) fib.odd <- c(fib.odd, as.raw(8))
  if (length(fib.odd) > length(fib.even)) fib.even <- c(fib.even, as.raw(8))
  result <- c(result,
              .amigaIntToRaw(.rawToAmigaInt(fib.odd)*0x10) | fib.even)
  return(result)
}

#' @rdname deltaFibonacciCompress
#' @name deltaFibonacciDecompress
#' @export
deltaFibonacciDecompress <- function(x, ...) {
  ## from http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/node02D6.html
  ## Unpack Fibonacci-delta encoded data from n byte source buffer into
  ## 2*(n-2) byte dest buffer. Source buffer has a pad byte, an 8-bit
  ## initial value, followed by n-2 bytes comprising 2*(n-2) 4-bit
  ## encoded samples.
  ## second byte indicates the base value:
  base.val <- x[2]
  ## first byte is a padding byte; second is already stored; skip them:
  x <- x[-1:-2]
  fibonacci <- c(-34,-21,-13,-8,-5,-3,-2,-1,0,1,2,3,5,8,13,21)
  result <- c(rbind(.hiNybble(x), .loNybble(x)))
  result <- fibonacci[result + 1]
  result <- .rawToAmigaInt(base.val, 8, TRUE) + cumsum(result)
  result <- ((result + 128) %% 256) - 128
  return(.amigaIntToRaw(result, 8, TRUE))
}

.is.colour <- function(x)
{
  unlist(lapply(x, function(y) {
    res <- try(col2rgb(y), silent = TRUE)
    return(!inherits(res, "try-error"))
  }))
}

.inverseViewPort <- function(display.mode, monitor) {
  adm <- AmigaFFH::amiga_display_modes
  camg <- adm$DISPLAY_MODE_ID[adm$DISPLAY_MODE == display.mode][[1]] |
    AmigaFFH::amiga_monitors$CODE[AmigaFFH::amiga_monitors$MONITOR_ID == monitor][[1]]
  new("IFFChunk", chunk.type = "CAMG", chunk.data = list(camg))
}

.amigaViewPortModes <- function(x) {
  MONITOR_ID_MASK <- as.raw(c(0xff, 0xff, 0x10, 0x00))
  UPPER_MASK      <- as.raw(c(0xff, 0xff, 0x00, 0x00))
  EXTENDED_MODE   <- as.raw(c(0x00, 0x00, 0x10, 0x00))
  SPRITES	        <- as.raw(c(0x00, 0x00, 0x40, 0x00))
  GENLOCK_AUDIO   <- as.raw(c(0x00, 0x00, 0x01, 0x00))
  GENLOCK_VIDEO   <- as.raw(c(0x00, 0x00, 0x00, 0x02))
  VP_HIDE         <- as.raw(c(0x00, 0x00, 0x20, 0x00))
  
  # Knock bad bits out of old-style CAMGs modes before checking
  # availability. (some ILBM CAMG's have these bits set in old 1.3 modes,
  # and should not) If not an extended monitor ID, or if marked as
  # extended but missing upper 16 bits, screen out inappropriate bits now.
  ## see: http://amigadev.elowar.com/read/ADCD_2.1/AmigaMail_Vol2_guide/node00FD.html
  
  if (!any(as.logical(x & MONITOR_ID_MASK)) ||
      (any(as.logical(x & EXTENDED_MODE)) && !any(as.logical(x & UPPER_MASK)))) {
    if (any(as.logical(x & (EXTENDED_MODE|SPRITES|GENLOCK_AUDIO|GENLOCK_VIDEO|VP_HIDE)))) {
      warning("CAMG / display mode contains old style bad bits, I will knock them out...")
      x <- x & !(EXTENDED_MODE|SPRITES|GENLOCK_AUDIO|GENLOCK_VIDEO|VP_HIDE)
    }
  }
  
  # Check for bogus CAMG like some brushes have, with junk in
  # upper word and extended bit NOT set not set in lower word.
  if (any(as.logical(x & UPPER_MASK)) && !(any(as.logical(x & EXTENDED_MODE)))) {
    warning("CAMG / display mode contains bogus bits, I will use the simplest display mode possible...")
    x <- as.raw(c(0x00, 0x00, 0x00, 0x00))
  }
  
  monitors <- AmigaFFH::amiga_monitors
  display_modes <- AmigaFFH::amiga_display_modes
  display_modes <- AmigaFFH::amiga_display_modes
  mon <- unlist(lapply(monitors$CODE, function(y) all(y == (x & MONITOR_ID_MASK))))
  mon <- monitors$MONITOR_ID[mon]
  if (length(mon) > 0 && mon %in% c("STANDARD", "DEFAULT_MONITOR_ID", "NTSC_MONITOR_ID", "PAL_MONITOR_ID")) {
    x <- x & !MONITOR_ID_MASK
  }
  disp <- unlist(lapply(display_modes$DISPLAY_MODE_ID, function(y) all(y == x)))
  disp <- display_modes$DISPLAY_MODE[disp]
  
  if (length(mon) > 0 && length(disp) == 0 && mon %in% c("EURO36_MONITOR_ID", "SUPER72_MONITOR_ID")) {
    x <- x & !MONITOR_ID_MASK
    disp <- unlist(lapply(display_modes$DISPLAY_MODE_ID, function(y) all(y == x)))
    disp <- display_modes$DISPLAY_MODE[disp]
  }
  
  return(list(monitor = mon, display.mode = disp))
}

.display.properties <- function(display.mode, monitor) {
  attribs <- list(
    is.lace         = grepl("LACE", display.mode),
    is.super        = grepl("SUPER", display.mode),
    is.hires        = grepl("HIRES", display.mode),
    is.HAM          = grepl("HAM", display.mode),
    is.extralores   = grepl("EXTRALORES", display.mode),
    is.noflicker    = grepl("FF", display.mode),
    is.scan.doubled = grepl("DBL", display.mode),
    is.productivity = grepl("PRODUCT", display.mode),
    is.halfbright   = grepl("HB|HALFBRITE", display.mode)
  )
  
  aspect.x <- ifelse(monitor == "A2024_MONITOR_ID", 14,
                     ifelse(monitor == "SUPER72_MONITOR_ID", 34, 44))
  width <- ifelse(monitor == "A2024_MONITOR_ID", 1024,
                  ifelse(monitor == "SUPER72_MONITOR_ID", 200, 320))
  if (length(attribs$is.extralores) > 0 && attribs$is.extralores) {
    width <- width/2
    aspect.x <- aspect.x*2
  }
  if (length(attribs$is.hires) > 0 && attribs$is.hires) {
    width <- width*2
    aspect.x <- aspect.x/2
  }
  if (length(attribs$is.super) > 0 && attribs$is.super) {
    width <- width*4
    aspect.x <- aspect.x/4
  }
  if (length(attribs$is.productivity) > 0 && attribs$is.productivity) {
    width <- width*2
    aspect.x <- aspect.x/2
  }
  
  height <- ifelse(monitor == "A2024_MONITOR_ID", 800,
                   ifelse(monitor == "SUPER72_MONITOR_ID", 300,
                          ifelse(monitor == "VGA_MONITOR_ID", 480,
                                 ifelse(monitor == "EURO72_MONITOR_ID", 400,
                                        ifelse(grepl("PAL", monitor), 256, 200)))))
  aspect.y <- ifelse(monitor == "A2024_MONITOR_ID", 11,
                     ifelse(monitor == "SUPER72_MONITOR_ID", 40,
                            ifelse(monitor == "VGA_MONITOR_ID", 22,
                                   ifelse(monitor == "EURO72_MONITOR_ID", 22,
                                          ifelse(grepl("NTSC", monitor), 52, 44)))))
  
  if (length(attribs$is.lace) > 0 && attribs$is.lace) {
    height <- height*2
    aspect.y <- aspect.y/2
  }
  if (length(attribs$is.scan.doubled) > 0 && attribs$is.scan.doubled) {
    height <- height*2
    aspect.y <- aspect.y/2
  }
  if (length(attribs$is.no.flicker) > 0 && attribs$is.no.flicker) {
    height <- height/2
    aspect.y <- aspect.y*2
  }
  attribs[["screenwidth"]] <- width
  attribs[["screenheight"]] <- height
  attribs[["aspect.x"]] <- aspect.x
  attribs[["aspect.y"]] <- aspect.y
  return(attribs)
}

.indexToBitmap <- function(x, depth, interleaved) {
  ## x should be a matrix of palette indices
  x <- cbind(x, matrix(1, ncol = -ncol(x)%%16, nrow = nrow(x)))
  x.dim <- dim(x)
  x <- .rawToBitmap(.amigaIntToRaw(c(x) - 1, 32, FALSE), TRUE, FALSE)
  sq <- c(outer(31:(32 - depth), seq(1, length(x), by = 32), "+"))
  x <- as.logical(x[sq])
  rm(sq)
  # dimensions are bitplane, height, width
  x <- array(x, dim = c(depth, x.dim))
  if (interleaved == TRUE) {
    ## rearrange dimensions to height, bitplane, width (non-interleaved.)
    x <- c(aperm(x, c(3, 1, 2)))
  } else {
    ## rearrange dimensions to bitplane, height, width (interleaved.)
    x <- c(aperm(x, c(3, 2, 1)))
  }
}

#' Get an Amiga timeval struct value from raw data
#'
#' Some Amiga applications use a timeval struct (see references) to represent a
#' time span in seconds. This function coerces raw data to such a numeric time span.
#'
#' Timeval is a structure (struct) as specified in device/timer.h on the Amiga (see
#' references). It represents a timespan in seconds. This function retrieves the
#' numeric value from `raw` data. Amongst others, the timeval struct was used
#' in the system-configuration file (see [SysConfig]) to specify key repeat speed,
#' key repeat delay and mouse double click speed. Use `as.raw` for the inverse
#' of this function and get the original raw data.
#' @rdname timeval
#' @name timeval
#' @param x a `vector` of `raw` data that need to be converted into
#' Amiga timeval structs.
#' @returns Returns a `numeric` `vector` of a timespan in seconds. It is
#' represented as an S3 AmigaTimeVal class.
#' @examples
#' ## First four raw values represent seconds, the latter four microseconds:
#' temp <- timeval(as.raw(c(0, 0, 0, 1, 0, 0, 0, 1)))
#' print(temp)
#' 
#' ## You can use 'as.raw' to get the original raw data again:
#' as.raw(temp)
#' @author Pepijn de Vries
#' @references <http://amigadev.elowar.com/read/ADCD_2.1/Includes_and_Autodocs_2._guide/node0053.html>
#' @export
timeval <- function(x) {
  ## get timeval struct from raw data
  if ((length(x) %% 8) != 0) stop("The length of x should be a multiple of 8.")
  if (typeof(x) != "raw") stop("x should be of type 'raw'.")
  x <- matrix(.rawToAmigaInt(x, 32, FALSE), ncol = 2, byrow = TRUE)
  result <- apply(x, 1, function(y) y[[1]] + y[[2]]/1e6)
  class(result) <- "AmigaTimeVal"
  return(result)
}

#' @name as.raw
#' @rdname as.raw
#' @export
as.raw.AmigaTimeVal <- function(x, ...) {
  ## convert a timval (time interval in seconds) to raw timeval struct
  if (!inherits(x, "AmigaTimeVal")) stop("x should be of S3 class AmigaTimeVal.")
  secs   <- floor(x)
  micros <- round((x - secs)*1e6)
  secs[secs >= 2^32] <-  (2^32) - 1
  micros[micros >= 2^32] <-  (2^32) - 1
  .amigaIntToRaw(c(rbind(secs, micros)), 32, FALSE)
}

#' @export
print.AmigaTimeVal <- function(x, ...) {
  invisible(lapply(x, function(y) cat(sprintf("%f [s] Amiga timeval struct\n", y, ...))))
}

.read.amigaData <- function(dat, n.bytes, signed, par.names) {
  ## read numeric and raw data from amiga raw input
  ## dat = vector of raw data
  ## n.bytes = vector of lengths of bytes to be read from input data. Negative values are negated and indicate that raw data should be read as is
  ## signed = vector of logicals. Indicate whether the values read from data are signed (TRUE) or unsigned (FALSE)
  ## par.names = parameter names for the data read from the input data
  n.bytes <- round(n.bytes)
  offset <- 0
  result <- mapply(function(n.b, sgnd) {
    res <- NULL
    if (n.b > 0) {
      res <- .rawToAmigaInt(dat[offset + (1:n.b)], n.b*8, sgnd)
    } else if (n.b < 0){
      n.b <- -n.b
      res <- dat[offset + (1:n.b)]
    }
    offset <<- offset + n.b
    return(res)
  }, n.b = n.bytes, sgnd = signed, SIMPLIFY = FALSE)
  names(result) <- par.names
  result
}

.write.amigaData <- function(lst, n.bytes, signed, par.names) {
  ## inverse function for .read.amigaData
  ## first make sure list is in correct order:
  lst <- lst[par.names]
  if (any(n.bytes > 0)) {
    lst[n.bytes > 0] <- lapply(1:sum(n.bytes > 0), function(y) {
      .amigaIntToRaw(lst[n.bytes > 0][[y]],
                                 8*n.bytes[n.bytes > 0][y],
                                 signed[n.bytes > 0][y])
    })
  }
  lst <- unlist(lst)
  names(lst) <- NULL
  return(lst)
}

.match.factor <- function(lst, element.name, vals, levs) {
  result <- match(lst[[element.name]], vals)
  if (is.na(result)) stop(sprintf("Unknown %s.", element.name))
  result <- factor(levs[result], levs)
  return(result)
}

.match.factor.inv <- function(lst, element.name, vals, levs) {
  result <- vals[which(levs %in% lst[[element.name]])]
  if (length(result) == 0) stop(sprintf("Unknown level for %s.", element.name))
  if (length(result) > 1) stop(sprintf("Only a single value for %s is allowed.", element.name))
  return(result)
}

.bitwOrAll <- function(x) {
  while (length(x) > 1) {
    x <- c(bitwOr(x[1], x[2]), x[-1:-2])
  }
  return(x)
}

.match.multi.factor <- function(lst, element.name, vals, levs) {
  result <- levs[bitwAnd(lst[[element.name]], vals) == vals]
  result <- factor(result, levs)
}

.match.multi.factor.inv <- function(lst, element.name, vals, levs) {
  result <- vals[which(levs %in% lst[[element.name]])]
  while (length(result) > 1) {
    result <- c(bitwOr(result[1], result[2]), result[-1:-2])
  }
  return(result)
}

.loNybble <-
  function(raw_dat)
    ## function that gets the value [0,16] of the 4 low bits of a raw byte
  {
    if (!inherits(raw_dat, "raw")) stop ("Only raw data is accepted as input")
    return(as.integer(raw_dat)%%16)
  }

.hiNybble <-
  function(raw_dat)
    ## function that gets the value [0,16] of the 4 high bits of a raw byte
  {
    if (!inherits(raw_dat, "raw")) stop ("Only raw data is accepted as input")
    return(as.integer(as.integer(raw_dat)/16))
  }

.rawToCharNull <- function(raw_dat) {
  result <- ""
  if (length(raw_dat) < 3) try(result <- (rawToChar(raw_dat)), silent = TRUE) else
  {
    result    <- raw_dat
    runlength <- rle(result)$lengths
    if (length(runlength) > 2)
    {
      rel_range <- (runlength[1] + 1):(length(result) - runlength[length(runlength)])
      if (result[[1]] != raw(1)) rel_range <- 1:rel_range[length(rel_range)]
      if (result[[length(result)]] != raw(1)) rel_range <- rel_range[1]:length(result)
      result[rel_range][result[rel_range] == as.raw(0x00)] <- as.raw(0x20)
      result <- result[result != raw(1)]
    }
    try(result <- rawToChar(result), silent = TRUE)
    if (inherits(result, "raw")) result <- ""
  }
  return(result)
}

.read.generic <- function(file) {
  ## If the file size can be determined from 'file', that size
  ## will be read. Other wise, the file will be read in 5 kB chunks.
  size <- 5*1024
  close_later <- FALSE
  if (inherits(file, "virtual_path")) {
    if (requireNamespace("adfExplorer")) {
      file <- adfExplorer::adf_file_con(file)
      close_later <- TRUE
    } else {
      stop("Package `adfExplorer` is required to read from virtual disks.")
    }
  }
  if (inherits(file, "character")) {
    close_later <- TRUE
    size <- file.size(file)
    file <- file(file, "rb")
  }
  if (inherits(file, "connection")) {
    con_info <- summary(file)
    if (con_info$`can read` != "yes" || con_info$text != "binary") stop("file is not a connection from which binary data can be read...")
  }
  result <- NULL
  repeat {
    l1 <- length(result)
    if (inherits(file, "adf_file_con")) {
      if (requireNamespace("adfExplorer")) {
        result <- c(result, adfExplorer::readBin(file, "raw", size))
      } else {
        stop("Package `adfExplorer` is required to read from virtual disks.")
      }
    } else {
      result <- c(result, readBin(file, "raw", size))
    }
    l2 <- length(result)
    if ((l2 - l1) < size) break
  }
  if (close_later) close(file)
  return(result)
}

.write.generic <- function(x, file, ...) {
  raw.dat <- as.raw(x, ...)
  if (inherits(file, "virtual_path")) {
    if (requireNamespace("adfExplorer")) {
      con <- adfExplorer::adf_file_con(file, writable = TRUE)
    } else {
      stop("Please install package `adfExplorer` to write to virtual disks.")
    }
  }
  if (inherits(file, "character")) con <- file(file, "wb")
  if (inherits(file, "connection")) {
    con_info <- summary(con)
    if (con_info$`can write` != "yes" || con_info$text != "binary") stop("file is not a connection to which binary data can be written...")
    con <- file
  }
  if (inherits(con, "adf_file_con")) {
    if (requireNamespace("adfExplorer")) {
      adfExplorer::writeBin(raw.dat, con, endian = "big")
    } else {
      stop("Please install package `adfExplorer` to write to virtual disks.")
    }
  } else {
    writeBin(raw.dat, con, endian = "big")
  }
  if (inherits(file, c("character", "virtual_path"))) return(close(con))
}
