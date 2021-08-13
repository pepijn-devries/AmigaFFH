## Desired Y position + vertical-offset of the display window = 25 + 44 = 69
## Desired X position + horizontal-offset of display window ~= 94 + 64 = 158

.validity.HWSprite <- function(object)
{
  if (length(object@VStart)       != 1) stop("VStart should have a length of 1")
  if (length(object@HStart)       != 1) stop("HStart should have a length of 1")
  if (length(object@VStop)        != 1) stop("VStop should have a length of 1")
  if (length(object@control.bits) != 8) stop("control.bits should have a length of 8")
  if (length(object@end.of.data)  != 4) stop("end.of.data should have a length of 4")
  if (length(object@colours)      != 3) stop("colours should have a length of 3")
  if (object@VStop  < object@VStart)    stop("VStop should be equal to or greater than VStart.")
  if (object@HStart < 0)                stop("HStart can't be negative")
  if (object@VStart < 0)                stop("VStart can't be negative")
  if (object@VStop  < 0)                stop("VStop can't be negative")
  if (!all(.is.colour(object@colours))) stop("colours should represent colours")
  if (!all(object@end.of.data == raw(4))) warning("Extended sprites are currently not supported")
  if (length(object@bitmap) != (object@VStop - object@VStart)*4) stop("bitmap should have a length of (VStop-VStart)*4")
  return(T)
}

#' The hardwareSprite class
#'
#' An S4 class object that represent graphical objects known as
#' hardware sprites on the Commodore Amiga.
#'
#' Amiga hardware supported sprites, which are graphical objects that
#' could be moved around the display and indepently from eachother.
#' Basic sprites were 16 pixels wide and any number of pixels high and
#' were composed of four colours, of which one is transparent.
#' 
#' More complex sprites could be formed by linking separate sprites
#' together. That way, sprites could become wider, or be composed of
#' more colours. Such extended sprites are currently not supported
#' by this package.
#' 
#' A well known example of hardware sprite on the Commodore Amiga is
#' the mouse pointer.
#' 
#' This object simply holds the basic information belonging to
#' hardware sprite. Use \code{\link{as.raster}} to convert it to
#' a raster which is a more usuful graphical element in R.
#'
#' @slot VStart The vertical starting position of a sprite.
#' @slot HStart The horizontal starting position of a sprite.
#' @slot VStop The vertical stopping position of a sprite. The
#' height of a sprite should be given by \code{VStart} - \code{VStop}.
#' @slot control.bits 8 \code{logical} values used for extending
#' sprites. The values are stored in this objects but extending
#' sprites is not (yet) supported.
#' @slot bitmap Interleaved bitmap data containing information on
#' the pixel colour numbers of the sprite.
#' @slot colours A vector of the 3 colours used for the sprite.
#' @slot end.of.data Sprite data can be followed by another sprite. It is terminated
#' with two WORDS equalling zero (\code{raw(4)}). Repeated sprite data is currently not
#' supported.
#' @references \url{http://amigadev.elowar.com/read/ADCD_2.1/Hardware_Manual_guide/node00AE.html}
#' @name hardwareSprite-class
#' @rdname hardwareSprite-class
#' @aliases hardwareSprite
#' @examples
#' ## This generates a sprite of a single line (16x1 pixels) with an empty bitmap:
#' new("hardwareSprite")
#' 
#' ## This generates a sprite of a single line (16x1 pixels) where
#' ## the bitmap contains some coloured pixels:
#' new("hardwareSprite", bitmap = as.raw(c(0x01,0x02,0x03,0x04)))
#' 
#' ## This generates a sprite of 16x16 pixels:
#' new("hardwareSprite",
#'     VStop = 16,
#'     bitmap = as.raw(sample.int(255, 64, replace = TRUE)))
#' @exportClass hardwareSprite
#' @author Pepijn de Vries
setClass("hardwareSprite",
         representation(VStart        = "numeric",
                        HStart        = "numeric",
                        VStop         = "numeric",
                        control.bits  = "logical",
                        bitmap        = "raw",
                        end.of.data   = "raw",
                        colours       = "character"),
         prototype(VStart        = 0,
                   HStart        = 0,
                   VStop         = 1,
                   control.bits  = rep(F, 8),
                   bitmap        = raw(4),
                   end.of.data   = raw(4),
                   colours       = c("#000000", "#888888", "#FFFFFF")),
         validity = .validity.HWSprite)

setGeneric("rawToHWSprite", function(x, col) standardGeneric("rawToHWSprite"))

#' Convert raw data into an Amiga hardware sprite
#'
#' Convert \code{raw} data structured conform a Commodore Amiga hardware
#' sprite (see references) into a \code{\link{hardwareSprite}} object.
#'
#' Information to set up a hardware sprite is stored as \code{raw} data
#' on Commodore Amigas. This method can be used to convert this data
#' into a \code{\link{hardwareSprite}} object. This object can in turn
#' be converted with \code{\link{as.raster}} such that it can be plotted in R.
#'
#' @docType methods
#' @rdname rawToHWSprite
#' @name rawToHWSprite
#' @aliases rawToHWSprite,raw,missing-method
#' @param x \code{raw} data structured as an Amiga hardware sprite
#' (see references).
#' @param col A \code{vector} of colours (\code{character}) to be used
#' for the hardware sprite. Specify the three visible colours for the
#' sprite. When missing some default colours (grayscale) will be used.
#' The colours have to be provided separately as they are usually not stored
#' together with the hardware sprite data.
#' @return Returns a \code{\link{hardwareSprite}} object based on the provided raw data
#' @references \url{http://amigadev.elowar.com/read/ADCD_2.1/Hardware_Manual_guide/node00B9.html}
#' @examples
#' ## Let's generate a 16x16 sprite with a random bitmap:
#' dat <- as.raw(c(0x00, 0x00, 0x10, 0x00,
#'               sample.int(255, 64, replace = TRUE),
#'               0x00, 0x00, 0x00, 0x00))
#' ## make it a hardware sprite object:
#' spr <- rawToHWSprite(dat)
#' ## and plot it:
#' plot(spr, interpolate = FALSE)
#' 
#' ## with some imagination when can make
#' ## a more structured image:
#' dat <- as.raw(c(0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xff, 0xf8,
#'                 0x7f, 0x80, 0x80, 0x70, 0x7f, 0x00, 0xbe, 0xe0,
#'                 0x7e, 0x00, 0x85, 0xc0, 0x7d, 0x80, 0x82, 0x40,
#'                 0x6b, 0xc0, 0x95, 0xa0, 0x57, 0xe0, 0xa8, 0xd0,
#'                 0x2f, 0xf0, 0xd1, 0x68, 0x4f, 0xf8, 0xb0, 0x34,
#'                 0x07, 0xfc, 0xf8, 0x5a, 0x03, 0xfe, 0xe4, 0x0d,
#'                 0x01, 0xfc, 0xc2, 0x12, 0x00, 0xf8, 0x81, 0x04,
#'                 0x00, 0x70, 0x00, 0x88, 0x00, 0x20, 0x00, 0x50,
#'                 0x00, 0x00, 0x00, 0x20, 0x00, 0x00, 0x00, 0x00))
#' spr <- rawToHWSprite(dat, c("#EE4444", "#000000", "#EEEECC"))
#' plot(spr, interpolate = FALSE)
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
setMethod("rawToHWSprite", c("raw", "missing"), function(x, col) {
  result <- methods::new("hardwareSprite")
  result@HStart <- adfExplorer::rawToAmigaInt(x[1], 8, F)
  result@VStart <- adfExplorer::rawToAmigaInt(x[2], 8, F)
  result@VStop  <- adfExplorer::rawToAmigaInt(x[3], 8, F)
  if (result@VStop == 0) result@VStop <- 16 ## This appears to be the case for the mouse pointer. Check if this is always the case
  result@control.bits <- as.logical(adfExplorer::rawToBitmap(x[4], invert.longs = F))
  vlen <- result@VStop - result@VStart
  result@bitmap <- x[4 + 1:(vlen*4)]
  offset <- vlen*4 + 4
  eod    <- x[offset + 1:4]
  result@end.of.data <- eod
  while (!all(eod == raw(4))) {
    result@end.of.data <- c(result@end.of.data, eod)
    offset <- offset + 4
    eod    <- x[offset + 1:4]
  }
  return(result)
})

#' @rdname rawToHWSprite
#' @aliases rawToHWSprite,raw,character-method
#' @export
setMethod("rawToHWSprite", c("raw", "character"), function(x, col) {
  result <- rawToHWSprite(x)
  result@colours <- col
  return(result)
})

#' @rdname as.raster
#' @name as.raster
#' @aliases as.raster,hardwareSprite-method
#' @export
setMethod("as.raster", c("hardwareSprite"), function(x, background = "#AAAAAA", ...) {
  cols <- c(background, x@colours)
  return(bitmapToRaster(x@bitmap, 16, length(x@bitmap)*8/(2*16), 2, cols)) # assume 2 bitplanes
})

#' @rdname plot
#' @name plot
#' @export
plot.hardwareSprite <- function(x, y, ...) {
  graphics::plot(as.raster(x), ...)
}

#' @rdname as.raw
#' @name as.raw
#' @aliases as.raw,hardwareSprite-method
#' @export
setMethod("as.raw", "hardwareSprite", function(x) {
  result <- c(
    adfExplorer::amigaIntToRaw(c(x@HStart, x@VStart, x@VStop), 8, F),
    adfExplorer::bitmapToRaw(x@control.bits, invert.longs = F, invert.bytes = F),
    x@bitmap,
    x@end.of.data
  )

  return(result)
})

#' @export
print.hardwareSprite <- function(x, ...) {
  cat(sprintf("A %i row high hardware sprite sprite", x@VStop - x@VStart))
}

setMethod("show", "hardwareSprite", function(object){
  print(object)
})
