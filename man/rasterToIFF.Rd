% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iff.ilbm.r
\name{rasterToIFF}
\alias{rasterToIFF}
\title{Convert a grDevices raster image into an IFF formated bitmap image}
\usage{
rasterToIFF(
  x,
  display.mode = as.character(AmigaFFH::amiga_display_modes$DISPLAY_MODE),
  monitor = as.character(AmigaFFH::amiga_monitors$MONITOR_ID),
  anim.options,
  ...
)
}
\arguments{
\item{x}{A raster object created with \code{\link[grDevices:as.raster]{grDevices::as.raster()}} which
needs to be converted into an IFF formated bitmap image. It is also possible to let \code{x} be
a matrix of \code{character}s, representing colours.}

\item{display.mode}{Specify the Amiga display mode that should be used.
See \code{\link[=amiga_display_modes]{amiga_display_modes()}} for all possible options.
"\code{LORES_KEY}" is used by default, this is the lowest resolution
possible on the Amiga.}

\item{monitor}{The Amiga monitor on which the needs to be displayed.
See \code{\link[=amiga_monitors]{amiga_monitors()}} for more details and posible options.
By default "\code{DEFAULT_MONITOR_ID}" is used.}

\item{anim.options}{Currently ignored. This argument will potentitally be implemented
in future versions of this package. Currently, animations are always encoded
with the "ByteVerticalCompression" in this package (when \code{x} is a list of
\code{raster} objects).}

\item{...}{Arguments passed on to \code{\link[=rasterToBitmap]{rasterToBitmap()}}.}
}
\value{
Returns an \code{\link[=IFFChunk]{IFFChunk()}} object holding an Interleaved
Bitmap (ILBM) image based on \code{x}.
}
\description{
Convert grDevices raster images (\code{\link[grDevices:as.raster]{grDevices::as.raster()}})
into a formal \code{\link[=IFFChunk]{IFFChunk()}} object, as an interleaved bitmap (ILBM)
image.
}
\details{
Convert any modern image into a interleaved bitmap (image) conform
Interchange File Format (IFF) specifications. If your original image
is in true colour (i.e., a 24 bit colour depth) it will be converted
into a bitmap image with an indexed palette.
}
\examples{
\dontrun{
## first: Let's make a raster out of the 'volcano' data, which we can use in the example:
volcano.raster <- as.raster(t(matrix(terrain.colors(1 + diff(range(volcano)))[volcano -
  min(volcano) + 1], nrow(volcano))))

## Turning the raster into an IFFChunk object is easy:
volcano.iff <- rasterToIFF(volcano.raster)

## This object can be saved as an IFF file using write.iff

## in special modes HAM6 and HAM 8 higher quality images
## can be obtained. See 'rasterToBitmap' for more info on the
## special HAM modes.
volcano.ham <- rasterToIFF(volcano.raster, "HAM_KEY", depth = "HAM8")

## The result can be further improved by applying dithering
volcano.ham.dither <- rasterToIFF(volcano.raster, "HAM_KEY", depth = "HAM8",
  indexing = function(x, length.out) {
    index.colours(x, length.out, dither = "JJN", iter.max = 20)
  })
}
}
\seealso{
Other iff.operations: 
\code{\link{IFFChunk-class}},
\code{\link{WaveToIFF}()},
\code{\link{as.raster.AmigaBasicShape}()},
\code{\link{getIFFChunk}()},
\code{\link{interpretIFFChunk}()},
\code{\link{rawToIFFChunk}()},
\code{\link{read.iff}()},
\code{\link{write.iff}()}

Other raster.operations: 
\code{\link{AmigaBitmapFont}},
\code{\link{as.raster.AmigaBasicShape}()},
\code{\link{bitmapToRaster}()},
\code{\link{dither}()},
\code{\link{index.colours}()},
\code{\link{rasterToAmigaBasicShape}()},
\code{\link{rasterToAmigaBitmapFont}()},
\code{\link{rasterToBitmap}()},
\code{\link{rasterToHWSprite}()}
}
\author{
Pepijn de Vries
}
\concept{iff.operations}
\concept{raster.operations}
