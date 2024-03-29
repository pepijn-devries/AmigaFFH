% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hardwareSprites.r
\name{rasterToHWSprite}
\alias{rasterToHWSprite}
\title{Convert a raster object into an hardwareSprite object}
\usage{
rasterToHWSprite(x, indexing = index.colours)
}
\arguments{
\item{x}{A \code{\link[=grDevices]{grDevices()}} raster object (\code{\link[grDevices:as.raster]{grDevices::as.raster()}})
that needs to be converted into a \code{\link[=hardwareSprite]{hardwareSprite()}} class object.
Note that a \code{\link[=hardwareSprite]{hardwareSprite()}} has a maximum width of 16 pixels.
When \code{x} is wider, it will be cropped.}

\item{indexing}{A function that accepts two arguments: \code{x} (a grDevices
\code{raster} object); \code{length.out}, a numeric value indicating the
desired size of the palette (i.e., the number of colours). It should return
a matrix with numeric palette indices (ranging from 1 up to the number of
colours in the palette). The result should have an attribute named \verb{palette' that contains the colours that correspond with the index numbers. The result should also carry an attribute with the name }transparent', with a single numeric value
representing which colour in the palette should be treated as transparent (or
\code{NA} when no transparency is required). By default the
function \code{\link[=index.colours]{index.colours()}} is used.}
}
\value{
Returns a \code{\link[=hardwareSprite]{hardwareSprite()}} class object based on \code{x}
}
\description{
Convert a grDevices raster object into an Amiga hardwareSprite class object.
}
\details{
A \code{\link[=grDevices]{grDevices()}} raster image can be converted into a
\code{\link[=hardwareSprite]{hardwareSprite()}} class object with this function. For this purpose
the any true-colour image will be converted to an indexed palette with 4 colours.
The Amiga hardware sprite will reserve one of the colours as transparent. Thos function
will use fully transparent colours in the original image (i.e., the alpha level equals 0)
for this purpose. Or when the image has no fully transparent colours, it will use the
most frequently occuring colour (at least when the default \code{indexing} function
is used).
}
\examples{
\dontrun{
## first create a raster object that can be used as input
## (making sure that the background is transparent):
rst <- as.raster(simpleSysConfig()$PointerMatrix, "#AAAAAA00")

## now turn it into a hardware sprite:
spr <- rasterToHWSprite(rst)

## and plot it as a check:
plot(spr)
}
}
\seealso{
Other raster.operations: 
\code{\link{AmigaBitmapFont}},
\code{\link{as.raster.AmigaBasicShape}()},
\code{\link{bitmapToRaster}()},
\code{\link{dither}()},
\code{\link{index.colours}()},
\code{\link{rasterToAmigaBasicShape}()},
\code{\link{rasterToAmigaBitmapFont}()},
\code{\link{rasterToBitmap}()},
\code{\link{rasterToIFF}()}

Other HWSprite.operations: 
\code{\link{rawToHWSprite}()}
}
\author{
Pepijn de Vries
}
\concept{HWSprite.operations}
\concept{raster.operations}
