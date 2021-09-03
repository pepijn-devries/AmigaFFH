#' A table of display modes on the Amiga and corresponding \code{raw} codes
#'
#' A table of display modes on the Amiga and corresponding \code{raw} codes
#' representing these modes.
#'
#' This table contains most display modes that were available on the Amiga.
#' It also contains \code{raw} codes that were used to represent these modes.
#' The table also contains the hardware monitors that could display the specific
#' modes, and the minimal chip set that was required for the display mode.
#' This data is used to interpret \code{\link{IFFChunk}} objects of type
#' `CAMG'. It is also used to interpret ILBM images and creating IFF files from
#' raster images.
#'
#' @docType data
#' @name amiga_display_modes
#' @format A \code{data.frame} with 4 columns:
#' \itemize{
#'   \item{The column named `DISPLAY_MODE': a \code{factor} reflecting
#'   the display mode}
#'   \item{The column named `DISPLAY_MODE_ID': A \code{list} containing a \code{vector}
#'   of 4 \code{raw} values as used by the Amiga to reflect specific display modes.
#'   These raw values are usually also stored with bitmap images in the Interchange
#'   File Format in a \code{\link{IFFChunk}} called `CAMG'.}
#'   \item{The column named `MONITOR_ID': A \code{character} string identifying
#'   the monitor that could display the specific mode.}
#'   \item{The column named `CHIPSET': a \code{factor} identifying the minimal
#'   chip set that was required to display the specific mode. OCS is the original
#'   chip set; ECS is the Enhanced Chip Set. AGA is the Advanced Graphics Architecture
#'   chip set (in some countries known as just Advanced Architecture). AGA could
#'   also display OCS and ECS modes, ECS could also display OCS modes, OCS could only
#'   display OCS modes.}
#' }
#' @references \url{https://wiki.amigaos.net/wiki/Display_Database#ModeID_Identifiers}
#' @references \url{http://amigadev.elowar.com/read/ADCD_2.1/AmigaMail_Vol2_guide/node00FD.html}
#' @examples
#' data("amiga_display_modes")
NULL

#' A list of special display modes
#'
#' A list of special display modes on the Amiga and corresponding \code{raw} keys.
#'
#' This table show specific special display modes and to which Amiga monitors they
#' relate. The \code{raw} codes can be used to interpret specific display modes
#' as listed in \code{\link{amiga_display_modes}}. This information is used to
#' interpret \code{\link{IFFChunk}} objects of type `CAMG'. It is also used to
#' interpret ILBM images and creating IFF files from raster images.
#'
#' @docType data
#' @name amiga_display_keys
#' @format a \code{data.frame} with 2 columns:
#' \itemize{
#'   \item{The column named `mode': a \code{factor} reflecting a display mode, monitor or bitwise mask}
#'   \item{The column named `code': vector of 4 \code{raw} values as used by the Amiga to reflect specific display modes}
#' }
#' @references \url{https://wiki.amigaos.net/wiki/Display_Database#ModeID_Identifiers}
#' @references \url{http://amigadev.elowar.com/read/ADCD_2.1/AmigaMail_Vol2_guide/node00FD.html}
#' @examples
#' data("amiga_display_keys")
NULL

#' A list of Amiga monitors
#'
#' This table lists Amiga monitors and corresponding \code{raw} codes that represent
#' these monitors.
#'
#' This table contains monitors that were compatible with the Amiga.
#' It also contains \code{raw} codes that were used to represent them.
#' This data is used to interpret \code{\link{IFFChunk}} objects of type
#' `CAMG'. It is also used to interpret ILBM images and creating IFF files from
#' raster images.
#'
#' @docType data
#' @name amiga_monitors
#' @format A \code{data.frame} with 2 columns:
#' \itemize{
#'   \item{The column named `MONITOR_ID': a \code{factor} representing an Amiga monitor}
#'   \item{The column named `CODE': A \code{list} containing a \code{vector} of 4
#'   \code{raw} values as used by the Amiga to represent a specific monitor.}
#' }
#' @references \url{https://wiki.amigaos.net/wiki/Display_Database#ModeID_Identifiers}
#' @examples
#' data("amiga_monitors")
NULL

#' An example file of a bitmap image stored in the Interchange File Format
#'
#' This file is provided to demonstrate the structure of an Interchange File
#' Format and is used in several examples throughout this package.
#'
#' The Interchange File Format stores information compartmentally in separate
#' containers called `chunks'. This file demonstrates how a bitmap image
#' is stored in this format. In addition to the raw bitmap data, the file
#' also contains meta-information on the bitmap dimensions, its colour palette and
#' the display mode that should be used on an Amiga. See also
#' \code{\link{interpretIFFChunk}}, \code{\link{IFFChunk-class}}
#' and the example for \code{\link{bitmapToRaster}}.
#'
#' @docType data
#' @name ilbm8lores.iff
#' @format See \code{\link{IFFChunk-class}} and references for more information
#' about the Interchange File Format.
#' @examples
#' \dontrun{
#' filename <- system.file("ilbm8lores.iff", package = "AmigaFFH")
#' example.iff <- read.iff(filename)
#' 
#' ## show the structure of the IFF file:
#' print(example.iff)
#' }
#' @references \url{https://en.wikipedia.org/wiki/Interchange_File_Format}
#' @references \url{https://wiki.amigaos.net/wiki/A_Quick_Introduction_to_IFF}
NULL

#' 'demo.bas', 'r_logo.shp' and 'ball.shp' as example files for AmigaBasic and AmigaBasicShape objects
#'
#' `demo.bas', `r_logo.shp' and `ball.shp' as example files for \code{\link{AmigaBasic}} and
#' \code{\link{AmigaBasicShape}} objects
#'
#' The `r_logo.shp' and `ball.shp' files are formatted such that they can be read with
#' \code{\link{read.AmigaBasicShape}}. They serve as an example of the \code{\link{AmigaBasicShape}} class, where
#' the first represents a blitter object, and the latter a sprite.
#' 
#' The `demo.bas' file is an example of a binary encoded \link[AmigaFFH:AmigaBasic]{Amiga Basic} script. It can be read with
#' \code{\link{read.AmigaBasic}}. The script demonstrates how the shape files could be used in Amiga Basic.
#'
#' @docType data
#' @aliases demo.bas r_logo.shp ball.shp
#' @name AmigaBasic-files
#' @rdname AmigaBasic-files
#' @format See \code{\link{AmigaBasic}} and \code{\link{AmigaBasicShape}} for more information
#' about the format.
#' @examples
#' \dontrun{
#' read.AmigaBasic(system.file("demo.bas", package = "AmigaFFH"))
#' read.AmigaBasicShape(system.file("ball.shp", package = "AmigaFFH"))
#' read.AmigaBasicShape(system.file("r_logo.shp", package = "AmigaFFH"))
#' }
NULL

#' Commonly used palettes on the Commodore Amiga
#'
#' \code{amiga_palettes} is a named list, where each element represents a commonly
#' used palette on the Commodore Amiga.
#'
#' Some files that contain bitmap images with an indexed palette did not store the
#' palette in the same file. Amiga Workbench icons (\code{\link{AmigaIcon}}) for instance
#' only store the index values of the palette, but not the palette itself.
#' \code{amiga_palettes} therefore provides some commonly used palettes on the Amiga,
#' such that these files can be interpreted.
#'
#' @docType data
#' @name amiga_palettes
#' @format A named list with the following elements:
#' \itemize{
#'   \item{\code{wb.os1}: A \code{vector} of 4 colours that were used as the default
#'   palette of the Workbench on Amiga OS 1.x.}
#'   \item{\code{wb.os2}: A \code{vector} of 8 colours. The first 4 colours are the default
#'   colours of a standard Workbench on Amiga OS 2.x. The latter 4 are additional
#'   colours used by the Workbench expansion MagicWB.}
#'   \item{\code{spr.os1}: A \code{vector} of 3 colours that were used by default
#'   for a mouse pointer sprite on Amiga OS 1.x.}
#'   \item{\code{spr.os2}: A \code{vector} of 3 colours that were used by default
#'   for a mouse pointer sprite on Amiga OS 2.x.}
#' }
#' @examples
#' data("amiga_palettes")
NULL

#' An example object for the AmigaBitmapFontSet class
#'
#' An example object for the \code{\link{AmigaBitmapFontSet}} class used in
#' examples throughout this package. It also contains a nested
#' \code{\link{AmigaBitmapFont}} class objects, which can be obtain by
#' using \code{getAmigaBitmapFont(font_example, 9)}.
#'
#' The \code{font_example} contains a font that was designed as an example
#' for this package. It holds bitmap glyphs for 8 and 9 pixels tall
#' characters.
#'
#' @docType data
#' @name font_example
#' @format \code{font_example} is an \code{\link{AmigaBitmapFontSet}}
#' object. For details see the object class documentation.
#' @family AmigaBitmapFont.operations
#' @examples
#' data("font_example")
NULL
