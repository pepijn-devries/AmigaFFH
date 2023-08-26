#' The Amiga File Format Handler package
#'
#' The Amiga File Format Handler package (AmigaFFH) is designed to interpret file formats that were native
#' to Commodore Amiga machines.
#'
#' In combination with the adfExplorer package this package can be used to interpret older file formats that
#' were native to the Commodore Amiga. The focus of this package will be on the older system (Amiga OS <= 3.0).
#' This will allow you to analyse and interpret these files in the scripted environment of R.
#' 
#' Note that all functions and methods in this package are implemented as scripted source code and may not run
#' very fast.
#' 
#' @section Supported File Formats:
#' This package supports a number of file formats. The ProTracker module file format requires sophisticated interpretation
#' and a dedicated package (\link[ProTrackR:ProTrackR-package]{ProTrackR}) is developed for that purpose.
#' 
#' The following formats are supported by this package (to some extend):
#' 
#' \itemize{
#'   \item{
#'     \link[AmigaFFH:AmigaBasic]{Amiga Basic} binary encode scripts and \link[AmigaFFH:AmigaBasicShape]{Amiga Basic shapes} which
#'     were used by such scripts to display specific graphics.
#'   }
#'   \item{
#'     Bitmap Font (*.font). Originally fonts were stored in separate files on
#'     the Amiga. An overarching *.font file contained generic information,
#'     amongst others the specific pixel heights that were available for a font.
#'     The actual font bitmap images were stored in separate files. There
#'     was a file available for each individual font height. For more details
#'     see \code{\link{AmigaBitmapFont}} and \code{\link{AmigaBitmapFontSet}}.
#'   }
#'   \item{
#'     Interchange File Format (IFF). This file format is actually a container for a wide variety of data flavours. Of which the
#'     following are supported:
#'     \itemize{
#'       \item{
#'         8SVX (8-bit sampled voices (i.e., audio)). There are no major restrictions in this package's implementation.
#'       }
#'       \item{
#'         ANIM (animations). Not all display modes are supported as per ILBM. Furthermore, the vertical byterun
#'         encoding for the animation frames is the only encoding currently supported.
#'       }
#'       \item{
#'         ILBM (InterLeaved BitMap images). Specific display modes (such as \sQuote{extra halfbrite}) can in
#'         some cases be decoded, but encoding for these modes may not (yet) be supported.
#'       }
#'     }
#'     For more details see \code{\link{IFFChunk}}, \code{\link{interpretIFFChunk}}, \code{\link{read.iff}} and
#'     \code{\link{write.iff}}.
#'   }
#'   \item{
#'     Hardware sprites. This format follows the hardware structure for displaying sprites on the screen. It is usually not used
#'     as a file format as such, but it can be found embedded in some files (for instance the mouse pointer is embeded as a
#'     hardware sprite in the \sQuote{system-configuration} file). For more details see \code{\link{hardwareSprite}}.
#'   }
#'   \item{
#'     System-configuration. A file that was stored in the \sQuote{devs} directory of a system disk.
#'     As the file name suggests, it holds many of the systems configurations. See \link{SysConfig}
#'     for more details.
#'   }
#'   \item{
#'     Workbench icons (*.info). Icons (i.e., graphical representation of files and directories
#'     on the Amiga) were stored as separate files with the extension *.info. See
#'     \code{\link{AmigaIcon}} for more details.
#'   }
#' }
#' 
#' In future versions of this package more file types may be added to this list.
#' 
#' @section In Addition...:
#' Several helper functions are also exported by this package. This will give you access
#' to older compression techniques, such as the run length encoding (\code{\link{packBitmap}})
#' and delta Fibonacci compression (\code{\link{deltaFibonacciCompress}}). But also other
#' techniques that will help in converting modern files into classic file formats and vice versa.
#' Such as for instance the function to \code{\link{dither}} full colour images to a limited
#' colour palette.
#' @references
#' Documentation on several Amiga File types:
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/}
#' @keywords internal
"_PACKAGE"
NULL