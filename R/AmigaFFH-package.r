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
#' @section Supported File Formats:
#' This package supports a number of file formats. The ProTracker module file format requires sophisticated interpretation
#' and a dedicated package ([ProTrackR][ProTrackR::ProTrackR-package]) is developed for that purpose.
#' 
#' The following formats are supported by this package (to some extend):
#' 
#'  * [Amiga Basic][AmigaFFH::AmigaBasic] binary encode scripts and
#'    [Amiga Basic shapes][AmigaFFH::AmigaBasicShape] which
#'    were used by such scripts to display specific graphics.
#'  * Bitmap Font (.font). Originally fonts were stored in separate files on
#'    the Amiga. An overarching '.font' file contained generic information,
#'    amongst others the specific pixel heights that were available for a font.
#'    The actual font bitmap images were stored in separate files. There
#'    was a file available for each individual font height. For more details
#'    see [AmigaBitmapFont()] and [AmigaBitmapFontSet()].
#'  * Interchange File Format (IFF). This file format is actually a container for a wide variety of
#'    data flavours. Of which the following are supported:
#'     * 8SVX (8-bit sampled voices (i.e., audio)). There are no major restrictions in this package's
#'       implementation.
#'     * ANIM (animations). Not all display modes are supported as per ILBM. Furthermore, the vertical byterun
#'       encoding for the animation frames is the only encoding currently supported.
#'     * ILBM (InterLeaved BitMap images). Specific display modes (such as 'extra halfbrite') can in
#'       some cases be decoded, but encoding for these modes may not (yet) be supported.
#'  * For more details see [IFFChunk()], [interpretIFFChunk()], [read.iff()] and
#'    [write.iff()].
#'  * Hardware sprites. This format follows the hardware structure for displaying sprites on the screen.
#'    It is usually not used
#'    as a file format as such, but it can be found embedded in some files (for instance the mouse pointer
#'    is embedded as a
#'    hardware sprite in the 'system-configuration' file). For more details see [hardwareSprite()].
#'  * System-configuration. A file that was stored in the 'devs' directory of a system disk.
#'    As the file name suggests, it holds many of the systems configurations. See [SysConfig]
#'    for more details.
#'  * Workbench icons (.info). Icons (i.e., graphical representation of files and directories
#'    on the Amiga) were stored as separate files with the extension '.info'. See
#'    [AmigaIcon()] for more details.
#' 
#' In future versions of this package more file types may be added to this list.
#' @section In Addition...:
#' Several helper functions are also exported by this package. This will give you access
#' to older compression techniques, such as the run length encoding ([packBitmap()])
#' and delta Fibonacci compression ([deltaFibonacciCompress()]). But also other
#' techniques that will help in converting modern files into classic file formats and vice versa.
#' Such as for instance the function to [dither()] full colour images to a limited
#' colour palette.
#' @references
#' Documentation on several Amiga File types:
#' <http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/>
#' @keywords internal
"_PACKAGE"
NULL
