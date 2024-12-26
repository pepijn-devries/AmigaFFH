#' The S3 SysConfig class
#' 
#' A comprehensive representation of an Amiga system-configuration file.
#' 
#' The system-configuration file is a binary file stored in the `devs' folder
#' of the root of a bootable Amiga DOS device, containing system preferences.
#' It was used in Amiga OS 1.x. Although it could be used in later versions, it was
#' gradually phased out and some settings may not be usable in the later versions.
#' See references below for more details.
#' 
#' Definitions of the system-configuration have file been revised at some points.
#' Revisions are minor and usually targeted at backward compatibility. Here
#' revision V38.2 (released on 16 September 1992) is implemented, which is the
#' latest documented version.
#' 
#' The sytem-configuration file contains settings for the serial and parallel
#' port and the printer. It also contains some settings for the `workbench'
#' which was the Amiga equivalent of what is now mostly known as the computers desktop.
#' Colours for the workbench and the shape of the mouse pointer are also stored
#' in this file. Settings for the mouse and basic screen resolution are also part
#' of the file.
#' 
#' The SysConfig object is a comprehensive representation of the binary
#' system-configuration file. It is a a `list` where the elements have identical
#' names as listed in the documents provided the references. The names are usually
#' self-explanatory, but the referred documents can also be
#' consulted to obtain more detailed information with respect to each of
#' these elements. The mouse pointer is included as a [hardwareSprite()] object
#' in the list. The pointer image can be replaced by a different [hardwareSprite()],
#' but make sure it has an height of 16 pixels.
#' 
#' It is possible to change the values of the list, but not all values may be valid.
#' Note that they will not be fully checked for validity. Invalid values may result in errors
#' when writing to a binary file using [write.SysConfig()], or may simply not
#' work properly on an Amiga or in an emulator.
#' 
#' Use [simpleSysConfig()] for creating a simple SysConfig object which can
#' be modified. Use [read.SysConfig()] to read, and [write.SysConfig()]
#' to write system-configuration files. With [rawToSysConfig()] and
#' [AmigaFFH::as.raw()] SysConfig can be coerced back and forth from and to
#' its raw (binary) form.
#' 
#' @docType class
#' @name SysConfig
#' @rdname SysConfig
#' @family SysConfig.operations
#' @author Pepijn de Vries
#' @references
#' <https://wiki.amigaos.net/wiki/Preferences#Preferences_in_1.3_and_Older_Versions_of_the_OS>
#' <http://amigadev.elowar.com/read/ADCD_2.1/Includes_and_Autodocs_2._guide/node00D5.html>
#' <http://amigadev.elowar.com/read/ADCD_2.1/Includes_and_Autodocs_3._guide/node063B.html>
NULL

#' @rdname plot
#' @name plot
#' @export
plot.SysConfig <- function(x, y, ...) {
  plot(as.raster(x$PointerMatrix, background = x$WBColours[[1]]), ...)
}

#' @export
print.SysConfig <- function(x, ...) {
  cat(sprintf("Amiga system-configuration\nFontHeight:\t%i\nPrinter:\t%s %s\nInterlaced:\t%s",
              x$FontHeight,
              tolower(strsplit(as.character(x$PrinterPort), "_")[[1]][[1]]),
              x$PrinterFilename,
              as.character(x$LaceWB == "LACE")), ...)
}

#' @rdname as.raw
#' @name as.raw
#' @export
as.raw.SysConfig <- function(x, ...) {
  class(x) <- NULL
  x$KeyRptSpeed           <- as.raw(x$KeyRptSpeed)
  x$KeyRptDelay           <- as.raw(x$KeyRptDelay)
  x$DoubleClick           <- as.raw(x$DoubleClick)
  x$PointerMatrix         <- as.raw(x$PointerMatrix)
  x$PointerMatrix[3]      <- raw(1) ## set VStop to 0
  x$WBColours             <- colourToAmigaRaw(x$WBColours, colour.depth = "12 bit", n.bytes = "2")
  x$spriteColours         <- colourToAmigaRaw(x$spriteColours, colour.depth = "12 bit", n.bytes = "2")
  x$PrinterFilename       <- charToRaw(x$PrinterFilename)[1:30]
  x$PrtDevName            <- charToRaw(x$PrtDevName)[1:16]
  x[names(.SysConfigFactors)] <-
    lapply(names(.SysConfigFactors), function(y) {
      .match.factor.inv(x, y, .SysConfigFactors[[y]]$vals, .SysConfigFactors[[y]]$levs)
    })
  x$PrintFlags <- .match.multi.factor.inv(x,
                                          "PrintFlags",
                                          .SysConfigMultiFactors[["PrintFlags"]]$vals,
                                          .SysConfigMultiFactors[["PrintFlags"]]$levs)
  x$SerRWBits <- .bitmapToRaw(x$SerRWBits, FALSE, FALSE)
  x$SerParShk <- .amigaIntToRaw(
    16*.match.factor.inv(x$SerParShk, "SerialParity", 0:4, c("SPARITY_NONE", "SPARITY_EVEN", "SPARITY_ODD",
                                                             "SPARITY_MARK", "SPARITY_SPACE")) +
      .match.factor.inv(x$SerParShk,
                        "HandshakeMode", 0:2, c("SHSHAKE_XON", "SHSHAKE_RTS", "SHSHAKE_NONE")),
    8, FALSE)
  x$SerStopBuf <- .amigaIntToRaw(
    16*(x$SerStopBuf$N.StopBits - 1) +
      .match.factor.inv(x$SerStopBuf,
                        "BufSize", 0:5, c("SBUF_512", "SBUF_1024", "SBUF_2048", "SBUF_4096", "SBUF_8000", "SBUF_16000")),
    8, FALSE)
  x <- .write.amigaData(x, .SysConfigData$byte, .SysConfigData$signed, .SysConfigData$par.names)
  return(x)
}

.SysConfigMultiFactors <- list(
  PrintFlags= data.frame(
    vals = c(0x0001, 0x0002, 0x0004, 0x0008, 0x0000, 0x0010, 0x0020, 0x0040,
             0x0080, 0x0100, 0x0000, 0x0200, 0x0400, 0x0800, 0x1000),
    levs = c("CORRECT_RED", "CORRECT_GREEN", "CORRECT_BLUE", "CENTER_IMAGE",
             "IGNORE_DIMENSIONS", "BOUNDED_DIMENSIONS", "ABSOLUTE_DIMENSIONS",
             "PIXEL_DIMENSIONS", "MULTIPLY_DIMENSIONS", "INTEGER_SCALING",
             "ORDERED_DITHERING", "HALFTONE_DITHERING", "FLOYD_DITHERING",
             "ANTI_ALIAS", "GREY_SCALE2"),
    stringsAsFactors = FALSE
  )
)

.SysConfigFactors <- list(
  PrinterPort  = data.frame(vals = c(0x00, 0x01),
                            levs = c("PARALLEL_PRINTER", "SERIAL_PRINTER"),
                            stringsAsFactors = FALSE),
  BaudRate     = data.frame(vals = 0:7,
                            levs = c("BAUD_110", "BAUD_300", "BAUD_1200", "BAUD_2400",
                                     "BAUD_4800", "BAUD_9600", "BAUD_19200", "BAUD_MIDI"),
                            stringsAsFactors = FALSE),
  PaperType    = data.frame(vals = c(0x00, 0x80),
                            levs = c("FANFOLD", "SINGLE"),
                            stringsAsFactors = FALSE),
  PrintPitch   = data.frame(vals = c(0x000, 0x400, 0x800),
                            levs = c("PICA", "ELITE", "FINE"),
                            stringsAsFactors = FALSE),
  PrintQuality = data.frame(vals = c(0x000, 0x100),
                            levs = c("DRAFT", "LETTER")),
  PrintSpacing = data.frame(vals = c(0x000, 0x200),
                            levs = c("SIX_LPI", "EIGHT_LPI"),
                            stringsAsFactors = FALSE),
  PrintImage   = data.frame(vals = c(0x00, 0x01),
                            levs = c("IMAGE_POSITIVE", "IMAGE_NEGATIVE"),
                            stringsAsFactors = FALSE),
  PrintAspect  = data.frame(vals = c(0x00, 0x01),
                            levs = c("ASPECT_HORIZ", "ASPECT_VERT"),
                            stringsAsFactors = FALSE),
  PrintShade   = data.frame(vals = c(0x00, 0x01, 0x02),
                            levs = c("SHADE_BW", "SHADE_GREYSCALE", "SHADE_COLOR"),
                            stringsAsFactors = FALSE),
  PaperSize    = data.frame(vals = (0:13)*16,
                            levs = c("US_LETTER", "US_LEGAL", "N_TRACTOR", "W_TRACTOR", "CUSTOM", paste0("EURO_A", 0:8)),
                            stringsAsFactors = FALSE),
  PrinterType  = data.frame(vals = 0:12,
                            levs = c("CUSTOM_NAME", "ALPHA_P_101", "BROTHER_15XL", "CBM_MPS1000",
                                     "DIAB_630", "DIAB_ADV_D25", "DIAB_C_150", "EPSON", "EPSON_JX_80",
                                     "OKIMATE_20", "QUME_LP_20", "HP_LASERJET", "HP_LASERJET_PLUS"),
                            stringsAsFactors = FALSE),
  LaceWB       = data.frame(vals = 0:1,
                            levs = c("NO_LACE", "LACE"),
                            stringsAsFactors = FALSE)
)

.SysConfigData <- data.frame(
  byte      = c(1, 1, 2, -8, -8, -8, -72, 1, 1, -6, 2, -8, 1, 1, 2, 2, -2, 2, -30, rep(2, 12),
                rep(-1, 3), 1, -12, -16, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1),
  signed    = c(TRUE, FALSE, FALSE,  FALSE,  FALSE,  FALSE,   FALSE, TRUE, TRUE,  FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,  FALSE, TRUE,   FALSE, rep(FALSE, 8), TRUE, FALSE, FALSE, FALSE,
                rep( FALSE, 3), FALSE,   FALSE,   FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  par.names = c("FontHeight", "PrinterPort", "BaudRate", "KeyRptSpeed",
                "KeyRptDelay", "DoubleClick", "PointerMatrix", "XOffset", "YOffset",
                "spriteColours", "PointerTicks", "WBColours", "ViewXOffset", "ViewYOffset", "ViewInitX",
                "ViewInitY", "EnableCLI", "PrinterType", "PrinterFilename", "PrintPitch",
                "PrintQuality", "PrintSpacing", "PrintLeftMargin", "PrintRightMargin",
                "PrintImage", "PrintAspect", "PrintShade", "PrintThreshold", "PaperSize",
                "PaperLength", "PaperType", "SerRWBits", "SerStopBuf", "SerParShk",
                "LaceWB", "Pad", "PrtDevName", "DefaultPrtUnit", "DefaultSerUnit",
                "RowSizeChange", "ColumnSizeChange", "PrintFlags", "PrintMaxWidth",
                "PrintMaxHeight", "PrintDensity", "PrintXOffset", "wb_Width", "wb_Height",
                "wb_Depth", "ext_size"),
  stringsAsFactors = FALSE
)

#' Read an Amiga system-configuration file
#'
#' Read a binary Amiga system-configuration file and return as [SysConfig] object.
#'
#' Amiga OS 1.x stored system preferences in a binary system-configuration file. This
#' function returns the file in a comprehensive format (a [SysConfig] object).
#'
#' @rdname read.SysConfig
#' @name read.SysConfig
#' @param file The file name of a system-configuration file to be read.
#' Can also be a connection that allows reading binary data.
#' @return Returns an S3 [SysConfig] class object based on the file that is read.
#' @examples
#' ## Put a simple SysConfig object into the tempdir:
#' write.SysConfig(simpleSysConfig(), file.path(tempdir(), "system-configuration"))
#' 
#' ## Now read the same file:
#' sc <- read.SysConfig(file.path(tempdir(), "system-configuration"))
#' 
#' ## and plot it
#' plot(sc)
#' @family SysConfig.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.SysConfig <- function(file) {
  dat <- .read.generic(file)
  rawToSysConfig(dat)
}

#' Write an Amiga system-configuration file
#'
#' Write a [SysConfig] class object to an Amiga binary system-configuration file.
#'
#' Amiga OS 1.x stored system preferences in a binary system-configuration file. This
#' function writes a [SysConfig] class object as such a binary file. This file
#' can be used on an Amiga or in an emulator.
#'
#' @rdname write.SysConfig
#' @name write.SysConfig
#' @param x An S3 [SysConfig] class object.
#' @param file A file name to which the binary file should be written.
#' @return Returns `NULL` or an `integer` status passed on by the
#' [close()] function, that is used to close the file connection.
#' It is returned invisibly.
#' 
#' @examples
#' ## First generate a simple SysConfig object to write to a file:
#' sc <- simpleSysConfig()
#' 
#' ## And write to the tempdir:
#' write.SysConfig(sc, file.path(tempdir(), "system-configuration"))
#' @family SysConfig.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.SysConfig <- function(x, file) {
  if (!inherits(x, "SysConfig")) stop("x should be of class SysConfig.")
  .write.generic(x, file)
}

#' Coerce raw data into a SysConfig class object
#'
#' [SysConfig] objects are comprehensive representations of binary Amiga
#' system-configuration files. Use this function to convert `raw` data from
#' such a file to a [SysConfig] object.
#'
#' The Amiga used the system-configuration file to store certain system preferences
#' in a binary file. With this function such `raw` data can be converted into
#' a more comprehensive [SysConfig] object. Use [AmigaFFH::as.raw()]
#' to achieve the inverse.
#'
#' @rdname rawToSysConfig
#' @name rawToSysConfig
#' @param x A vector of `raw` data that needs to be converted into an S3
#' [SysConfig] class object. It should have a length of at least 232. Although
#' system-configurations can be extended, such extended files are not supported here.
#' @return Returns a [SysConfig] class object based on `x`.
#' @examples
#' if (requireNamespace("adfExplorer", quietly = TRUE)) {
#'   library(adfExplorer)
#'   
#'   ## get the system-configuration from the adfExplorer example disk:
#'   disk <- connect_adf(
#'     system.file("example.adz", package = "adfExplorer")
#'   )
#'   virtual_file_con <- adf_file_con(disk, "devs/system-configuration")
#'   sc <- readBin(virtual_file_con, "raw", 1024)
#'   close(disk)
#' 
#'   ## This will get you the raw data from the file:
#'   typeof(sc)
#' 
#'   ## Convert the raw data to a more comprehensive named list (and S3 SysConfig class):
#'   sc <- rawToSysConfig(sc)
#' }
#' @family SysConfig.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToSysConfig <- function(x) {
  ## assuming pref. version 38.2 (slightly different from 37.4) see urls
  system.configuration <- .read.amigaData(x, .SysConfigData$byte, .SysConfigData$signed, .SysConfigData$par.names)
  ## EnableCLI is obsolete
  system.configuration$KeyRptSpeed           <- timeval(system.configuration$KeyRptSpeed)
  system.configuration$KeyRptDelay           <- timeval(system.configuration$KeyRptDelay)
  system.configuration$DoubleClick           <- timeval(system.configuration$DoubleClick)
  ## When converting back to raw, we need to adjust the third byte (VSTOP). It is set to 16 in the call
  system.configuration$PointerMatrix         <- rawToHWSprite(system.configuration$PointerMatrix)
  system.configuration$WBColours             <- amigaRawToColour(system.configuration$WBColours, colour.depth = "12 bit", n.bytes = "2")
  system.configuration$spriteColours         <- amigaRawToColour(system.configuration$spriteColours, colour.depth = "12 bit", n.bytes = "2")
  system.configuration$PointerMatrix@colours <- system.configuration$spriteColours
  system.configuration$PrinterFilename       <- .rawToCharNull(system.configuration$PrinterFilename)
  system.configuration$PrtDevName            <- .rawToCharNull(system.configuration$PrtDevName)
  system.configuration[names(.SysConfigFactors)] <-
    lapply(names(.SysConfigFactors), function(y) {
      .match.factor(system.configuration, y, .SysConfigFactors[[y]]$vals, .SysConfigFactors[[y]]$levs)
    })
  system.configuration$PrintFlags <- .match.multi.factor(system.configuration,
                                                         "PrintFlags",
                                                         .SysConfigMultiFactors[["PrintFlags"]]$vals,
                                                         .SysConfigMultiFactors[["PrintFlags"]]$levs)
  if (sum(grepl("DIMENSIONS", system.configuration$PrintFlags)) > 1)
    system.configuration$PrintFlags <- system.configuration$PrintFlags[!grepl("IGNORE_DIMENSIONS",
                                                                              system.configuration$PrintFlags)]
  if (sum(grepl("DITHERING", system.configuration$PrintFlags)) > 1)
    system.configuration$PrintFlags <- system.configuration$PrintFlags[!grepl("ORDERED_DITHERING",
                                                                              system.configuration$PrintFlags)]
  system.configuration$SerRWBits <- as.logical(.rawToBitmap(system.configuration$SerRWBits, FALSE, FALSE))
  names(system.configuration$SerRWBits) <- c(t(outer(c("write.bit", "read.bit"), 0:3, paste0)))
  system.configuration$SerParShk <- list(
    SerialParity  = .match.factor(list(SerialParity = ProTrackR::hiNybble(system.configuration$SerParShk)),
                                  "SerialParity", 0:4, c("SPARITY_NONE", "SPARITY_EVEN", "SPARITY_ODD",
                                                         "SPARITY_MARK", "SPARITY_SPACE")),
    HandshakeMode = .match.factor(list(HandshakeMode = ProTrackR::loNybble(system.configuration$SerParShk)),
                                  "HandshakeMode", 0:2, c("SHSHAKE_XON", "SHSHAKE_RTS", "SHSHAKE_NONE"))
  )
  system.configuration$SerStopBuf <- list(
    N.StopBits = ProTrackR::hiNybble(system.configuration$SerStopBuf) + 1L,
    BufSize    = .match.factor(list(BufSize = ProTrackR::loNybble(system.configuration$SerStopBuf)),
                               "BufSize", 0:5, c("SBUF_512", "SBUF_1024", "SBUF_2048", "SBUF_4096", "SBUF_8000", "SBUF_16000"))
  )
  class(system.configuration) <- "SysConfig"
  return(system.configuration)
}

#' Function to generate a simple Amiga system-configuration representation
#'
#' [SysConfig] objects are comprehensive representations of binary Amiga
#' system-configuration files. Use this function to create a simple [SysConfig] object.
#'
#' The Amiga used the system-configuration file to store certain system preferences
#' in a binary file. In the AmigaFFH package such files can be represented by the more
#' comprehensive [SysConfig] class object. Use this function to create such an object
#' with basic settings (which can be modified).
#'
#' @rdname simpleSysConfig
#' @name simpleSysConfig
#' @param options A named `list` with elements of the target
#' [SysConfig()] object that need to be modified.
#' @return Returns a comprehensive representation of a system-configuration file in the
#' for of a [SysConfig] class object.
#' @examples
#' ## Create a simple system-configuration (S3 SysConfigClass)
#' sc <- simpleSysConfig()
#' 
#' ## And modify it as you wish.
#' ## in this case change the setting for the printer
#' ## from the parallel port to the serial port:
#' sc$PrinterPort <- factor("SERIAL_PRINTER", levels(sc$PrinterPort))
#' 
#' ## It is also to provide modifications to the configuration
#' ## via the 'options' argument:
#' sc <- simpleSysConfig(options = list(FontHeight = 9))
#' @family SysConfig.operations
#' @author Pepijn de Vries
#' @export
simpleSysConfig <- function(options) {
  result <- paste0("789ce3606060650083c301608a53fd0090646460f0538088373430341c3c",
                   "e078e0b142f2836b02e51fce70d4ff38c152ffe70053fdbf038c60fac701",
                   "96fa0f0738ea1f1c1000d320fe9f06260686ffff4126fcffc8e7c2c0c0f7",
                   "8681916b150303ff7fb66ea06023830e482e3d352fb5283399012f6065f0",
                   "86b2d81914189cf02b46058c501a004f782dc6")
  pos <- 1:(nchar(result)/2)
  pos <- pos*2 - 1
  result <- as.raw(as.numeric(
    paste0("0x",c(sapply(result, substring, first = pos, last = pos + 1)))))
  result <- memDecompress(result, "gzip")
  result <- rawToSysConfig(result)
  if (!missing(options)) {
    for (opt in names(options)) {
      result[[opt]] <- options[[opt]]
    }
  }
  return(result)
}

#' @export
`$<-.SysConfig` <- function(x, i, value) {
  x[[i]] <- value
  x
}

#' @export
`[[<-.SysConfig` <- function(x, i, value) {
  if (!inherits(i, "character")) stop("Refer to elements by name, not by index number, when replacing them.")
  # make sure that x is in the correct order:
  cl <- class(x)
  x <- x[.SysConfigData$par.names]
  if (!(i %in% .SysConfigData$par.names)) stop(sprintf("Element \"%s\" is not part of SysConfig and cannot be assigned.", i))
  class(x) <- NULL
  if (i %in% c("WBColours", "spriteColours") && !all(.is.colour(value))) stop(sprintf("Can only assign colours to %s.", i))
  if (i == "WBColours" && length(value) != 4) stop("WBColours needs a vector of 4 colours.")
  if (i == "spriteColours" && length(value) != 3) stop("spriteColours needs a vector of 3 colours.")
  if (i == "Pad") {
    value <- as.raw(value)
    if (length(value) != 12) stop("'Pad' should be a vector of 12 raw values.")
  }
  if (i == "SerRWBits") {
    value <- as.logical(value)
    if (length(value) != 12) stop("'SerRWBits' should be a vector of 8 logical values.")
    names(sc$SerRWBits) <- c(paste0("write.bit", 0:3), paste0("read.bit", 0:3))
  }
  if (i == "SerStopBuf") {
    if (typeof(value) == "list" && all(names(value) == c("N.StopBits", "BufSize"))) {
      value$N.StopBits <- as.numeric(value$N.StopBits)
      if (value$N.StopBits < 0 || value$N.StopBits > 15) stop("value is out of range.")
      bfs <- c("SBUF_512", "SBUF_1024", "SBUF_2048", "SBUF_4096", "SBUF_8000", "SBUF_16000")
      if (is.numeric(value$BufSize)) value$BufSize <- bfs[match(value$BufSize, 0:5)]
      if (is.factor(value$BufSize))  value$BufSize <- as.character(value$BufSize)
      value$BufSize <- factor(value$BufSize[1], bfs)
      if (is.na(value$BufSize)) stop("Illegal value for SerStopBuf.")
    } else {
      stop("SerStopBuf should be a list with elements N.StopBits and BufSize")
    }
  }
  if (i == "SerParShk") {
    if (typeof(value) == "list" && all(names(value) == c("SerialParity", "HandshakeMode"))) {
      sp <- c("SPARITY_NONE", "SPARITY_EVEN", "SPARITY_ODD", "SPARITY_MARK", "SPARITY_SPACE")
      if (is.numeric(value$SerialParity)) value$SerialParity <- sp[match(value$SerialParity, 0:4)]
      if (is.factor(value$SerialParity))  value$SerialParity <- as.character(value$SerialParity)
      value$SerialParity <- factor(value$SerialParity[1], sp)
      hs <- c("SHSHAKE_XON", "SHSHAKE_RTS", "SHSHAKE_NONE")
      if (is.numeric(value$HandshakeMode)) value$HandshakeMode <- hs[match(value$HandshakeMode, 0:2)]
      if (is.factor(value$HandshakeMode))  value$HandshakeMode <- as.character(value$HandshakeMode)
      value$HandshakeMode <- factor(value$HandshakeMode[1], hs)
      if (is.na(value$HandshakeMode) || is.na(value$SerialParity)) stop("Illegal value for SerParShk")
    } else {
      stop("SerStopBuf should be a list with elements N.StopBits and BufSize")
    }
  }
  if (i == "EnableCLI") {
    value <- as.raw(value)
    if (length(value) != 2) stop("'EnableCLI' should be a vector of 2 raw values.")
  }
  if (i %in% c("PrtDevName", "PrinterFilename")) value <- as.character(value)
  if (i %in% c("KeyRptSpeed", "KeyRptDelay", "DoubleClick")) {
    if (is.numeric(value)) {
      if (value < 0) stop(sprintf("Negative numbers are not allowed for %s", i))
      class(value) <- "AmigaTimeVal"
    }
    if (!inherits(value, "AmigaTimeVal")) stop("Value cannot be cast to 'timeval' class object.")
  }
  if (i == "PointerMatrix") {
    if (inherits(value, "raster")) value <- rasterToHWSprite(value)
    if (!inherits(value, "hardwareSprite"))
      stop ("PointerMatrix element and its replacement should be an S4 class hardwareSprite object.")
    if (!all(dim(value) == 16)) stop("The pointer sprite should be 16 pixels wide and 16 pixels high.")
    x[["spriteColours"]] <- value@colours
  }
  fct <- .SysConfigFactors[[i]]
  if (!is.null(fct)) {
    if (is.factor(value)) {
      if (!all(levels(value) == fct$levs)) stop(sprintf("Illegal levels for factor %s.", i))
      value <- value[[1]]
    }
    if (is.numeric(value)) {
      if (!(value[[1]] %in% fct$vals)) stop (sprintf("Illegal value for %s.", i))
      value <- factor(fct$levs[fct$vals == value[[1]]], fct$levs)
    }
    if (is.character(value)) {
      if (!(value[[1]] %in% fct$levs)) stop(sprintf("Illegal level for factor %s.", i))
      value <- factor(value[[1]], fct$levs)
    }
  } else {
    fct <- .SysConfigMultiFactors[[i]]
    if (!is.null(fct)) {
      if (is.factor(value)) value <- as.character(value)
      if (is.character(value)) value <- match(value, .SysConfigMultiFactors[[i]]$levs)
      if (is.numeric(value)) {
        if (any(is.na(value))) stop(sprintf("Illegal value for %s.", i))
        value <- .bitwOrAll(value)
        if (value < 0) stop(sprintf("Illegal value for %s.", i))
        temp <- eval(parse(text = sprintf("list(%s = %i)", i, value)))
        value <- .match.multi.factor(temp, i,
                                    .SysConfigMultiFactors[[i]]$vals,
                                    .SysConfigMultiFactors[[i]]$levs)
      }
    } else {
      bt <- .SysConfigData$byte[.SysConfigData$par.names == i]
      sn <- .SysConfigData$signed[.SysConfigData$par.names == i]
      if (bt > 0) {
        value <- value[[1]]
        rn <- c(0, 2^(bt*8) - 1)
        if (sn) rn <- rn - ceiling(rn[2]/2)
        if (!sn && value < 0) stop("Negative values are not allowed for %s", i)
        if (value < rn[1] || value > rn[2]) stop("Value is out of range.")
      }
    }
  }
  x[[i]] <- value
  class(x) <- cl
  return(x)
}
