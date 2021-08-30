#' The S3 AmigaBitmapFont and AmigaBitmapFontSet classes
#' 
#' A comprehensive representation of monochromous Amiga bitmap fonts.
#' 
#' Nowadays fonts are represented by vector graphics an computer systems.
#' On the original Commodore Amiga, the screen resolution, system memory
#' and cpu speed were limited. On those systems, it was more efficient
#' to use bitmap images to represent the glyphs in fonts. The
#' \code{AmigaBitmapFontSet} and \code{AmigaBitmapFont} classes can be used
#' to represent Amiga bitmap fonts.
#' 
#' The Commodore Amiga had a directory named 'FONTS' located in the
#' root, where (bitmap) fonts were stored. Font sets were stored
#' under the font name with a *.font extension. Files with the *.font
#' extension did not contain the bitmap images of the font. Rather
#' the *.font file contained information on which font heights (in
#' pixels) are available, in addition to some other meta-information.
#' 
#' The bitmap images were stored in separate files for each individual
#' height. The \code{AmigaBitmapFontSet} is an S3 class that forms
#' a comprehensive format (named \code{list}) to represent the *.font
#' files. The \code{AmigaBitmapFont} is an S3 class is a comprehensive
#' format (named \code{list}) that represent each font bitmap and glyph
#' information. The \code{AmigaBitmapFontSet} objects will hold one or more
#' \code{AmigaBitmapFont} objects.
#' 
#' The \code{AmigaBitmapFont} and \code{AmigaBitmapFontSet} objects are
#' essentially named \code{list}s. Their structure and most important
#' elements are described below. Although it is possible to replace
#' elements manually, it is only advisable when you know what you
#' are doing as it may break the validity of the font.
#' 
#' @section AmigaBitmapFontSet:
#' \describe{
#'   \item{\code{fch_FileID}}{
#'     A \code{factor} with levels 'FontContents', 'TFontContents' and
#'     'ScalableOutline'. It specifies the type of font.
#'     Currently only the first level is supported.
#'   }
#'   \item{\code{fch_NumEntries}}{
#'     number of font heights available for this font. It should
#'     match with the length of \code{FontContents}. Do not change
#'     this value manually.
#'   }
#'   \item{\code{FontContents}}{
#'     This is a \code{list} with bitmap entries for each specific font
#'     height (in pixels). The name of each element in this list is
#'     'pt' followed by the height. Each element in this list holds
#'     the elements:
#'     \describe{
#'       \item{Miscellaneous}{Miscellaneous information from the *.font file}
#'       \item{\code{fc_FileName}}{
#'         This element represents the filename of the
#'         nested font bitmap images. Note that it should be a valid
#'         Commodore Amiga filename. It is best to modify this name
#'         using \code{\link{fontName}}. Note that this field could cause
#'         problems as Commodore Amiga filenames can contain characters
#'         that most modern platforms would not allow (such as the
#'         question mark).
#'       }
#'       \item{\code{BitmapFont}}{
#'         This element is of type \code{AmigaBitmapFont} and is structured
#'         as described in the following section. The information in this
#'         element is no longer part of the *.font file.
#'       }
#'     }
#'   }
#' }
#' 
#' @section AmigaBitmapFont:
#' Information represented by a \code{AmigaBitmapFont} is not stored
#' in *.font files. Rather it is stored in sub-directories of the font
#' in separate files. It has the following structure:
#' \describe{
#'   \item{Miscellaneous}{Elements with information on the font
#'   properties and style, and also relative file pointers.}
#'   \item{\code{glyph.info}}{
#'     A \code{data.frame} containing glyph info with information
#'     for specific glyphs on each row. Each row matches with a specific
#'     ASCII code, ranging from \code{tf_LoChar} up to \code{tf_HiChar}. There is an additional
#'     row that contains information for the default glyph that is
#'     out of the range of the \code{tf_LoChar} and \code{tf_HiChar}. The \code{data.frame}
#'     thus has \code{2 + tf_HiChar - tf_LoChar} rows. This
#'     table is used to extract and plot a glyph from the
#'     \code{bitmap} image correctly.
#'   }
#'   \item{\code{bitmap}}{
#'     Is a monochromous bitmap image of all the font's glyphs in a
#'     single line. It is a simple \code{raster} object
#'     (see \code{\link[grDevices]{as.raster}}) with an additional
#'     attribute 'palette', which lists the two colours in the image. In
#'     this palette, the first colour is the background colour and the
#'     second colour is interpreted as the foregroundcolour.
#'   }
#' }
#' 
#' @section Useful functions:
#' For importing and exporting the following functions are useful:
#' \code{\link{read.AmigaBitmapFont}}, \code{\link{read.AmigaBitmapFontSet}},
#' \code{\link{write.AmigaBitmapFont}} and \code{\link{write.AmigaBitmapFontSet}}.
#' 
#' The following generic functions are implemented for these objects:
#' \code{\link[AmigaFFH]{plot}}, \code{print},
#' \code{\link[AmigaFFH]{as.raster}} and \code{\link[AmigaFFH]{as.raw}}.
#' 
#' Use \code{\link[AmigaFFH]{c}} to combine one or more
#' \code{AmigaBitmapFont} objects into a \code{AmigaBitmapFontSet}.
#' 
#' @docType class
#' @aliases AmigaBitmapFontSet
#' @name AmigaBitmapFont
#' @rdname AmigaBitmapFont
#' @family AmigaBitmapFont.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @references
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Libraries_Manual_guide/node03E0.html}
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Libraries_Manual_guide/node03DE.html}
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Libraries_Manual_guide/node05BA.html}
#' @examples
#' \dontrun{
#' ## 'font_example' is an example of the AmigaBitmapFontSet object:
#' data(font_example)
#' 
#' ## An AmigaBitmapFont object can also be extracted from this object:
#' font_example_9 <- getAmigaBitmapFont(font_example, 9)
#' 
#' ## the objects can be printed, plotted, converted to raw data or a raster:
#' print(font_example)
#' plot(font_example)
#' font_example_raw    <- as.raw(font_example)
#' font_example_raster <- as.raster(font_example)
#' 
#' ## You can also format text using the font:
#' formated_raster     <- as.raster(font_example, text = "Foo bar", style = "bold")
#' plot(font_example, text = "Foo bar", style = "underlined", interpolate = F)
#' }
NULL

.print_to_raster <- function(text, font, style = NULL, palette = NULL) {
  if (inherits(font, "AmigaBitmapFontSet")) {
    h <- availableFontSizes(font)
    font <- getAmigaBitmapFont(font, h[length(h)])
  }
  if (!inherits(font, "AmigaBitmapFont")) stop("'font' should be of class AmigaBitmapFont or AmigaBitmapFontSet.")
  if (!is.null(style)) {
    style <- match.arg(style, c("bold", "italic", "underlined"), T)
    ## You can't apply styles to a font that is already styled:
    style <- style[!font$tf_Style[toupper(style)]]
  }
  if (length(text) > 1) {
    warning("'text' has multiple elements, only using the first.")
    text <- text[1]
  }
  # split text along lines:
  text <- strsplit(text, "\n")[[1]]
  pal <- attr(font$bitmap, "palette")
  bm  <- as.matrix(font$bitmap)
  bm  <- apply(bm, 2, function(y) y == pal[[2]])
  if (!is.null(palette)) pal <- palette
  result <- lapply(text, function(y) {
    ## get ascii codes:
    y <- utf8ToInt(enc2utf8(y))
    ## the final glyph is the default glyph when it is out of range
    y[y < font$tf_LoChar] <- font$tf_HiChar + 1
    y[y > font$tf_HiChar] <- font$tf_HiChar + 1
    y <- 1 + y - font$tf_LoChar
    positions <- apply(
      font$glyph.info[y,][,names(font$glyph.info) %in% c("glyphWidth", "charSpace"), drop = F],
      1, max
    )
    if (!is.null(font$glyph.info$charKern)) positions <- positions + c(font$glyph.info$charKern[y][-1], 0)
    positions <- cumsum(c(0, positions))[1:length(positions)]
    positions <- 1 + positions - min(positions)
    h <- font$tf_YSize
    glyphs <- mapply(function(loc, w, sp, kern) {
      result <- matrix(F, h + 2,
                       max(c(sp, w)) +
                         ifelse("bold" %in% style, font$tf_BoldSmear, 0) +
                         ifelse(kern > 0, kern, 0) +
                         ifelse("italic" %in% style, ceiling((h + 1)/2), 0))
      if (w > 0) {
        for (i in 0:ifelse("bold" %in% style, font$tf_BoldSmear, 0)) {
          result[1:h, i + (1:w)] <- result[1:h, i + (1:w)] | bm[,(loc + 1):(loc + w)]
        }
        if ("italic" %in% style) {
          for (j in 1:(h - ifelse(h > 15, 2, 1))) { ## The displacement seems to shift for h > 15
            displacement <- floor((h - j + ifelse(h > 15, 0, 1))/2)
            result[j,] <- c(rep(F, displacement), utils::head(result[j,], -displacement))
          }
        }
      }
      if (ncol(result) > 0) {
        if ("underlined" %in% style) {
          result[font$tf_Baseline + 2,] <- result[font$tf_Baseline + 2,] |
            (c(!result[font$tf_Baseline + 2,][-1], T) & c(T, !result[font$tf_Baseline + 2,][-ncol(result)]))
        }
      }
      result
    },
    loc = font$glyph.info$glyphLocation[y],
    w   = font$glyph.info$glyphWidth[y],
    sp  = if(is.null(font$glyph.info$charSpace)) rep(0, length(y)) else
      font$glyph.info$charSpace[y],
    kern  = if(is.null(font$glyph.info$charKern)) rep(0, length(y)) else
      c(font$glyph.info$charKern[y][-1], 0),
    SIMPLIFY = F)
    widths <- unlist(lapply(glyphs, ncol))
    result <- matrix(F, nrow = h + 2, ncol = max(positions + widths))
    lapply(1:length(glyphs), function(i) {
      result[,positions[i]:(positions[i] + widths[i] - 1)] <<-
        result[,positions[i]:(positions[i] + widths[i] - 1)] |
        glyphs[[i]]
    })
    ## make background transparent if no palette is specified
    if (is.null(palette)) pal[1] <<- grDevices::adjustcolor(pal[1], alpha.f = 0)
    result <- apply(result, 2, function(y) pal[1 + as.numeric(y)])
    return(result)
  })
  result2 <- matrix(pal[1],
                    sum(unlist(lapply(result, nrow))),
                    sum(unlist(lapply(result, ncol))))
  lapply(1:length(result), function(i) {
    result2[(i - 1)*(font$tf_YSize + 2) + 1:(font$tf_YSize + 2),][,1:ncol(result[[i]])] <<- result[[i]]
    NULL
  })
  result2 <- as.raster(result2)
  return(result2)
}

.amiga.node <- data.frame(
  byte      = c(4, 4, 1, 1, 4),
  signed    = c(F, F, F, T, F),
  par.names = c("ln_Succ", "ln_Pred", "ln_Type", "ln_Pri", "ln_Name"),
  stringsAsFactors = F
)

.amiga.node.types <- c("NT_UNKNOWN", "NT_TASK", "NT_INTERRUPT", "NT_DEVICE", "NT_MSGPORT", "NT_MESSAGE",
                       "NT_FREEMSG", "NT_REPLYMSG", "NT_RESOURCE", "NT_LIBRARY", "NT_MEMORY", "NT_SOFTINT",
                       "NT_FONT", "NT_PROCESS", "NT_SEMAPHORE", "NT_SIGNALSEM", "NT_BOOTNODE", "NT_KICKMEM",
                       "NT_GRAPHICS", "NT_DEATHMESSAGE", "NT_USER", "NT_EXTENDED")

.amiga.font.sets <- data.frame(
  byte      = c(2, 2),
  signed    = c(F, F),
  par.names = c("fch_FileID", "fch_NumEntries"),
  stringsAsFactors = F
)

.amiga.font.types <- c("FontContents", "TFontContents", "ScalableOutline")

.read.amiga.node <- function(dat) {
  ## Currently only works for simple nodes as included in font files
  result <- with(.amiga.node, .read.amigaData(dat, byte, signed, par.names))
  if (result$ln_Type < 22 || result$ln_Type >= 254) {
    result$ln_Type <- result$ln_Type + 1
    if (result$ln_Type > 254) result$ln_Type <- result$ln_Type - 234
    result$ln_Type <- .match.factor(result, "ln_Type", 1:length(.amiga.node.types),
                                    .amiga.node.types)
    }
  return(result)
}

.amiga.node.to.raw <- function(node) {
  node <- node[.amiga.node$par.names]
  node$ln_Type <- match(as.character(node$ln_Type), .amiga.node.types)
  if (node$ln_Type >= 22) node$ln_Type <- node$ln_Type + 234
  node$ln_Type <- node$ln_Type - 1
  node <- with(.amiga.node, .write.amigaData(node, byte, signed, par.names))
  node
}

.amiga.font.header <- data.frame(
  byte      = c(-36, -14, 2, 2, 4, -32, -14, 4, 2, 2, -1, -1, 2, 2, 2, 2, 1, 1, 4, 2, 4, 4, 4),
  signed    = c(  F,   F, F, F, F,   F,   F, F, F, F,  F,  F, F, F, F, F, F, F, F, F, F, F, F),
  par.names = c("leadingHunks", "node.disklink", "dfh_FileID", "dfh_Revision", "dfh_Segment",
                 "fontName", "node.message", "mn_ReplyPort", "mn_Length", "tf_YSize",
                 "tf_Style", "tf_Flags", "tf_XSize", "tf_Baseline", "tf_BoldSmear",
                 "tf_Accessors", "tf_LoChar", "tf_HiChar", "tf_CharData", "tf_Modulo",
                 "tf_CharLoc", "tf_CharSpace", "tf_CharKern"),
  stringsAsFactors = F
)

#' Read AmigaBitmapFontSet from *.font file
#'
#' Reads \code{\link{AmigaBitmapFontSet}} from *.font file including
#' all nested bitmap images for all font heights.
#'
#' The *.font file only holds meta-information. The bitmap images for
#' each font height are stored in separate files, which are listed
#' in the *.font file. The function reads the *.font file, including
#' all nested bitmap files and returns it as a
#' \code{\link{AmigaBitmapFontSet}}.
#'
#' It can also read *.font files
#' from \code{\link[adfExplorer]{amigaDisk-class}} objects,
#' but that requires the adfExplorer package to be installed.
#' @rdname read.AmigaBitmapFontSet
#' @name read.AmigaBitmapFontSet
#' @param file A \code{character} string of the filename of the *.font file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @param ... Currently ignored.
#' @return Returns an \code{\link{AmigaBitmapFontSet}} object read from the specified file.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## in order to read, we first need to write a file"
#' write.AmigaBitmapFontSet(font_example, tempdir())
#' 
#' ## The font is written as 'AmigaFFH.font' as that name
#' ## is embedded in the AmigaBitmapFontSet object 'font_example'.
#' ## We can read it as follows:
#' font.read <- read.AmigaBitmapFontSet(file.path(tempdir(), "AmigaFFH.font"))
#' 
#' ## similarly, the file can also be written and read from and to
#' ## a virtual amiga disk. The following codes requires the 'adfExplorer'
#' ## package:
#' adf <- adfExplorer::blank.amigaDOSDisk("font.disk")
#' adf <- adfExplorer::dir.create.adf(adf, "FONTS")
#' adf <- write.AmigaBitmapFontSet(font_example, "DF0:FONTS", adf)
#' font.read <- read.AmigaBitmapFontSet("DF0:FONTS/AmigaFFH.font", adf)
#' }
#' @family AmigaBitmapFont.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaBitmapFontSet <- function(file, disk = NULL, ...) {
  dat <- .read.generic(file, disk)
  rawToAmigaBitmapFontSet(dat, file, disk)
}

#' Coerce raw data into an AmigaBitmapFontSet class object
#'
#' \code{\link{AmigaBitmapFontSet}} objects are comprehensive representations of binary Amiga
#' font files (*.font). Use this function to convert \code{raw} data from
#' such a file to an \link{AmigaBitmapFontSet} object.
#'
#' This function converts \code{raw} data as stored in *.font
#' files. The function also needs the file location, in order
#' to load the nested bitmap images for each font height.
#' This function is effectively the inverse of \code{\link[AmigaFFH]{as.raw}}.
#'
#' @rdname rawToAmigaBitmapFontSet
#' @name rawToAmigaBitmapFontSet
#' @param x A \code{vector} of \code{raw} data that needs to be
#' converted into an \code{\link{AmigaBitmapFontSet}}.
#' @param file The \code{raw} version of the \code{\link{AmigaBitmapFontSet}}
#' does not contain the nested font bitmap images. In order to correctly
#' construct an \code{\link{AmigaBitmapFontSet}} the file location of the
#' original *.font file is required in order to read and include the
#' font bitmap image information. \code{file} should thus be a \code{character}
#' string specifying the file location of the *.font file.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Returns an \code{\link{AmigaBitmapFontSet}} object.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## First create raw font set data. Note that this raw data
#' ## does not include the nested font bitmap images.
#' fontset.raw <- as.raw(font_example)
#' 
#' ## Therefore it is necesary to have the entire font stored as files:
#' write.AmigaBitmapFontSet(font_example, tempdir())
#' 
#' font.restored <- rawToAmigaBitmapFontSet(fontset.raw, file.path(tempdir(), "AmigaFFH.font"))
#' }
#' @family AmigaBitmapFont.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBitmapFontSet <- function(x, file, disk = NULL) {
  result <- with(.amiga.font.sets, .read.amigaData(x[1:4], byte, signed, par.names))

  result$fch_FileID <- .match.factor(result, "fch_FileID", c(0x0f00, 0x0f02, 0x0f03),
                                  .amiga.font.types)
  
  if (result$fch_FileID != "FontContents") stop(sprintf("%s font type is not (yet) supported.", as.character(result$fch_FileID)))

  result[["FontContents"]] <- lapply(result$fch_NumEntries:1, function(i) {
    if (result$fch_FileID == "FontContents") {
      offset <- 260*(i - 1)
      result <- list(
        fc_FileName = .rawToCharNull(x[-1:-4][offset + 1:256]),
        fc_YSize    = .rawToAmigaInt(x[-1:-4][offset + 257:258], 16, F),
        fc_Style    = x[-1:-4][offset + 259],
        fc_Flags    = x[-1:-4][offset + 260]
      )
      result$fc_Style        <- as.logical(.rawToBitmap(result$fc_Style, F, F))
      names(result$fc_Style) <- c("UNDERLINED", "BOLD", "ITALIC", "EXTENDED", "RESERVED1", "RESERVED2", "COLORFONT", "TAGGED")
      result$fc_Flags        <- as.logical(.rawToBitmap(result$fc_Flags, F, F))
      names(result$fc_Flags) <- c("ROMFONT", "DISKFONT", "REVPATH", "TALLDOT", "WIDEDOT", "PROPORTIONAL", "DESIGNED", "REMOVED")
      
      if (is.null(disk)) {
        result[["BitmapFont"]] <-
          ## replace the file amiga file separator with the platform dependent file separator,
          ## and read font from file:
          read.AmigaBitmapFont(file.path(dirname(file), gsub("[/]", .Platform$file.sep, result$fc_FileName)))
      } else {
        result[["BitmapFont"]] <-
          read.AmigaBitmapFont(paste(c(utils::head(strsplit(file, "/", T)[[1]], -1), result$fc_FileName), collapse = "/"), disk)
      }
      if (any(result$fc_Style != result$BitmapFont$tf_Style)) warning(sprintf("Styles defined in main font (*.font) and bitmap file (%s) do not match.", result$fc_FileName))
      if (any(result$fc_Flags != result$BitmapFont$tf_Flags)) warning(sprintf("Flags defined in main font (*.font) and bitmap file (%s) do not match.", result$fc_FileName))
      result
    } else {
      stop("This bitmap font type is not (yet) supported")
    }
  })

  pt.size <- as.numeric(unlist(lapply(strsplit(unlist(lapply(result$FontContents, function(x) x$fc_FileName)), "/"),
                                      function(x) x[[2]])))
  result$FontContents <- result$FontContents[order(pt.size)]
  names(result$FontContents) <- paste0("pt", sort(pt.size))
  class(result) <- "AmigaBitmapFontSet"
  return(result)
}

#' Read an AmigaBitmapFont class object from a file
#'
#' Amiga Font Bitmaps of distinctive font heights are stored in separate
#' files, which in combination form a font collection or set. This
#' function can be used to read a specific bitmap from a set and returns
#' it as an \code{\link{AmigaBitmapFont}} class object.
#'
#' Individual font bitmaps are stored in a font's subdirectory where
#' the file name is usually equal to the font height in pixels. This
#' function will read such a font bitmap file and return it as an
#' \code{\link{AmigaBitmapFont}} class object. It can also read such
#' files from \code{\link[adfExplorer]{amigaDisk-class}} objects,
#' but that requires the adfExplorer package to be installed.
#'
#' @rdname read.AmigaBitmapFont
#' @name read.AmigaBitmapFont
#' @param file The file name of a font subset is usually simply a numeric number
#' indicating the font height in pixels. Use \code{file} as a \code{character}
#' string representing that file location.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @param ... Arguments passed on to \code{\link{rawToAmigaBitmapFont}}.
#' @return Returns an \code{\link{AmigaBitmapFont}} object read from the specified file.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## Let's store the example font first:
#' write.AmigaBitmapFontSet(font_example, tempdir())
#' 
#' ## Now read a specific subset from the font files:
#' font.sub <- read.AmigaBitmapFont(file.path(tempdir(), "AmigaFFH", "9"))
#' 
#' ## The same can be done with a virtual Amiga disk. The following
#' ## examples require the 'adfExplorer' package.
#' font.disk <- adfExplorer::blank.amigaDOSDisk("font.disk")
#' font.disk <- adfExplorer::dir.create.adf(font.disk, "FONTS")
#' font.disk <- write.AmigaBitmapFontSet(font_example, "DF0:FONTS", font.disk)
#' font.sub <- read.AmigaBitmapFont("DF0:FONTS/AmigaFFH/9", font.disk)
#' }
#' @family AmigaBitmapFont.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaBitmapFont <- function(file, disk = NULL, ...) {
  dat <- .read.generic(file, disk)
  rawToAmigaBitmapFont(dat, file, ...)
}

#' Write an AmigaBitmapFont(set) file
#'
#' Functions to write \code{\link{AmigaBitmapFont}} and \code{\link{AmigaBitmapFontSet}}
#' class objects to files.
#'
#' \code{\link{AmigaBitmapFontSet}} class objects are written to a *.font
#' file. The filename used for this purpose is obtained from the object
#' itself using \code{\link{fontName}}. In addition, a subdirectory is
#' created automatically (when it doesn't already exist)
#' to which al the separate bitmap images for each font height are written
#' to individual files.
#'
#' \code{\link{AmigaBitmapFont}} class objects can also be written to a
#' file. In order to use it on a Commodore Amiga or emulator, it is better
#' to embed the font bitmap in a font set (using \code{\link[AmigaFFH]{c}})
#' and write the set to corresponding files.
#' @rdname write.AmigaBitmapFont
#' @name write.AmigaBitmapFont
#' @param x Respectively an \code{\link{AmigaBitmapFont}} or a
#' \code{\link{AmigaBitmapFontSet}} object depending on which of the
#' write-functions is called. This is the object that will be written
#' to the specified file.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBitmapFont}} object) needs to be written.
#' It is common practice on the Amiga to use the font height in pixels as
#' file name.
#' @param path A \code{character} string specifying the path where
#' \code{x} (an \code{\link{AmigaBitmapFontSet}} object) needs to be stored.
#' The filename for the font set will be extracted from \code{x} using
#' \code{\link{fontName}} followed by the *.font extension. A subdirectory
#' will be created with the same name (without the extension) if it doesn't
#' already exists. In this subdirectory all the nested \code{\link{AmigaBitmapFont}}
#' objects are stored.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Invisibly returns the result of the call of \code{close} to the
#' file connection. Or, when \code{disk} is specified, a copy of
#' \code{disk} is returned to which the file(s) is/are written.
#' @examples
#' \dontrun{
#' ## obtain a bitmap font set:
#' data(font_example)
#' 
#' ## write the font set to their files. The file name
#' ## is extracted from the font object, so you only have
#' ## to provide the path:
#' write.AmigaBitmapFont(font_example, temp.dir())
#' 
#' ## extract a font bitmap:
#' font <- getAmigaBitmapFont(font_example, 9)
#' 
#' ## and write it to the temp dir:
#' write.AmigaBitmapFont(font, file.path(temp.dir(), "9"))
#' 
#' ## The following examples require the 'adfExplorer' package:
#' font.disk <- adfExplorer::blank.amigaDOSDisk("font.disk")
#' font.disk <- adfExplorer::dir.create.adf(font.disk, "FONTS")
#' font.disk <- write.AmigaBitmapFontSet(font_example, "DF0:FONTS", font.disk)
#' }
#' @family AmigaBitmapFont.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.AmigaBitmapFont <- function(x, file, disk = NULL) {
  if (!inherits(x, "AmigaBitmapFont")) stop("x should be of class AmigaBitmapFont.")
  .write.generic(x, file, disk)
}

#' @rdname write.AmigaBitmapFont
#' @name write.AmigaBitmapFontSet
#' @export
write.AmigaBitmapFontSet <- function(x, path = getwd(), disk = NULL) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  filenames <- unlist(lapply(x$FontContents, function(y) y$fc_FileName))
  filenames <- do.call(rbind, strsplit(filenames, "/"))
  if (ncol(filenames) != 2) stop("Unexpected file structure.")
  if (length(unique(filenames[,1])) != 1) stop("Not a single base name for the font.")
  fn <- sprintf("%s.font", filenames[1, 1])
  if (is.null(disk)) {
    if (path != "") fn <- file.path(path, fn)
  } else {
    if (path != "") fn <- paste(path, fn, sep = ifelse(substr(path, nchar(path), nchar(path)) == "/", "", "/"))
  }
  disk <- .write.generic(x, fn, disk)
  dr <- filenames[1, 1]
  if (is.null(disk)) {
    if (path != "") dr <- file.path(path, dr)
    if (!dir.exists(dr))
      dir.create(dr)
    result <- lapply(1:nrow(filenames), function(y) {
      fn <- file.path(filenames[y, 1], filenames[y, 2])
      if (path != "") fn <- file.path(path, fn)
      write.AmigaBitmapFont(x$FontContents[[y]]$BitmapFont, fn)
    })
    return(invisible(result[[length(result)]]))
  } else {
    if (path != "") dr <- paste(path, dr, sep = ifelse(substr(path, nchar(path), nchar(path)) == "/", "", "/"))
    if (!adfExplorer::adf.file.exists(disk, dr)) ## better to use 'dir.exists', this needs to be implemented in adfExplorer
      disk <- adfExplorer::dir.create.adf(disk, dr)
    lapply(1:nrow(filenames), function(y) {
      fn <- paste(filenames[y, 1], filenames[y, 2], sep = "/")
      if (path != "") fn <- paste(path, fn, sep = ifelse(substr(path, nchar(path), nchar(path)) == "/", "", "/"))
      disk <<- write.AmigaBitmapFont(x$FontContents[[y]]$BitmapFont, fn, disk = disk)
    })
    return(disk)
  }
}

#' Coerce raw data into an AmigaBitmapFont class object
#'
#' \code{\link{AmigaBitmapFont}} objects are comprehensive representations of binary Amiga
#' font subset files. The file name is usually simply a numeric number
#' indicating the font height in pixels. Use this function to convert
#' \code{raw} content from such a file to an \code{\link{AmigaBitmapFont}} object.
#'
#' This function converts \code{raw} data as stored in font bitmap
#' files. These files are stored in subdirectories with the font's
#' name and usually have the font height in pixels as file name.
#' This function is effectively the inverse of \code{\link[AmigaFFH]{as.raw}}.
#'
#' @rdname rawToAmigaBitmapFont
#' @name rawToAmigaBitmapFont
#' @param x An \code{\link{AmigaBitmapFont}} object which needs to be converted
#' into \code{raw} data.
#' @param ... Currently ignored.
#' @return A \code{vector} of \code{raw} data representing \code{x}.
#' @examples
#' \dontrun{
#' ## first create raw data that can be converted into a AmigaBitmapFont
#' data(font_example)
#' font.raw <- as.raw(getAmigaBitmapFont(font_example, 9))
#' 
#' ## Convert it back into an AmigaBitmapFont object:
#' font <- rawToAmigaBitmapFont(font.raw)
#' }
#' @family AmigaBitmapFont.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBitmapFont <- function(x, ...) {
  result <- with(.amiga.font.header, .read.amigaData(x, byte, signed, par.names))
  index.trailing.hunks     <- (1 + 4*(.rawToAmigaInt(result$leadingHunks[21:24], 32, F) + 8))
  ## remove first part as it is not required for interpreting the data
  
  ## Check the leading hunks
  result$leadingHunks[c(-21:-24, -29:-32)]
  idx.check <- c(3,4,12,27,28,33,35,36)

  if (any(result$leadingHunks[idx.check] != as.raw(c(0x03, 0xf3, 0x01, 0x03, 0xe9, 0x70, 0x4e, 0x75))) ||
      any(result$leadingHunks[c(-21:-24, -29:-32, -34, -idx.check)] != raw(1)) ||
      !(result$leadingHunks[34] %in% as.raw(c(0x00, 0x64, 0xff)))) {
    warning("Unexpected file header. This file may not be a font bitmap.")
  }
  ## remove leading hunks after checks. They are no longer needed...
  result$leadingHunks   <- NULL
  
  if (index.trailing.hunks > length(x)) stop("Unexpected end of file")
  ## remove trailing hunks from data
  trailing.hunks         <- x[index.trailing.hunks:length(x)]
  x                      <- x[1:(index.trailing.hunks - 1)]
  result$dfh_FileID      <- .match.factor(result, "dfh_FileID", 0xf80, "DFH_ID") ## Disk Font Header
  result$tf_Style        <- as.logical(.rawToBitmap(result$tf_Style, F, F))
  names(result$tf_Style) <- c("UNDERLINED", "BOLD", "ITALIC", "EXTENDED", "RESERVED1", "RESERVED2", "COLORFONT", "TAGGED")
  result$tf_Flags        <- as.logical(.rawToBitmap(result$tf_Flags, F, F))
  names(result$tf_Flags) <- c("ROMFONT", "DISKFONT", "REVPATH", "TALLDOT", "WIDEDOT", "PROPORTIONAL", "DESIGNED", "REMOVED")
  
  result$node.disklink     <- .read.amiga.node(result$node.disklink)
  result$node.message      <- .read.amiga.node(result$node.message)
  result$fontName          <- .rawToCharNull(result$fontName)
  
  ## mn_ReplyPort points to FontExtension. Current implementation ignores these extensions
  n.glyphs <- 2 + result$tf_HiChar - result$tf_LoChar # +1 for difference in index base; another +1 for the default character
  
  glyph.info <- .rawToAmigaInt(x[result$tf_CharLoc + 32 + (1:(4*n.glyphs))],
                               16, F)
  # glyph.info <- as.data.frame(matrix(glyph.info, ncol = 2, byrow = T))
  glyph.info <- as.data.frame(matrix(glyph.info, ncol = 2, byrow = T))
  names(glyph.info) <- c("glyphLocation", "glyphWidth")

  if (result$tf_CharSpace != 0) {
    glyph.info$charSpace <- .rawToAmigaInt(x[result$tf_CharSpace + 32 + (1:(2*n.glyphs))],
                                           16, T)
  }
  if (result$tf_CharKern != 0) {
    glyph.info$charKern <- .rawToAmigaInt(x[result$tf_CharKern + 32 + (1:(2*n.glyphs))],
                                          16, T)
  }
  result[["glyph.info"]] <- glyph.info

  ## trailing.hunks
  trailing.hunks <- .rawToAmigaInt(trailing.hunks, 32, F)
  hunk.dat <- c(2, 7, 19, 21, if(result$tf_CharSpace != 0) 22, if(result$tf_CharKern != 0) 23)
  hunk.dat <- cumsum(abs(.amiga.font.header$byte[-1]))[hunk.dat - 1]
  if (any(trailing.hunks[c(1, length(trailing.hunks))] != c(1004, 1010)) ||
      trailing.hunks[2] != (length(trailing.hunks) - 5) ||
      any(trailing.hunks[c(3, length(trailing.hunks) - 1)] != 0) ||
      !all(trailing.hunks[4:(length(trailing.hunks) - 2)] %in% hunk.dat))
    warning("Unexpected trailing file hunks.")
  
  font.bitmap.data <- x[result$tf_CharData + 32 + (1:(result$tf_Modulo*result$tf_YSize))]
  
  result[["bitmap"]] <- bitmapToRaster(font.bitmap.data,
                                       w = result$tf_Modulo*8,
                                       h = result$tf_YSize,
                                       depth = 1, palette = c("white", "black"))
  attr(result[["bitmap"]], "palette") <- c("white", "black")
  
  class(result) <- "AmigaBitmapFont"
  return(result)
}

#' @rdname plot
#' @name plot
#' @export
plot.AmigaBitmapFont <- function(x, y, ...) {
  if (!inherits(x, "AmigaBitmapFont")) stop("x should be of class AmigaBitmapFont.")
  args <- list(...)
  raster.args <- list(x = x)
  for (elm in c("text", "palette", "style")) {
    raster.args[[elm]] <- args[[elm]]
    args[[elm]] <- NULL
  }
  args$x <- do.call(as.raster, raster.args)
  if (is.null(args$asp)) {
    args$asp <- 1
    if (x$tf_Flags["TALLDOT"]) args$asp <- args$asp*2
    if (x$tf_Flags["WIDEDOT"]) args$asp <- args$asp/2
  }
  do.call(plot, args)
}

#' @rdname plot
#' @name plot
#' @export
plot.AmigaBitmapFontSet <- function(x, y, ...) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  if (missing(y)) {
    args <- list(...)
    raster.args <- list(x = x)
    for (elm in c("text", "palette", "style")) {
      raster.args[[elm]] <- args[[elm]]
      args[[elm]] <- NULL
    }
    args$x <- do.call(as.raster, raster.args)

    if (is.null(args$asp)) {
      dot <- do.call(rbind, lapply(x$FontContents, function(z) z$fc_Flags[c("TALLDOT", "WIDEDOT")]))
      dot <- as.data.frame(table(as.data.frame(dot)))
      dot <- dot[which(dot$Freq == max(dot$Freq))[[1]],]
      args$asp <- 1
      if (dot$TALLDOT == "TRUE") args$asp <- args$asp*2
      if (dot$WIDEDOT == "TRUE") args$asp <- args$asp/2
    }
    do.call(plot, args)
  } else {
    plot(getAmigaBitmapFont(x, y), ...)
  }
}

#' @export
print.AmigaBitmapFont <- function(x, ...) {
  if (!inherits(x, "AmigaBitmapFont")) stop("x should be of class AmigaBitmapFont.")
  cat(sprintf("  y-size %i, %s",
              x$tf_YSize,
              paste(tolower(c(names(x$tf_Flags)[x$tf_Flags],
                              names(x$tf_Style)[x$tf_Style])), collapse = ", ")))
  cat("\n")
  invisible(NULL)
}

#' @export
print.AmigaBitmapFontSet <- function(x, ...) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  cat(fontName(x))
  cat("\n")
  lapply(x$FontContents, function(y) print(y$BitmapFont))
  invisible(NULL)
}

#' @rdname as.raw
#' @name as.raw
#' @export
as.raw.AmigaBitmapFont <- function(x, ...) {
  withCallingHandlers({ ## TODO remove handlers when replace functions are implemented
    ## initial checks. Throw errors when checks are unsuccessful
    if (!inherits(x, "AmigaBitmapFont")) stop("x should be of class AmigaBitmapFont.")
    max.loc <- max(x$glyph.info$glyphLocation)
    if (max.loc + max(x$glyph.info$glyphWidth[x$glyph.info$glyphLocation == max.loc]) > dim(x$bitmap)[2]) {
      stop("Glyph information exceeds bitmap dimensions!")
    }
    rm(max.loc)
    if (x$tf_YSize != dim(x$bitmap)[1]) stop ("tf_YSize does not bitmap height.")
    if (x$dfh_FileID != "DFH_ID") stop("Unexpected file ID...")
    if ((x$tf_HiChar - x$tf_LoChar + 2) != nrow(x$glyph.info)) stop("Glyph information does not match the number of characters...")
    if (x$tf_Modulo*8 != dim(x$bitmap)[2]) stop("tf_Modulo*8 does not equal the bitmap width.")
    
    x$leadingHunks    <- raw(36)
    x$leadingHunks[c(3:4, 12, 27:28, 33:36)] <- as.raw(c(0x03, 0xF3, 0x01, 0x03, 0xE9, 0x70, 0x64, 0x4E, 0x75))
    
    header               <- x[.amiga.font.header$par.names]
    header$node.disklink <- .amiga.node.to.raw(header$node.disklink)
    header$node.message  <- .amiga.node.to.raw(header$node.message)
    header$dfh_FileID    <- 0xF80
    header$fontName      <- charToRaw(header$fontName)[1:32]
    header$tf_Style      <- .bitmapToRaw(header$tf_Style, F, F)
    header$tf_Flags      <- .bitmapToRaw(header$tf_Flags, F, F)
    
    header$tf_CharLoc <- 110  # This is where the first data always start
    addToPointer      <- 2*prod(dim(x$glyph.info[,c("glyphLocation", "glyphWidth"), drop = F]))
    if (is.null(x$glyph.info$charSpace)) {
      header$tf_CharSpace  <- 0
    } else {
      header$tf_CharSpace <- header$tf_CharLoc + addToPointer
      addToPointer        <- 2*prod(dim(x$glyph.info[,"charSpace", drop = F]))
    }
    if (is.null(x$glyph.info$charKern)) {
      header$tf_CharKern <- 0
    } else {
      header$tf_CharKern <- max(c(header$tf_CharLoc, header$tf_CharSpace)) + addToPointer
      addToPointer       <- 2*prod(dim(x$glyph.info[,"charKern", drop = F]))
    }
    header$tf_CharData <- max(with(header, c(tf_CharLoc, tf_CharSpace, tf_CharKern))) + addToPointer
    
    ## create HUNK_RELOC32 trailing relocator hunk. Points to relative
    ## addresses that contain relative pointers and should be reallocated when loaded in memory
    ## they are the nodes, the tf_CharData, tf_CharLoc, tf_CharSpace and tf_CharKern
    ## (the latter two are optional, and should only be included when not equal to 0)
    trailing.hunks <- c(2, 7, 19, 21, if(header$tf_CharSpace != 0) 22, if(header$tf_CharKern != 0) 23)
    trailing.hunks <- cumsum(abs(.amiga.font.header$byte[-1]))[trailing.hunks - 1]
    trailing.hunks <- c(1004,                   # HUNK_RELOC32 id
                        length(trailing.hunks), # number of addresses to relocate
                        0,                      # hunk id number
                        trailing.hunks,         # adresses
                        0,                      # terminator (no more data follows)
                        1010)                   # stop loading HUNKS
    trailing.hunks <- .amigaIntToRaw(trailing.hunks, 32, F)
    
    header <- with(.amiga.font.header,
                   .write.amigaData(header, byte, signed, par.names))
    
    font.data <- .amigaIntToRaw(unlist(c(t(x$glyph.info[,c("glyphLocation", "glyphWidth")])), use.names = F), 16, F)
    
    if (!is.null(x$glyph.info$charSpace))
      font.data <- c(font.data, .amigaIntToRaw(x$glyph.info$charSpace, 16, T))
    if (!is.null(x$glyph.info$charKern))
      font.data <- c(font.data, .amigaIntToRaw(x$glyph.info$charKern, 16, T))
    
    palette   <- attr(x$bitmap, "palette")
    bm        <- apply(as.matrix(x$bitmap), 1, function(y) c(F, T)[match(y, palette)])
    bm        <- .bitmapToRaw(bm, invert.bytes = T, invert.longs = F)
    font.data <- c(header, font.data, bm)
    
    ## Add padding bytes to align the data along 32 bit.
    font.data <- font.data[1:(4*ceiling(length(font.data)/4))]
    
    ## specify in the leading hunks where the trailing hunks start:
    font.data[c(21:24, 29:32)] <- .amigaIntToRaw(ceiling((length(font.data) - 32)/4), 32, F)
    
    return(c(font.data, trailing.hunks))
  },
  warning=function(w) {
    if (startsWith(conditionMessage(w), "Replacement operator for AmigaBitmapFont"))
      invokeRestart("muffleWarning")
  })
  
}

#' @rdname as.raw
#' @name as.raw
#' @export
as.raw.AmigaBitmapFontSet <- function(x, ...) {
  withCallingHandlers({ ## TODO remove handlers when replace functions are implemented
    if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
    if (x$fch_FileID != "FontContents") stop("Sorry, currently only 'FontContents' font sets are supported.")
    ## put list in correct order...
    x <- x[c("fch_FileID", "fch_NumEntries", "FontContents")]
    .as.raw.FontContents <- function(y) {
      y <- y[c("fc_FileName", "fc_YSize", "fc_Style", "fc_Flags")]
      y$fc_FileName <- charToRaw(y$fc_FileName)[1:256]
      y$fc_YSize    <- .amigaIntToRaw(y$fc_YSize, 16, F)
      y$fc_Style    <-  .bitmapToRaw(y$fc_Style, F, F)
      y$fc_Flags    <-  .bitmapToRaw(y$fc_Flags, F, F)
      unlist(y, use.names = F)
    }
    x$FontContents   <- unlist(lapply(x$FontContents, .as.raw.FontContents), use.names = F)
    x$fch_FileID     <- c(0x0f00, 0x0f02, 0x0f03)[match(as.character(x$fch_FileID), .amiga.font.types)]
    x$fch_FileID     <- .amigaIntToRaw(x$fch_FileID, 16, F)
    x$fch_NumEntries <- .amigaIntToRaw(x$fch_NumEntries, 16, F)
    return(unlist(x, use.names = F))
  },
  warning=function(w) {
    if (startsWith(conditionMessage(w), "Replacement operator for AmigaBitmapFont"))
      invokeRestart("muffleWarning")
  })
}

#' @param text Text (a \code{character} string) to be formated
#' with \code{x} (when \code{x} is an \code{\link{AmigaBitmapFont}}
#' or an \code{\link{AmigaBitmapFontSet}}.
#' @param style Argument is only valid when \code{x} is an \code{\link{AmigaBitmapFont}}
#' or an \code{\link{AmigaBitmapFontSet}}. No styling is applied
#' when missing or \code{NULL}. One or more of the following styles
#' can be used '\code{bold}', '\code{italic} or '\code{underlined}'.
#' @param palette Argument is only valid when \code{x} is an \code{\link{AmigaBitmapFont}}
#' or an \code{\link{AmigaBitmapFontSet}}. Should be a \code{vector} of
#' two colours. The first is element is used as background colour, the
#' second as foreground. When missing, transparent white and black are used.
#' @family raster.operations
#' @rdname as.raster
#' @name as.raster
#' @export
as.raster.AmigaBitmapFont <- function(x, text, style, palette, ...) {
  if (!inherits(x, "AmigaBitmapFont")) stop("x should be of class AmigaBitmapFont.")
  if (missing(text)) return(x$bitmap)
  if (missing(style)) style <- NULL
  if (missing(palette)) palette <- NULL
  .print_to_raster(text, x, style, palette)
}

#' @family raster.operations
#' @rdname as.raster
#' @name as.raster
#' @export
as.raster.AmigaBitmapFontSet <- function(x, text, style, palette, ...) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  if (missing(text)) {
    dims <- do.call(rbind, lapply(x$FontContents, function(y) dim(y$BitmapFont$bitmap)))
    pals <- do.call(rbind, lapply(x$FontContents, function(y) attr(y$BitmapFont$bitmap, "palette")))
    if (!all(apply(pals, 2, function(y) all(y == y[1])))) stop("No consitent palettes used for bitmaps in the font set")
    result <- matrix(pals[1,1], sum(dims[,1]), max(dims[,2]))
    y.offset <- c(0, cumsum(dims[,1]))
    lapply(1:length(x$FontContents), function(y) {
      ys <- x$FontContents[[y]]$BitmapFont$tf_YSize
      xs <- ncol(x$FontContents[[y]]$BitmapFont$bitmap)
      result[y.offset[y] + (1:ys), 1:xs] <<- as.matrix(x$FontContents[[y]]$BitmapFont$bitmap)
    })
    result <- grDevices::as.raster(result)
    return(result)
  }
  if (missing(style)) style <- NULL
  if (missing(palette)) palette <- NULL
  .print_to_raster(text, x, style, palette)
}

#' Combine multiple AmigaFFH objects
#'
#' Use this function to correctly combine one or more \code{\link{AmigaBitmapFont}}
#' class objects into a single \code{\link{AmigaBitmapFontSet}} class
#' object, or to combine multiple \code{\link{AmigaBasic}} class objects.
#'
#' In case \code{...} are one or more \code{\link{AmigaBasic}} class objects:
#' 
#' \code{\link{AmigaBasic}} class objects are combined into a single
#' \code{\link{AmigaBasic}} class object in the same order as they
#' are given as argument to this function. for this purpose the lines of
#' Amiga Basic codes are simply concatenated.
#' 
#' In case \code{...} are one or more \code{\link{AmigaBitmapFont}} class objects:
#' 
#' \code{\link{AmigaBitmapFontSet}} class objects can hold multiple
#' \code{\link{AmigaBitmapFont}} class objects. Use this method to
#' combine font bitmaps into such a font set. Make sure each bitmap
#' represents a unique font height (in pixels). When heights are duplicated
#' an error will be thrown.
#' 
#' You can also specify a \code{name} for the font, that will be embeded
#' in the object. As this name will also be used as a file name when
#' writing the font to a file, make sure that it is a valid filename.
#'
#' @rdname c
#' @name c
#' @param ... Either \code{\link{AmigaBasic}} or \code{\link{AmigaBitmapFont}}
#' class objects. In case of \code{\link{AmigaBitmapFont}} objects:
#' Each \code{\link{AmigaBitmapFont}} object should have a
#' unique Y-size.
#' @param name This argument is only valid when \code{...} are one or more
#' \code{\link{AmigaBitmapFont}} class objects.
#' 
#' A \code{character} string specifying the name that needs to be
#' applied to the font set. When unspecified, the default name 'font' is
#' used. Note that this name will also be used as a file name when writing
#' the font to a file. So make sure the name is also a valid file name. This
#' will not be checked for you and may thus result in errors.
#' @return Returns an \code{\link{AmigaBitmapFontSet}} in which the
#' \code{\link{AmigaBitmapFont}} objects are combined. Or when \code{\link{AmigaBasic}}
#' objects are combined, an \code{\link{AmigaBasic}} object is returned
#' in which the lines of Amiga Basic code are combined.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## first get some AmigaBitmapFont objects:
#' font8 <- getAmigaBitmapFont(font_example, 8)
#' font9 <- getAmigaBitmapFont(font_example, 9)
#' 
#' ## now bind these bitmaps again in a single set
#' font.set <- c(font8, font9, name = "my_font_name")
#' 
#' ## Amiga Basic codes can also be combined:
#' bas1 <- as.AmigaBasic("LET a = 1")
#' bas2 <- as.AmigaBasic("PRINT a")
#' bas  <- c(bas1, bas2)
#' }
#' @family AmigaBitmapFont.operations
#' @author Pepijn de Vries
#' @export
c.AmigaBitmapFont <- function(..., name = "font") {
  fonts <- list(...)
  lapply(fonts, function(f) {
    if (!inherits(f, "AmigaBitmapFont")) stop("'...' should be all of class AmigaBitmapFont.")
  })
  sz <- unlist(lapply(fonts, function(x) x$tf_YSize))
  if (any(duplicated(sz))) stop("The font Y-sizes are not unique.")
  fonts <- fonts[order(sz)]
  result <- list(
    fch_FileID = factor(.amiga.font.types[1],
                        .amiga.font.types),
    fch_NumEntries = length(fonts),
    FontContents = lapply(fonts, function(x) {
      list(
        fc_FileName = paste(name, x$tf_YSize, sep = "/"),
        fc_YSize = x$tf_YSize,
        fc_Style = x$tf_Style,
        fc_Flags = x$tf_Flags,
        BitmapFont = x
      )
    })
  )
  names(result$FontContents) <- sprintf("pt%i", sort(sz))
  class(result) <- "AmigaBitmapFontSet"
  result
}

#' Extract or replace a font name
#'
#' Extract or replace a font name from an \code{\link{AmigaBitmapFontSet}}
#' object.
#'
#' The name of a font is embeded at multiple locations of an \code{\link{AmigaBitmapFontSet}}
#' object. This function can be used to extract or replace the font name
#' correctly. This is also the name that will be used when writing the
#' font to a file with \code{\link{write.AmigaBitmapFontSet}}.
#' @rdname fontName
#' @name fontName
#' @param x An \code{\link{AmigaBitmapFontSet}} for which the font name
#' needs to be changed.
#' @param value A \code{character} string specifying the name you
#' wish to use for the font.
#' @return Returns the font name. In case of the replace function, a copy
#' of \code{x} is returned with the name replaced by '\code{value}'.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## show the name of the example font:
#' fontName(font_example)
#' 
#' ## This is how you change the name into "foo"
#' fontName(font_example) <- "foo"
#' 
#' ## see it worked:
#' fontName(font_example)
#' }
#' @family AmigaBitmapFont.operations
#' @author Pepijn de Vries
#' @export
fontName <- function(x) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet")
  filenames <- unlist(lapply(x$FontContents, function(y) y$fc_FileName))
  filenames <- do.call(rbind, strsplit(filenames, "/"))
  if (ncol(filenames) != 2) stop("Unexpected file structure.")
  if (length(unique(filenames[,1])) != 1) stop("Not a single base name for the font.")
  as.character(filenames[1,1])
}

#' @rdname fontName
#' @name fontName<-
#' @export
`fontName<-` <- function(x, value) {
  withCallingHandlers({ ## TODO remove handlers when replace functions are implemented
    lapply(1:length(x$FontContents), function(i) {
      x$FontContents[[i]]$fc_FileName <<-
        gsub(".+?([/])", paste0(value, "/"),
             x$FontContents[[i]]$fc_FileName)
    })
    return(x)
  },
  warning=function(w) {
    if (startsWith(conditionMessage(w), "Replacement operator for AmigaBitmapFont"))
      invokeRestart("muffleWarning")
  })
}

#' Extract a specific AmigaBitmapFont from a AmigaBitmapFontSet
#'
#' Extract a specific \code{\link{AmigaBitmapFont}} from a
#' \code{\link{AmigaBitmapFontSet}}.
#'
#' An \code{\link{AmigaBitmapFontSet}} object can hold one or more
#' bitmaps for specific font sizes (heights). Use this function to
#' obtain such a specific \code{\link{AmigaBitmapFont}}.
#' @rdname getAmigaBitmapFont
#' @name getAmigaBitmapFont
#' @param x An \code{\link{AmigaBitmapFontSet}} object, from which the
#' specific \code{\link{AmigaBitmapFont}} object needs to be extracted.
#' @param size A single \code{numeric} value specifying the desired font
#' size in pixels. Use \code{\link{availableFontSizes}} to get available
#' sizes.
#' @return Returns an \code{\link{AmigaBitmapFont}} of the requested size.
#' An error is thrown when the requested size is not available.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## get the font object for the first available size:
#' font <- getAmigaBitmapFont(font_example,
#'                            availableFontSizes(font_example)[1])
#' }
#' @family AmigaBitmapFont.operations
#' @author Pepijn de Vries
#' @export
getAmigaBitmapFont <- function(x, size) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  if (length(size) > 1) {
    warning("Multiple sizes specified, using only the first value.")
    size <- size[1]
  }
  if (!(size %in% availableFontSizes(x))) stop(paste0("Requested font size (", size, ") not available."))
  return (x$FontContents[[paste0("pt", size)]]$BitmapFont)
}

#' Get available font sizes from an AmigaBitmapFontSet
#'
#' Get available font sizes (height) from an \code{\link{AmigaBitmapFontSet}} in pixels.
#'
#' An \code{\link{AmigaBitmapFontSet}} can hold bitmaps of multiple font
#' sizes. Use this function to obtain the available size from such a set.
#' @rdname availableFontSizes
#' @name availableFontSizes
#' @param x An \code{\link{AmigaBitmapFontSet}} for which the available
#' font sizes (height) in number of pixels need to be obtained.
#' @return Returns a \code{vector} of \code{numeric} values specifying
#' the available font sizes (height in pixels) for \code{x}.
#' @examples
#' \dontrun{
#' data(font_example)
#' 
#' ## The example font holds two font sizes (8 and 9):
#' availableFontSizes(font_example)
#' }
#' @family AmigaBitmapFont.operations
#' @author Pepijn de Vries
#' @export
availableFontSizes <- function(x) {
  if (!inherits(x, "AmigaBitmapFontSet")) stop("x should be of class AmigaBitmapFontSet.")
  as.numeric(unlist(lapply(x$FontContents, function(y) y$fc_YSize)))
}

#' Convert a raster image into an AmigaBitmapFont
#'
#' Convert a two-coloured \code{\link[grDevices]{as.raster}} image into
#' an \code{\link{AmigaBitmapFont}} class object.
#'
#' Create an \code{\link{AmigaBitmapFont}} class object by providing
#' a two-coloured raster image and specifying which characters are
#' depicted by the image.
#' @rdname rasterToAmigaBitmapFont
#' @name rasterToAmigaBitmapFont
#' @param x A \code{raster} (see grDevices package) object composed of
#' two colours only. Make sure that all glyphs (graphical representation
#' of characters) are next to eachother on a single line. The height
#' of this raster (in pixels) is taken automatically as font height.
#' @param glyphs Specify which glyphs are included in the image
#' \code{x} from left to right. It can be specified in one of the
#' following ways:
#' 
#' A single \code{character} string, where the length of the string
#' (\code{nchar}) equals the number of displayed glyphs in \code{x}.
#' 
#' A \code{vector} of \code{numeric} ASCII codes. The length of
#' the \code{vector} should equal the number of displayed glyphs
#' in \code{x}.
#' 
#' A \code{list} of either \code{character} strings or \code{vector} of
#' \code{numeric}s.  The length of the \code{list} should equal the
#' number of displayed glyphs in \code{x}. Each element can represent
#' multiple characters, meaning that the nth element of the list
#' uses the nth glyph shown in \code{x} to represent all the characters
#' included in that element.
#' 
#' Note that Amiga bitmap fonts represent ASCII characters and may
#' not include all special characters or symbols.
#' @param baseline The baseline of the font, specified in number of
#' pixels from the top (\code{numeric}). Should be a whole number
#' between 0 and the font height (height of \code{x}) minus 1.
#' @param default_glyph A single \code{character} or ASCII code
#' (\code{numeric}) that should be used by default. This means
#' that all characters that are not specified by \code{glyphs} will
#' be represented by this \code{default_glyph}. \code{default_glyph} should
#' be included in \code{glyphs}.
#' @param glyph_width A \code{numeric} \code{vector} with the same number
#' of elements or characters as used for \code{glyphs}. It specifies
#' the width in pixels for each glyph reserved in the raster image \code{x}.
#' They should be whole numbers greater or equal to 0.
#' @param glyph_space A \code{numeric} \code{vector} with the same number
#' of elements or characters as used for \code{glyphs}. It specifies
#' the width in pixels for each glyph that should be used when formatting.
#' text. Note that these values can be smaller or larger than the values
#' specified for \code{glyph_width}.
#' They should be whole numbers greater or equal to 0.
#' @param glyph_kern Not that in Amiga bitmap fonts not the formal
#' definition from typography is used for kerning. Here, kerning is
#' used as the number of pixels the cursor should be moved forward or
#' backward after typesetting a character. It should be a
#' \code{numeric} \code{vector} with the same number of elements or
#' characters as used for \code{glyphs}. It can hold both positive
#' and negative values.
#' @param palette A \code{vector} of two colours. Both colours should
#' be in \code{x}. The first colour is used as background colour,
#' the second as foreground colour.
#' 
#' When missing, it will be checked whether \code{x} has a palette
#' as attribute, and uses that. If that attribute is also missing,
#' the palette will be guessed from \code{x}, where the most
#' frequently occurring colour is assumed to be the background
#' colour.
#' @param ... Currently ignored.
#' @return Returns a \code{\link{AmigaBitmapFont}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' data("font_example")
#' 
#' ## make a raster that we can use to create a bitmap font
#' font9.rast <- as.raster(getAmigaBitmapFont(font_example, 9))
#' 
#' ## note the glyphs and the order in which they are included in
#' ## the raster image:
#' plot(font9.rast)
#' 
#' ## let's build a simple font, using only the first few glyphs
#' ## in the raster:
#' font9 <- rasterToAmigaBitmapFont(
#'   ## 'x' needs the raster image:
#'   x             = font9.rast,
#'   
#'   ## 'glyphs' are the graphical representation of the characters
#'   ## that we will include in our font. We will only use the
#'   ## first 7 characters in the raster image:
#'   glyphs        = " !\"#$%&",
#'   
#'   ## We will use the '&' glyph to represent all characters that
#'   ## are not specified in the font:
#'   default_glyph = "&",
#'   
#'   ## The raster image is 9 pixels tall, as will be the font.
#'   ## Let's use 7 as the base (it needs to be less than the height)
#'   baseline      = 7,
#'   
#'   ## Let's define the width in pixels for each of the 7
#'   ## characters. This is their width in the raster image:
#'   glyph_width   = c(0, 1, 3, 6, 5, 5, 5),
#'   
#'   ## Let's define the space the character should take in pixels
#'   ## when it is used to format text:
#'   glyph_space   = c(4, 2, 4, 7, 6, 6, 6),
#'   
#'   ## the raster uses white as background colour and black as
#'   ## foreground:
#'   palette       = c("white", "black")
#' )
#' 
#' ## note that for all characters that are not specified,
#' ## the default glyph ('&') is used:
#' plot(font9, text = "!@#$%ABCD")
#' 
#' ## Let's take a subset from the font's bitmap (rasteer):
#' font9abc.rast <- font9.rast[,263:282]
#' 
#' ## as you can see this bitmap only contains the lowercase
#' ## characters 'a', 'b', 'c', 'd' and 'e':
#' plot(font9abc.rast)
#' 
#' font9.abc <- rasterToAmigaBitmapFont(
#'   x             = font9abc.rast,
#'   ## Each glyph in the image can be represented by a single
#'   ## element in a list. By specifying multiple characters in
#'   ## each element, you can recycle a glyph to represent different
#'   ## characters. So in this case, the glyph 'a' is used for
#'   ## all the accented variants of the character 'a'.
#'   glyphs        = list("a\xE0\xE1\xE2\xE3\xE4\xE5",
#'                        "b",
#'                        "c\xA2\xE7",
#'                        "d",
#'                        "e\xE8\xE9\xEA\xEB"),
#'   default_glyph = "c", ## 'c' is used as default glyph for all other characters
#'   baseline      = 7,
#'   glyph_width   = c(4, 4, 4, 4, 4),
#'   glyph_space   = c(5, 5, 5, 5, 5),
#'   palette       = c("white", "black")
#' )
#' 
#' ## see what happens when you format text using the font we just created:
#' plot(font9.abc, text = "a\xE0\xE1\xE2\xE3\xE4\xE5bc\xA2\xE7de\xE8\xE9\xEA\xEB, foo bar")
#' }
#' @family AmigaBitmapFont.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
rasterToAmigaBitmapFont <- function(x, glyphs, default_glyph, baseline, glyph_width, glyph_space, glyph_kern, palette, ...) {
  glyph_width <- round(glyph_width)
  if (glyph_width < 0 || glyph_width > 65535) stop("'glyph_width' out of range (0, 65535).")
  if (baseline < 0 || baseline > (nrow(x) - 1) || baseline != round(baseline)) stop("'baseline' should be whole number between 0 and tf_YSize - 1.")
  if (is.character(default_glyph)) default_glyph <- utf8ToInt(enc2utf8(default_glyph))
  if (length(default_glyph) != 1) stop("'default_glyph' should have a length of 1.")
  if (is.list(glyphs)) {
    test.default <- F
    glyphs <- lapply(1:length(glyphs), function(i) {
      if (is.character(glyphs[[i]])) glyphs[[i]] <- utf8ToInt(enc2utf8(glyphs[[i]]))
      if (default_glyph %in% glyphs[[i]]) {
        default_glyph <<- i
        test.default <<- T
      }
      data.frame(
        idx = i,
        glyphs = glyphs[[i]]
      )
    })
    glyphs <- do.call(rbind, glyphs)
    if (!test.default) stop("'default_glyph' should be in 'glyphs'.")
  }
  if (is.character(glyphs)) glyphs <- utf8ToInt(enc2utf8(glyphs))
  if (is.numeric(glyphs)) {
    default_glyph <- which(glyphs == default_glyph)
    if (length(default_glyph) == 0) stop("'default_glyph' should be in 'glyphs'.")
    glyphs <- data.frame(idx = 1:length(glyphs), glyphs)
  }
  if (any(duplicated(glyphs$glyphs))) stop("Can't handle duplicated characters or ascii codes in 'glyphs'.")
  char_lo <- min(glyphs$glyphs)
  char_hi <- max(glyphs$glyphs)
  if (char_hi == glyphs$glyphs[default_glyph]) char_hi <- char_hi - 1
  if (!inherits(x, "raster")) stop("'x' should be of class raster.")
  if (char_lo < 0 || char_lo > 255 || char_hi < 0 || char_hi > 255) stop("ASCII codes for 'glyphs' are out of range (0-255).")
  if (sum(glyph_width) > dim(x)[2]) stop("Sum of char width is wider than the provided raster image.")
  if (baseline > (dim(x)[1] - 1)) stop("'baseline' should not be greater then the height of 'x' minus 1.")
  if (missing(palette)) {
    if (is.null(attr(x, "palette"))) {
      ## If a palette is missing, take a guess based on the raster
      ## assume that the most frequent colour is the background colour
      palette <- table(x)
      palette <- names(palette)[order(-palette)]
      if (length(palette) != 2) stop("'x' does not contain 2 unique values/colours.")
    } else {
      palette <- attr(x, "palette")
    }
  }
  if (any(!(unique(x) %in% palette)) || length(palette) != 2) stop ("'palette' doesn't specify two colours, or 'x' contains different colours.")
  attr(x, "palette") <- palette
  font.result              <- as.list(rep(0, length(.amiga.font.header$par.names)))
  names(font.result)       <- .amiga.font.header$par.names
  font.result$leadingHunks <- NULL
  font.result$node.disklink <- font.result$node.message <- list(
    ln_Succ = 0,
    ln_Pred = 0,
    ln_Type = factor("NT_FONT", .amiga.node.types),
    ln_Pri = 0,
    ln_Name = 26
  )
  font.result$fontName     <- ""
  font.result$dfh_FileID   <- factor("DFH_ID", "DFH_ID")
  font.result$dfh_Revision <- 1
  font.result$tf_LoChar    <- char_lo
  font.result$tf_HiChar    <- char_hi
  font.result$tf_Modulo    <- ceiling(dim(x)[2]/8)
  if (dim(x)[2] != 8*font.result$tf_Modulo) {
    x <- cbind(as.matrix(x), matrix(palette[1],
                                    dim(x)[1],
                                    8*font.result$tf_Modulo - dim(x)[2]))
    x <- as.raster(x)
    attr(x, "palette") <- palette
  }
  font.result$tf_YSize     <- dim(x)[1]
  font.result$tf_Baseline  <- baseline
  font.result$tf_XSize     <- stats::median(glyph_width)
  font.result$tf_BoldSmear <- 1
  font.result$tf_Style     <- rep(F, 8)
  font.result$tf_Flags     <- c(F, T, F, F, F, T, T, F)
  names(font.result$tf_Style) <- c("UNDERLINED", "BOLD", "ITALIC", "EXTENDED", "RESERVED1",
                                   "RESERVED2", "COLORFONT", "TAGGED")
  names(font.result$tf_Flags) <- c("ROMFONT", "DISKFONT", "REVPATH", "TALLDOT",
                                   "WIDEDOT", "PROPORTIONAL", "DESIGNED", "REMOVED")
  glyphs <- merge(glyphs, data.frame(glyphs = char_lo:(char_hi + 1)), all.x = T, all.y = T)
  glyphs$idx[is.na(glyphs$idx)] <- default_glyph
  loc <- cumsum(c(0, glyph_width))
  font.result$glyph.info <- data.frame (
    glyphLocation = loc[glyphs$idx],
    glyphWidth    = glyph_width[glyphs$idx]
  )
  offs <- 0
  font.result$tf_CharLoc   <- 110
  if (!missing(glyph_space)) {
    glyph_space <- round(glyph_space)
    if (glyph_space < 0 || glyph_space > 65535) stop("'glyph_space' out of range (0, 65535).")
    font.result$glyph.info$charSpace <- glyph_space[glyphs$idx]
    font.result$tf_CharSpace <- 110 + 2*2*nrow(font.result$glyph.info)
    offs <- offs + 1
  }
  if (!missing(glyph_kern)) {
    glyph_kern  <- round(glyph_kern)
    if (glyph_kern < -32768 || glyph_kern > 32767) stop("'glyph_kern' out of range (-32768, 32767).")
    font.result$glyph.info$charKern <- glyph_kern[glyphs$idx]
    font.result$tf_CharKern <- 110 + (2 + offs)*2*nrow(font.result$glyph.info)
    offs <- offs + 1
  } 
  font.result$tf_CharData  <- 110 + (2 + offs)*2*nrow(font.result$glyph.info)
  font.result$bitmap       <- x
  class(font.result) <- "AmigaBitmapFont"
  return(font.result)
}

#' @export
`$<-.AmigaBitmapFont` <- function(x, i, value) {
  x[[i]] <- value
  x
}

#' @export
`[[<-.AmigaBitmapFont` <- function(x, i, value) {
  cl <- class(x)
  class(x) <- NULL
  x[[i]] <- value
  class(x) <- cl
  ## TODO update this replacement function and remove warning
  warning(paste0("Replacement operator for AmigaBitmapFont objects ",
                 "will be modified in future versions of this package. ",
                 "Note that not all replacement operations may be ",
                 "allowed in future versions of this package."))
  x
}

#' @export
`$<-.AmigaBitmapFontSet` <- function(x, i, value) {
  x[[i]] <- value
  x
}

#' @export
`[[<-.AmigaBitmapFontSet` <- function(x, i, value) {
  cl <- class(x)
  class(x) <- NULL
  x[[i]] <- value
  class(x) <- cl
  ## TODO update this replacement function and remove warning
  warning(paste0("Replacement operator for AmigaBitmapFontSet objects ",
                 "will be modified in future versions of this package. ",
                 "Note that not all replacement operations may be ",
                 "allowed in future versions of this package."))
  x
}
