% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitmapfont.r
\name{write.AmigaBitmapFont}
\alias{write.AmigaBitmapFont}
\alias{write.AmigaBitmapFontSet}
\title{Write an AmigaBitmapFont(set) file}
\usage{
write.AmigaBitmapFont(x, file)

write.AmigaBitmapFontSet(x, path = getwd())
}
\arguments{
\item{x}{Respectively an \code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}} or a
\code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}} object depending on which of the
write-functions is called. This is the object that will be written
to the specified file.}

\item{file}{A \code{character} string specifying the file location
to which \code{x} (an \code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}} object) needs to be written.
It is common practice on the Amiga to use the font height in pixels as
file name.}

\item{path}{A \code{character} string specifying the path where
\code{x} (an \code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}} object) needs to be stored.
The filename for the font set will be extracted from \code{x} using
\code{\link[=fontName]{fontName()}} followed by the *.font extension. A subdirectory
will be created with the same name (without the extension) if it doesn't
already exists. In this subdirectory all the nested \code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}}
objects are stored.}
}
\value{
Invisibly returns the result of the call of \code{close} to the
file connection.
}
\description{
Functions to write \code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}} and \code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}}
class objects to files.
}
\details{
\code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}} class objects are written to a *.font
file. The filename used for this purpose is obtained from the object
itself using \code{\link[=fontName]{fontName()}}. In addition, a subdirectory is
created automatically (when it doesn't already exist)
to which al the separate bitmap images for each font height are written
to individual files.

\code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}} class objects can also be written to a
file. In order to use it on a Commodore Amiga or emulator, it is better
to embed the font bitmap in a font set (using \code{\link[=c]{c()}})
and write the set to corresponding files.
}
\examples{
## obtain a bitmap font set:
data(font_example)

## write the font set to their files. The file name
## is extracted from the font object, so you only have
## to provide the path:
write.AmigaBitmapFontSet(font_example, tempdir())

## extract a font bitmap:
font <- getAmigaBitmapFont(font_example, 9)

## and write it to the temp dir:
write.AmigaBitmapFont(font, file.path(tempdir(), "9"))

## The following examples require the 'adfExplorer' package:
if (requireNamespace("adfExplorer")) {
  library("adfExplorer")
  virtual_disk_file <- tempfile(fileext = ".adf") |>
    create_adf_device(write_protected = FALSE) |>
    prepare_adf_device("font_disk") |>
    make_adf_dir("FONTS")
  
  dest <- virtual_path(virtual_disk_file, "DF0:FONTS")
  write.AmigaBitmapFontSet(font_example, dest)
  close(virtual_disk_file)
}
}
\seealso{
Other AmigaBitmapFont.operations: 
\code{\link{AmigaBitmapFont}},
\code{\link{availableFontSizes}()},
\code{\link{c}()},
\code{\link{fontName}()},
\code{\link{font_example}},
\code{\link{getAmigaBitmapFont}()},
\code{\link{rasterToAmigaBitmapFont}()},
\code{\link{rawToAmigaBitmapFontSet}()},
\code{\link{rawToAmigaBitmapFont}()},
\code{\link{read.AmigaBitmapFontSet}()},
\code{\link{read.AmigaBitmapFont}()}

Other io.operations: 
\code{\link{read.AmigaBasicBMAP}()},
\code{\link{read.AmigaBasicShape}()},
\code{\link{read.AmigaBasic}()},
\code{\link{read.AmigaBitmapFontSet}()},
\code{\link{read.AmigaBitmapFont}()},
\code{\link{read.AmigaIcon}()},
\code{\link{read.SysConfig}()},
\code{\link{read.iff}()},
\code{\link{write.AmigaBasicShape}()},
\code{\link{write.AmigaBasic}()},
\code{\link{write.AmigaIcon}()},
\code{\link{write.SysConfig}()},
\code{\link{write.iff}()}
}
\author{
Pepijn de Vries
}
\concept{AmigaBitmapFont.operations}
\concept{io.operations}
