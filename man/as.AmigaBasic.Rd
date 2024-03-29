% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic.r
\name{as.AmigaBasic}
\alias{as.AmigaBasic}
\title{Coerce raw or character data to an AmigaBasic class object}
\usage{
as.AmigaBasic(x, ...)
}
\arguments{
\item{x}{\code{x} should be a \code{vector} of \code{raw} data or
\code{character} strings. When \code{x} is \code{raw} data, it
is interpreted as if it where from an Amiga Basic binary encoded file.

When \code{x} is a \code{vector} of \code{character} strings,
each element of the vector should represent one line of Basic code.
Each line should not contain line break or other special characters,
as this will result in errors. The text should represent valid
Amiga Basic syntax. The syntax is only checked to a limited extent as
this package does not implement an interpreter for the code.}

\item{...}{Currently ignored.}
}
\value{
Returns an \code{\link[=AmigaBasic]{AmigaBasic()}} class object based on \code{x}.
}
\description{
Coerce raw or character data to an \code{\link[=AmigaBasic]{AmigaBasic()}} S3 class object
}
\details{
Convert text to an \code{\link[=AmigaBasic]{AmigaBasic()}} S3 class object. The text should
consist of valid Amiga BASIC syntaxis. This function does not perform a
full check of the syntaxis, but will break on some fundamental syntaxis malformations
}
\examples{
\dontrun{
## An AmigaBasic object can be created from text.
## Note that each line of code is a seperate element
## in the vector:
bas <- as.AmigaBasic(c(
  "CLS ' Clear the screen",
  "PRINT \"Hello world!\" ' Print a message on the screen"
))

## Let's make it raw data:
bas.raw <- as.raw(bas)

## We can also use the raw data to create an Amiga Basic object:
## Note that this effectively the same as calling 'rawToAmigaBasic'
bas <- as.AmigaBasic(bas.raw)
}
}
\references{
\url{https://en.wikipedia.org/wiki/AmigaBASIC}
}
\seealso{
Other AmigaBasic.operations: 
\code{\link{AmigaBasic.reserved}()},
\code{\link{AmigaBasicBMAP}},
\code{\link{AmigaBasic}},
\code{\link{[.AmigaBasic}()},
\code{\link{as.AmigaBasicBMAP}()},
\code{\link{as.character}()},
\code{\link{check.names.AmigaBasic}()},
\code{\link{names.AmigaBasic}()},
\code{\link{rawToAmigaBasicBMAP}()},
\code{\link{rawToAmigaBasic}()},
\code{\link{read.AmigaBasicBMAP}()},
\code{\link{read.AmigaBasic}()},
\code{\link{write.AmigaBasic}()}

Other raw.operations: 
\code{\link{as.raw.AmigaBasic}()},
\code{\link{colourToAmigaRaw}()},
\code{\link{packBitmap}()},
\code{\link{rawToAmigaBasicBMAP}()},
\code{\link{rawToAmigaBasicShape}()},
\code{\link{rawToAmigaBasic}()},
\code{\link{rawToAmigaBitmapFontSet}()},
\code{\link{rawToAmigaBitmapFont}()},
\code{\link{rawToAmigaIcon}()},
\code{\link{rawToHWSprite}()},
\code{\link{rawToIFFChunk}()},
\code{\link{rawToSysConfig}()},
\code{\link{simpleAmigaIcon}()}
}
\author{
Pepijn de Vries
}
\concept{AmigaBasic.operations}
\concept{raw.operations}
