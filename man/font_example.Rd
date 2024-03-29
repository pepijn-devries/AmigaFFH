% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{font_example}
\alias{font_example}
\title{An example object for the AmigaBitmapFontSet class}
\format{
{
\code{font_example} is an \code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}}
object. For details see the object class documentation.}
}
\description{
An example object for the \code{\link[=AmigaBitmapFontSet]{AmigaBitmapFontSet()}} class used in
examples throughout this package. It also contains a nested
\code{\link[=AmigaBitmapFont]{AmigaBitmapFont()}} class objects, which can be obtain by
using \code{getAmigaBitmapFont(font_example, 9)}.
}
\details{
The \code{font_example} contains a font that was designed as an example
for this package. It holds bitmap glyphs for 8 and 9 pixels tall
characters.
}
\examples{
data("font_example")
}
\seealso{
Other AmigaBitmapFont.operations: 
\code{\link{AmigaBitmapFont}},
\code{\link{availableFontSizes}()},
\code{\link{c}()},
\code{\link{fontName}()},
\code{\link{getAmigaBitmapFont}()},
\code{\link{rasterToAmigaBitmapFont}()},
\code{\link{rawToAmigaBitmapFontSet}()},
\code{\link{rawToAmigaBitmapFont}()},
\code{\link{read.AmigaBitmapFontSet}()},
\code{\link{read.AmigaBitmapFont}()},
\code{\link{write.AmigaBitmapFont}()}
}
\concept{AmigaBitmapFont.operations}
