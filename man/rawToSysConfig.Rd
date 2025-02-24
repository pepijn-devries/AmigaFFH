% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/system_configuration.r
\name{rawToSysConfig}
\alias{rawToSysConfig}
\title{Coerce raw data into a SysConfig class object}
\usage{
rawToSysConfig(x)
}
\arguments{
\item{x}{A vector of \code{raw} data that needs to be converted into an S3
\link{SysConfig} class object. It should have a length of at least 232. Although
system-configurations can be extended, such extended files are not supported here.}
}
\value{
Returns a \link{SysConfig} class object based on \code{x}.
}
\description{
\link{SysConfig} objects are comprehensive representations of binary Amiga
system-configuration files. Use this function to convert \code{raw} data from
such a file to a \link{SysConfig} object.
}
\details{
The Amiga used the system-configuration file to store certain system preferences
in a binary file. With this function such \code{raw} data can be converted into
a more comprehensive \link{SysConfig} object. Use \code{\link[=as.raw]{as.raw()}}
to achieve the inverse.
}
\examples{
if (requireNamespace("adfExplorer", quietly = TRUE)) {
  library(adfExplorer)
  
  ## get the system-configuration from the adfExplorer example disk:
  disk <- connect_adf(
    system.file("example.adz", package = "adfExplorer")
  )
  virtual_file_con <- adf_file_con(disk, "devs/system-configuration")
  sc <- readBin(virtual_file_con, "raw", 1024)
  close(disk)

  ## This will get you the raw data from the file:
  typeof(sc)

  ## Convert the raw data to a more comprehensive named list (and S3 SysConfig class):
  sc <- rawToSysConfig(sc)
}
}
\seealso{
Other SysConfig.operations: 
\code{\link{SysConfig}},
\code{\link{read.SysConfig}()},
\code{\link{simpleSysConfig}()},
\code{\link{write.SysConfig}()}

Other raw.operations: 
\code{\link{as.AmigaBasic}()},
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
\code{\link{simpleAmigaIcon}()}
}
\author{
Pepijn de Vries
}
\concept{SysConfig.operations}
\concept{raw.operations}
