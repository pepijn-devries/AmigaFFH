#' The S3 AmigaIcon class
#' 
#' A comprehensive representation of an Amiga Workbench icon file.
#' 
#' Files, directories and other similar objects were depicted as icons on the
#' Amiga Workbench (the Amiga's equivalent of what is now mostly known as the compmuter's
#' desktop). Icons were actually separate files with the exact same name as the
#' file or directory it represents, except for an additional `.info' extension.
#' 
#' In addition of being a graphical representation of files or directories, icon
#' files also contained additional information about the file. It could for instance
#' indicate which tool would be required to open the file.
#' 
#' The classic Amiga Workbench icon file has a rather complex structure as it is
#' basically a dump of how it is stored in memory. As a result it contains many
#' memory pointers that are really not necassary to store in a file.
#' 
#' The S3 AmigaIcon class is used to represent these complex files as a named
#' \code{list}. The elements in that \code{list} have mostly identical
#' names as listed in the document at the top referenced below. The names are usually
#' self-explenatory, but the referred documents can also be
#' consulted to obtain more detailed information with respect to each of
#' these elements. As pointed out ealier, not all elements will have a meaningful
#' use.
#' 
#' It is possible to change the values of the list, but not all values may be valid.
#' Note that they will not be fully checked for validity. Invalid values may result in errors
#' when writing to a binary file using \code{\link{write.AmigaIcon}}, or may simply not
#' work properly on an Amiga or in an emulator.
#' 
#' The original `.info' file could be extended with NewIcon or with an OS3.5
#' \code{\link{IFFChunk}} data, that allowed for icons with larger colour depths.
#' These extensions are currently not implemented.
#' 
#' Use \code{\link{simpleAmigaIcon}} for creating a simple \code{AmigaIcon} object which can
#' be modified. Use \code{\link{read.AmigaIcon}} to read, and \code{\link{write.AmigaIcon}}
#' to write workbench icon files (*.info). With \code{\link{rawToAmigaIcon}} and
#' \code{\link[AmigaFFH]{as.raw}} \code{AmigaIcon} can be coerced back and forth from
#' and to its raw (binary) form.
#' @docType class
#' @name AmigaIcon
#' @rdname AmigaIcon
#' @family AmigaIcon.operations
#' @author Pepijn de Vries
#' @references
#' \url{http://www.evillabs.net/wiki/index.php/Amiga_Icon_Formats}
#' \url{http://fileformats.archiveteam.org/wiki/Amiga_Workbench_icon}
#' \url{http://krashan.ppa.pl/articles/amigaicons/}
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Libraries_Manual_guide/node0241.html}
#' \url{http://amigadev.elowar.com/read/ADCD_2.1/Includes_and_Autodocs_3._guide/node05D6.html}
NULL

.icon.data.head <- data.frame(
  byte      = c(2, 2, -44, 1, 1, 4, 4, 4, 4, 4, 4, 4),
  signed    = c(F, F,   F, F, 1, F, F, T, T, F, F, T),
  par.names = c("ic_Magic", "ic_Version", "ic_Gadget", "ic_Type", "ic_Pad",
                "ic_DefaultTool", "ic_ToolTypes", "ic_CurrentX", "ic_CurrentY",
                "ic_DrawerData", "ic_ToolWindow", "ic_StackSize"),
  stringsAsFactors = F
)

.icon.gadget.data <- data.frame(
  byte      = c(4, 2, 2, 2, 2, -2, 2, 2, 4, 4, 4, 4, 4, 2, 4),
  signed    = c(F, T, T, T, T,  F, F, F, F, F, F, T, F, F, F),
  par.names = c("ga_Next", "ga_LeftEdge", "ga_TopEdge", "ga_Width", "ga_Height", "ga_Flags",
                "ga_Activation", "ga_GadgetType", "ga_GadgetRender", "ga_SelectRender", 
                "ga_GadgetText", "ga_MutualExclude", "ga_SpecialInfo", "ga_GadgetID",
                "ga_UserData"),
  stringsAsFactors = F
)

.icon.drawer.data <- data.frame(
  byte      = c(-48, 4, 4),
  signed    = c(  F, T, T),
  par.names = c("NewWindow", "dd_CurrentX", "dd_CurrentY"),
  stringsAsFactors = F
)

.icon.new.window.data <- data.frame(
  byte      = c(2, 2, 2, 2, 1, 1, -4, -4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2),
  signed    = c(T, T, T, T, F, F,  F,  F, F, F, F, F, F, T, T, F, F, F),
  par.names = c("nw_LeftEdge", "nw_TopEdge", "nw_Width", "nw_Height", "nw_DetailPen",
                "nw_BlockPen", "nw_IDCMPFlags", "nw_Flags", "nw_FirstGadget",
                "nw_CheckMark", "nw_Title", "nw_Screen", "nw_BitMap", "nw_MinWidth",
                "nw_MinHeight", "nw_MaxWidth", "nw_MaxHeight", "nw_Type"),
  stringsAsFactors = F
)

.icon.image.data <- data.frame(
  byte      = c(2, 2, 2, 2, 2, 4, 1, 1, 4),
  signed    = c(T, T, T, T, T, F, F, F, F),
  par.names = c("im_LeftEdge", "im_TopEdge", "im_Width", "im_Height", "im_Depth",
                "im_Data", "im_PlanePick", "im_PlaneOnOff", "im_Next"),
  stringsAsFactors = F
)

#' Create simple AmigaIcon objects
#'
#' Graphical representation of files and directories (icons) are stored as
#' separate files (with the .info extension) on the Amiga. This function writes
#' \code{\link{AmigaIcon}} class objects to such files.
#'
#' This function creates basic \code{\link{AmigaIcon}} objects which
#' can be modified afterwards. It uses simple generic images to represent
#' different types of files or directories.
#'
#' @rdname simpleAmigaIcon
#' @name simpleAmigaIcon
#' @param version A \code{character} string indicating the Amiga OS version
#' with which the icon should be compatible. "\code{OS2.x}" indicates
#' >=OS2.0 and "\code{OS1.x}" indicates <OS2.0.
#' @param type A \code{character} string indicating the type of object (file, disk, directory, etc.)
#' the icon should represent. See the `Usage' section for all posible options.
#' @param two.images A single \code{logical} value, indicating whether
#' the selected icon is depicted as a second image (in which case the
#' icon contains two images). The default value is \code{TRUE}.
#' @param back.fill A single \code{logical} value, indicating whether
#' the selected image of the icon should use the `back fill' mode (default).
#' If set to \code{FALSE} `complement' mode is used. Note that
#' back fill is not compatible when the icon holds two images. In the
#' `complement' mode, the image colours are inverted when selected.
#' In the `back fill' exterior first colour is not inverted.
#' @param ... Reserverd for additional arguments. Currently ignored.
#' @return Returns a simple S3 object of class \code{\link{AmigaIcon}}.
#' @examples
#' \dontrun{
#' ## Create an AmigaIcon object using the default arguments:
#' icon <- simpleAmigaIcon()
#' }
#' @family AmigaIcon.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
simpleAmigaIcon <- function(version    = c("OS1.x", "OS2.x"),
                            type       = c("WBDISK", "WBDRAWER", "WBTOOL", "WBPROJECT", "WBGARBAGE", "WBDEVICE", "WBKICK", "WBAPPICON"),
                            two.images = TRUE,
                            back.fill  = FALSE,
                            ...) {
  two.images <- as.logical(two.images[[1]])
  back.fill  <- as.logical(back.fill[[1]])
  if (back.fill && two.images) {
    warning("Backfill is not possible when a second image is used for the selected icon. Backfill is set to FALSE.")
    back.fill <- F
  }
  version          <- match.arg(version)
  type             <- match.arg(type)
  result           <- as.list(c(0xE310, 1, rep(0, nrow(.icon.data.head) - 2)))
  names(result)    <- .icon.data.head$par.names
  result$ic_Type   <- factor(type, c("WBDISK", "WBDRAWER", "WBTOOL", "WBPROJECT", "WBGARBAGE", "WBDEVICE", "WBKICK", "WBAPPICON"))
  result$ic_Gadget <- as.list(rep(0, nrow(.icon.gadget.data)))
  names(result$ic_Gadget)          <- .icon.gadget.data$par.names
  result$ic_Gadget$ga_UserData     <- factor(version, c("OS1.x", "OS2.x"))
  result$ic_Gadget$ga_GadgetRender <- 1
  result$ic_Gadget$ga_SelectRender <- as.numeric(two.images)
  result$ic_Gadget$ga_Flags        <- c(back.fill, two.images, T, rep(F, 13))
  names(result$ic_Gadget$ga_Flags) <- c("BACKFILL_MODE", "TWO_IMAGE_MODE", "IMAGE_MODE", paste0("UNDEFINED", 1:13))
  if (type %in% c("WBDISK", "WBDRAWER", "WBGARBAGE")) {
    result$ic_DrawerData <- 1
    result$drawer <- list(
      NewWindow = NULL,
      dd_CurrentX = 0,
      dd_CurrentY = 0
    )
    result$drawer$NewWindow               <- as.list(c(0, 0, 400, 100, 255, 255, 0, 0, 0, 0, 0, 0, 0, 90, 65, 640, 200, 1))
    names(result$drawer$NewWindow)        <- .icon.new.window.data$par.names
    result$drawer$NewWindow$nw_IDCMPFlags <- rep(F, 32)
    result$drawer$NewWindow$nw_Flags      <- rep(F, 32)
  } else {
    result[["drawer"]] <- list(NewWindow = list())
  }
  make_img <- function(type, idx) {
    img <- as.list(rep(0, nrow(.icon.image.data)))
    names(img) <- .icon.image.data$par.names
    img$im_Bitmap <- icon.images[[paste0("project", idx)]]
    if (type == "WBDISK") {
      img$im_Bitmap <- icon.images[[paste0("disk", idx)]]
    } else if (type == "WBDRAWER") {
      img$im_Bitmap <- icon.images[[paste0("drawer", idx)]]
    } else if (type == "WBTOOL") {
      img$im_Bitmap <- icon.images[[paste0("tool", idx)]]
    } else if (type == "WBGARBAGE") {
      img$im_Bitmap <- icon.images[[paste0("garbage", idx)]]
    }
    ## use 4 colour palette for default icons
    attributes(img$im_Bitmap)$palette <- attributes(img$im_Bitmap)$palette[1:4]
    img$im_Depth  <- log2(length(attributes(img$im_Bitmap)$palette))
    img$im_Width  <- dim(img$im_Bitmap)[2]
    img$im_Height <- dim(img$im_Bitmap)[1]
    img$im_Data   <- 1
    attributes(img$im_Bitmap)$bitmap.size <-
      sum(abs(.icon.image.data$byte)) + 16*ceiling(img$im_Width/16)*img$im_Height*img$im_Depth/8
    attribs <- attributes(img$im_Bitmap)
    if (version == "OS1.x") {
      img$im_Bitmap[img$im_Bitmap == AmigaFFH::amiga_palettes$wb.os2[4]] <-
        AmigaFFH::amiga_palettes$wb.os1[[4]]
      img$im_Bitmap[img$im_Bitmap == AmigaFFH::amiga_palettes$wb.os2[2]] <-
        AmigaFFH::amiga_palettes$wb.os1[[3]]
      img$im_Bitmap[img$im_Bitmap == AmigaFFH::amiga_palettes$wb.os2[3]] <-
        AmigaFFH::amiga_palettes$wb.os1[[2]]
      img$im_Bitmap[img$im_Bitmap == AmigaFFH::amiga_palettes$wb.os2[1]] <-
        AmigaFFH::amiga_palettes$wb.os1[[1]]
      attribs$palette <- AmigaFFH::amiga_palettes$wb.os1
    }
    attributes(img$im_Bitmap) <- attribs
    img$im_PlanePick <- length(attributes(img$im_Bitmap)$palette) - 1
    img
  }
  result$firstImage <- make_img(type, 1)
  if (two.images) result$secondImage <- make_img(type, 2)
  result$ic_Gadget$ga_Width   <- result$firstImage$im_Width
  result$ic_Gadget$ga_Height  <- result$firstImage$im_Height
  if (two.images) {
    if (result$ic_Gadget$ga_Width < result$secondImage$im_Width)
      result$ic_Gadget$ga_Width   <- result$secondImage$im_Width
    if (result$ic_Gadget$ga_Height < result$secondImage$im_Height)
      result$ic_Gadget$ga_Height   <- result$secondImage$im_Height
  }
  result$defaultTool          <- ""
  result$toolTypes            <- ""
  result$toolWindow           <- ""
  result$dd_Flags             <- factor(NULL, c("DDFLAGS_SHOWDEFAULT", "DDFLAGS_SHOWICONS",
                                                "DDFLAGS_SHOWALL"))
  result$dd_ViewModes         <- factor(NULL, c("DDVM_BYDEFAULT", "DDVM_BYICON", "DDVM_BYNAME",
                                                "DDVM_BYDATE", "DDVM_BYSIZE", "DDVM_BYTYPE"))
  class(result)               <- "AmigaIcon"
  return(result)
}

#' Coerce raw data into an AmigaIcon class object
#'
#' \code{\link{AmigaIcon}} objects are comprehensive representations of binary Amiga
#' Workbench icon files (*.info). Use this function to convert \code{raw} data from
#' such a file to an \code{\link{AmigaIcon}} object.
#'
#' Icons files (*.info) were used as a grapical representations of files and
#' directories on the Commodore Amiga. This function will convert the raw data from such files
#' into a more comprehensive names list (see \code{\link{AmigaIcon}}). Use
#' \code{\link[AmigaFFH]{as.raw}} to achieve the inverse.
#'
#' @rdname rawToAmigaIcon
#' @name rawToAmigaIcon
#' @param x A vector of \code{raw} data that needs to be converted into an S3
#' \code{\link{AmigaIcon}} class object.
#' @param palette Provide a palette (\code{vector} of colours) for the icon bitmap image.
#' When set to \code{NULL} (default) the standard Amiga Workbench palette will be used.
#' @return Returns an \code{\link{AmigaIcon}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' ## generate a simple AmigaIcon object:
#' icon <- simpleAmigaIcon()
#' 
#' ## convert it into raw data:
#' icon.raw <- as.raw(icon)
#' 
#' ## convert the raw data back into an icon:
#' icon.restored <- rawToAmigaIcon(icon.raw)
#' }
#' @family AmigaIcon.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaIcon <- function(x, palette = NULL) {
  if (!all(is.raw(x))) stop("x should be a vector of raw data.")
  if (!is.null(palette) && (!all(.is.colour(palette)) || length(palette) < 4))
    stop ("The palette should consist of at least 4 colours.")
  if (any(duplicated(palette))) warning("Cannot convert this icon correctly back to raw data due to duplicated colours in the palette.")
  result <- with(.icon.data.head, .read.amigaData(x, byte, signed, par.names))
  if (result$ic_Magic != 0xe310) stop("This is not Amiga icon data")
  result$ic_Type <- .match.factor(result,
                                  "ic_Type", 1:8,
                                  c("WBDISK", "WBDRAWER", "WBTOOL", "WBPROJECT", "WBGARBAGE", "WBDEVICE", "WBKICK", "WBAPPICON"))
  result$ic_Gadget <- with(.icon.gadget.data, .read.amigaData(result$ic_Gadget, byte, signed, par.names))
  result$ic_Gadget$ga_Flags <- rev(as.logical(.rawToBitmap(result$ic_Gadget$ga_Flags, T, F)))
  names(result$ic_Gadget$ga_Flags) <- c("BACKFILL_MODE", "TWO_IMAGE_MODE", "IMAGE_MODE", paste0("UNDEFINED", 1:13))
  result$ic_Gadget$ga_UserData <- .match.factor(result$ic_Gadget,
                                                "ga_UserData", 0:1,
                                                c("OS1.x", "OS2.x"))
  
  
  ## get remainder of x:
  x <- x[-1:-sum(abs(.icon.data.head$byte))]
  result$drawer <- list(NewWindow = list())
  
  if (result$ic_DrawerData != 0) {
    result$drawer <- with(.icon.drawer.data, .read.amigaData(x, byte, signed, par.names))
    result$drawer$NewWindow <- with(.icon.new.window.data, .read.amigaData(result$drawer$NewWindow, byte, signed, par.names))
    result$drawer$NewWindow$nw_IDCMPFlags <- as.logical(.rawToBitmap(
      result$drawer$NewWindow$nw_IDCMPFlags, invert.bytes = T, T
    ))
    result$drawer$NewWindow$nw_Flags <- as.logical(.rawToBitmap(
      result$drawer$NewWindow$nw_Flags, invert.bytes = T, T
    ))
    x <- x[-1:-sum(abs(.icon.drawer.data$byte))]
  }
  
  if (is.null(palette)) {
    palette <- AmigaFFH::amiga_palettes[["wb.os2"]]
    if (result$ic_Gadget$ga_UserData == "OS1.x") {
      palette <- AmigaFFH::amiga_palettes[["wb.os1"]]
    }
  }
  
  .get.icon.image <- function(y, p = palette) {
    img <- with(.icon.image.data, .read.amigaData(y, byte, signed, par.names))
    y <- y[-1:-sum(abs(.icon.image.data$byte))]
    w <- 16*ceiling(img$im_Width/16)
    h <- img$im_Height
    bm <- bitmapToRaster(y[1:((w*h*img$im_Depth)/8)],
                         img$im_Width,
                         h, img$im_Depth,
                         p[1:(2^img$im_Depth)],
                         interleaved = F)
    attributes(bm) <- c(list(palette     = p,
                             bitmap.size = sum(abs(.icon.image.data$byte)) + w*h*img$im_Depth/8
    ),
    attributes(bm))
    img$im_Bitmap <- bm
    return(img)
  }
  
  result$firstImage <- .get.icon.image(x)
  x <- x[-1:-attributes(result$firstImage$im_Bitmap)[["bitmap.size"]]]
  
  result$secondImage <- list()
  
  if (result$ic_Gadget$ga_SelectRender != 0) {
    result$secondImage <- .get.icon.image(x)
    x <- x[-1:-attributes(result$secondImage$im_Bitmap)[["bitmap.size"]]]
  }
  
  result$defaultTool <- ""
  if (result$ic_DefaultTool != 0) {
    t.len <- .rawToAmigaInt(x[1:4], 32, F)
    result$defaultTool <- .rawToCharNull(x[5:(4 + t.len)])
    x <- x[-1:-(4 + t.len)]
  }
  
  result$toolTypes <- ""
  if (result$ic_ToolTypes != 0) {
    entries <- .rawToAmigaInt(x[1:4], 32, F)
    result$toolTypes <- NULL
    x <- x[-1:-4]
    for (i in 1:entries) {
      t.len <- .rawToAmigaInt(x[1:4], 32, F)
      result$toolTypes <- c(result$toolTypes,
                            ProTrackR::rawToCharNull(x[5:(4 + t.len)]))
      x <- x[-1:-(4 + t.len)]
    }
  }  
  
  result$toolWindow = ""
  if (result$ic_ToolWindow != 0) {
    t.len <- .rawToAmigaInt(x[1:4], 32, F)
    result$toolWindow <- .rawToCharNull(x[5:(4 + t.len)])
    x <- x[-1:-(4 + t.len)]
  }  
  
  result$dd_Flags     <- factor(NULL, c("DDFLAGS_SHOWDEFAULT", "DDFLAGS_SHOWICONS", "DDFLAGS_SHOWALL"))
  result$dd_ViewModes <- factor(NULL, c("DDVM_BYDEFAULT", "DDVM_BYICON", "DDVM_BYNAME",
                                        "DDVM_BYDATE", "DDVM_BYSIZE", "DDVM_BYTYPE"))
  if (result$ic_DrawerData != 0 && result$ic_Gadget$ga_UserData == "OS2.x") {
    result$dd_Flags <- .rawToAmigaInt(x[1:4], 32, F)
    result$dd_Flags <- .match.factor(result, "dd_Flags", 0:2,
                                     c("DDFLAGS_SHOWDEFAULT", "DDFLAGS_SHOWICONS", "DDFLAGS_SHOWALL"))
    result$dd_ViewModes <- .rawToAmigaInt(x[5:6], 16, F)
    result$dd_ViewModes <- .match.factor(result, "dd_ViewModes", 0:5,
                                         c("DDVM_BYDEFAULT", "DDVM_BYICON", "DDVM_BYNAME",
                                           "DDVM_BYDATE", "DDVM_BYSIZE", "DDVM_BYTYPE"))
    x <- x[-1:-6]
  }
  class(result) <- "AmigaIcon"
  return(result)
}

#' @rdname plot
#' @name plot
#' @export
plot.AmigaIcon <- function(x, y, asp = 2, ...) {
  if (missing(y)) y <- 1
  ## Note that the aspect ratio is not set and is probably a bit off...
  plot(as.raster(x, selected = y), asp = asp, ...)
}

#' @rdname as.raster
#' @name as.raster
#' @export
as.raster.AmigaIcon <- function(x, selected = F, ...) {
  y <- selected[[1]]
  if (is.logical(selected)) y <- as.numeric(selected[[1]]) + 1
  img <- x$firstImage$im_Bitmap
  if (x$ic_Gadget$ga_Flags[["TWO_IMAGE_MODE"]] && y == 2) img <- x$secondImage$im_Bitmap
  pal <- attributes(img)[["palette"]]
  img <- apply(img, 2, function(y) match(y, pal))
  if (x$ic_Gadget$ga_Flags[["BACKFILL_MODE"]]) {
    ## set all pixels at the edge that are equal to the background colour to NA
    img[1,][img[1,] == 1] <- NA
    img[,1][img[,1] == 1] <- NA
    img[nrow(img),][img[nrow(img),] == 1] <- NA
    img[,ncol(img)][img[,ncol(img)] == 1] <- NA
    ## Then flood fill the NA values to the remainder of the icon
    center.sel    <- as.matrix(expand.grid(2:(nrow(img) - 1), 2:(ncol(img) - 1)))
    center.left   <- center.sel
    center.right  <- center.sel
    center.top    <- center.sel
    center.bottom <- center.sel
    center.left[,1]   <- center.left[,1] - 1
    center.right[,1]  <- center.right[,1] + 1
    center.top[,2]    <- center.top[,2] - 1
    center.bottom[,2] <- center.bottom[,2] + 1
    img.center <- img[center.sel]
    while (T) {
      img.center[img.center == 1 & is.na(img[center.left])]   <- NA
      img.center[img.center == 1 & is.na(img[center.right])]  <- NA
      img.center[img.center == 1 & is.na(img[center.top])]    <- NA
      img.center[img.center == 1 & is.na(img[center.bottom])] <- NA
      if (all(is.na(img.center) == is.na(img[center.sel]))) break
      img[center.sel] <- img.center
    }
  }
  
  if (!x$ic_Gadget$ga_Flags[["TWO_IMAGE_MODE"]] && !x$ic_Gadget$ga_Flags[["BACKFILL_MODE"]] && y == 2) img <- 1 + 2^x$firstImage$im_Depth - img
  img <- apply(img, 2, function(y) pal[y])
  img <- grDevices::as.raster(img)
  attributes(img)$palette <- pal
  return(img)
}

#' @export
print.AmigaIcon <- function(x, ...) {
  print(sprintf("A %s type Amiga Icon with %s in %s mode.",
                substring(tolower(x$ic_Type), 3),
                ifelse(x$ic_Gadget$ga_Flags["TWO_IMAGE_MODE"], "two images", "one image"),
                ifelse(x$ic_Gadget$ga_Flags["BACKFILL_MODE"], "back fill", "complement")),
        ...)
}

#' @rdname as.raw
#' @name as.raw
#' @export
as.raw.AmigaIcon <- function(x, ...) {
  x$ic_Gadget$ga_Flags <- .bitmapToRaw(rev(x$ic_Gadget$ga_Flags), T, F)
  x$ic_Gadget$ga_UserData <- .match.factor.inv(x$ic_Gadget,
                                               "ga_UserData", 0:1,
                                               c("OS1.x", "OS2.x"))
  sec.img <- x$ic_Gadget$ga_SelectRender != 0
  x$ic_Gadget <- .write.amigaData(x$ic_Gadget,
                                  .icon.gadget.data$byte,
                                  .icon.gadget.data$signed,
                                  .icon.gadget.data$par.names)
  x$ic_Type <- .match.factor.inv(x,
                                 "ic_Type", 1:8,
                                 c("WBDISK", "WBDRAWER", "WBTOOL", "WBPROJECT", "WBGARBAGE", "WBDEVICE", "WBKICK", "WBAPPICON"))
  if (x$ic_DrawerData != 0) {
    x$drawer$NewWindow$nw_IDCMPFlags <- .bitmapToRaw(x$drawer$NewWindow$nw_IDCMPFlags, F, T)
    x$drawer$NewWindow$nw_Flags      <- .bitmapToRaw(x$drawer$NewWindow$nw_Flags, F, T)
    x$drawer$NewWindow <- with(.icon.new.window.data, .write.amigaData(x$drawer$NewWindow, byte, signed, par.names))
    x$drawer           <- with(.icon.drawer.data,     .write.amigaData(x$drawer, byte, signed, par.names))
  } else {
    x$drawer <- NULL
  }
  iconImgToRaw <- function(y) {
    pal <- attributes(y$im_Bitmap)[["palette"]][1:(2^y$im_Depth)]
    list(
      bmhead = .write.amigaData(y,
                                .icon.image.data$byte,
                                .icon.image.data$signed,
                                .icon.image.data$par.names),
      bm     = .bitmapToRaw(rasterToBitmap(
        y$im_Bitmap,
        depth = y$im_Depth,
        interleaved = F,
        indexing = function(x, length.out) index.colours(x, length.out,
                                                         palette = pal)),
        T, F)
    )
  }
  x$firstImage <- iconImgToRaw(x$firstImage)
  if (sec.img) {
    x$secondImage <- iconImgToRaw(x$secondImage)
  } else {
    x$secondImage <- NULL
  }
  if (x$defaultTool != "") {
    x$defaultTool <- c(
      .amigaIntToRaw(nchar(x$defaultTool) + 1, 32, F),
      charToRaw(x$defaultTool),
      raw(1))
  } else {
    x$defaultTool <- NULL
  }
  if (length(x$toolTypes) == 1 && x$toolTypes == "") {
    x$toolTypes <- NULL
  } else {
    x$toolTypes <- c(
      .amigaIntToRaw(length(x$toolTypes), 32, F),
      unlist(lapply(x$toolTypes, function(y){
        nc <- nchar(y)
        if (nc == 0) yc <- raw(0) else yc <- charToRaw(y)
        c(.amigaIntToRaw(nc, 32, F),
          yc,
          raw(1))
      }))
    )
  }
  if (x$toolWindow != "") {
    x$toolWindow <- c(
      .amigaIntToRaw(nchar(x$toolWindow) + 1, 32, F),
      charToRaw(x$toolWindow),
      raw(1))
  } else {
    x$toolWindow <- NULL
  }
  x[.icon.data.head$par.names] <- lapply(1:nrow(.icon.data.head), function(y) {
    .write.amigaData(x[.icon.data.head$par.names[y]],
                     .icon.data.head$byte[y],
                     .icon.data.head$signed[y],
                     .icon.data.head$par.names[y])
  })
  if (length(x$dd_Flags) == 1) {
    x$dd_Flags <- .match.factor.inv(x, "dd_Flags", 0:2,
                                    c("DDFLAGS_SHOWDEFAULT", "DDFLAGS_SHOWICONS", "DDFLAGS_SHOWALL"))
    x$dd_Flags <- .amigaIntToRaw(x$dd_Flags, 32, F)
  } else {
    x$dd_Flags <- NULL
  }
  if (length(x$dd_ViewModes) == 1) {
    x$dd_ViewModes <- .match.factor.inv(x, "dd_ViewModes", 0:5,
                                        c("DDVM_BYDEFAULT", "DDVM_BYICON", "DDVM_BYNAME",
                                          "DDVM_BYDATE", "DDVM_BYSIZE", "DDVM_BYTYPE"))
    x$dd_ViewModes <- .amigaIntToRaw(x$dd_ViewModes, 16, F)
  } else {
    x$dd_ViewModes <- NULL
  }
  x <- unlist(x)
  names(x) <- NULL
  return(x)
}

#' Write an Amiga Workbench icon (info) file
#'
#' Graphical representation of files and directories (icons) are stored as
#' separate files (with the .info extension) on the Amiga. This function writes
#' \code{\link{AmigaIcon}} class objects to such files.
#'
#' The \code{\link{AmigaIcon}} S3 object provides a comprihensive format
#' for Amiga icons, which are used as a graphical representation of files
#' and directories on the Amiga. The \code{\link{AmigaIcon}} is a named
#' list containing all information of an icon. Use this function to
#' write this object to a file which can be used on the Commodore Amiga
#' or emulator.
#'
#' @rdname write.AmigaIcon
#' @name write.AmigaIcon
#' @param x An \code{\link{AmigaIcon}} class object.
#' @param file A \code{character} string representing the file name to which the
#' icon data should be written.
#' @return Returns \code{NULL} or an \code{integer} status passed on by the
#' \code{\link{close}} function, that is used to close the file connection.
#' It is returned invisibly.
#' @examples
#' \dontrun{
#' ## create a simple AmigaIcon:
#' icon <- simpleAmigaIcon()
#' 
#' ## write the icon to the temp dir:
#' write.AmigaIcon(icon, file.path(tempdir(), "icon.info"))
#' }
#' @family AmigaIcon.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.AmigaIcon <- function(x, file) {
  if (class(x) != "AmigaIcon") stop("x should be of S3 class AmigaIcon.")
  con <- file(file, "wb")
  writeBin(as.raw(x), con)
  close(con)
}

#' Read an Amiga Workbench icon (info) file
#'
#' Graphical representation of files and directories (icons) are stored as
#' separate files (with the .info extension) on the Amiga. This function reads such files
#' and imports them as \code{\link{AmigaIcon}} class objects.
#'
#' The \code{\link{AmigaIcon}} S3 object provides a comprihensive format
#' for Amiga icons, which are used as a graphical representation of files
#' and directories on the Amiga. The \code{\link{AmigaIcon}} is a named
#' list containing all information of an icon. Use this function to
#' read an Amiga icon (with the .info extension) from a file and convert
#' it into an \code{\link{AmigaIcon}} object.
#'
#' @rdname read.AmigaIcon
#' @name read.AmigaIcon
#' @param file A \code{character} string representing the file name from which the
#' icon data should be read.
#' @param ... Arguments passed on to \code{\link{rawToAmigaIcon}}.
#' @return Returns an \code{\link{AmigaIcon}} class object as read from the \code{file}.
#' @examples
#' \dontrun{
#' ## create a simple AmigaIcon:
#' icon <- simpleAmigaIcon()
#' 
#' ## write the icon to the temp dir:
#' write.AmigaIcon(icon, file.path(tempdir(), "icon.info"))
#' 
#' ## read the same file:
#' icon2 <- read.AmigaIcon(file.path(tempdir(), "icon.info"))
#' }
#' @family AmigaIcon.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaIcon <- function(file, ...) {
  fz <- file.size(file)
  con <- file(file, "rb")
  dat <- readBin(con, "raw", fz)
  close(con)
  rawToAmigaIcon(dat, ...)
}

#' @export
`$<-.AmigaIcon` <- function(x, i, value) {
  x[[i]] <- value
  x
}

#' @export
`[[<-.AmigaIcon` <- function(x, i, value) {
  cl <- class(x)
  class(x) <- NULL
  x[[i]] <- value
  class(x) <- cl
  ## XXX update this replacement function and remove warning
  warning(paste0("Replacement operator for AmigaIcon objects ",
          "will be modified in future versions of this package. ",
          "Note that not all replacement operations may be ",
          "allowed in future versions of this package."))
  x
}