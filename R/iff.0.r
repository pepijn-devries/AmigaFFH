.validity.IFFChunk <- function(object) {
  if (length(object@chunk.type) !=1) stop("Slot 'chunk.type' should be of length 1")
  if (nchar(object@chunk.type) != 4) stop("Slot 'chunk.type' should consist of 4 characters.")
  data.types <- unlist(lapply(object@chunk.data, class))
  if (length(object@chunk.data) == 0) stop("Chunk should have at least one element in it's data list.")
  if (!(all(data.types == "IFFChunk") || all(data.types == "raw"))) stop("Chunk data should be a list of IFFChunk objects or a list with a single element with raw data.")
  if (all(data.types == "raw") && length(object@chunk.data) > 1) stop ("Chunk data can only hold one element of raw data.")
  ## note that the chunk.size can be odd, but in that case a padding byte should be read from
  ## an iff file...
  if (all(data.types == "raw") && length(object@chunk.data[[1]]) > (2^32 - 1)) stop("Chunk data cannot be larger than 4,294,967,295 bytes")
  if (object@chunk.type %in% c("FORM", "LIST", "CAT ", "PROP") && data.types == "raw") stop("IFF containers should contain IFF chunks, not raw data.")
  ## TODO note that the validity of all the IFFChunk objects in the data list are not checked
  ## TODO Doing that will make the object fool-proof, but also a lot slower...
  return(TRUE)
}

#' A class structure to represent IFF files
#'
#' An S4 class structure to represent data stored in the Interchange File
#' Format (IFF).
#'
#' The Interchange File Format (IFF) was introduced in 1985 by Electronic Arts.
#' This format stores files in standardised modular objects, called `chunks'.
#' At the start of each chunk it is specified what type of data can be expected
#' and what the size of this data is. This was a very forward thinking way of
#' storing data, similar structures are still used in modern file formats (such
#' as PNG images and XML files).
#' 
#' Although the IFF format is still in use, and new standardised chunk types can
#' still be registered, this package will focus on the older chunk types that
#' were primarily used on the Commodore Amiga (OS <= 3.0). IFF files could
#' contain any kind of information. It could contain bitmap images, but also
#' audio clips or (formatted) texts.
#' 
#' The `IFFChunk` class is designed such that it theoretically can hold
#' any type of IFF data. This package will mostly focus on the early IFF file types
#' (i.e., IFF chunks as originally registered by Electronic Arts). IFF files are
#' read by this package in a none lossy way ([read.iff()]), such that all
#' information is preserved (even if it is of an unknown type, as long as the chunk
#' identifier is 4 characters long).
#' 
#' This means that the object needs to be interpreted in order to make sense out of
#' it ([interpretIFFChunk()]). This interpretation returns simplified
#' interpretations of class `IFF.ANY` when it is supported (see
#' [IFFChunk-method()] for supported chunk types). Note that in the
#' interpretation process (meta-)information may get lost. converting
#' `IFF.ANY` objects back into [IFFChunk()] objects (if possible)
#' could therefore result in an object that is different from then one stored in the
#' original file and could even destroy the correct interpretation of IFF objects.
#' IFF files should thus be handled with care.
#'
#' @slot chunk.type A four `character` long code reflecting the type of
#' information represented by this chunk.
#' @slot chunk.data A `list` that holds either one or more valid
#' `IFFChunk`s or a single `vector` of `raw` data. This data
#' can only be interpreted in context of the specified type or in some cases
#' information from other `IFFChunk`s.
#' @references <https://wiki.amigaos.net/wiki/IFF_Standard>
#' @references <https://wiki.amigaos.net/wiki/IFF_FORM_and_Chunk_Registry>
#' @references <https://en.wikipedia.org/wiki/Interchange_File_Format>
#' @name IFFChunk-class
#' @rdname IFFChunk-class
#' @examples
#' ## load an IFF file
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## 'example.iff' is of class IFFChunk:
#' class(example.iff)
#' 
#' ## let's plot it:
#' plot(example.iff)
#' 
#' ## The default constructor will create an empty FORM:
#' new("IFFChunk")
#' 
#' ## The constructor can also be used to create simple chunks:
#' new("IFFChunk",
#'     chunk.type = "TEXT",
#'     chunk.data = list(charToRaw("A simple chunk")))
#' @family iff.operations
#' @exportClass IFFChunk
#' @author Pepijn de Vries
setClass("IFFChunk",
         representation(
           chunk.type = "character",
           chunk.data = "list" ## of either a single element of type raw, or a list of iffchunks
         ),
         prototype(
           chunk.type = "FORM",
           chunk.data = list(raw(0))
         ),
         validity = .validity.IFFChunk)

#' Read Interchange File Format (IFF)
#'
#' Read the Interchange File Format (IFF) as an [IFFChunk()] object.
#'
#' Information is stored as `chunks' in IFF files (see [IFFChunk()]).
#' Each chunk should at least contain a label of the type of chunk and the data
#' for that chunk. This function reads all chunks from a valid IFF file, including
#' all nested chunks and stores them in an [IFFChunk()] object. IFF
#' files can hold any kind of data (e.g. images or audio), this read function
#' does not interpret the file. Use [interpretIFFChunk()] for that
#' purpose.
#'
#' @rdname read.iff
#' @name read.iff
#' @param file A filename of an IFF file to be read, or a connection from which
#' binary data can be read.
#' @returns Returns a [IFFChunk()] object read from the specified file.
#' @examples
#' ## let's read a bitmap image stored in IFF as provided with this package:
#' filename <- system.file("ilbm8lores.iff", package = "AmigaFFH")
#' example.iff <- read.iff(filename)
#' 
#' ## And plot it:
#' plot(example.iff)
#' @family io.operations
#' @family iff.operations
#' @author Pepijn de Vries
#' @export
read.iff <- function(file) {
  result <- .read.generic(file)
  result <- rawToIFFChunk(result)
  if (result@chunk.type != "FORM") stop("FORM is currently the only supported IFF container. LIST, CAT and others are not...")
  return (result)
}

#' Convert AmigaFFH objects into raw data
#'
#' Convert AmigaFFH objects into raw data, as they would be stored in the Commodore
#' Amiga's memory or files.
#'
#' Objects originating from this package can in some cases be converted into
#' raw data, as they would be stored on an original Amiga. See the usage section
#' for the currently supported objects.
#' 
#' Not all information from `x` may be included in the `raw`
#' data that is returned, so handle with care.
#' 
#' As this package grows additional objects can be converted with this method.
#'
#' @docType methods
#' @rdname as.raw
#' @name as.raw
#' @aliases as.raw,IFFChunk-method
#' @param x An AmigaFFH object that needs to be converted into raw data.
#' See usage section for all supported objects.
#' @returns Returns a `vector` of `raw` data based on `x`.
#' @examples
#' ## read an IFF file as an IFFChunk object:
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## This will recreate the exact raw data as it was read from the file:
#' example.raw <- as.raw(example.iff)
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
setMethod("as.raw", "IFFChunk", function(x) {
  get.data <- function(y, parent.is.container = FALSE) {
    result <- charToRaw(y@chunk.type)
    if (inherits(y@chunk.data[[1]], "raw")) {
      ## only store chunk size if parent is not a container (i.e., FORM, LIST or CAT)
      if (!parent.is.container) result <- c(result, .amigaIntToRaw(length(y@chunk.data[[1]]), 32, FALSE))
      result <- c(result, y@chunk.data[[1]])
      ## Chunks should always be WORD aligned (pad with zeros if
      ## it is not):
      if ((length(y@chunk.data[[1]]) %% 2) != 0) result <- c(result, raw(1))
      return(result)
    } else if (inherits(y@chunk.data[[1]], "IFFChunk")) {
      container <- y@chunk.type %in% c("FORM", "LIST", "CAT ", "PROP")
      dat <- unlist(lapply(y@chunk.data, get.data,
                           parent.is.container = container))
      ## only store chunk size if parent is not a container (i.e., FORM, LIST, CAT or PROP)
      if (!parent.is.container) result <- c(result, .amigaIntToRaw(length(dat), 32, FALSE))
      result <- c(result, dat)
    } else {
      stop("IFFChunk contains invalid data")
    }
  }
  return(get.data(x, FALSE))
})

#' @rdname as.raw
#' @method as.raw IFF.ANY
#' @param ... Arguments passed on to [IFFChunk-method()] when `x` is
#' of class `IFF.ANY`.
#' @export
as.raw.IFF.ANY <- function(x, ...) {
  as.raw(IFFChunk(x, ...))
}

#' Write Interchange File Format (IFF)
#'
#' Write an [IFFChunk()] object conform the Interchange File Format (IFF).
#'
#' Writes an [IFFChunk()] object (including all nested chunks) to the
#' specified file. Only the structure of the object needs to be valid, however,
#' a correctly structured file does not necessarily result in an interpretable file
#' (see examples).
#'
#' @rdname write.iff
#' @name write.iff
#' @param x An [IFFChunk()] object that needs to be written to a file.
#' @param file A filename for the IFF file to which the [IFFChunk()] needs
#' to be saved, or a connection to which the data should be written.
#' @returns Returns either `NULL` or an `integer` status invisibly as passed
#' by the [`close()`][base::connections] statement used to close the file connection.
#' 
#' @references <https://en.wikipedia.org/wiki/Interchange_File_Format>
#' @examples
#' ## read an IFF file as an IFFChunk object:
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## This will write the IFF file (in this case a bitmap image)
#' ## to the temp directory:
#' write.iff(example.iff, file.path(tempdir(), "image.iff"))
#' 
#' @family io.operations
#' @family iff.operations
#' @author Pepijn de Vries
#' @export
write.iff <- function(x, file) {
  if (!inherits(x, "IFFChunk")) stop("x should be of class IFFChunk.")
  .write.generic(x, file)
}

setGeneric("getIFFChunk", function(x, chunk.path, chunk.number) standardGeneric("getIFFChunk"))

#' Get a specific IFFChunk nested inside other IFFChunks
#'
#' [IFFChunk()]s can be nested in a tree-like structure. Use this method
#' to get a specific chunk with a specific label.
#'
#' `IFFChunk` objects have 4 `character` identifiers, indicating what type
#' of chunk you are dealing with. These chunks can be nested inside of each other.
#' Use this method to extract specific chunks by referring to there respective
#' identifiers. The identifiers are shown when calling `print` on an
#' [IFFChunk()]. If a specified path doesn't exist, this method throws a
#' `subscript out of range' error.
#'
#' @docType methods
#' @rdname getIFFChunk
#' @name getIFFChunk
#' @aliases getIFFChunk,IFFChunk,character,integer-method
#' @param x An [IFFChunk()] object from which the nested
#' [IFFChunk()] should be extracted an returned.
#' @param chunk.path A `vector` of 4 `character` long strings
#' of IFF chunk labels, specifying the path of the target IFF chunk.
#' For example: `c("ILBM", "BODY")` means, get the "BODY" chunk
#' from inside the "ILBM" chunk.
#' @param chunk.number A `vector` of the same length as `chunk.path`,
#' with `integer` index numbers. Sometimes a chunk can contain a list of
#' chunks with the same label. With this argument you can specify which element
#' should be returned. By default (when missing), the first element is always
#' returned.
#' @returns Returns an [IFFChunk()] object nested inside `x` at the
#' specified path. Or in case of the replace method the original chunk `x` is
#' returned with the target chunk replaced by `value`.
#' @examples
#' ## load an IFF file
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## Get the BMHD (bitmap header) from the ILBM (interleaved bitmap) chunk:
#' bmhd <- getIFFChunk(example.iff, c("ILBM", "BMHD"))
#' 
#' ## This is essentially doing the same thing, but we now explicitly
#' ## tell the method to get the first element for each specified label:
#' bmhd <- getIFFChunk(example.iff, c("ILBM", "BMHD"), c(1L, 1L))
#' 
#' ## Let's modify the bitmap header and replace it in the parent IFF chunk.
#' bmhd.itpt <- interpretIFFChunk(bmhd)
#' 
#' ## Let's disable the masking, the bitmap will no longer be transparent:
#' bmhd.itpt$Masking <- "mskNone"
#' bmhd <- IFFChunk(bmhd.itpt)
#' 
#' ## Now replace the header from the original iff with the modified header:
#' getIFFChunk(example.iff, c("ILBM", "BMHD")) <- bmhd
#' @family iff.operations
#' @author Pepijn de Vries
#' @export
setMethod("getIFFChunk", c("IFFChunk", "character", "integer"), function(x, chunk.path, chunk.number) {
  chunk <- x
  id <- unlist(lapply(1:length(chunk.path), function(y){
    id <- which(unlist(lapply(chunk@chunk.data, function(z) z@chunk.type)) == chunk.path[[y]])
    chunk <<- chunk@chunk.data[[id[[chunk.number[[y]]]]]]
    return(id)
  }))
  return(chunk)
})

#' @rdname getIFFChunk
#' @name getIFFChunk
#' @aliases getIFFChunk,IFFChunk,character,missing-method
#' @export
setMethod("getIFFChunk", c("IFFChunk", "character", "missing"), function(x, chunk.path, chunk.number) {
  getIFFChunk(x, chunk.path, rep(1L, length(chunk.path)))
})

setGeneric("getIFFChunk<-", function(x, chunk.path, chunk.number, value) standardGeneric("getIFFChunk<-"))

#' @rdname getIFFChunk
#' @name getIFFChunk<-
#' @aliases getIFFChunk<-,IFFChunk,character,missing,IFFChunk-method
#' @param value An [IFFChunk()] with which the target chunk should be
#' replaced. Make sure that `value` is of the same `chunk.type` as the last
#' chunk specified in the `chunk.path`.
#' @export
setReplaceMethod("getIFFChunk", c("IFFChunk", "character", "missing", "IFFChunk"), function(x, chunk.path, chunk.number = NULL, value) {
  getIFFChunk(x, chunk.path, rep(1L, length(chunk.path))) <- value
  return(x)
})

#' @rdname getIFFChunk
#' @name getIFFChunk<-
#' @aliases getIFFChunk<-,IFFChunk,character,integer,IFFChunk-method
#' @export
setReplaceMethod("getIFFChunk", c("IFFChunk", "character", "integer", "IFFChunk"), function(x, chunk.path, chunk.number = NULL, value) {
  if (value@chunk.type != chunk.path[[length(chunk.path)]]) stop("'value' should be of the same IFF type as the last type in the 'chunk.path'")
  iff.list <- as.list(x)
  val.list <- as.list(value)
  eval(parse(text = paste0("iff.list",
                           paste0(sprintf("[[\"%s_%i\"]]", chunk.path, chunk.number), collapse = ""),
                           "<-val.list")))
  
  ## code below doesn't seem to work correctly:
  list.to.chunk <- function(x, nam) {
    if (typeof(x) == "raw") {
      return(new("IFFChunk",
                 chunk.type = substr(nam, 1, 4),
                 chunk.data = list(x)))
    } else if (typeof(x) == "list") {
      result <- lapply(seq_along(x), function(y) {
        list.to.chunk(x[[y]], names(x)[y])
      })
      return(new("IFFChunk",
                 chunk.type = substr(nam, 1, 4),
                 chunk.data = result))
    }
  }
  return(list.to.chunk(iff.list, x@chunk.type))
})

setGeneric("interpretIFFChunk", function(x, ...) standardGeneric("interpretIFFChunk"))

#' Interpret an IFFChunk object
#'
#' [IFFChunk()]s represent the structure of the Interchange File Format well,
#' but the iformation is stored as `raw` data. This method tries to interpret and
#' translate the information into a more comprehensive format.
#'
#' Interchange File Format chunks can hold any kind of information (images, audio,
#' (formatted) text, etc.). This method will try to convert this information into
#' something useful. Information may get lost in the translation, so be careful when
#' converting back to an [IFFChunk-class()] object using
#' [IFFChunk-method()].
#' 
#' An error is thrown when the [IFFChunk()] object is currently not
#' interpretable by this package. See [IFFChunk-method()] for an overview
#' of currently supported IFF chunks. This list may increase while this package
#' matures.
#'
#' @docType methods
#' @rdname interpretIFFChunk
#' @name interpretIFFChunk
#' @aliases interpretIFFChunk,IFFChunk-method
#' @param x An [IFFChunk()] object which needs to be interpreted.
#' @param ... Currently ignored.
#' @returns If `x` is interpretable by this package an S3 class object of
#' `IFF.ANY` is returned. The content of the returned object will depend
#' on the type of [IFFChunk()] provided for `x`. The result can
#' for instance be a `raster` image ([grDevices::as.raster()]),
#' a list of audio [tuneR::Wave()]s, a `character` string or a named
#' `list`.
#' @examples
#' ## load an IFF file
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## in this case, the file is a FORM container with a bitmap image, and a
#' ## list with a raster object is returned when interpreted:
#' example.itpt <- interpretIFFChunk(example.iff)
#' class(example.itpt)
#' typeof(example.itpt)
#' class(example.itpt[[1]])
#' 
#' ## Let's extraxt the bitmap header from the main chunk:
#' bmhd <- getIFFChunk(example.iff, c("ILBM", "BMHD"))
#' 
#' ## When interpreted, a named list is returned with (meta-)information
#' ## on the bitmap image:
#' bmhd.itpt <- interpretIFFChunk(bmhd)
#' class(bmhd.itpt)
#' typeof(bmhd.itpt)
#' print(bmhd.itpt)
#' @family iff.operations
#' @author Pepijn de Vries
#' @export
setMethod("interpretIFFChunk", "IFFChunk", function(x, ...) {
  type <- x@chunk.type
  if (!inherits(x@chunk.data[[1]], "raw")) sub.types <- unlist(lapply(x@chunk.data, function(x) x@chunk.type))
  dat  <- x@chunk.data[[1]]
  if (type == "FORM") {
    ## FORM can hold multiple objects, return the interpreation of these objects
    ## as a list:
    result <- lapply(x@chunk.data, interpretIFFChunk, ...)
    class(result) <- c("IFF.FORM", "IFF.ANY")
    return(result)
  } else if (type == "ILBM") {
    ## ILBM = interleaved bitmap; return as raster
    ## If subtype contains ANHD and not a BMHD, it should also contain a DLTA chunk rather than a BODY
    ## note that DLTA chunks are actualy not interleaved, in contrast to what the parent ILBM chunk
    ## suggests.
    if ("ANHD" %in% sub.types && !("BODY" %in% sub.types)) {
      anhd <- interpretIFFChunk(x@chunk.data[sub.types %in% "ANHD"][[1]], ...)
      dlta <- x@chunk.data[sub.types %in% "DLTA"][[1]]@chunk.data[[1]]
      wbm <- with(anhd, w + (-w %% 16))
      if (anhd$operation == "ByteVerticalCompression") {
        result <- .byteVerticalDecompression(dlta, wbm, anhd$h, anhd$interleave, anhd$flags[2], list(...)$hidden)
        cl <- class(result)
        result <- result[,0:anhd$w]
        class(result) <- cl
        return(result)
      } else {
        stop("Sorry this animation format is not (yet) supported by this package.")
      }
    }
    if ("hidden" %in% names(list(...))) {
      result <- grDevices::as.raster(x, palette = NULL)
    } else {
      result <- grDevices::as.raster(x)
    }
    class(result) <- c("IFF.ILBM", "IFF.ANY", class(result))
    return(result)
  } else if (type == "CMAP") {
    ## COLOUR MAP
    ## when all low bits are 0 assume 12 bit, otherwise asume 24 bit
    colour.depth  <- ifelse(all(ProTrackR::loNybble(x@chunk.data[[1]]) == 0), "12 bit", "24 bit")
    result        <- amigaRawToColour(dat, colour.depth, "3")
    class(result) <- c("IFF.CMAP", "IFF.ANY")
    return(result)
  } else if(type == "BMHD") {
    ## BITMAP HEADER
    result <- list(
      w                  = .rawToAmigaInt(dat[1:2],   16, FALSE),
      h                  = .rawToAmigaInt(dat[3:4],   16, FALSE),
      x                  = .rawToAmigaInt(dat[5:6],   16, TRUE),
      y                  = .rawToAmigaInt(dat[7:8],   16, TRUE),
      nPlanes            = .rawToAmigaInt(dat[9],     8, FALSE),
      Masking            = .rawToAmigaInt(dat[10],    8, FALSE),
      Compression        = .rawToAmigaInt(dat[11],    8, FALSE),
      pad                = dat[12],
      transparentColour  = .rawToAmigaInt(dat[13:14], 16, FALSE),
      xAspect            = .rawToAmigaInt(dat[15],    8, FALSE),
      yAspect            = .rawToAmigaInt(dat[16],    8, FALSE),
      pageWidth          = .rawToAmigaInt(dat[17:18], 16, TRUE),
      pageHeight         = .rawToAmigaInt(dat[19:20], 16, TRUE)
    )
    if (result$Masking > 3) result$Masking <- 4
    result$Masking <- c("mskNone", "mskHasMask", "mskHasTransparentColour", "mskLasso", "mskUnknown")[result$Masking + 1]
    if (result$Compression > 1) result$Compression <- 2
    result$Compression <- c("cmpNone", "cmpByteRun1", "cmpUnknown")[result$Compression + 1]
    class(result) <- c("IFF.BMHD", "IFF.ANY")
    return(result)
  } else if (type == "CAMG") {
    ## Amiga Viewport Mode
    result <- .amigaViewPortModes(dat)
    class(result) <- c("IFF.CAMG", "IFF.ANY")
    return(result)
  } else if (type == "CRNG") {
    ## DPaint colour range (used for colour cycling)
    result <- list(
      padding = dat[1:2],
      rate    = .rawToAmigaInt(dat[3:4], 16, FALSE)*60/(2^14), ## steps per second
      flags   = .rawToAmigaInt(dat[5:6], 16, FALSE),
      low     = .rawToAmigaInt(dat[7], 8, FALSE),
      high    = .rawToAmigaInt(dat[8], 8, FALSE)
    )
    result$flags[result$flags > 2] <- 2
    result$flags <- c("RNG_OFF", "RNG_ACTIVE", "RNG_REVERSE", "RNG_UNKNOWN")[result$flags + 1]
    class(result) <- c("IFF.CRNG", "IFF.ANY")
    return (result)
  } else if (type == "ANIM") {
    result <- list()
    dpan <- NULL
    suppressWarnings(try(dpan <- interpretIFFChunk(getIFFChunk(x, c("FORM", "ILBM", "DPAN"))), TRUE))
    interleave <- interpretIFFChunk(getIFFChunk(x, c("FORM", "ILBM", "ANHD"), c(2L, 1L, 1L)))$interleave
    interleave[interleave == 0] <- 2
    if (length(x@chunk.data) < (2 + interleave)) stop("Animation contains insufficient information for at least 2 frames.")
    for (i in 1:length(x@chunk.data)) {
      ## Assuming that each sub-form contains a single frame (hence [[1]]):
      ## if first frame contains DLTA, it should be the difference from blank screen
      ## in that case result == list(), and is handled by .byteVerticalDecompression
      result[[i]] <- interpretIFFChunk(x@chunk.data[[i]], hidden = result, ...)[[1]]
    }
    for (k in 1:interleave) {
      if (!all(result[[length(result) + k - interleave]] == result[[k]])) stop("Could not interpret animation DLTA chunks correctly")
    }
    result <- result[1:(length(result) - interleave)]

    palettes <- lapply(result, function(z) attributes(z)$palette)
    asps     <- lapply(result, function(z) attributes(z)$asp)
    modes    <- lapply(result, function(z) attributes(z)$mode)
    
    for (i in 1:length(palettes)) {
      if (is.null(palettes[[i]])) palettes[[i]] <- palettes[[i - 1]]
      if (is.null(asps[[i]])) asps[[i]] <- asps[[i - 1]]
      if (i > 1 && !is.null(modes[[1]]) && is.null(modes[[i]])) modes[[i]] <- modes[[i - 1]]

      if (!is.null(modes[[i]]) && modes[[i]] %in% c("HAM6", "HAM8")) {
        result[[i]] <- .indexToHAMraster(result[[i]],
                                         ifelse(modes[[i]] == "HAM8", 8, 6),
                                         palettes[[i]], 0) ## assume transparent colour is zero        
      } else {
        result[[i]] <- grDevices::as.raster(apply(result[[i]], 2, function(z) palettes[[i]][z + 1]))
      }
      class(result[[i]]) <- c("IFF.ILBM", "IFF.ANY", class(result[[i]]))
      attributes(result[[i]])[["asp"]] <- asps[[i]]
    }

    if (!is.null(dpan) && dpan$nframes != length(result)) warning("Number of frames does not match the number specified in the DPAN chunk.")
    class(result) <- c("IFF.ANIM", "IFF.ANY")
    return(result)
  } else if (type == "ANHD") {
    result <- list(
      operation       = .rawToAmigaInt(dat[1], 8, FALSE),
      mask            = as.logical(.rawToBitmap(dat[2], FALSE, FALSE)),
      w               = .rawToAmigaInt(dat[3:4], 16, FALSE),
      h               = .rawToAmigaInt(dat[5:6], 16, FALSE),
      x               = .rawToAmigaInt(dat[7:8], 16, TRUE),
      y               = .rawToAmigaInt(dat[9:10], 16, TRUE),
      abstime         = .rawToAmigaInt(dat[11:14], 32, TRUE),
      reltime         = .rawToAmigaInt(dat[15:18], 32, TRUE),
      interleave      = .rawToAmigaInt(dat[19], 8, FALSE),
      pad0            = dat[20],
      flags           = as.logical(.rawToBitmap(dat[21:24], FALSE, TRUE)),
      pad1            = dat[25:40]
    )
    if (result$operation > 7) result$operation <- "UnknownMode" else
      result$operation <- c("standard", "XOR", "LongDeltaMode", "ShortDeltaMode", "GeneralDeltamode", "ByteVerticalCompression", "StereoOp5", "ShortLongVerticalDeltaMode")[1 + result$operation]
    class(result) <- c("IFF.ANHD", "IFF.ANY")
    return(result)
  } else if (type == "DLTA") {
    ## without context (its parent ANIM chunk or neighbouring ANHD chunk),
    ## DLTA chunks can't be interpreted. Just return the raw data
    class(dat) <- c("IFF.DLTA", "IFF.ANY")
    return(dat)
  } else if (type == "DPAN") {
    ## DPaint Animation chunk, is only used to determine number of animation frames
    result <- list(
      version         = .rawToAmigaInt(dat[1:2], 16, FALSE),
      nframes         = .rawToAmigaInt(dat[3:4], 16, FALSE),
      flags           = as.logical(.rawToBitmap(dat[5:8], FALSE, TRUE))
    )
    return(result)
  } else if (type == "VHDR") {
    result <- as.list(.rawToAmigaInt(dat[1:(3*4)], 32, FALSE))
    result <- c(result,
                .rawToAmigaInt(dat[13:14], 16, FALSE),
                .rawToAmigaInt(dat[15:16], 8, FALSE),
                .rawToAmigaInt(dat[17:20], 32, FALSE))
    names(result) <- c("oneShotHiSamples",
                       "repeatHiSamples",
                       "samplesPerHiCycle",
                       "samplesPerSec",
                       "ctOctave",
                       "sCompression",
                       "volume")
    result$sCompression[result$sCompression > 2] <- 2
    result$sCompression <- c("sCmpNone", "sCmpFibDelta", "sCmpUnknown")[result$sCompression + 1]
    class(result) <- c("IFF.VHDR", "IFF.ANY")
    return(result)
  } else if (type == "CHAN") {
    result <- which(c(2, 4, 6) %in% .rawToAmigaInt(dat, 32, FALSE))[[1]]
    if (length(result) == 0) result <- list(channel = "UNKOWN") else
      result <- list(channel = c("LEFT", "RIGHT", "STEREO")[result])
    class(result) <- c("IFF.CHAN", "IFF.ANY")
    return(result)
  } else if (type == "8SVX") {
    body <- getIFFChunk(x, "BODY")@chunk.data[[1]]
    if ("CHAN" %in% sub.types) {
      chan <- interpretIFFChunk(getIFFChunk(x, "CHAN"), ...)
    } else {
      chan <- interpretIFFChunk(IFFChunk("CHAN"), ...)
    }
    if ("VHDR" %in% sub.types) {
      vhdr <- interpretIFFChunk(getIFFChunk(x, "VHDR"), ...)
    } else {
      warning("Voice header is missing, going to make some assumptions here...")
      vhdr <- interpretIFFChunk(IFFChunk("VHDR"), ...)
      vhdr$oneShotHiSamples <- length(body)
      if (chan == "STEREO") vhdr$oneShotHiSamples <- floor(vhdr$oneShotHiSamples/2)
    }
    if (vhdr$sCompression == "sCmpUnknown") warning("An unknown form of compression is applied to the wave. Trying to continue anyway.")
    if (vhdr$sCompression == "sCmpFibDelta") body <- deltaFibonacciDecompress(body)
    samp.offset <- 0
    wav <- lapply(1:ifelse(chan == "STEREO", 2, 1), function(z){
      lapply(1:vhdr$ctOctave, function(y) {
        l <- ((vhdr$oneShotHiSamples + vhdr$repeatHiSamples)*(2^(y - 1)))
        result <- body[samp.offset + (1:l)]
        samp.offset <<- samp.offset + l
        .rawToAmigaInt(result, 8, TRUE) + 128
      })
    })
    ## wav will be a list of one or more Wave objects.
    ## We have let wav be a list. We can't downgrade the S4 class
    ## to an S3 IFF.ANY class. Ortherwise the Wave-methods are
    ## no longer applicable.
    wav <- lapply(1:vhdr$ctOctave, function(y) {
      ## Wave objects are always specified for left channel when audio is mono
      ## eventhough amiga iff is specified for right channel only.
      right <- numeric(0)
      if (chan == "STEREO") right <- wav[[2]][[y]]
      Wave(left = wav[[1]][[y]],
           right = right,
           bit = 8,
           pcm = TRUE,
           samp.rate = vhdr$samplesPerSec)
    })
    class(wav) <- c("IFF.8SVX", "IFF.ANY", class(wav))
    return(wav)
  } else if (type == "BODY") {
    ## without context, the body cannot be converted to anything usefull
    ## return as unprocessed raw data
    class(dat) <- c("IFF.BODY", "IFF.ANY")
    return(dat)
  } else if (type %in% c("ANNO", "AUTH", "CHRS", "NAME", "TEXT", "(c) ")) {
    ## These are simply just ASCII texts, return them as such...
    result <- .rawToCharNull(dat)
    type[type == "(c) "] <- "copyright"
    class(result) <- c(paste0("IFF.", type), "IFF.ANY")
    return(result)
  } else {
    stop("Can't handle this chunk type (yet).")
  }
})

#' Coerce to and create IFFChunk objects
#' 
#' Convert `IFF.ANY` objects (created with [interpretIFFChunk()]) into
#' [IFFChunk()] objects. A basic [IFFChunk()] can also be
#' created with this method by providing the chunk type name.
#'
#' IFF data is stored in a [IFFChunk-class()] object when read from an
#' IFF file ([read.iff()]). These objects reflect the file structure
#' well, but the data is stored as `raw` information. IFF files can contain
#' a wide variety of information types, ranging from bitmap images to audio
#' clips. The raw information stored in [IFFChunk()] objects can
#' be interpreted into more meaningful representations that can be handled in
#' R. This is achieved with the [interpretIFFChunk()] method, which
#' returns `IFF.ANY` objects.
#' 
#' These `IFF.ANY` objects are a less strict representation of the
#' IFF Chunk, but are easier to handle in R. The interpretation method is lossy
#' and may not preserve all information in the `IFF.ANY` object.
#' The [IFFChunk-method()] can coerce `IFF.ANY` back
#' to the more strictly defined [IFFChunk-class()] objects.
#' Be careful with conversions between [IFFChunk-class()] and
#' `IFF.ANY` objects and vice versa, as information may get lost.
#' 
#' More detailed information about IFF chunks can be found in the IFF chunk registry
#' (see references).
#' 
#'  * `IFF.FORM` represents a FORM chunk, which is a container that can hold any kind of chunk.
#'    When interpreted, it is represented as a `list`, where each element is an interpreted chunk
#'    nested inside the FORM.
#'  * `IFF.BODY` represents the actual data in an IFF file. However, without context
#'    this chunk cannot be interpreted and is therefore interpreted as a vector of `raw` data.
#'  * `IFF.ANIM` represents an animation (ANIM) chunk. When interpreted, it will return a `list` where each
#'    element is an animation frame represented as an `IFF.ILBM` object. Each animation frame should be
#'    nested inside an ILBM chunk nested inside a FORM chunk, nested inside an ANIM chunk.
#'     * `IFF.ANHD` represents an ANimation HeaDer (ANHD) chunk. When interpreted,
#'       it returns a named `list` containing the
#'       following information:
#'        * `operation` is a `character` string indicating how the bitmap
#'          data for the animation frame is encoded. Can be one of the following:
#'          "`standard`", "`XOR`", "`LongDeltaMode`",
#'          "`ShortDeltaMode`", "`GeneralDeltamode`",
#'          "`ByteVerticalCompression`", "`StereoOp5`", or
#'          "`ShortLongVerticalDeltaMode`". Currently, only the
#'          ByteVerticalCompression is implemented in this package.
#'        * `mask` is a `vector` of 8 `logical` values. It is currently
#'          ignored.
#'        * `w` and `h` are positive `numeric` values, specifying
#'          the width and height of the frame (should be identical for all frames).
#'        * `x` and `y` are `numeric` values, specifying the plotting
#'          position for the frame.
#'        * `abstime` is a positive `numeric` value - currently unused - used for
#'          timing the frame relative to the time the first frame was displayed. In
#'          jiffies (1/60 sec).
#'        * `reltime` is a positive `numeric` value for timing the frame
#'          relative to time previous frame was displayed. In jiffies (1/60 sec).
#'        * `interleave` is currently unused. It should be set to 0.
#'        * `pad0` is a padding byte (`raw`) for future use.
#'        * `flags` is a `vector` of 32 `logical` values. They contain
#'          information on how the bitmap data is stored.
#'        * `pad1` are 16 padding bytes (`raw`) for future use.
#'     * `IFF.DPAN` represents an DPaint ANimation (DPAN) chunk. Some software will
#'       require this chunk to correctly derive the total number of frames in the animation.
#'       When interpreted, it will return a named `list` with the following elements:
#'        * `version` a `numeric` version number.
#'        * `nframes` a positive `numeric` value, indicating the number
#'          of frames in the animation.
#'        * `flags` a `vector` of 32 `logical` values. Ignored in
#'          this package as it was intended for future implementations.
#'     * `IFF.DLTA` represents a delta mode data chunk (DLTA). The first animation
#'       frame is stored as a normal InterLeaved BitMap (ILBM) image as described below.
#'       The following frames only store differences in bitmap data compared to the
#'       previous frames but is not interleaved. They are thus incorrectly embedded in
#'       an ILBM chunk (but is kept so for backward compatibility). When interpreted,
#'       a `grDevices` raster object is returned only showing the differences. It
#'       is not very meaningful to interpret these chunks on their own, but rather the
#'       entire parent ANIM chunk.
#'  * `IFF.ILBM` represents InterLeaved BitMap (ILBM) chunks. It is interpreted here as a
#'    raster image (see [grDevices::as.raster()]). ILBM chunks are usually nested inside
#'    a FORM container.
#'     * `IFF.BMHD` represents the header chunk of a bitmap (BMHD), and should always be present
#'       (nested inside) an ILBM chunk. It is interpreted as a named list containing the following elements:
#'        * `w` and `h` are positive `numeric` values specifying
#'          the bitmap width and height in pixels. Note that the width
#'          can be any positive whole number, whereas the bitmap data always
#'          has a width divisible by 16.
#'        * `x` and `y` are `numeric` values specifying the plotting
#'          position relative to the top left position of the screen.
#'          Although required in the bitmap header. It is ignored in the
#'          interpretation of bitmap images.
#'        * `nPlanes` is a positive value indicating the number of
#'          bitplanes in the image. The number of colours in an image
#'          can be calculated as `2^nPlanes`.
#'        * `Masking` indicates whether there are bitplanes that should
#'          be masked (i.e. are treated as transparent). It is a `character`
#'          string equalling any of the following: "`mskNone`",
#'          "`mskHasMask`", "`mskHasTransparentColour`",
#'          "`mskLasso`" or "`mskUnknown`". Only the first (no transparency)
#'          and third (one of the colours should be treated as transparent)
#'          id is currently interpreted correctly. The others are ignored.
#'          "`mskUnknown`" means that an undocumented mask is applied
#'          to the image.
#'        * `Compression` indicates whether the bitmap data is
#'          compressed. It is a `character` string that can equal any
#'          of the following: "`cmpNone`", "`cmpByteRun1`" or
#'          "`cmpUnknown`". The latter means an undocumented form of
#'          compression is applied and is currently ignored. In most cases
#'          bitmap data is compressed with the `cmpByteRun1` algorithm
#'          ([packBitmap()]). In some cases, bitmap data is not
#'          compressed (`cmpNone`).
#'        * `pad` is a `raw` byte that is only used to
#'          align data. It is ignored in the interpretation.
#'        * `transparentColour` is a `numeric` value that indicates
#'          which colour number in the palette should be treated as fully
#'          transparent (when `Masking` equals
#'          "`mskHasTransparentColour`").
#'        * `xAspect` and `yAspect` or positive `numeric`
#'          values that indicate the aspect ratio of
#'          the pixels in the image. Amiga screen modes allowed for some
#'          extreme pixel aspect ratios. These values are used to stretch
#'          the image to their intended display mode.
#'        * `pageWidth` and `pageHeight` are positive
#'          `numeric` values indicating the size of the screen in which
#'          the image should be displayed. They are ignored in the
#'          interpretation of the image.
#'     * `IFF.CMAP` represents the colour map (CMAP) or palette of a bitmap image. Although common,
#'       the chunk is optional and can be omitted from the parent ILBM chunk. It is interpreted as a
#'       vector of colours (i.e., a `character` string formatted as '#RRGGBB' or named colours such as
#'       'blue').
#'     * `IFF.CAMG` represents a chunk with information with respect
#'       to the display mode in which the bitmap image should be displayed.
#'       This information can be used to determine the correct pixel aspect
#'       ratio, or is sometimes required to correctly interpret the bitmap
#'       information. The `IFF.CAMG` chunk is interpreted as a named list
#'       containing the following elements:
#'        * `monitor`: a `factor` indicating the hardware monitor
#'          on which the image was created and should be displayed (see 
#'          [amiga_monitors()]).
#'        * `display.mode`: a `factor` indicating the display
#'          mode in which the image should be displayed (see
#'          [amiga_display_modes()]).
#'     * `IFF.CRNG` is an optional chunk nested in an ILBM chunk.
#'       It represents a `colour range' and is used to cycle through
#'       colours in the bitmap's palette in order to achieve
#'       animation effects. It is interpreted as a named list with the
#'       following elements. This chunk is currently not used with
#'       the interpretation of ILBM images.
#'        * `padding` are two `raw` padding bytes and are
#'          ignored when interpreted.
#'        * `rate` is a `numeric` value specifying the rate
#'          at which the colours are cycled. The rate is in steps per
#'          second.
#'        * `flags` is a flag that indicates how colours should
#'          be cycled. It is a `character` string that can equal
#'          any of the following: "`RNG_OFF`", "`RNG_ACTIVE`",
#'          "`RNG_REVERSE`" or "`RNG_UNKNOWN`". When equal to the
#'          first, colours are not cycled. When equal to the second, colours
#'          are cycled. When equal to the third, colours are cycled in
#'          reverse direction. When equal to the latter, an undocumented
#'          form of cycling is applied.
#'        * `low` and `high` are `numeric` indices of
#'          colours between which should be cycled. Only colour from
#'          index `low` up to index `high` are affected.
#'  * `IFF.8SVX` represents 8-bit sampled voice chunks (8SVX). The original
#'    Amiga supported 8-bit audio which could be stored using the IFF. 8SVX chunks
#'    can contain separate audio samples for each octave. 8SVX chunks are usually
#'    stored inside a FORM container. Its body chunk contains 8-bit PCM wave data that
#'    could be compressed. When the 8SVX chunk is
#'    interpreted with this package, a `list` is returned where each element
#'    represents an octave given as a [tuneR::Wave()] object. Possible
#'    chunks nested in 8SVX chunks and currently supported by this package are
#'    as follows.
#'     * `IFF.VHDR` represents voice header chunks (VHDR). It contains (meta-)information about
#'       the audio stored in the body of the parent 8SVX chunk. When interpreted, a named `list` is
#'       returned with the following elements:
#'        * `oneShotHiSamples` is a `numeric` value indicating how many samples there are in the
#'          audio wave of the first octave in the file, that should not be looped (repeated).
#'        * `repeatHiSamples` is a `numeric` value indicating how many samples there are in the
#'          audio wave of the first octave in the file, that should be looped (repeated).
#'        * `samplesPerHiCycle` is a `numeric` value specifying the
#'          number of samples per repeat cycle in the first octave, or 0 when unknown.
#'          The number of `repeatHiSamples` should be an exact multiple of
#'          `samplesPerHiCycle`.
#'        * `samplesPerSec` is a `numeric` value specifying the data
#'          sampling rate.
#'        * `ctOctave` a positive whole `numeric` value indicating how many octaves are included.
#'          In 8SVX files the audio wave is resampled for each octave. The wave data in the body starts with
#'          the audio sample in the highest octave (least number of samples). The data is then followed by
#'          each subsequent octave, where the number of samples increase by a factor of 2 for each octave.
#'        * `sCompression` is a `character` string indicating whether and how the wave data in the body
#'          is compressed. It can have one of the following values: "`sCmpNone`" (no compression),
#'          "`sCmpFibDelta`" ([deltaFibonacciCompress()]ion is applied), "`sCmpUnknown`" (an
#'          undocumented and unknown form of compression is applied).
#'        * `volume` is a numeric value between `0` (minimum) and `0x10000` (maximum) playback volume.
#'     * `IFF.CHAN` represents the channel chunk (CHAN). When interpreted it returns a named list
#'       with 1 named element:
#'       "`channel`". It's value can be one of the following `character` strings "`LEFT`", "`RIGHT`" or
#'       "`STEREO`". This indicates for how many (one or two) audio channels data is available in the
#'       body of the parent
#'       8SVX chunk. It also indicates two which channels the audio should be played back.
#'  * `IFF.ANNO`, `IFF.AUTH`, `IFF.CHRS`, `IFF.NAME`, `IFF.TEXT` and `IFF.copyright`
#'    are all unformatted text chunks that can be included optionally in any of the chunk types.
#'    Respectively, they
#'    represent an annotation, the author's name, a generic character string, the name of the work,
#'    generic unformatted text,
#'    and copyright text. They are interpreted as a `character` string.
#' @param x An S3 class `IFF.ANY` object that needs to be coerced into an
#' [IFFChunk-class()] object. `IFF.ANY` objects are created with the
#' [interpretIFFChunk()] method. `x` can also be a `character` string
#' of a IFF chunk type (e.g., "`FORM`" or "`BMHD`"). In that case an
#' [IFFChunk()] object of that type is created with some basic content.
#' @param ... Arguments passed onto methods underlying the interpretation of the
#' specific IFF chunks. Allowed arguments depend on the specific type of IFF chunk that
#' `x` represents.
#' @returns Returns an [IFFChunk-class()] representation of `x`.
#' @examples
#' ## load an IFF file
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## interpret the IFF file (in some cases information
#' ## will get lost in this step):
#' example.itpt <- interpretIFFChunk(example.iff)
#' 
#' ## now coerce back to a formal IFFChunk class object.
#' ## Only information in the interpreted object is used
#' ## The coerced object may therefore depart from the
#' ## original read from the file.
#' example.coerce <- IFFChunk(example.itpt)
#' 
#' ## and indeed the objects are not identical, as shown below.
#' ## In this case the difference is not disastrous, the order
#' ## of the colours in the palette have shifted. But be careful
#' ## with switching between formal IFFChunk objects and
#' ## interpreted IFF.ANY objects.
#' identical(example.iff, example.coerce)
#' 
#' ## It is also possible to create simple IFFChunk objects
#' ## by providing the desired chunk type identifier as a
#' ## character string.
#' 
#' ## This creates a basic bitmap header:
#' bmhd <- IFFChunk("BMHD")
#' 
#' ## This creates a basic colour palette:
#' cmap <- IFFChunk("CMAP")
#' @references <https://wiki.amigaos.net/wiki/IFF_FORM_and_Chunk_Registry>
#' @name IFFChunk-method
#' @rdname IFFChunk
#' @export
IFFChunk <- function (x, ...) {
  UseMethod("IFFChunk", x)
}

#' @rdname IFFChunk
#' @method IFFChunk character
#' @export
IFFChunk.character <- function(x, ...) {
  if (length(x) > 1) warning("x has more than 1 element, only the first will be used...")
  x <- x[[1]]
  if (x == "FORM") {
    return(new("IFFChunk"))
  } else if (x == "BODY") {
    return(new("IFFChunk", chunk.type = "BODY"))
  } else if (x == "CMAP") {
    return(new("IFFChunk", chunk.type = "CMAP",
               chunk.data = list(colourToAmigaRaw(grDevices::gray(round(seq(0, 15, length.out = 16))/15),
                                                  "12 bit", "3"))))
  } else if (x == "BMHD") {
    BMHD <- list(
      w                 = 16,
      h                 = 16,
      x                 = 0,
      y                 = 0,
      nPlanes           = 4,
      Masking           = "mskNone",
      Compression       = "cmpByteRun1",
      pad               = raw(1),
      transparentColour = 0,
      xAspect           = 44,
      yAspect           = 44,
      pageWidth         = 320,
      pageHeight        = 200
    )
    class(BMHD) <- "IFF.BMHD"
    return(IFFChunk(BMHD))
  } else if(x == "CAMG") {
    return(.inverseViewPort("LORES_KEY", "DEFAULT_MONITOR_ID"))
  } else if (x == "CRNG") {
    dat <- list(
      padding = raw(2),
      rate = 0,
      flag = "RNG_OFF",
      low = 0,
      high = 0
    )
    class(dat) <- "IFF.CRNG"
    return(IFFChunk(dat))
  } else if (x == "ILBM") {
    bm.dep <- 1
    bmhd <- IFFChunk("BMHD")
    bmhd.int <- interpretIFFChunk(bmhd)
    bmhd.int$nPlanes <- 1
    bmhd <- IFFChunk(bmhd.int)
    pal <- c("black", "white")
    class(pal) <- "IFF.CMAP"
    pal <- IFFChunk(pal)
    body <- new("IFFChunk", chunk.type = "BODY",
                chunk.data = list(
                  as.raw(c(0x01, 0x03, 0xc0, 0x01, 0x0c, 0x30, 0x01, 0x10, 0x08, 0x01,
                           0x20, 0x04, 0x01, 0x40, 0x02, 0x01, 0x44, 0x22, 0x01, 0x8e,
                           0x71, 0x01, 0x84, 0x21, 0x01, 0x80, 0x01, 0x01, 0x8e, 0x71,
                           0x01, 0x44, 0x22, 0x01, 0x44, 0x22, 0x01, 0x23, 0xc4, 0x01,
                           0x10, 0x08, 0x01, 0x0c, 0x30, 0x01, 0x03, 0xc0))
                ))
    return(new("IFFChunk", chunk.type = "ILBM",
               chunk.data = list(
                 bmhd,
                 pal,
                 IFFChunk("CAMG"),
                 body
               )
    ))
  } else if (x == "ANIM") {
    anhd   <- lapply(1:4, function(y) {
      anhdy <- interpretIFFChunk(IFFChunk("ANHD"))
      anhdy$abstime <- y*2
      IFFChunk(anhdy)
    })
    frame1 <- IFFChunk("ILBM")@chunk.data
    frame1 <- c(frame1[1:(length(frame1) - 1)], anhd[[1]], frame1[length(frame1)])
    frame1 <- new("IFFChunk", chunk.type = "ILBM", chunk.data = frame1)
    dlta2  <- new("IFFChunk", chunk.type = "DLTA", chunk.data = list(unPackBitmap(as.raw(c(0x03, 0x00, 0x00, 0x00, 0x40,
                                                                                           0xc4, 0x00, 0x06, 0x03, 0x05,
                                                                                           0x83, 0x02, 0x11, 0x61, 0x08)))))
    dlta3  <- new("IFFChunk", chunk.type = "DLTA", chunk.data = list(unPackBitmap(as.raw(c(0xc1, 0x00)))))
    result <- new("IFFChunk", chunk.type = "ANIM", chunk.data = list(
      new("IFFChunk", chunk.type = "FORM", chunk.data = list(frame1)),
      new("IFFChunk", chunk.type = "FORM", chunk.data = list(
        new("IFFChunk", chunk.type = "ILBM", chunk.data = list(anhd[[2]], dlta2))
      )),
      new("IFFChunk", chunk.type = "FORM", chunk.data = list(
        new("IFFChunk", chunk.type = "ILBM", chunk.data = list(anhd[[3]], dlta3))
      )),
      new("IFFChunk", chunk.type = "FORM", chunk.data = list(
        new("IFFChunk", chunk.type = "ILBM", chunk.data = list(anhd[[4]], dlta3))
      ))
    ))
    return(result)
  } else if (x == "ANHD") {
    result <- list(
      operation       = "ByteVerticalCompression",
      mask            = rep(FALSE, 8),
      w               = 16,
      h               = 16,
      x               = 0,
      y               = 0,
      abstime         = 2,
      reltime         = 2,
      interleave      = 0,
      pad0            = raw(1),
      flags           = rep(FALSE, 32),
      pad1            = raw(16)
    )
    class(result) <- "IFF.ANHD"
    return(IFFChunk(result))
  } else if (x == "DPAN") {
    result <- list(
      version = 4,
      nframes = 2,
      flags   = rep(FALSE, 32)
    )
    class(result) <- "IFF.DPAN"
    return(IFFChunk(result))
  } else if (x == "CHAN") {
    chan <- list(channel = "LEFT")
    class(chan) <- "IFF.CHAN"
    return(IFFChunk(chan))
  } else if(x == "VHDR") {
    vhdr <- list(
      oneShotHiSamples = 0,
      repeatHiSamples = 0,
      samplesPerHiCycle = 0,
      samplesPerSec = 16574.28,
      ctOctave = 1,
      sCompression = "sCmpNone",
      volume = 0x10000L
    )
    class(vhdr) <- "IFF.VHDR"
    return(IFFChunk(vhdr))
  } else if(x == "8SVX") {
    return(WaveToIFF(Wave(round(128*sin((0:(10*62))*2*pi/62) - 0.5) + 128,
                          bit = 8,
                          samp.rate = 16574.28),
                     loop.start = 0))
  } else if (x %in% c("ANNO", "AUTH", "CHRS", "NAME", "TEXT", "(c) ")) {
    return (new("IFFChunk", chunk.type = x, chunk.data = raw(0)))
  }
  stop("Cannot generate this chunk type...")
}

#' @rdname IFFChunk
#' @method IFFChunk IFF.FORM
#' @export
IFFChunk.IFF.FORM <- function(x, ...) {
  new("IFFChunk", chunk.type = "FORM", chunk.data = lapply(x, IFFChunk, ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.BODY <- function(x, ...) {
  x <- as.raw(x, ...)
  return(new("IFFChunk", chunk.type = "BODY", chunk.data = list(x)))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.ANNO <- function(x, ...) {
  return(.text.chunk(x, "ANNO", ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.AUTH <- function(x, ...) {
  return(.text.chunk(x, "AUTH", ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.CHRS <- function(x, ...) {
  return(.text.chunk(x, "CHRS", ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.NAME <- function(x, ...) {
  return(.text.chunk(x, "NAME", ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.TEXT <- function(x, ...) {
  return(.text.chunk(x, "TEXT", ...))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.copyright <- function(x, ...) {
  return(.text.chunk(x, "(c) ", ...))
}

#' @rdname plot
#' @name plot
#' @export
plot.IFFChunk <- function(x, y, ...) {
  if (missing(y)) y <- NULL
  graphics::plot(interpretIFFChunk(x), ...)
}

#' @rdname plot
#' @name plot
#' @export
plot.IFF.FORM <- function(x, y, ...) {
  invisible(lapply(x, graphics::plot, ...))
}

#' @export
print.IFFChunk <- function(x, ...) {
  ctypes <- function(z) {
    result <- paste0("- ", z@chunk.type)
    temp <- unlist(lapply(z@chunk.data, function(a) {
      if (inherits(a, "IFFChunk")) return(ctypes(a)) else return(NULL)
    }))
    if (length(temp) > 0) result <- c(result, paste0("  ", temp))
    return(result)
  }
  typ <- paste0(ctypes(x), collapse = "\n")
  cat(typ)
}

#' @export
as.list.IFFChunk <- function(x, ...) {
  slotList <- function(x) {
    if (inherits(x@chunk.data[[1]], "raw")) {
      return(x@chunk.data[[1]])
    } else {
      result <- lapply(x@chunk.data, slotList)
      names(result) <- unlist(lapply(x@chunk.data, function(x) x@chunk.type))
      id.num <- stats::ave(seq_along(names(result)), names(result), FUN = seq_along)
      names(result) <- sprintf("%s_%i", names(result), id.num)
      return(result)
    }
  }
  return(slotList(x))
}

.rawToIFFChunk <- function(x, skip.size = FALSE) {
  offset <- 0
  chunks <- list()
  while (offset < length(x)) {
    chunk.type <- rawToChar(x[offset + 1:4])
    if (skip.size) {
      offset <- offset + 4
      return(list(methods::new("IFFChunk",
                               chunk.type = chunk.type,
                               chunk.data = .rawToIFFChunk(x[(offset + 1):length(x)]))))
    }
    chunk.size <- .rawToAmigaInt(x[offset + 5:8], 32, FALSE)
    offset <- offset + 8
    chunk.data <- x[offset + 1:chunk.size]
    offset <- offset + chunk.size
    ## Data should be word aligned:
    if ((offset %% 2) == 1) offset <- offset + 1
    if (chunk.type %in% "FORM") { ## contains nested chunks, the container directly following 'FORM' will not specify the chunk length (hence the skip)
      chunk.data <- .rawToIFFChunk(chunk.data, skip.size = TRUE)
    } else if (chunk.type %in% c("ILBM", "8SVX", "ANIM")) {
      chunk.data <- .rawToIFFChunk(chunk.data, skip.size = FALSE)
    }
    if (inherits(chunk.data, "raw")) chunk.data <- list(chunk.data)
    chunks[[length(chunks) + 1]] <- methods::new("IFFChunk",
                                                 chunk.type = chunk.type,
                                                 chunk.data = chunk.data)
  }
  return(chunks)
}

.text.chunk <- function(x, type, ...) {
  return(new("IFFChunk", chunk.type = type, chunk.data = list(charToRaw(x))))
}

setGeneric("rawToIFFChunk", function(x) standardGeneric("rawToIFFChunk"))

#' Coerce raw data to an IFFChunk class object
#'
#' Coerce raw data, as it would be stored in the Interchange File Format (IFF), and
#' convert it into an [IFFChunk()] class object.
#'
#' This method should work for all IFF chunk types that are implemented in this package
#' (see [IFFChunk-method()] for details). For non-implemented chunks this method
#' may work properly as long as the chunks are nested inside a FORM type container chunk.
#' This method is provided for your convenience, but it is recommended to import IFFChunk
#' methods using the [read.iff()] function. Use [AmigaFFH::as.raw()]
#' to achieve the inverse of this method.
#'
#' @docType methods
#' @rdname rawToIFFChunk
#' @name rawToIFFChunk
#' @aliases rawToIFFChunk,raw-method
#' @param x A vector of raw data that needs to be converted into a [IFFChunk()]
#' class object.
#' @returns Returns an [IFFChunk()] class object based on `x`.
#' @examples
#' ## Get an IFFChunk object:
#' example.iff <- read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH"))
#' 
#' ## Coerce it to raw data:
#' example.raw <- as.raw(example.iff)
#' 
#' ## Coerce raw data to IFF chunk:
#' example.iff.new <- rawToIFFChunk(example.raw)
#' 
#' ## These conversions were non-destructive:
#' identical(example.iff, example.iff.new)
#' @family iff.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
setMethod("rawToIFFChunk", "raw", function(x) {
  result <- .rawToIFFChunk(x)
  if (length(result) == 1) result <- result[[1]]
  return(result)
})
