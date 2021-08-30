#' @rdname IFFChunk
#' @export
IFFChunk.IFF.CHAN <- function(x, ...) {
  x <- c(2, 4, 6)[c("LEFT", "RIGHT", "STEREO") == x$channel[[1]]]
  x <- list(.amigaIntToRaw(x, 32, F))
  return(new("IFFChunk", chunk.type = "CHAN", chunk.data = x))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.VHDR <- function(x, ...) {
  compr <- x$sCompression
  compr <- which(c("sCmpNone", "sCmpFibDelta") == x$sCompression) - 1
  ## if the compression type is unknown, set to 0
  if (length(compr) == 0) compr <- 0
  result <- c(
    .amigaIntToRaw(c(x$oneShotHiSamples,
                                 x$repeatHiSamples,
                                 x$samplesPerHiCycle), 32, F),
    .amigaIntToRaw(x$samplesPerSec, 16, F),
    .amigaIntToRaw(c(x$ctOctave, compr), 8, F),
    .amigaIntToRaw(x$volume, 32, F)
  )
  return(new("IFFChunk", chunk.type = "VHDR", chunk.data = list(result)))
}

#' @rdname IFFChunk
#' @export
IFFChunk.IFF.8SVX <- function(x, ...) {
  ## check lengths of the waveforms in the list:
  testlen <- unlist(lapply(x, length))/length(x[[1]])
  if (!all(testlen == 2^((1:length(x)) - 1))) {
    warning("Each subsequent waveform in the list should have twice the length of the previous waveform. Continuing with the first element only.")
    x <- x[1]
  }
  if (length(unique(unlist(lapply(x, function(y) y@stereo)))) != 1) {
    warning("x contains both mono and stereo samples, can't combine these. Continuing with the first element of x only.")
    x <- x[1]
  }
  if (!all(unlist(lapply(x, function(y) y@bit)) == 8)) {
    stop("All waves in x should be 8-bit.")
  }
  if (!all(unlist(lapply(x, function(y) y@pcm)))) {
    stop("All waves in x should be pcm formatted.")
  }
  wav  <- .amigaIntToRaw(c(do.call(c, lapply(x, function(y) y@left)) - 128,
                                       do.call(c, lapply(x, function(y) y@right)) - 128), 8, T)
  if ((length(wav) %% 2) != 0) wav <- c(wav, raw(1))
  vhdr <- list(
    oneShotHiSamples = length(x[[1]]),
    repeatHiSamples = 0,
    samplesPerHiCycle = 0,
    samplesPerSec = x[[1]]@samp.rate,
    ctOctave = length(x),
    sCompression = "sCmpNone",
    volume = 0x10000L
  )
  class(vhdr) <- "IFF.VHDR"
  vhdr <- IFFChunk(vhdr)
  chan <- list(channel = ifelse(x[[1]]@stereo, "STEREO", "LEFT"))
  class(chan) <- "IFF.CHAN"
  chan <- IFFChunk(chan)
  return(new("IFFChunk", chunk.type = "8SVX", chunk.data = list(
    vhdr,
    chan,
    new("IFFChunk", chunk.type = "BODY", chunk.data = list(wav))
  )))
}

#' @rdname plot
#' @name plot
#' @export
plot.IFF.8SVX <- function(x, y, ...) {
  invisible(lapply(x, tuneR::plot, ...))
}

#' Playing Amiga audio data
#' 
#' A wrapper for \code{\link{tuneR}}-package's \code{\link[tuneR]{play}} routine. Allowing it to play
#' Amiga audio (for instance stored in an 8SVX Interchange File Format).
#' 
#' A wrapper for \code{\link{tuneR}}-package's \code{\link[tuneR]{play}} routine. It will try to play
#' audio using an external audio player. When 8SVX audio is played, each octave is played separately.
#' When a FORM container contains multiple 8SVX samples, they are also played successively.
#' 
#' Note that a separate package is developed to interpret and play ProTracker modules and samples
#' (\code{\link[ProTrackR]{ProTrackR}}).
#' @rdname play
#' @name play
#' @aliases play,ANY-method
#' @param object An \code{\link{IFFChunk-class}} object that needs to be played. The \code{\link{IFFChunk}}
#' should be of type FORM, containing an 8SVX chunk, or an 8SVX itself. \code{object} can also be of class
#' \code{IFF.FORM} or \code{IFF.8SVX}. See \code{\link[tuneR]{play}} for other objects that can be played.
#' @param player Path to the external audio player. See \code{\link[tuneR]{play}} for more details.
#' @param ... Arguments passed onto the tuneR \code{\link{play}} routine.
#' @return Returns a list of data returned by tuneR's \code{\link[tuneR]{play}}, for which the output
#' is undocumented.
#' @examples
#' \dontrun{
#' ## First get an audio sample from the ProTrackR package
#' snare.samp <- ProTrackR::PTSample(ProTrackR::mod.intro, 2)
#' 
#' ## Coerce it into an IFFChunk object:
#' snare.iff <- WaveToIFF(snare.samp)
#' 
#' ## Play the 8SVX sample:
#' play(snare.iff)
#' }
#' @author Pepijn de Vries
#' @export
setMethod("play", "ANY", function(object, player = NULL, ...) {
  if ("IFF.FORM" %in% class(object)) {
    invisible(lapply(object, function(x) {
      lapply(x, function(y) {
        tuneR::play(y, ...)
      })
    }))
  } else if ("IFF.8SVX" %in% class(object)) {
    invisible(lapply(object, function(x) {
      tuneR::play(x, ...)
    }))
  } else {
    stop(sprintf("Sorry, can't play %s object", class(object)))
  }
})

#' @rdname play
#' @name play
#' @aliases play,IFFChunk-method
#' @export
setMethod("play", "IFFChunk", function(object, player = NULL, ...) {
  play(interpretIFFChunk(object), player, ...)
})

#' Convert WaveMC objects into an Interchange File Format object
#'
#' Convert \code{\link[tuneR]{WaveMC}} objects (or objects that can be coerced to
#' \code{WaveMC} objects) into an \code{\link{IFFChunk-class}} object which
#' can be stored as a valid Iterchange File Format (\code{\link{write.iff}}).
#'
#' \code{\link[tuneR]{WaveMC}} objects can be read from contemporary file containers
#' with \code{\link[tuneR]{readWave}} or \code{\link[tuneR]{readMP3}}. With this
#' function such objects can be converted into an \code{\link{IFFChunk-class}} object
#' which can be stored conform the Interchange File Format (\code{\link{write.iff}}).
#' 
#' When \code{x} is not a pcm formatted 8-bit sample, \code{x} will first be
#' normalised and scaled to a pcm-formatted 8-bit sample using
#' \code{\link[tuneR]{normalize}}. If you don't like the result you need to convert
#' the sample to 8-bit pcm yourself before calling this function.
#'
#' @rdname WaveToIFF
#' @name WaveToIFF
#' @param x A \code{\link[tuneR]{WaveMC}} object that needs to be converted into an \code{\link{IFFChunk}} object. \code{x}
#' can also be any other class object that can be coerced into a \code{\link[tuneR]{WaveMC}} object. \code{\link[tuneR]{Wave}}
#' and \code{\link[ProTrackR:PTSample-class]{PTSample}} objects are therefore also allowed.
#' @param loop.start If the sample should be looped from a specific position to the
#' end of the sample, this argument specifies the starting position in samples (with
#' a base of 0) for looping. \code{loop.start} therefore should be a whole non-negative
#' number. When set to \code{NA} or negative values, the sample will not be looped.
#' @param octaves A whole positive \code{numeric} value indicating the number of octaves
#' that should be stored in the resulting IFF chunk. The original wave will be resampled
#' for each value larger than 1. Each subsequent octave will contain precisely twice
#' as many samples as the previous octave.
#' @param compress A \code{character} string indicating whether compression should be applied to the waveform. "\code{sCmpNone}"
#' (default) applies no compression, "\code{sCmpFibDelta}" applies the lossy \code{\link{deltaFibonacciCompress}}ion.
#' @param ... Currently ignored.
#' @return Returns an \code{\link{IFFChunk-class}} object with a FORM container that
#' contains an 8SVX waveform based on \code{x}.
#' @examples
#' \dontrun{
#' ## First get an audio sample from the ProTrackR package
#' snare.samp <- ProTrackR::PTSample(ProTrackR::mod.intro, 2)
#' 
#' ## The sample can easily be converted into an IFFChunk:
#' snare.iff <- WaveToIFF(snare.samp)
#' 
#' ## You could also first convert the sample into a Wave object:
#' snare.wav <- as(snare.samp, "Wave")
#' 
#' ## And then convert into an IFFChunk. The result is the same:
#' snare.iff <- WaveToIFF(snare.wav)
#' 
#' ## You could also use a sine wave as input (although you will get some warnings).
#' ## This will work because the vector of numeric data can be coerced to
#' ## a WaveMC object
#' sine.iff <- WaveToIFF(sin((0:2000)/20))
#' }
#' @family iff.operations
#' @references \url{https://en.wikipedia.org/wiki/8SVX}
#' @author Pepijn de Vries
#' @export
WaveToIFF <- function(x, loop.start = NA, octaves = 1, compress = c("sCmpNone", "sCmpFibDelta"), ...) {
  octaves <- round(octaves[[1]])
  loop.start <- round(loop.start[[1]])
  x <- methods::as(x, "WaveMC")
  compress <- match.arg(compress)
  if (x@bit != 8 || !x@pcm) {
    warning(sprintf("Original %i-bit wave is normalized to an 8-bit wave.", x@bit))
    ## Note tuneR's implementation let's 8-bit audio range from 0-254 instead of 0-255
    ## writeWave from the same package uses the range of 0-255. Maybe contact package
    ## maintainer to check whether this discrapency is intentional
    x <- tuneR::normalize(x, "8", pcm = T)
  }
  if (is.null(colnames(x@.Data))) colnames(x@.Data) <- MCnames$name[1:ncol(x@.Data)]
  x <- tuneR::WaveMC(data  = cbind(FL = rowMeans(x@.Data[,grepl("L", colnames(x@.Data)), drop = F]),
                                   FR = rowMeans(x@.Data[,grepl("R", colnames(x@.Data)), drop = F])),
                     bit   = x@bit,
                     samp.rate = x@samp.rate,
                     pcm   = x@pcm)
  if (any(is.nan(x@.Data[,"FL"]))) x@.Data <- x@.Data[,colnames(x@.Data) != "FL", drop = F]
  if (any(is.nan(x@.Data[,"FR"]))) x@.Data <- x@.Data[,colnames(x@.Data) != "FR", drop = F]
  if (octaves > 1) {
    temp <- lapply(2:octaves, function(y) {
      tuneR::WaveMC(data = apply(x@.Data, 2, function(z) {
        ProTrackR::resample(z, 1, y)
      }),
      bit = x@bit,
      samp.rate = x@samp.rate,
      pcm = x@pcm)
    })
  }
  x <- list(x)
  if (octaves > 1) {
    x<- c(x, temp)
    rm(temp)
  }
  oneshot <- length(x[[1]]@.Data[,1])
  repeatsamp <- 0
  if (!is.na(loop.start) && loop.start >= 0) {
    if (loop.start >= oneshot) stop("'loop.start' cannot be equal or larger than the sample length.")
    repeatsamp <- oneshot - loop.start
    oneshot <- loop.start
  }
  vhdr <- list(
    oneShotHiSamples = oneshot,
    repeatHiSamples = repeatsamp,
    samplesPerHiCycle = 0,
    samplesPerSec = x[[1]]@samp.rate,
    ctOctave = octaves,
    sCompression = compress,
    volume = 0x10000L
  )
  class(vhdr) <- "IFF.VHDR"
  vhdr <- IFFChunk(vhdr)
  chan <- list(channel = ifelse(all(c("FL", "FR") %in% colnames(x[[1]]@.Data)),
                                "STEREO",
                                ifelse("FL" %in% colnames(x[[1]]@.Data), "LEFT", "RIGHT")))
  class(chan) <- "IFF.CHAN"
  if (chan$channel == "STEREO" && octaves > 1) warning(paste0("Multiple octaves in stereo may not be handled correctly by Amiga software.\n",
                                                              "Convert your sample to mono first or only use 1 octave to avoid problems."))
  chan <- IFFChunk(chan)

  wav <- unlist(lapply(1:ncol(x[[1]]@.Data), function(z) {
    unlist(lapply(x, function(y) .amigaIntToRaw(y@.Data[,z] - 128, 8, T)))
  }))
  
  if (compress == "sCmpFibDelta") {
    wav <- deltaFibonacciCompress(wav)
  }

  ## Pad with zero if the body has an odd length:
  if ((length(wav) %% 2) !=0) wav <- c(wav, raw(1))
  wav <- list(wav)
  wav <- list(vhdr, chan, new("IFFChunk", chunk.type = "BODY", chunk.data = wav))
  wav <- list(new("IFFChunk", chunk.type = "8SVX", chunk.data = wav))
  return(new("IFFChunk", chunk.type = "FORM", chunk.data = wav))
}