.amigaIntToRaw <- function(x, bits = 8, signed = FALSE) {
  x <- round(x)
  if (!signed && any(x < 0)) stop("negative values not allowed for unsigned values.")
  val.range <- c(0, 2^bits - 1)
  if (signed) val.range <- c(-(2^bits)/2,(2^bits)/2 - 1)
  if (any(x < val.range[1]) || any(x > val.range[2])) {
    warning("One or more values are out of the specified bit-range. They will be clipped...")
    x[x < val.range[1]] <- val.range[1]
    x[x > val.range[2]] <- val.range[2]
  }
  if (signed) x[x < 0] <- (2^bits) + x[x < 0]
  ## used later on to reorder bits for the little-endian bytes
  idx <- sort(rep(((1:(bits/8)) - 1)*8, 8), TRUE) + rep(1:8, bits/8)
  result <- unlist(lapply(x, function(y) {
    bitlist <- NULL
    while (y > 0) {
      bitlist <- c(bitlist, y %% 2)
      y <- floor(y/2)
    }
    bitlist <- c(bitlist, numeric(bits - length(bitlist)))
    res <- packBits(as.logical(bitlist)[idx], "raw")
    return(res)
  }))
  return(result)
}

.bitmapToRaw <- function(x, invert.bytes = TRUE, invert.longs = TRUE) {
  # 'x' should be anything that is accepted by packBits
  if (!all("logical" %in% c(typeof(invert.bytes), typeof(invert.longs)))) stop ("Both 'invert.bytes' and 'invert.longs' should be a logical value.")
  if (length(invert.bytes) != 1 || length(invert.longs) != 1) stop("Both 'invert.bytes' and 'invert.longs' should have a length of 1.")
  true.len <- length(x)
  ## pad with zeros
  x <- c(x, raw(32 - (true.len %% 32)))
  len <- length(x)/8
  if (invert.bytes) {
    ord <- 1 + sort(rep((0:(len - 1))*8, 8)) + (7:0)
  } else {
    ord <- 1:(8*len)
  }
  if (invert.longs) {
    l2 <- ceiling(8*len/32)
    ord2 <- 1 + sort(rep((0:(l2 - 1))*32, 32)) + (31:0)
    ord2 <- ord2[1:(8*len)]
    x <- x[ord2]
  }
  ## order results and trim length to correspond with input
  x <- packBits(x[ord])[1:ceiling(true.len/8)]
  return(x)
}

.rawToAmigaInt <- function(x, bits = 8, signed = FALSE) {
  # Convert raw values into Amiga integers (BYTE (8 bit signed), UBYTE (8 bit unsigned),
  # WORD (16 bit signed), UWORD (16 bit unsigned), LONG (32 bit signed), ULONG (32 bit unsigned))
  if ((bits %% 8) != 0 || bits < 8) stop("Bits should be positive, it should also be a multitude of 8 (or 8 itself).")
  # pad x with zeros when it does not consist of a multitude of specified bits
  x <- c(x, raw(length(x) %% (bits/8)))
  i.start <- 1:floor(length(x)/(bits/8))
  i.stop  <- i.start*(bits/8)
  i.start <- (i.start - 1)*(bits/8) + 1
  result <- mapply(function(start, stop) {
    y <- x[start:stop]
    result <- as.numeric(unlist(lapply(y, function(z) rev(rawToBits(z)))))
    result <- sum(2^(which(rev(result) == as.raw(0x01)) - 1))
    return(result)
  }, start = i.start, stop = i.stop)
  if (signed) {
    result[result >= (2^bits)/2] <- result[result >= (2^bits)/2] - (2^bits)
    return(result)
  } else {
    return(result)
  }
}

.rawToBitmap <- function(x, invert.bytes = FALSE, invert.longs = TRUE) {
  if (typeof(x) != "raw") stop("Argument 'x' should be a vector of raw data.")
  if (!all("logical" %in% c(typeof(invert.bytes), typeof(invert.longs)))) stop ("Both 'invert.bytes' and 'invert.longs' should be a logical value.")
  if (length(invert.bytes) != 1 || length(invert.longs) != 1) stop("Both 'invert.bytes' and 'invert.longs' should have a length of 1.")
  ## pad data with zeros and trim at the end
  true.len <- length(x)
  x <- c(x, raw(4 - (true.len %% 4)))
  len <- length(x)
  if (invert.longs) {
    l2 <- ceiling(len/4)
    ord2 <- 1 + sort(rep((0:(l2 - 1))*4, 4)) + (3:0)
    ord2 <- ord2[1:len]
    x <- x[ord2]
  }
  if (invert.bytes) {
    ord <- 1 + sort(rep((0:(len - 1))*8, 8)) + (7:0)
  } else {
    ord <- 1:(8*len)
  }
  ## trim the result to correspond with the input length (data might get lost!)
  rawToBits(x)[ord][1:(true.len*8)]
}
