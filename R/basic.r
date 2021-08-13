#' The S3 AmigaBasic class
#' 
#' XXX
#' 
#' XXX add to package documentation
#' 
#' @note No documentation XXX all the product of reverse engineering from my part
#' @docType class
#' @aliases AmigaBasic
#' @name AmigaBasic
#' @rdname AmigaBasic
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @examples
#' \dontrun{
#' ## XXX
#' }
NULL

#' The S3 AmigaBasicShape class
#' 
#' XXX implement plot+as.raster function
#' 
#' XXX add to package documentation
#' 
#' @docType class
#' @aliases AmigaBasicShape
#' @name AmigaBasicShape
#' @rdname AmigaBasicShape
#' @family AmigaBasicShape.operations
#' @author Pepijn de Vries
#' @examples
#' \dontrun{
#' ## XXX
#' }
NULL

.amigabasic_commands <- read.table(text =
"code1,code2,command
80,00,ABS
81,00,ASC
82,00,ATN
83,00,CALL
84,00,CDBL
85,00,CHR$
86,00,CINT
87,00,CLOSE
88,00,COMMON
89,00,COS
8a,00,CVD
8b,00,CVI
8c,00,CVS
8d,00,DATA
8e,00,ELSE
8f,00,EOF
90,00,EXP
91,00,FIELD
92,00,FIX
93,00,FN
94,00,FOR
95,00,GET
96,00,GOSUB
97,00,GOTO
98,00,IF
99,00,INKEY$
9a,00,INPUT
9b,00,INT
9c,00,LEFT$
9d,00,LEN
9e,00,LET
9f,00,LINE
a1,00,LOC
a2,00,LOF
a3,00,LOG
a4,00,LSET
a5,00,MID$
a6,00,MKD$
a7,00,MKI$
a8,00,MKS$
a9,00,NEXT
aa,00,ON
ab,00,OPEN
ac,00,PRINT
ad,00,PUT
ae,00,READ
af,00,REM
af,e8,'
b0,00,RETURN
b1,00,RIGHT$
b2,00,RND
b3,00,RSET
b4,00,SGN
b5,00,SIN
b6,00,SPACE$
b7,00,SQR
b8,00,STR$
b9,00,STRING$
ba,00,TAN
bc,00,VAL
bd,00,WEND
be,ec,WHILE
bf,00,WRITE
c0,00,ELSEIF
c1,00,CLNG
c2,00,CVL
c3,00,MKL$
c4,00,AREA
e3,00,STATIC
e4,00,USING
e5,00,TO
e6,00,THEN
e7,00,NOT
e9,00,>
ea,00,=
eb,00,<
ec,00,+
ed,00,-
ee,00,*
ef,00,/
f0,00,^
f1,00,AND
f2,00,OR
f3,00,XOR
f4,00,EQV
f5,00,IMP
f6,00,MOD
f7,00,\\
f8,81,CHAIN
f8,82,CLEAR
f8,83,CLS
f8,84,CONT
f8,85,CSNG
f8,86,DATE$
f8,87,DEFINT
f8,88,DEFSNG
f8,89,DEFDBL
f8,8a,DEFSTR
f8,8b,DEF
f8,8c,DELETE
f8,8d,DIM
f8,8f,END
f8,90,ERASE
f8,91,ERL
f8,92,ERROR
f8,93,ERR
f8,94,FILES
f8,95,FRE
f8,96,HEX$
f8,97,INSTR
f8,98,KILL
f8,9a,LLIST
f8,9b,LOAD
f8,9c,LPOS
f8,9d,LPRINT
f8,9e,MERGE
f8,9f,NAME
f8,a0,NEW
f8,a1,OCT$
f8,a2,OPTION
f8,a3,PEEK
f8,a4,POKE
f8,a5,POS
f8,a6,RANDOMIZE
f8,a8,RESTORE
f8,a9,RESUME
f8,aa,RUN
f8,ab,SAVE
f8,ad,STOP
f8,ae,SWAP
f8,af,SYSTEM
f8,b0,TIME$
f8,b1,TRON
f8,b2,TROFF
f8,b3,VARPTR
f8,b4,WIDTH
f8,b5,BEEP
f8,b6,CIRCLE
f8,b8,MOUSE
f8,b9,POINT
f8,ba,PRESET
f8,bb,PSET
f8,bd,TIMER
f8,be,SUB
f8,bf,EXIT
f8,c0,SOUND
f8,c2,MENU
f8,c3,WINDOW
f8,c5,LOCATE
f8,c6,CSRLIN
f8,c7,LBOUND
f8,c8,UBOUND
f8,c9,SHARED
f8,ca,UCASE$
f8,cb,SCROLL
f8,cc,LIBRARY
f8,d2,PAINT
f8,d3,SCREEN
f8,d4,DECLARE
f8,d5,FUNCTION
f8,d6,DEFLNG
f8,d7,SADD
f8,d8,AREAFILL
f8,d9,COLOR
f8,da,PATTERN
f8,db,PALETTE
f8,dc,SLEEP
f8,dd,CHDIR
f8,de,STRIG
f8,df,STICK
f9,f4,OFF
f9,f5,BREAK
f9,f6,WAIT
f9,f8,TAB
f9,f9,STEP
f9,fa,SPC
f9,fb,OUTPUT
f9,fc,BASE
f9,fd,AS
f9,fe,APPEND
f9,ff,ALL
fa,80,WAVE
fa,81,POKEW
fa,82,POKEL
fa,83,PEEKW
fa,84,PEEKL
fa,85,SAY
fa,86,TRANSLATE$
fa,87,OBJECT.SHAPE
fa,88,OBJECT.PRIORITY
fa,89,OBJECT.X
fa,8a,OBJECT.Y
fa,8b,OBJECT.VX
fa,8c,OBJECT.VY
fa,8d,OBJECT.AX
fa,8e,OBJECT.AY
fa,8f,OBJECT.CLIP
fa,90,OBJECT.PLANES
fa,91,OBJECT.HIT
fa,92,OBJECT.ON
fa,93,OBJECT.OFF
fa,94,OBJECT.START
fa,95,OBJECT.STOP
fa,96,OBJECT.CLOSE
fa,97,COLLISION
fb,ff,PTAB", header = T, sep = ",", quote = "", as.is = T)
.amigabasic_commands$code1 <- as.raw(paste0("0x", .amigabasic_commands$code1))
.amigabasic_commands$code2 <- as.raw(paste0("0x", .amigabasic_commands$code2))

#' XXX
#'
#' XXX
#'
#' XXX
#' 
#' @rdname rawToAmigaBasic
#' @name rawToAmigaBasic
#' @param x A \code{vector} of \code{raw} data that is to be converted
#' into an \code{\link{AmigaBasic}} class object.
#' @param ... Currently ignored.
#' @return An \code{\link{AmigaBasic}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' ## First create an AmigaBAsic object:
#' bas <- as.AmigaBasic("PRINT \"Hello world!\"")
#' 
#' ## Make it raw:
#' bas.raw <- as.raw(bas)
#' 
#' ## Now convert it back to an AmigaBasic object:
#' bas <- rawToAmigaBasic(bas.raw)
#' }
#' @family AmigaBasic.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBasic <- function(x, ...) {
  cursor <- 3
  result <- list()
  attr(result, "basic_header") <- x[1:2] ## Seems to be an identifier for basic scripts
  codelines <- T
  while (cursor < length(x)) {
    prev <- cursor
    cursor <- cursor + .rawToAmigaInt(x[cursor], 8, F)
    if (cursor == prev) {
      ## encountered a terminator. From here on no more code lines
      codelines <- F
      cursor <- cursor + 1
      if (x[cursor] == raw(1)) {
        cursor <- cursor + 1
      }
      if ((cursor %% 2) == 1) {
        if (x[cursor] != raw(1)) warning("Non-zero padding data encountered")
        cursor <- cursor + 1
      }
    } else {
      r <- x[(prev + 1):(cursor - ifelse(codelines, 1, 0))]
      if (codelines) {
        result[[length(result) + 1]] <- r
      } else {
        attr(result, "basic_names") <- c(attr(result, "basic_names"), rawToChar(r))
        cursor <- cursor + 1
      }
    }
  }
  class(result) <- "AmigaBasic"
  return(result)
}

#' @rdname as.raw
#' @name as.raw.AmigaBasic
#' @export
as.raw.AmigaBasic <- function(x, ...) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  nms <- attr(x, "basic_names")
  basic_header <- attr(x, "basic_header")
  lngths <- 1 + unlist(lapply(x, length))
  if (any(lngths > 255)) stop(sprintf("Lines %s are to long to encode", paste(which(lngths > 255), collapse = ", ")))
  x <- lapply(1:length(x), function(i) c(as.raw(lngths[i]), x[i]))
  c(basic_header,                            ## file header
    unlist(x),                               ## code lines
    raw(2 + ((length(unlist(x)) + 1) %% 2)), ## terminator and padding
    if(length(nms) > 0) {                    ## append variable/label/etc. names
      unlist(lapply(1:length(nms), function(i) c(as.raw(nchar(nms[i])), charToRaw(nms[i]))))
    } else {
      raw()
    })
}

#' Read Amiga Basic files
#'
#' Read an \code{\link{AmigaBasic}} script from its binary format.
#'
#' Normally Amiga Basic code is stored encoded in a binary format
#' (\code{\link{rawToAmigaBasic}}).
#' This function reads the binary data from a file (which can be
#' stored on a virtual disk (\code{\link[adfExplorer]{amigaDisk}}))
#' and converts in into an \code{\link{AmigaBasic}} class objec.
#' @rdname read.AmigaBasic
#' @name read.AmigaBasic
#' @param file A \code{character} string of the filename of the Amiga Basic file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @param ... Currently ignored
#' @return Returns an \code{\link{AmigaBasic}} class object read from the \code{file}.
#' @examples
#' \dontrun{
#' ## First create an AmigaBasic file
#' write.AmigaBasic(as.AmigaBasic("PRINT \"Hello world\""),
#'                  file.path(tempdir(), "helloworld.bas"))
#' 
#' ## Now let's read the same file:
#' bas <- read.AmigaBasic(file.path(tempdir(), "helloworld.bas"))
#' }
#' @family AmigaBasic.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaBasic <- function(file, disk = NULL, ...) {
  dat <- .read.generic(file, disk)
  rawToAmigaBasic(dat, ...)
}

#' Write an AmigaBasic object to a file
#'
#' Write an \code{\link{AmigaBasic}} class object to a file in its binary format.
#'
#' This function encodes the Amiga Basic code in its binary format
#' (using \code{\link{as.raw}}) and writes it to a file. The file
#' can also be stored onto a virtual Amiga disk
#' (\code{\link[adfExplorer]{amigaDisk}}).
#' 
#' @rdname write.AmigaBasic
#' @name write.AmigaBasic
#' @param x The \code{\link{AmigaBasic}} class object that needs to be
#' stored.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBasic}} object) needs to be written.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Invisibly returns the result of the call of \code{close} to the
#' file connection. Or, when \code{disk} is specified, a copy of
#' \code{disk} is returned to which the file(s) is/are written.
#' @examples
#' \dontrun{
#' ## First create an AmigaBasic object:
#' bas <- as.AmigaBasic("PRINT \"hello world!\"")
#' 
#' ## write to tempdir:
#' write.AmigaBasic(bas, file.path(tempdir(), "helloworld.bas"))
#' }
#' @family AmigaBasic.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.AmigaBasic <- function(x, file, disk = NULL) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  .write.generic(x, file, disk)
}

#' Coerce an AmigaBasic class object to its character representation
#'
#' XXX
#'
#' Amiga Basic files are encoded in a binary format. XXX
#' 
#' @rdname as.character
#' @name as.character
#' @param x An \code{\link{AmigaBasic}} class object that needs to be
#' coerced to its \code{character} representation.
#' @param ... Currently ignored.
#' @return A \code{vector} of \code{character} strings, where
#' each element of the \code{vector} is a \code{character} representation
#' of a line of Amiga Basic code stored in \code{x}.
#' @examples
#' \dontrun{
#' ## First create an Amiga Basic object:
#' bas <- as.AmigaBasic("PRINT \"Hello world!\"")
#' 
#' ## now convert the object back into text:
#' bas.txt <- as.character(bas)
#' }
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
as.character.AmigaBasic <- function(x, ...) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  nms <- attr(x, "basic_names")
  class(x) <- NULL
  i <- 0 #XXX
  code8 <- data.frame(where = character(), byte1=raw(), byte2 = raw(), byte3 = raw(), value = numeric()) #XXX
  x <- lapply(x, function(ln) {
    cmdln <- strrep(" ", .rawToAmigaInt(ln[1], 8, F))
    ln <- ln[-1]
    i <<- i + 1 ###XXX
    while (length(ln) > 0) {
      # if (i >= 4) browser()##XXX
      ln1 <- ln[1]
      ln <- ln[-1]
      #if (grepl("ELSEIF", cmdln)) browser()#XXX
      m1 <- which(.amigabasic_commands$code1 %in% ln1)
      if ((length(m1) > 1 && length(ln) > 0) ||
          (length(m1) == 1 && .amigabasic_commands$code2[m1] != raw(1))) {
        if (ln1 == as.raw(0xaf) && ln[1] != as.raw(0xe8)) { ## if command is "REM"
          m1 <- which(.amigabasic_commands$code1 %in% ln1 &
                        .amigabasic_commands$code2 %in% raw(1))
        } else {
          m1 <- which(.amigabasic_commands$code1 %in% ln1 &
                        .amigabasic_commands$code2 %in% ln[1])
          ln <- ln[-1]
        }
      }
      if (length(m1) == 0) {
        if (ln1 %in% as.raw(0x01:0x03)) {
          ## 0x01: variable, function or static sub names. Or an option name, such as 'P' in {SAVE "helloworld.bas",P}
          ## 0x02: label definition (using colon)
          ## 0x03: label reference
          if (ln1 == as.raw(0x03)) {
            if (ln[1] != raw(1)) warning(sprintf("Encountered non-zero padding byte (%02x).", as.numeric(ln[1])))
            ln <- ln[-1]
          }
          ## XXX check if this long is signed.
          ## XXX this will only cause problems at very large numbers...
          idx <- readBin(ln[1:2], "integer", size = 2, endian = "big") + 1
          cmdln <- paste0(cmdln, nms[idx])
          ln <- ln[-1:-2]
        } else if(ln1 == as.raw(0x08)) {
          #if (grepl("ELSEIF", cmdln)) browser()#XXX
          ## XXX what the heck is this? ignore ln[1] bytes of code
          ## XXX this always occurs after if then statement, or within elseif line
          code8 <<- rbind(code8,
                          data.frame(where = cmdln,
                                     byte1 = ln[1],
                                     byte2 = ln[2],
                                     byte3 = ln[3],
                                     value = readBin(ln[2:3], "integer", size = 2, signed = F, endian = "big"))) #XXX remove code8 later on
          ln <- ln[-1:-3] ## byte1 says nothing about the length of the data. byte1 is mostly 2, but can have different values. the amount of data is always fixed
        } else if (ln1 %in% as.raw(0x11:0x1a)) {
          cmdln <- paste0(cmdln, as.numeric(ln1) - 0x11)
        } else if (ln1 == as.raw(0x0f)) {
          cmdln <- paste0(cmdln, as.numeric(ln[1]))
          ln <- ln[-1]
        } else if (ln1 == as.raw(0x1c)) { ## short signed integer
          cmdln <- paste0(cmdln, readBin(ln[1:2], "integer", size = 2, endian = "big", signed = T))
          ln <- ln[-1:-2]
        } else if (ln1 == as.raw(0x1d)) { ## single precision float
          num <- readBin(ln[1:4], "numeric", size = 4, endian = "big")
          numopt1 <- gsub("0[.]", ".", toupper(format(num, digits = 7, scientific = F)))
          if (nchar(numopt1) > 8) {
            num <- toupper(format(num, digits = 7, scientific = T))
          } else num <- numopt1
          if (!grepl("[.]|E", num)) num <- paste0(num, "!")
          cmdln <- paste0(cmdln, num)
          ln <- ln[-1:-4]
        } else if (ln1 == as.raw(0x1e)) { ## long signed integer
          cmdln <- paste0(cmdln, readBin(ln[1:4], "integer", size = 4, endian = "big"), "&")
          ln <- ln[-1:-4]
        } else if (ln1 == as.raw(0x1f)) { ## double precision float
          num <- readBin(ln[1:8], "numeric", size = 8, endian = "big")
          numopt1 <- gsub("0[.]", ".", toupper(format(num, digits = 16, scientific = F)))
          if (nchar(numopt1) > 17) {
            num <- toupper(format(num, digits = 16, scientific = T))
          } else num <- numopt1
          num <- gsub("E", "D", num)
          if (!grepl("D", num)) num <- paste0(num, "#")
          cmdln <- paste0(cmdln, num)
          ln <- ln[-1:-8]
        } else if (ln1 == as.raw(0x22)) {
          quo <- ln1
          ln1 <- raw(1)
          ## If command is double quote, keep reading text
          ## until the end of the line or until another
          ## double quote is encountered
          while (length(ln) > 0 && ln1 != as.raw(0x22)) {
            ln1 <- ln[1]
            ln <- ln[-1]
            quo <- c(quo, ln1)
          }
          quo <- quo[quo != raw(1)] ## remove any possible null characters in case a closing quote is missing
          cmdln <- paste0(cmdln, rawToChar(quo))
        } else {
          try({ ## These are probably white space, brackets, &, %, #, etc.
            cmdln <- paste0(cmdln, rawToChar(ln1))
          })
        }
      } else {
        cmd <- .amigabasic_commands$command[m1]
        if (cmd %in% c("ELSE", "REM", "'")) {
          ## remove preceeding ":" if it is there
          if (nchar(cmdln) > 0) {
            if (substr(cmdln, nchar(cmdln), nchar(cmdln)) == ":")
              cmdln <- substr(cmdln, 0, nchar(cmdln) - 1)
          }
          if (cmd %in% c("REM", "'")) {
            ## all remaining data on the line should be treated as text when the command is REM or '
            cmd <- paste0(cmd, rawToChar(ln))
            ln <- raw(0)
          }
        }
        cmdln <- paste0(cmdln, cmd)
      }
    }
    return(cmdln)
  })
  x <- unlist(x)
  if (!is.null(x)) attr(x, "code8") <- code8 #XXX
  return(x)
}


#' List Amiga Basic reserved words.
#'
#' Obtain a list of reserved Amiga Basic words. These words are not
#' allowed as names of variables or labels in Amiga Basic.
#'
#' This function will return a full list of reserved Amiga Basic
#' words. This list does not serve as a manual for basic (for
#' that purpose consult external resources). This list is meant to
#' consult when chosing label names in Amiga Basic code. These
#' reserved words are not allowed as names.
#' 
#' @rdname AmigaBasic.reserved
#' @name AmigaBasic.reserved
#' @return Returns a \code{vecor} of \code{character} strings of
#' reserved Amiga Basic words.
#' @examples
#' AmigaBasic.reserved()
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
AmigaBasic.reserved <- function() {
  sort(.amigabasic_commands$command)
}

#' Coerce raw or character data to an AmigaBasic class object
#'
#' XXX
#'
#' XXX
#' 
#' @rdname as.AmigaBasic
#' @name as.AmigaBasic
#' @param x \code{x} should be a \code{vector} of \code{raw} data or
#' \code{character} strings. When \code{x} is \code{raw} data, it
#' is interpreted as if it where from an Amiga Basic binary encoded file.
#' 
#' When \code{x} is a \code{vector} of \code{character} strings,
#' each element of the vector should represent one line of Basic code.
#' Each line should not contain line break or other special characters,
#' as this will result in errors. The text should represent valid
#' Amiga Basic syntax. The syntax is only checked to a limited extent as
#' this package does not implement an interpreter for the code.
#' @param ... Currently ignored.
#' @return Returns an \code{\link{AmigaBasic}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' ## An AmigaBasic object can be created from text.
#' ## Note that each line of code is a seperate element
#' ## in the vector:
#' bas <- as.AmigaBasic(c(
#'   "CLS ' Clear the screen",
#'   "PRINT \"Hello world!\" ' Print a message on the screen"
#' ))
#' 
#' ## Let's make it raw data:
#' bas.raw <- as.raw(bas)
#' 
#' ## We can also use the raw data to create an Amiga Basic object:
#' ## Note that this effectively the same as calling 'rawToAmigaBasic'
#' bas <- as.AmigaBasic(bas.raw)
#' }
#' @family AmigaBasic.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
as.AmigaBasic <- function(x, ...) {
  if (inherits(x, "character")) {
    ## remember the number of leading spaces, then remove them
    leading.spaces <- unlist(lapply(gregexpr("^ +", x), function(y) attr(y, "match.length")))
    leading.spaces[leading.spaces < 0] <- 0
    leading.spaces <- as.raw(leading.spaces)
    x <- gsub("^ +", "", x)
    
    ## split the lines at special characters:
    x <- strsplit(x, "(?=[ !#$%\\^&*()\\-+=/?,<>:;\"'])", perl = T)
    
    nms <- NULL
    x <- lapply(x, function(y) {
      result <- raw(0)
      while (length(y) > 0) {
        ## This cannot be moved outside the loop, as the extra split modifies the text.
        ## And that should not happen when text is between double quotes or follows a REM statement
        if (grepl("^[0-9.]", y[1], perl = T)) {
          exponent <- gregexpr("[d-e][-]?\\d", tolower(paste(y, collapse = "")), perl = T)[[1]]
          ## Find first non-numeric/period that is not !#%&:
          extrasplit <- gregexpr("[^\\d.!#%&]", y[1], perl = T)[[1]]
          ## except for when it is an exponent
          extrasplit <- extrasplit[extrasplit != exponent[1]]
          ## Find first occurence of !#%&
          extrasplit <- c(extrasplit, 1 + gregexpr("[!#%&]", y[1], perl = T)[[1]])
          ## second occurence of period:
          extrasplit <- c(extrasplit, gregexpr("(?:.*?\\K[.]){2}", y[1], perl = T)[[1]])
          extrasplit[extrasplit < 1] <- 1 + nchar(y[1])
          extrasplit <- extrasplit[extrasplit == min(extrasplit)]
          if (extrasplit > 1 && extrasplit <= nchar(y[1])) {
            extrasplit <- c(substr(y[1], 1, extrasplit - 1), " ", substr(y[1], extrasplit, nchar(y[1])))
          } else extrasplit <- y[1]
          extrasplit[grepl("^[.]", extrasplit)] <- paste0("0", extrasplit[grepl("^[.]", extrasplit)])
          y <- c(extrasplit,
                 y[-1])
          y <- y[nchar(y) > 0]
        }
        if (y[1] == "?") y[1] <- "PRINT "
        cmd <- match(toupper(y[1]), .amigabasic_commands$command)
        if (!is.na(cmd)) {
          cmd <- unlist(.amigabasic_commands[cmd, c("code1", "code2")], use.names = F)
          cmd <- cmd[cmd != raw(1)]
          if (toupper(y[1]) %in% c("'", "REM")) {
            cmd <- c(if(length(result) > 0) charToRaw(":"),
                     cmd,
                     charToRaw(paste(y[-1], collapse = "")))
            y <- character()
          }
          result <- c(result, cmd)
        } else if (grepl("[ !#$%&()/,:;]", y[1])) {
          result <- c(result, charToRaw(y[1]))
        } else if (y[1] =="\"") { ## if it's between double quotes, it's a string
          repeat {
            result <- c(result, charToRaw(y[1]))
            y <- y[-1]
            if (length(y) == 0 || y[1] == "\"") break
          }
          if (length(y)> 0) result <- c(result, charToRaw(y[1]))
        } else if (suppressWarnings(!is.na(as.numeric(gsub("D", "E", toupper(y[1])))))) { ## if it is a numeric
          dot <- grepl("[.]", y[1])
          tp <- NA
          if (grepl("D", toupper(y[[1]]))) tp = "#" ## if The exponent is noted with a 'D', it is a double precision float
          if (grepl("E", toupper(y[[1]]))) tp = "!" ## if The exponent is noted with a 'E', it is a single precision float
          if (grepl("E|D", toupper(y[[1]])) && (length(y) > 1) && y[2] %in% c("+", "-")) {
            y[1] <- paste0(y[1], y[2])
            y <- y[-2]

            if (length(y) > 1 && grepl("^[0-9]", y[2])) {
              nonnum <- gregexpr("\\D", y[2], perl = T)[[1]][1]
              if (nonnum != -1) {
                y <- c(y[1],
                       substr(y[2], 1, nonnum - 1), " ",
                       substr(y[2], nonnum, nchar(y[2])),
                       y[-1:-2])
              }
              y[1] <- paste0(y[1], y[2])
              y <- y[-2]
            }
          }
          num <- as.numeric(gsub("D", "E", toupper(as.character(y[1]))))
          if (is.na(tp) && length(y) > 1 && y[2] %in% c("!", "#", "%", "&")) {
            if (length(y) > 2 && !grepl("^[ ()/,:;-]", y[3])) y <- c(y[1:2], " ", y[-1:-2])
            tp <- y[2]
            y <- y[-1]
          }
          ## if the number type is not specified, guess...
          if (is.na(tp)) {
            if (round(num) == num && !dot) {
              ## it's either a short or a long integer
              tp <- ifelse(num <= 32767 && num >= -32767, "%", "&")
            } else {
              ## it's either a single or a double float
              fm <- format(num, scientific = T)
              fm <- strsplit(fm, "e")[[1]]
              ## XXX need to test this more extensively
              tp <- ifelse(nchar(fm[1]) > 8 || abs(as.numeric(fm[2])) > 38,
                           "#", "!")
            }
          }
          if (tp == "%" && num < 255) {
            if (num < 10) {
              result <- c(
                result,
                as.raw(num + 0x11)
              )
            } else {
              result <- c(
                result,
                as.raw(c(0x0f, num))
              )
            }
          } else {
            result <- c(
              result,
              as.raw(c(0x1c, 0x1e, 0x1d, 0x1f))[match(tp, c("%", "&", "!", "#"))],
              writeBin(
                ifelse(tp %in% c("%", "&"), as.integer(num), num),
                raw(),
                size = ifelse(tp == "%", 2, ifelse(tp == "#", 8, 4)),
                endian = "big")
            )
          }
        } else { ## If it's none of the above, it must be a name/label
          ## If a name starts with FN, it should be a function defined with "DEF FN"
          if (grepl("^FN", toupper(y[1]))) {
            result <- c(result, as.raw(0x93))
            y[1] <- substr(y[1], 3, nchar(y[1]))
          }
          
          if (any(as.logical(check.AmigaBasic.names(y[1])))) stop(sprintf("Fatal syntax error Basic code at: '%s'", paste(y, collapse = "")))
          
          ## Check if the name was already used before, otherwise
          ## append it to the vector of names.
          nm <- match(toupper(y[1]), toupper(nms)) ## the first definition will deterimine the case of the characters in the name
          
          cd <- 1 ## variable or static sub label
          ## if the name is directly followed by a colon, it is a label definition
          if (length(y) > 1 && y[2] == ":") {
            cd <- 2 ## label definition
          } else {
            goto.gosub <- c(
              which(result %in% as.raw(c(0x96, 0x97))), # GOSUB and GOTO
              grepRaw(as.raw(c(0xf8, 0xa8)), result, fixed = T) # RESTORE
            )
            if (length(goto.gosub) > 0) {
              goto.gosub <- max(goto.gosub)
              ## assume the label is meant to go with the GOSUB, GOTO or RESTORE command when the line is not split by ":" in between
              if (goto.gosub == length(result) || all(result[(goto.gosub + 1):length(result)] != as.raw(0x3a))) cd <- 3
            }
          }
          if (is.na(nm)) {
            nms <<- c(nms, y[1])
            nm <- length(nms)
          }
          result <- c(
            result,
            as.raw(cd),
            writeBin(as.integer(nm - 1), raw(), size = 2, endian = "big")
          )
        }
        y <- y[-1]
      }
      return (result)
    })
    x <- lapply(1:length(x), function(i) {
      c(leading.spaces[[i]], x[[i]], raw(2))
    })
    ## XXX protected basic files have header 0xf4, 0xc2. U suspect that this is followed by a 5 byte encryption key, followed by encrypted data
    attr(x, "basic_header") <- as.raw(c(0xf5, 0x00))
    attr(x, "basic_names") <- nms
    class(x) <- "AmigaBasic"
    return (x)
  } else if (inherits(x, "raw")) {
    return (rawToAmigaBasic(x, ...))
  } else {
    print("Cannot convert 'x' to S3 AmigaBasic class object.")
  }
}

#' Extract or replace lines of Amiga Basic code
#'
#' XXX
#'
#' XXX
#' 
#' @rdname ExtractBasic
#' @name [.AmigaBasic
#' @param x An \code{AmigaBasic} class object from which specific lines
#' need to be extracted or replaced.
#' @param i A \code{vector} of \code{numeric} indices. This index
#' is used to select specific lines. Negative values will deselect lines.
#' @param value A \code{vector} of \code{character} strings or an
#' \code{\link{AmigaBasic}} class object that is used to replace
#' the selected indices \code{i}. \code{value} should represent the
#' same number of lines of code as the selected number of lines.
#' @return XXX
#' @examples
#' \dontrun{
#' ## First generate a few lines of Basic code:
#' bas <- as.AmigaBasic(c(
#'   "LET a = 1",
#'   "a = a + 1",
#'   "PRINT \"a now equals\";a",
#'   "INPUT \"clear screen (y/n)? \", b$",
#'   "IF UCASE$(b$) = \"Y\" THEN CLS"
#' ))
#' 
#' ## Select only lines 4 and 5:
#' bas[4:5]
#' 
#' ## use negative indices to deselect specific lines.
#' ## deselect line 2:
#' bas[-2]
#' 
#' ## replace line 2
#' bas[2] <- "a = a + 2"
#' 
#' ## You can also use AmigaBasic class object as replacement
#' bas[2] <- as.AmigaBasic("a = a + 3")
#' }
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
`[.AmigaBasic` <- function(x, i) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  a <- attributes(x)
  cl <- class(x)
  class(x) <- NULL
  x <- x[i]
  attributes(x) <- a
  class(x) <- cl
  x
}

#' @rdname ExtractBasic
#' @name [<-.AmigaBasic
#' @export
`[<-.AmigaBasic` <- function(x, i, value) {
  if (!inherits(x, "AmigaBasic")) stop("'x' should be of class AmigaBasic.")
  value <- as.character(value)
  x <- as.character(x)
  x[i] <- value
  return(as.AmigaBasic(x))
}

## XXX extract function is called recursively, even without selecting an element
# `[[.AmigaBasic` <- function(x, i) {
#   browser()#XXX
#   if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
#   sink();print("Maar ik kom ook hier [[")#XXX
#   a <- attributes(x)
#   cl <- class(x)
#   class(x) <- NULL
#   x <- list(x[[i]])
#   attributes(x) <- a
#   class(x) <- cl
#   x
# }

#' @export
`[[<-.AmigaBasic` <- function(x, i, value) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  browser()#XXX
  x
}

#' @export
print.AmigaBasic <- function(x, ...) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  cat(paste(as.character(x), collapse = "\n"))
}

#' Extract or replace variable and label names from Amiga Basic scripts
#'
#' In the binary Amiga Basic files, names for labels and variables
#' in the code are stored at the end of the file. In the encoded
#' there is only a reference to the index of the name in that list. Use
#' this function to list, select or replace names included in the code
#' 
#' XXX
#' 
#' @rdname AmigaBasic.names
#' @name AmigaBasic.names
#' @param x XXX
#' @param value XXX
#' @return XXX
#' @examples
#' ## Let's create some Basic code with labels and variables:
#' bas <- as.AmigaBasic(c(
#'   "REM - This will loop forever...",
#'   "my.label:",
#'   "  my.variable% = 0",
#'   "  WHILE my.variable% < 10",
#'   "    my.variable% = my.variable% + 1",
#'   "  WEND",
#'   "  GOTO my.label"
#' ))
#' 
#' ## list the names in the script above:
#' AmigaBasic.names(bas)
#' 
#' ## change the first name:
#' AmigaBasic.names(bas)[1] <- "better.label"
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
AmigaBasic.names <- function(x) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  result <- attr(x, "basic_names")
  if (length(result) == 0) character() else result
}

#' @rdname AmigaBasic.names
#' @name AmigaBasic.names<-
#' @export
`AmigaBasic.names<-` <- function(x, value) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  if (!is.character(value) || length(value) != length(attr(x, "basic_names")))
    stop("Replacement should be a vector of characters of the same length.")
  if (any(duplicated(toupper(value)))) stop("All names should be unique.")
  if (any(nchar(value) > 255 | nchar(value) < 1)) stop("All names should be more than 1 and less than 255 characters long.")
  if (any(toupper(value) %in% .amigabasic_commands$command)) stop("Names cannot be reserved AmigaBasic words.")
  if (any(grepl("[^a-zA-Z0-9.]", value, perl = TRUE))) stop("Names should consist of alphanumerics or periods.")
  if (any(grepl("[0-9]", substr(value, 1, 1), perl = T))) stop("Names should not start with numeric characters.")
  if (any(grepl("[ !#$%\\^&*()\\-+=/?,<>:;\"']", value, perl = T))) stop("Names should not contain special characters.")
  attr(x, "basic_names") <- value
  x
}

#' @export
check.AmigaBasic.names <- function(x, ...) {
  data.frame(
    duplicated = duplicated(toupper(x)),
    length     = nchar(x) > 255 | nchar(x) < 1,
    reserved   = toupper(x) %in% .amigabasic_commands$command,
    characters = grepl("[^a-zA-Z0-9.]", x, perl = TRUE),
    start      = grepl("^[0-9.]", x, perl = TRUE)
  )
}


#' @rdname c
#' @name c
#' @export
c.AmigaBasic <- function(...) {
  bas.codes <- list(...)
  if (!all(unlist(lapply(bas.codes, inherits, what = "AmigaBasic")))) stop ("All arguments should be of type AmigaBasic.")
  ## It's not fastest, but it is safest to convert all codes to text
  ## than back to AmigaBasic.
  bas.codes <- as.AmigaBasic(unlist(lapply(bas.codes, as.character)))
}

.basic.shape.header <- data.frame(
  byte      = c(4, 4, 4, 4, 4, -2, 2, -2),
  signed    = c(F, F, F, F, F,  F, F,  F),
  par.names = c("colorset", "dataset", "depth", "width", "height", "flags", "planePick", "planeOnOff"),
  stringsAsFactors = F
)

#' XXX
#'
#' XXX
#'
#' XXX
#' 
#' @rdname rawToAmigaBasicShape
#' @name rawToAmigaBasicShape
#' @param x A \code{vector} of \code{raw} data that is to be converted
#' into an \code{\link{AmigaBasicShape}} class object.
#' @param ... Arguments passed on to \code{\link{bitmapToRaster}}.
#' As the palette is not stored in the file. This can be used to
#' pass on a custom \code{palette} to \code{\link{bitmapToRaster}}.
#' @return XXX
#' @examples
#' \dontrun{
#' ## XXX
#' }
#' @family AmigaBasicShapes.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBasicShape <- function(x, ...) {
  ## colorset and dataset seem to be ignored in the basic object editor.
  ## They seem not required for interpretation of the file
  result <- with(.basic.shape.header, AmigaFFH:::.read.amigaData(x, byte, signed, par.names))
  x <- x[-1:-sum(abs(.basic.shape.header$byte))]
  bm <- with(result, AmigaFFH:::bitmapToRaster(x, width, height, depth, interleaved = F, ...))
  sz <- 2*ceiling(result$width/16)*result$height*result$depth
  x <- x[-1:-sz]
  browser()
  result$flags <- rev(as.logical(AmigaFFH:::.rawToBitmap(result$flags, T, F)))
  names(result$flags) <- c("fVSprite", "collisionPlaneIncluded", "imageShadowIncluded", "saveBack", "overlay", "saveBob",
                           sprintf("reserved%02i", 1:10))
  if (result$flags["fVSprite"] && result$depth != 2) warning("Unexpected bitmap depth for sprite mode.")
  if (result$flag["imageShadowIncluded"]) {
    #TODO read collision bitmap plane
  }
  if (result$flag["collisionPlaneIncluded"]) {
    #TODO read collision bitmap plane
  }
  if (result$flags["fVSprite"]) {
    #TODO read 3x2 bytes (color of sprite)
  }
  if (length(x) > 0) warning("Unexpected and unused trailing data!")
  ## XXX for now, just return the bitmap
  bm
}

#' Read Amiga Basic Shape files
#'
#' XXX
#'
#' XXX
#' @rdname read.AmigaBasicShape
#' @name read.AmigaBasicShape
#' @param file A \code{character} string of the filename of the Amiga Basic Shape file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @param ... XXX
#' @return Returns an \code{\link{AmigaBasicShape}} class object read from the \code{file}.
#' @examples
#' \dontrun{
#' ## XXX
#' }
#' @family AmigaBasicShape.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaBasicShape <- function(file, disk = NULL, ...) {
  dat <- .read.generic(file, disk)
  rawToAmigaBasicShape(dat, ...)
}

#' Write an AmigaBasicShape object to a file
#'
#' Write an \code{\link{AmigaBasicShape}} class object to a file in its binary format.
#'
#' This function coerces the Amiga Basic Shape into its binary format
#' (using \code{\link{as.raw}}) and writes it to a file. The file
#' can also be stored onto a virtual Amiga disk
#' (\code{\link[adfExplorer]{amigaDisk}}).
#' 
#' @rdname write.AmigaBasicShape
#' @name write.AmigaBasicShape
#' @param x The \code{\link{AmigaBasicShape}} class object that needs to be
#' stored.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBasicShape}} object) needs to be written.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Invisibly returns the result of the call of \code{close} to the
#' file connection. Or, when \code{disk} is specified, a copy of
#' \code{disk} is returned to which the file(s) is/are written.
#' @examples
#' \dontrun{
#' ## XXX
#' }
#' @family AmigaBasicShape.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.AmigaBasicShape <- function(x, file, disk = NULL) {
  if (!inherits(x, "AmigaBasicShape")) stop("x should be of class AmigaBasicShape.")
  .write.generic(x, file, disk)
}

#' @rdname as.raw
#' @name as.raw.AmigaBasicShape
#' @export
as.raw.AmigaBasicShape <- function(x, ...) {
  ## XXX update documentation
  if (!inherits(x, "AmigaBasicShape")) stop("x should be of class AmigaBasicShape.")
  ## XXX
}
