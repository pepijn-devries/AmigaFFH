#' The S3 AmigaBasic class
#' 
#' A class that represents the content of Amiga Basic files.
#' 
#' Amiga Basic is a \href{https://en.wikipedia.org/wiki/BASIC}{BASIC}-style programming language that was shipped
#' with early Commodore Amiga machines. It requires an interpreter to run an Amiga Basic script. The AmigaFFH
#' package does not interpret Amiga Basic scripts. It does allow for encoding and decoding scripts in the binary
#' format in which it was originally stored on the Amiga. Amiga Basic scripts were stored as encoded binaries instead
#' of ASCII text files in order to save (at the time precious) memory and disk space.
#' 
#' Amiga Basic binary files start with a file header (as an identifier) and is followed by each line of the script
#' as binary data. The \code{AmigaBasic}-class object stores each line of the script as a \code{list} item as a \code{vector}
#' of \code{raw} data. Use \code{\link{as.character}} and \code{\link{as.AmigaBasic}} to switch between
#' \code{character} data and \code{AmigaBasic}-class objects.
#' 
#' @note Although there is ample reference material on the Amiga BASIC language, there is no documentation
#' available on the script file storage format. The implementation in the AmigaFFH package is all the result of
#' painstaking reverse engineering on my part. Consequently the Amiga Basic file encoders and decoders implemented
#' here may not be infallible.
#' @docType class
#' @name AmigaBasic
#' @rdname AmigaBasic
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @examples
#' \dontrun{
#' ## This creates an AmigaBasic-class object:
#' bas <- as.AmigaBasic("PRINT \"hello world!\"")
#' 
#' ## This will decode the object as plain text:
#' as.character(bas)
#' }
#' @references \url{https://en.wikipedia.org/wiki/AmigaBASIC}
NULL

#' The S3 AmigaBasicShape class
#' 
#' A class that represents the file format used by Amiga Basic to store bitmap graphics: blitter objects and sprites.
#' 
#' Amiga Basic used a specific format to store bitmap images that could be displayed using Basic code. Both
#' sprites and blitter objects can be stored and used. This class is used to represent such files.
#' 
#' @docType class
#' @name AmigaBasicShape
#' @rdname AmigaBasicShape
#' @family AmigaBasicShape.operations
#' @author Pepijn de Vries
#' @examples
#' \dontrun{
#' ball   <- read.AmigaBasicShape(system.file("ball.shp", package = "AmigaFFH"))
#' r_logo <- read.AmigaBasicShape(system.file("r_logo.shp", package = "AmigaFFH"))
#' 
#' plot(ball)
#' plot(r_logo)
#' }
NULL

.amigabasicshape.flags <- c("fVSprite", "collisionPlaneIncluded", "imageShadowIncluded", "saveBack", "overlay", "saveBob",
                            sprintf("reserved%02i", 1:10))

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
f8,bc,RESET
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

.valid_code <- function(x) {
  apply(
    .amigabasic_commands[,c("code1", "code2")],
    1, function(y) {
      y <- as.raw(paste0("0x", y))
      (length(x) > 0 && y[[1]] == x[[1]]) &&
        (y[[2]] == 0x00 || (length(x) > 1 && x[[2]] == y[[2]]))
    })
}

#' Coerce raw data into an AmigaBasic class object
#'
#' \code{\link{AmigaBasic}} objects are comprehensive representations of binary-encode Amiga Basic scripts.
#' Use this function to convert raw content from encoded Amiga Basic scripts to an \code{\link{AmigaBasic}}
#' object.
#'
#' This function will convert raw data as stored in Amiga Basic files into its corresponding S3
#' \code{\link{AmigaBasic}}-class object.
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
        ## This padding byte only occurs align files to word (2byte) size. It does not seem to have any other function.
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
  x <- unclass(x)
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
#' stored on a virtual disk (\code{\link[adfExplorer:amigaDisk-class]{amigaDisk}}))
#' and converts in into an \code{\link{AmigaBasic}} class objec.
#' @rdname read.AmigaBasic
#' @name read.AmigaBasic
#' @param file A \code{character} string of the filename of the Amiga Basic file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
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
#' 
#' ## There's also a demo file included with the package
#' demo.bas <- read.AmigaBasic(system.file("demo.bas", package = "AmigaFFH"))
#' demo.bas
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
#' (\code{\link[adfExplorer:amigaDisk-class]{amigaDisk}}).
#' 
#' @rdname write.AmigaBasic
#' @name write.AmigaBasic
#' @param x The \code{\link{AmigaBasic}} class object that needs to be
#' stored.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBasic}} object) needs to be written.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
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
#' Coerce an \code{\link{AmigaBasic}}-class object to its character representation
#'
#' Amiga Basic files are encoded in a binary format and are also stored as such
#' in \code{\link{AmigaBasic}}-class objects. Use this function to convert
#' these objects into legible \code{character} data.
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
  x <- lapply(x, function(ln) {
    cmdln <- strrep(" ", .rawToAmigaInt(ln[1], 8 , F))
    ln <- ln[-1]
    while (length(ln) > 2) {
      ln1 <- ln[1]
      ln <- ln[-1]
      
      if (cmdln == "" && length(ln) > 1) { ## Check if line starts with numeric label
        if (((utils::tail(ln, 1) & as.raw(0x80)) != 0x00) || !(any(.valid_code(c(ln1, ln[1]))) || ln1 %in% as.raw(1:3)) &&
            (length(ln) < 3 || (any(.valid_code(ln[2:3])) || ln[2] %in% as.raw(1:3)))) {
          ln[length(ln)] <- xor(utils::tail(ln, 1), as.raw(0x80))
          cmdln <- paste0(cmdln, as.character(readBin(c(ln1, ln[1]), "integer", 2, 2, F, "big")))
          ln1 <- ln[2]
          ln <- ln[-1:-2]
          if (length(ln) > 2) cmdln <- paste0(cmdln, " ")
        }
      }
      #TODO 'THEN' and 'ELSEIF' seems to be followed by redundant binaries
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
          idx <- .rawToAmigaInt(ln[1:2], 16, F) + 1
          cmdln <- paste0(cmdln, nms[idx])
          ln <- ln[-1:-2]
        } else if(ln1 == as.raw(0x08)) {
          ## TODO This code appears to be redundant and follows after THEN or ELSEIF statements
          ln <- ln[-1:-3] ## byte1 says nothing about the length of the data. byte1 is mostly 2, but can have different values. the amount of data is always fixed. There seems to be a correlation with the line number in which this occurs
        } else if (ln1 %in% as.raw(0x11:0x1a)) {
          cmdln <- paste0(cmdln, as.numeric(ln1) - 0x11)
        } else if (ln1 == as.raw(0x0b)) { ## octal number
          cmdln <- sprintf("%s&O%o", cmdln, readBin(ln[1:2], "integer", size = 2, endian = "big", signed = F))
          ln <- ln[-1:-2]
        } else if (ln1 == as.raw(0x0c)) { ## hexadecimal short signed integer
          cmdln <- sprintf("%s&H%X", cmdln, readBin(ln[1:2], "integer", size = 2, endian = "big", signed = F))
          ln <- ln[-1:-2]
        } else if (ln1 == as.raw(0x0e)) { ## longish unsigned integer
          cmdln <- paste0(cmdln, readBin(c(raw(1), ln[1:3]), "integer", size = 4, endian = "big"))
          ln <- ln[-1:-3]
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
          ## remove preceding ":" if it is there
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
#' consult when choosing label names in Amiga Basic code. These
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
#' Coerce raw or character data to an \code{\link{AmigaBasic}} S3 class object
#'
#' Convert text to an \code{\link{AmigaBasic}} S3 class object. The text should
#' consist of valid Amiga BASIC syntaxis. This function does not perform a
#' full check of the syntaxis, but will break on some fundamental syntaxis malformations
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
#' @references \url{https://en.wikipedia.org/wiki/AmigaBASIC}
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
      trailing_marker <- raw(1)
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
          ## Find first occurrence of !#%&
          extrasplit <- c(extrasplit, 1 + gregexpr("[!#%&]", y[1], perl = T)[[1]])
          ## second occurrence of period:
          extrasplit <- c(extrasplit, gregexpr("(?:.*?\\K[.]){2}", y[1], perl = T)[[1]])
          extrasplit[extrasplit < 1] <- 1 + nchar(y[1])
          extrasplit <- extrasplit[which(extrasplit == min(extrasplit))[1]]
          if (extrasplit > 1 && extrasplit <= nchar(y[1])) {
            extrasplit <- c(substr(y[1], 1, extrasplit - 1), " ", substr(y[1], extrasplit, nchar(y[1])))
          } else extrasplit <- y[1]
          extrasplit[grepl("^[.]", extrasplit)] <- paste0("0", extrasplit[grepl("^[.]", extrasplit)])
          y <- c(extrasplit,
                 y[-1])
          y <- y[nchar(y) > 0]
        }
        if (y[1] == "?") y <- c("PRINT", " ", y[-1])
        cmd <- match(toupper(y[1]), .amigabasic_commands$command)
        if (is.na(cmd) && length(y) > 1 && y[2] == "$") {
          cmd <- match(toupper(paste0(y[1], y[2])), .amigabasic_commands$command)
          if (!is.na(cmd)) {
            y[1] <- paste0(y[1], y[2])
            y <- y[-2]
          }
        }
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
        } else if (suppressWarnings(!is.na(as.numeric(gsub("^O", "", gsub("^H", "0x", gsub("D", "E", toupper(y[1])))))))) { ## if it is a numeric
          if (any(y[y != "  "][1] == "$")) stop("Fatal syntax error, numeric cannot be followed by '$'.")
          dot <- grepl("[.]", y[1])
          tp <- NA
          if (grepl("D", toupper(y[[1]]))) tp <- "#" ## if The exponent is noted with a 'D', it is a double precision float
          if (grepl("E", toupper(y[[1]]))) tp <- "!" ## if The exponent is noted with a 'E', it is a single precision float
          if (startsWith(toupper(y[1]), "H") && length(result) > 0 && utils::tail(result, 1) == as.raw(0x26)) tp <- "&H" ## a signed hexadecimal short int (two bytes)
          if (startsWith(toupper(y[1]), "O") || (is.na(tp) && length(result) > 0 && utils::tail(result, 1) == as.raw(0x26))) tp <- "&O" ## an octal number
          if (length(result) == 0) tp <- "numeric_label"
          if (length(result) > 2) {
            comm_check <- .amigabasic_commands$command %in% c("GOTO", "GOSUB", "BREAK", "COLLISION", "ERROR", "MENU", "MOUSE", "TIMER")
            comm_check <- .amigabasic_commands[comm_check, c("code1", "code2")]
            check <- any(apply(comm_check, 1, function(cc) {
              code <- as.raw(paste0("0x", cc))
              code <- if (!any(code == 0x00)) c(as.raw(c(0xaa, 0x20)), code) else code[code != 0x00] # 0xaa20 == "on "
              code <- c(code, as.raw(0x20))
              return(length(result) >= length(code) && identical(code, utils::tail(result, length(code))))
            }))
            if (check) tp <- "numeric_ref"
          }
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
          num <- suppressWarnings(as.numeric(gsub("D", "E", toupper(as.character(y[1])))))
          if (is.na(num) && !(tp %in% c("&H", "&O"))) stop(sprintf("Fatal error in number format: %s%s.", tp, y[1]))
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
              ## TODO need to test this more extensively
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
          } else if(tp %in% c("numeric_label", "numeric_ref")) {
            if (as.integer(y[1]) < 0) stop("Fatal syntax error: numeric labels cannot be negative.")
            if (as.integer(y[1]) > 0xfff9) {
              y <- c(substr(y[1], 1, 4), substr(y[1], 5, nchar(y[1])), y[-1])
            }
            result <- c(result, if(tp == "numeric_ref") as.raw(c(0x0e, 0x00)) else raw(),
                        writeBin(as.integer(y[1]), raw(), 2, "big"))
            z <- y[y != " "]
            valid_current <- any(.valid_code(result)) || result[[1]] %in% as.raw(1:3)
            valid_next    <- length(z) > 1 &&
              (toupper(z[2]) %in% .amigabasic_commands$command || !is.na(suppressWarnings(try(as.numeric(z[2])))) ||
                                              (!is.null(nms) && toupper(z[2]) %in% toupper(nms)) ||
                 grepl("[ !#$%\\^&*()\\-+=/?,<>:;\"']", z[2], perl = T))
            ## TODO This check isn't 100% similar to original amigabasic. It will produce a different outcome on as.AmigaBasic("1 PRINT \"ja\":GOTO 9"), but it will still work
            if (length(result) == 2 &&
                (!(valid_current && (length(y) == 1 || valid_next)) || valid_current)) trailing_marker <- as.raw(0x80)
            if (length(y) > 2 && y[2] == " ") y <- y[-2]
          } else {
            if (tp %in% c("&H", "&O")) result <- utils::head(result, -1) ## remove the '&' previously written to the result
            result <- c(
              result,
              as.raw(c(0x0b, 0x0c, 0x1c, 0x1d, 0x1e, 0x1f))[match(tp, c("&O", "&H", "%", "!", "&", "#"))],
              writeBin(
                ifelse(tp %in% c("%", "&"), as.integer(num),
                       ifelse(tp == "&H", { ## two byte hexadecimal
                         num <- as.integer(gsub("^H", "0x", y[1]))
                         as.integer(ifelse(num > 0x8000, num - 0x10000, num))
                       },
                       ifelse(tp =="&O", { ## two byte octal
                         y[1] <- gsub("^O", "", toupper(y[1]))
                         num <- strsplit(y[1], "[^0-7]")[[1]]
                         y_insert <- substr(y[1], nchar(num) + 1, nchar(y[1]))
                         if (y_insert != "") y <- c(y[1], " ", y_insert, utils::tail(y, -1))
                         num <- as.numeric(strsplit(num, "")[[1]])
                         num <- as.integer(sum(num * (8 ^ (rev(seq_along(num) - 1)))))
                         as.integer(ifelse(num > 0x8000, num - 0x10000, num))
                       }, num))),
                raw(),
                size = ifelse(tp %in% c("&H", "&O", "%"), 2, ifelse(tp == "#", 8, 4)),
                endian = "big")
            )
          }
        } else { ## If it's none of the above, it must be a name/label
          ## If a name starts with FN, it should be a function defined with "DEF FN"
          if (grepl("^FN", toupper(y[1]))) {
            result <- c(result, as.raw(0x93))
            y[1] <- substr(y[1], 3, nchar(y[1]))
          }
          
          if (any(as.logical(check.names.AmigaBasic(y[1])))) stop(sprintf("Fatal syntax error Basic code at: '%s'", paste(y, collapse = "")))
          
          ## Check if the name was already used before, otherwise
          ## append it to the vector of names.
          nm <- match(toupper(y[1]), toupper(nms)) ## the first definition will determine the case of the characters in the name
          
          cd <- 1 ## variable or static sub label
          ## if the name is directly followed by a colon (and nothing else), it is a label definition
          if (length(y) == 2 && y[2] == ":") {
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
            as.raw(cd), if(cd == 3) raw(1) else raw(),
            writeBin(as.integer(nm - 1), raw(), size = 2, endian = "big")
          )
        }
        y <- y[-1]
      }
      attributes(result)$trailing_marker <- trailing_marker
      return (result)
    })
    x <- lapply(seq_along(x), function(i) {
      padding <- raw(2)
      if (!is.null(attr(x[[i]], "trailing_marker"))) {
        padding[2] <- attr(x[[i]], "trailing_marker")
        attr(x[[i]], "trailing_marker") <- NULL
      }
      c(leading.spaces[[i]], x[[i]], padding)
    })
    x <- lapply(x, function(y) {
      if (length(y) > 2 && identical(y[1:3], as.raw(c(0x00, 0xaf, 0xe8)))) y <- c(y[1], as.raw(0x3a), y[-1])
      y
    })
    ## TODO protected basic files have header 0xf4, 0xc2. U suspect that this is followed by a 5 byte encryption key, followed by encrypted data
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
#' Extract or replace lines of Amiga Basic code
#'
#' Extract or replace specific lines in an \code{\link{AmigaBasic}}-class object.
#' 
#' @rdname ExtractBasic
#' @name [.AmigaBasic
#' @param x An \code{AmigaBasic} class object from which specific lines
#' need to be extracted or replaced.
#' @param i In case of `[[', an integer index, representing the line-number of basic code to be selected.
#' In case of `[': a \code{vector} of \code{numeric} indices. This index
#' is used to select specific lines. Negative values will deselect lines.
#' @param value A \code{vector} of \code{character} strings or an
#' \code{\link{AmigaBasic}} class object that is used to replace
#' the selected indices \code{i}. \code{value} should represent the
#' same number of lines of code as the selected number of lines.
#' @return The extraction method returns an \code{\link{AmigaBasic}} object based in the lines selected with \code{i}.
#' The replacement method returns an \code{\link{AmigaBasic}} object with the selected lines replaced with \code{value}.
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
#' 
#' ## single lines can also be selected with '[['
#' bas[[2]]
#' }
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
`[.AmigaBasic` <- function(x, i) {
  vctrs::vec_restore(NextMethod(), x)
}

#' @rdname ExtractBasic
#' @name [<-.AmigaBasic
#' @export
`[<-.AmigaBasic` <- function(x, i, value) {
  if (!inherits(x, "AmigaBasic")) stop("'x' should be of class AmigaBasic.")
  x <- as.character(x)
  x[i] <- as.character(value)
  return(as.AmigaBasic(x))
}

#' @rdname ExtractBasic
#' @name `[[.AmigaBasic`
#' @export
`[[.AmigaBasic` <- function(x, i) {
  vctrs::vec_restore(list(NextMethod()), x)
}

#' @rdname ExtractBasic
#' @name `[[<-.AmigaBasic`
#' @export
`[[<-.AmigaBasic` <- function(x, i, value) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  x <- as.character(x)
  x[[i]] <- as.character(value)
  return(as.AmigaBasic(x))
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
#' there is only a pointer to the index of the name in that list. Use
#' this function to list, select or replace names included in the code
#' 
#' Make sure that variable and label names are valid for the basic script (see \link{check.names.AmigaBasic}).
#' 
#' @rdname names.AmigaBasic
#' @name names.AmigaBasic
#' @param x An \code{\link{AmigaBasic}}-class object for which to obtain or change variable and/or label names
#' @param value A (\code{vector} of) \code{character} string of desired replacement variable/label names.
#' @return A \code{vector} of \code{character} strings with label and variable names in the basic script.
#' In case of the replacement method a \code{\link{AmigaBasic}}-class with replaced names is returned.
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
#' names(bas)
#' 
#' ## change the first name:
#' names(bas)[1] <- "better.label"
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
names.AmigaBasic <- function(x) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  result <- attr(x, "basic_names")
  if (length(result) == 0) character() else result
}

#' @rdname names.AmigaBasic
#' @name names<-.AmigaBasic
#' @export
`names<-.AmigaBasic` <- function(x, value) {
  if (!inherits(x, "AmigaBasic")) stop("x should be of class AmigaBasic.")
  if (!is.character(value) || length(value) != length(attr(x, "basic_names")))
    stop("Replacement should be a vector of characters of the same length.")
  if (any(duplicated(toupper(value)))) stop("All names should be unique.")
  if (any(nchar(value) > 255 | nchar(value) < 1)) stop("All names should be one or more and less than 255 characters in length.")
  if (any(toupper(value) %in% .amigabasic_commands$command)) stop("Names cannot be reserved AmigaBasic words.")
  if (any(grepl("[^a-zA-Z0-9.]", value, perl = TRUE))) stop("Names should consist of alphanumerics or periods.")
  if (any(grepl("[0-9]", substr(value, 1, 1), perl = T))) stop("Names should not start with numeric characters.")
  if (any(grepl("[ !#$%\\^&*()\\-+=/?,<>:;\"']", value, perl = T))) stop("Names should not contain special characters.")
  attr(x, "basic_names") <- value
  x
}

#' Check Amiga Basic label/variable names for validity
#'
#' Check Amiga Basic label/variable names for validity
#' 
#' Names for variables and labels should adhere to the following rules in Amiga Basic:
#' 
#' \itemize{
#'   \item{Length of the names should be in the range of 1 up to 255 character}
#'   \item{Names cannot be \code{\link{AmigaBasic.reserved}} words}
#'   \item{Names should only contain alphanumeric characters or periods and
#'   should not contain special characters (i.e., reserved for type definition,
#'   such as dollar- or percentage sign)}
#'   \item{Names should not start with a numeric character}
#' }
#' 
#' This function tests names against each of these criteria.
#' 
#' @rdname check.names.AmigaBasic
#' @name check.names.AmigaBasic
#' @param x A \code{vector} of \code{character} strings that need to be checked
#' @param ... Currently ignored.
#' @return A \code{data.frame} with \code{logical} values with the same number of rows as the length of \code{x}.
#' Columns in the data.frame corresponds with the criteria listed in the details.
#' \code{FALSE} for invalid names.
#' @examples
#' \dontrun{
#' ## These are valid names in Amiga Basic:
#' check.names.AmigaBasic(c("Foo", "Bar"))
#' 
#' ## Reserved words and repeated names are not allowed:
#' 
#' check.names.AmigaBasic(c("Print", "Foo", "Foo"))
#' }
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
check.names.AmigaBasic <- function(x, ...) {
  nm <- if (inherits(x, "AmigaBasic")) names(x) else as.character(x)
  result <- data.frame(
    duplicated = duplicated(toupper(nm)),
    length     = nchar(nm) > 255 | nchar(nm) < 1,
    reserved   = toupper(nm) %in% .amigabasic_commands$command,
    characters = grepl("[^a-zA-Z0-9.]", nm, perl = TRUE),
    start      = grepl("^[0-9.]", nm, perl = TRUE)
  )
  row.names(result) <- nm
  result
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

#' Coerce raw data into an AmigaBasicShape class object
#'
#' Coerce raw data into an \code{\link{AmigaBasicShape}}-class object
#'
#' \code{\link{AmigaBasicShape}} objects are comprehensive representations of blitter
#' and sprite graphics that can be used in \code{\link{AmigaBasic}} scripts. Use this function
#' to convert \code{raw} content to an \code{\link{AmigaBasicShape}} object.
#'
#' @rdname rawToAmigaBasicShape
#' @name rawToAmigaBasicShape
#' @param x A \code{vector} of \code{raw} data that is to be converted
#' into an \code{\link{AmigaBasicShape}} class object.
#' @param palette A \code{vector} of \code{character} strings, where each element represents a colour in the palette.
#' This palette will be used to display the graphics (note that the raw format does not store the palette, but this
#' S3 class does). When this argument is omitted a grey scale palette will be generated.
#' @return returns an \code{\link{AmigaBasicShape}}-class object.
#' @examples
#' \dontrun{
#' filename <- system.file("ball.shp", package = "AmigaFFH")
#' 
#' ## read as binary:
#' con      <- file(filename, "rb")
#' ball.raw <- readBin(con, "raw", file.size(filename))
#' close(con)
#' 
#' ## convert raw data into something useful:
#' ball     <- rawToAmigaBasicShape(ball.raw)
#' 
#' ## A shortcut would be to call read.AmigaBasicShape
#' ball2    <- read.AmigaBasicShape(filename)
#' }
#' @family AmigaBasicShapes.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBasicShape <- function(x, palette) {
  ## colorset and dataset seem to be ignored in the basic object editor.
  ## They seem not required for interpretation of the file
  if (missing(palette)) palette <- NULL
  result        <- with(.basic.shape.header, .read.amigaData(x, byte, signed, par.names))
  if (is.null(palette)) {
    palette <- grDevices::grey(seq(0, 1, length.out = 2^result$depth))
  } else {
    if (length(palette) != 2^result$depth)
      stop(sprintf("This shape requires a palette of %i colours, but got %i.", 2^result$depth, length(palette)))
  }
  x             <- x[-1:-sum(abs(.basic.shape.header$byte))]
  result$bitmap <- with(result, bitmapToRaster(x, width, height, depth, palette, interleaved = F))
  attributes(result$bitmap)$palette <- palette
  sz_alt        <- 2*ceiling(result$width/16)*result$height
  sz            <- sz_alt*result$depth
  x             <- x[-1:-sz]
  result$flags <- rev(as.logical(.rawToBitmap(c(raw(2), result$flags), T, F)))[1:16]
  result$planeOnOff <- rev(as.logical(.rawToBitmap(c(raw(2), result$planeOnOff), T, F)))[1:16]
  names(result$flags) <- .amigabasicshape.flags
  if (result$flags["fVSprite"] && result$depth != 2) warning("Unexpected bitmap depth for sprite mode.")
  if (result$flag["imageShadowIncluded"]) {
    result$shadow <- with(result, bitmapToRaster(x, width, height, 1, interleaved = F, palette = c("black", "white")))
    x <- x[-1:-sz_alt]
  }
  if (result$flag["collisionPlaneIncluded"]) {
    result$collision <- with(result, bitmapToRaster(x, width, height, 1, interleaved = F, palette = c("black", "white")))
    x <- x[-1:-sz_alt]
  }
  if (result$flags["fVSprite"]) {
    result$sprite_palette <- amigaRawToColour(x[1:6], "12 bit", "2")
    x <- x[-1:-6]
  }
  if (length(x) > 0) warning("Unexpected and unused trailing data!")
  class(result) <- "AmigaBasicShape"
  result
}

#' Read Amiga Basic Shape files
#'
#' Read Amiga Basic Shape files
#'
#' AmigaBasic used the term 'shapes' for graphics (sprites and blitter objects) which it could display.
#' These graphics were stored in a specific binary format, which can be read with this function. See
#' \code{\link{AmigaBasicShape}} for more details. The file can also be read from a virtual Amiga disk
#' (\code{\link[adfExplorer:amigaDisk-class]{amigaDisk}}).
#' @rdname read.AmigaBasicShape
#' @name read.AmigaBasicShape
#' @param file A \code{character} string of the filename of the Amiga Basic Shape file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @param ... Arguments passed to \code{\link{rawToAmigaBasicShape}}.
#' @return Returns an \code{\link{AmigaBasicShape}} class object read from the \code{file}.
#' @examples
#' \dontrun{
#' filename <- system.file("ball.shp", package = "AmigaFFH")
#' ball     <- read.AmigaBasicShape(filename)
#' ## This is a sprite:
#' ball$flags[["fVSprite"]]
#' 
#' filename <- system.file("r_logo.shp", package = "AmigaFFH")
#' ## The palette is not stored with an Amiga Basic Shape, so let's provide one:
#' r_logo   <- read.AmigaBasicShape(filename,
#'                                  palette = c("#FFFFFF", "#2266BB", "#3366BB", "#4477AA",
#'                                  "#778899", "#999999", "#AAAAAA", "#BBBBBB"))
#' ## This is a blitter object:
#' r_logo$flags[["fVSprite"]]
#' 
#' ## Just for fun, plot it:
#' plot(r_logo)
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
#' (using \code{\link[AmigaFFH]{as.raw}}) and writes it to a file. The file
#' can also be stored onto a virtual Amiga disk
#' (\code{\link[adfExplorer:amigaDisk-class]{amigaDisk}}).
#' 
#' @rdname write.AmigaBasicShape
#' @name write.AmigaBasicShape
#' @param x The \code{\link{AmigaBasicShape}} class object that needs to be
#' stored.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBasicShape}} object) needs to be written.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Invisibly returns the result of the call of \code{close} to the
#' file connection. Or, when \code{disk} is specified, a copy of
#' \code{disk} is returned to which the file(s) is/are written.
#' @examples
#' \dontrun{
#' filename <- system.file("ball.shp", package = "AmigaFFH")
#' ball     <- read.AmigaBasicShape(filename)
#' write.AmigaBasicShape(ball, file.path(tempdir(), "ball.shp"))
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
  if (!inherits(x, "AmigaBasicShape")) stop("x should be of class AmigaBasicShape.")
  sprite        <- x$flags[["fVSprite"]]
  shadow        <- x$flags[["imageShadowIncluded"]]
  collision     <- x$flags[["collisionPlaneIncluded"]]
  x$flags       <- .bitmapToRaw(c(x$flags, rep(F, 16)), T, T)[3:4]
  x$planeOnOff  <- .bitmapToRaw(c(x$planeOnOff, rep(F, 16)), T, T)[3:4]
  result        <- with(.basic.shape.header, .write.amigaData(x[par.names], byte, signed, par.names))
  pal <- attributes(x$bitmap)$palette
  result        <- c(result,
                     .bitmapToRaw(rasterToBitmap(
                       x$bitmap,
                       depth = x$depth,
                       interleaved = F,
                       indexing = function(x, length.out) index.colours(x, length.out,
                                                                        palette = pal)),
                       T, F)
  )
  if (shadow) {
    if (is.null(x$shadow)) stop("Expected shadow layer, but found nothing") else {
      result        <- c(result,
                         .bitmapToRaw(rasterToBitmap(
                           x$shadow,
                           depth = 1,
                           interleaved = F,
                           indexing = function(x, length.out) index.colours(x, length.out,
                                                                            palette = c("black", "white"))),
                           T, F)
      )
                         
    }
  }
  if (collision) {
    if (is.null(x$collision)) stop("Expected shadow layer, but found nothing") else {
      result        <- c(result,
                         .bitmapToRaw(rasterToBitmap(
                           x$collision,
                           depth = 1,
                           interleaved = F,
                           indexing = function(x, length.out) index.colours(x, length.out,
                                                                            palette = c("black", "white"))),
                           T, F)
      )
    }
  }
  if (sprite) {
    if (is.null(x$sprite_palette)) stop("Expected sprite palette, but found nothing") else {
      result <- c(result, colourToAmigaRaw(x$sprite_palette, "12 bit", "2"))
    }
  }
  return(result)
}

#' @export
print.AmigaBasicShape <- function(x, ...) {
  print(sprintf("A %i x %i %s with %i colours to be used in Amiga Basic.",
                x$width, x$height,
                c("blitter object", "sprite")[as.numeric(x$flags[["fVSprite"]]) + 1], 2^x$depth), ...)
}

#' @rdname plot
#' @name plot.AmigaBasicShape
#' @export
plot.AmigaBasicShape <- function(x, y, ...) {
  if (!inherits(x, "AmigaBasicShape")) stop("x should be of class AmigaBasicShape.")
  if (missing(y)) y <- "bitmap"
  plot(as.raster(x, selected = y), ...)
}

#' @rdname as.raster
#' @name as.raster.AmigaBasicShape
#' @export
as.raster.AmigaBasicShape <- function(x, selected = c("bitmap", "shadow", "collision"), ...) {
  if (!inherits(x, "AmigaBasicShape")) stop("x should be of class AmigaBasicShape.")
  bm <- x[[match.arg(selected, c("bitmap", "shadow", "collision"))]]
  if (is.null(bm)) stop(sprintf("No %s layer available in this object!", selected))
  bm
}

#' Convert a grDevices raster object into an AmigaBasicShape class object.
#'
#' Convert a \code{\link[grDevices:as.raster]{raster}} object into an \code{\link{AmigaBasicShape}} class object.
#'
#' TODO
#' 
#' @rdname rasterToAmigaBasicShape
#' @name rasterToAmigaBasicShape
#' @param x A \code{\link[grDevices:as.raster]{raster}} class object to convert into a \code{\link{AmigaBasicShape}} class obejct.
#' @param type A \code{character} string indicating what type of graphic needs to be created: "\code{blitter object}" (default) or "\code{sprite}".
#' @param palette A \code{vector} of \code{character} strings, where each element represents a colour. This palette is used to quantize the
#' colours that occur in the \code{raster} \code{x}.
#' @param ... Arguments passed onto \code{\link{index.colours}}. Can be used, for instance, to achieve specific dithering effects.
#' @return Returns an \code{\link{AmigaBasicShape}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' ## get a raster image:
#' ilbm <- as.raster(read.iff(system.file("ilbm8lores.iff", package = "AmigaFFH")))
#' 
#' ## convert to an Amiga Basic blitter object:
#' bob <- rasterToAmigaBasicShape(ilbm, "blitter object")
#' }
#' @family AmigaBasicShape.operations
#' @family raster.operations
#' @author Pepijn de Vries
#' @export
rasterToAmigaBasicShape <- function(x, type = c("blitter object", "sprite"), palette, ...) {
  if (missing(palette)) {
    palette <- table(x)
    palette <- palette[order(-palette)]
    palette <- names(palette)
  }
  depth <- ceiling(log2(length(palette)))
  ## if not all pixel colours are in palette, the bitmap needs to be quantized
  if (!all((x %in% palette))) {
    x <- apply(
      index.colours(x, length.out = 2^depth, palette = palette, ...), 2,
      function(y) palette[y])
    x <- as.raster(x)
    attributes(x)$palette <- palette
  }
  type <- match.arg(type, c("blitter object", "sprite"))
  if (type == "sprite" && (length(palette) != 4 || attributes(x)$dim[[2]] != 16))
    stop("AmigaBasicShape sprites have to be 16 pixels wide and consist of 4 colours!")
  result <- sapply(.basic.shape.header$par.names, function(x) NULL)
  result$colorset     <- 0
  result$dataset      <- 0
  result$depth        <- depth
  result$width        <- attributes(x)$dim[[2]]
  result$height       <- attributes(x)$dim[[1]]
  result$flags        <- rep(F, 16)
  names(result$flags) <- .amigabasicshape.flags
  result$flags[c("saveBack", "overlay")] <- T
  result$flags["fVSprite"] <- type == "sprite"
  result$planePick    <- 2^depth - 1
  result$planeOnOff   <- rep(F, 16)
  result$bitmap       <- x
  attributes(result$bitmap)$palette <- palette
  # TODO add missing elements (shadow and collision layer)
  if (type == "sprite") result$sprite_palette <- palette[-1] ## background colour is not stored for sprite, hence -1
  class(result) <- "AmigaBasicShape"
  result
}

#' The S3 AmigaBasicBMAP class
#' 
#' A class that represents the content of Amiga Basic BMAP files.
#' 
#' The Amiga operating system made use of library files to execute specific (repetitive/routine) tasks. Amiga Basic
#' was also able to call such routines from library files. In order to do so, it required a 'bmap' file for each
#' library. This file contains a map of the library where it specifies: the name of routine; the `Library Vector Offset'
#' (explained below); and used CPU registers (explained below).
#' 
#' The `Library Vector Offset' is an offset to the base address of a library in memory. This offsets indicates where
#' a specific executable routine starts. The CPU registers are used to (temporary) store (pointers to) input data
#' used by the routine. The BMAP file thus lists which CPU registers are used by specified routines.
#' 
#' @docType class
#' @name AmigaBasicBMAP
#' @rdname AmigaBasicBMAP
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @references \url{https://en.wikipedia.org/wiki/AmigaOS#Libraries_and_devices}
NULL

#' Read Amiga Basic BMAP files
#'
#' Read \code{\link{AmigaBasicBMAP}} binary file format.
#'
#' TODO
#' @rdname read.AmigaBasicBMAP
#' @name read.AmigaBasicBMAP
#' @param file A \code{character} string of the filename of the Amiga Basic BMAP file to be read.
#' @param disk A virtual Commodore Amiga disk from which the \code{file} should be
#' read. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Returns an \code{\link{AmigaBasicBMAP}} class object read from the \code{file}.
#' @examples
#' \dontrun{
#' ## TODO
#' }
#' @family AmigaBasic.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
read.AmigaBasicBMAP <- function(file, disk = NULL) {
  dat <- .read.generic(file, disk)
  rawToAmigaBasicBMAP(dat)
}

#' Write an AmigaBasicBMAP object to a file
#'
#' Write an \code{\link{AmigaBasicBMAP}} class object to a file in its binary format.
#'
#' TODO details
#' 
#' @rdname write.AmigaBasicBMAP
#' @name write.AmigaBasicBMAP
#' @param x The \code{\link{AmigaBasicBMAP}} class object that needs to be
#' stored.
#' @param file A \code{character} string specifying the file location
#' to which \code{x} (an \code{\link{AmigaBasicBMAP}} object) needs to be written.
#' @param disk A virtual Commodore Amiga disk to which the \code{file} should be
#' written. This should be an \code{\link[adfExplorer:amigaDisk-class]{amigaDisk}} object. Using
#' this argument requires the adfExplorer package.
#' When set to \code{NULL}, this argument is ignored.
#' @return Invisibly returns the result of the call of \code{close} to the
#' file connection. Or, when \code{disk} is specified, a copy of
#' \code{disk} is returned to which the file(s) is/are written.
#' @examples
#' \dontrun{
#' ## TODO
#' }
#' @family AmigaBasic.operations
#' @family io.operations
#' @author Pepijn de Vries
#' @export
write.AmigaBasicBMAP <- function(x, file, disk = NULL) {
  .write.generic(x, file, disk)
}

#' @rdname as.raw
#' @name as.raw.AmigaBasicBMAP
#' @export
as.raw.AmigaBasicBMAP <- function(x) {
  .validate_AmigaBasicBMAP(x)
  unlist(lapply(seq_along(x), function(i) {
    c(charToRaw(names(x)[i]), raw(1),
      .amigaIntToRaw(x[[i]]$libraryVectorOffset, 16, T),
      x[[i]]$registers, raw(1))
  }))
}

#' @export
print.AmigaBasicBMAP <- function(x, ...) {
  print(sprintf("An AmigaBasicBMAP with %i references.", length(x)), ...)
}

#' Coerce raw data into an AmigaBasicBMAP class object
#'
#' Coerce raw data into an \code{\link{AmigaBasicBMAP}} class object
#'
#' TODO
#' 
#' @rdname rawToAmigaBasicBMAP
#' @name rawToAmigaBasicBMAP
#' @param x A \code{vector} of \code{raw} data that is to be converted
#' into an \code{\link{AmigaBasicBMAP}} class object.
#' @param ... Currently ignored.
#' @return An \code{\link{AmigaBasicBMAP}} class object based on \code{x}.
#' @examples
#' \dontrun{
#' ##TODO
#' }
#' @family AmigaBasic.operations
#' @family raw.operations
#' @author Pepijn de Vries
#' @export
rawToAmigaBasicBMAP <- function(x, ...) {
  x <- c(raw(1), x)
  terminator <- c(which(x == 0x00), length(x) + 1)
  
  nm <- utils::head(unlist(lapply(
    mapply(seq,
           from = terminator[seq.int(1, length(terminator), 2)] + 1,
           to   = terminator[seq.int(2, length(terminator), 2)] - 1),
    function(y) rawToChar(x[y]))), -1)
  
  lvo <- utils::head(unlist(lapply(
    mapply(seq,
           from = terminator[seq.int(2, length(terminator), 2)] + 1,
           to   = terminator[seq.int(2, length(terminator), 2)] + 2, SIMPLIFY = F),
    function(y) .rawToAmigaInt(as.raw(x[y]), 16, T))), -1)
  registers <- lapply(
    mapply(seq,
           from = utils::head(terminator[seq.int(2, length(terminator), 2)] + 3, -1),
           to   = utils::head(terminator[seq.int(2, length(terminator), 2) + 1] -1, -1), SIMPLIFY = F),
    function(y) as.raw(x[y]))
  no_reg <- diff(terminator) < 4
  no_reg <- which(no_reg[c(F, T)])
  registers[no_reg] <- lapply(seq_along(no_reg), function(i) raw())
  
  result <- lapply(seq_along(nm), function(i) {
    list(libraryVectorOffset = lvo[[i]],
         registers           = registers[[i]])
  })
  names(result) <- nm
  class(result) <- "AmigaBasicBMAP"
  .validate_AmigaBasicBMAP(result)
  return (result)
}

#' @export
as.list.AmigaBasicBMAP <- function(x, ...) {
  unclass(x)
}

#' Coerce raw or named list to an AmigaBasicBMAP class object
#'
#' Coerce raw or named list to an \code{\link{AmigaBasicBMAP}} class object
#'
#' TODO
#' 
#' @rdname as.AmigaBasicBMAP
#' @name as.AmigaBasicBMAP
#' @param x When \code{x} is a \code{vector} of \code{raw} data, it needs to be structured as it would be
#' when stored in a binary file (see \code{\link{read.AmigaBasicBMAP}}). \code{x} can also be a named \code{list},
#' where the name of each element corresponds with a routine in the library. Each element should than consist
#' of a \code{list} with 2 elements: The first should be named `libraryVectorOffset' and should hold the \code{numeric}
#' offset of the routine in the library (see details). The second element should be named `registers' and should
#' contain a \code{vector} of \code{raw} values refering to CPU registers used by the routine (see details).
#' @return Returns a \code{\link{AmigaBasicBMAP}} based on \code{x}
#' @examples
#' \dontrun{
#' ## For the dos.library, the start of the bmap list would look like:
#' dos.list <- list(
#'   xOpen = list(
#'     libraryVectorOffset = -30,
#'     registers = as.raw(2:3)
#'   ),
#'   xClose = list(
#'     libraryVectorOffset = -36,
#'     registers = as.raw(2)
#'   ),
#'   xRead = list(
#'     libraryVectorOffset = -42,
#'     registers = as.raw(2:4)
#'   )
#' )
#' 
#' ## Note that the list above is incomplete, the dos.library holds more routines than shown here.
#' ## This merely serves as an example.
#' ## This list can be converted to an S3 class as follows:
#' dos.bmap <- as.AmigaBasicBMAP(dos.list)
#' }
#' @family AmigaBasic.operations
#' @author Pepijn de Vries
#' @export
as.AmigaBasicBMAP <- function(x) {
  if (typeof(x) == "raw") return(rawToAmigaBasicBMAP(x))
  if (typeof(x) != "list") stop("No method available for converting this object into AmigaBasicBMAP")
  x <- as.list(x)
  class(x) <- "AmigaBasicBMAP"
  .validate_AmigaBasicBMAP(x)
  x
}

.validate_AmigaBasicBMAP <- function(x) {
  if (!inherits(x, "AmigaBasicBMAP")) stop("AmigaBasicBMAP should inherit AmigaBasicBMAP class")
  if (typeof(x) != "list") stop("AmigaBasicBMAP should be of type list")
  if (any(apply(check.names.AmigaBasic(names(x)), 1, any))) stop("AmigaBasicBMAP should contain valid basic names")
  registers_ok  <- unlist(lapply(x, function(y) {
    y$libraryVectorOffset >= -32768 && y$libraryVectorOffset < 0 &&
      (length(y$registers) == 0 || ((y$registers %in% as.raw(1:15)) && !any(duplicated(y$registers))))
  }))
  if (!any(registers_ok)) stop("Register numbers should be unique raw values ranging from 1 to 15, and library vector offsets should be in the range of -1 to -32768")
  return (T)
}
