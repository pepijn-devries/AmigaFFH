AmigaFFH v0.4.7
=============

 * Fixed broken `adfExplorer` dependency
 * Manual and examples updated conform CRAN policies

AmigaFFH v0.4.5
=============

 * Added README
 * Switched to markdown documentation with `roxygen2`
 * Added `pkgdown` website

AmigaFFH v0.4.3
=============

 * Fixed link in documentation in order to pass CRAN
   checks

AmigaFFH v0.4.2
=============

 * Updates in order to compy with latest CRAN policies
   and `roxygen2` standards

AmigaFFH v0.4.1
=============

 * Updated URLs in manual conform CRAN policy

AmigaFFH v0.4.0
=============

 * Fixed bug in Amiga Icon implementation
 * Added support for Amiga Basic files
 * Minor corrections and updates to the manual.

AmigaFFH v0.3.1
=============

 * Resubmission to CRAN due to offline URL references in documentation

AmigaFFH v0.3.0
=============

 * Added support for Amiga Bitmap fonts
 * Added support for interpretation and creation of
   HAM6 and HAM8 mode (special modes on the Amiga) bitmap images.
 * Added support for directly reading files from and writing files to
   amigaDisk class objects (from the adfExplorer package).
 * Modified 'simpleSysConfig' such that specific settings can be
   passed as arguments
 * Fixed bug in viewportmode checks
 * Corrected replace functions for SysConfig objects
 * Suppressed unnecessary warnings in the as.raw.AmigaIcon function
 * Minor corrections to the manual

AmigaFFH v0.2.0
=============

 * Added support for Amiga Workbench Icons (*.info)
 * Added support for the Amiga's system-configuration file
 * Added rasterToHWSprite function
 * Changed 'as.raster' from S4 method to S3 function
   for more flexibility.
 * 'rasterToBitmap' did not decode continuous bitmaps
   correctly. This has been fixed. Fortunately the Amiga
   mostly uses interleaved bitmaps anyway.
 * Added constrains to error distribution in dither method
   of 'index.colours'. Without these constrains errors could
   accumulate for certain images and cause undesirable results.
 * Fixed transparency detection in 'index.colour' function
 * Bug-fix in the interpretation of IFF ANIM frames
 * Moved ProTrackR and adfExplorer from 'Depends' to 'Suggests'
   to avoid circular dependancy
 * Several minor corrections to the manual

AmigaFFH v0.1.2
=============

 * A minor fix in the examples to pass CRAN checks.

AmigaFFH v0.1.1
=============

First release:

 * This is the first release that supports
   - reading and writing IFF files
   - interpretation of several common data flavours in IFF files
   - interpret hardware sprites