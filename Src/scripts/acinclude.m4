# $Id$
#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#	+                                                               +
#	+                       Copyright (c) 2004                      +
#	+          The Regents of the University of California          +
#	+                      All Rights Reserved                      +
#	+                                                               +
#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#	*****************************************************************
#	*								*
#	*	     Griz4 configure - File "configure.ac"    		*
#	*								*
#	*****************************************************************
#
#	Process this file with "autoconf" to produce "configure".
#
#	-----------------------------------------------------------------
#	Revision History:
#
#	20-Feb-04 IRC: Created
#                      This file contains the macros used by the Griz 
#                      configure script.
#
#	04-Oct-04 IRC: Updated
#                      Added -bmaxdata: option to AIX build.
#
#	03-Nov-04 IRC: Updated
#                      Constructed configure options string to pass to
#		       executable.
#
#	05-Nov-04 trt: Updated
#                      Renamed acloca.m4 to the standard acinclude.m4 plan
#                        on using aclocal generated aclocl.m4 file.
#                      Using AC_HELP_STRING for AC_ARG_WITH and 
#                        AC_ARG_ENABLE formating.
#		       Added AC_DEFINE_UNQUOTEDs for config.h define's.
#		       Cleaned up "underquoted call to AC_DEFUN" problems.
#
#       21-Jan-05 IRC: Updated
#                      Changes to support Intel Opteron processor.
#
#       01-Jun-05 JKD: added ability to have user define the directory that 
#		       configure creates for build. Also deleted some reduncies
#                      in include statement for GRIZ build.
#
#       10-Jul-06 JKD: added new parameter HOST_PREFIX set CONFIGURE_HOST
#                        which is used to set the directory prefix and will be used 
#                        in conjuction with the PNG files. 
#                        Fixed the search path for the number of tasks to only look
#                        for the psrinfo file only if it exists on the system. Defaulting 
#                        to 16 otherwise.
#                        Fixed the X11 library path for AIX to /usr/lib instead of 
#                        /usr/X11R6/lib, fixed include directory to match.
#
#       24-Aug-06 IRC: Updated
#                      Add hostname to the build directory path.
#
#       11-Nov-06 IRC: Updated
#                      Added a path for gmake - needed for local version
#                      of gmake.
#
#       21-Dec-06 IRC: Updated
#                      Fixed problems with building on new 64-bit platforms
#                      Changed to pathscale compilers.
#
#       15-Feb-07 IRC: Updated
#                      Fixed batch OSMesa build to use system Mesa Libraries.
#
#       16-Feb-07 IRC: Updated
#                      Add softlink for alpha and beta for testing griz
#		       script.
#
#       02-Oct-07 IRC: Updated
#                      Add support for KlockWorks build - Re: Bill Oliver.
#		       SCR # 429
#
#       25-Apr-08 IRC: Updated
#                      Add support for building on CHAOS4
#		       SCR # 429
#	-----------------------------------------------------------------

#
# The Regents of the University of California.  
# All rights reserved.
# 
# This work was produced at the University of California, 
# Lawrence Livermore National Laboratory (UC LLNL) under 
# contract no. W-7405-ENG-48 between the U.S. Department of 
# Energy (DOE) and The Regents of the University of California
# (University) for   the operation of UC LLNL. Copyright is 
# reserved to the University for purposes of controlled dis-
# semination, commercialization through formal licensing, or 
# other disposition under terms of Contract 48; DOE policies, 
# regulations and orders; and U.S. statutes.
# 
# DISCLAIMER OF LIABILITY
# 
# This document was prepared as an account of work sponsored 
# by an agency of the United States Government.  Neither the 
# United States Government nor the University of California 
# nor any of their employees, makes any warranty, express
# or implied, including the warranties of merchantability
# and fitness for a particular purpose, or assumes any legal
# liability or responsibility for the accuracy, completeness,
# or usefulness of any information, apparatus, product, or 
# process disclosed, or represents that its use would not 
# infringe privately-owned rights.
# 
# DISCLAIMER OF ENDORSEMENT
# 
# Reference herein to any specific commercial products,
# process, or service by trade name, trademark, manufacturer,
# or otherwise, does not necessarily constitute or imply its
# endorsement, recommendation, or favoring by the United States
# Government or the University of California.  The views and 
# opinions of authors expressed herein do not necessarily state
# or reflect those of the United States Government or the 
# University of California, and shall not be used for advertising
# or product endorsement purposes.
#################################################################

AC_DEFUN([CONFIGURE_INIT],
  [
        #
        # Init Configure options string
        #
        CONFIG_OPTIONS="CONFIGURE=configure "
        AC_SUBST(CONFIG_OPTIONS)
  ])

AC_DEFUN([CONFIGURE_VERSION],
  [
        #
        # Set the codes major version
        #
#       GRIZ_VERSION="V4_10"
#       AC_ARG_WITH(griz_version,
#    AC_HELP_STRING([--with-griz_version=[PATH]],[Set Griz Version]),
#            GRIZ_VERSION="${withval}" &&
#            AC_MSG_RESULT("Griz Version String : $withval")
#           )
#       AC_MSG_RESULT(Code Version=$GRIZ_VERSION)
#       AC_SUBST(GRIZ_VERSION)
  ])

AC_DEFUN([CONFIGURE_OS],
  [
        #
        # Determine the name and version of the operating system.
        #
        AC_MSG_CHECKING(for name and version of operating system)
        OS_NAME="`(uname -s) 2> /dev/null`"
        OS_NAME_VERSION="`(uname -s -r) 2> /dev/null`"
        
        PROCESSOR_TYPE="`(uname -p) 2> /dev/null`"

        AC_DEFINE_UNQUOTED(OS_NAME,        $OS_NAME,        [Operating system name])
        AC_DEFINE_UNQUOTED(PROCESSOR_TYPE, $PROCESSOR_TYPE, [Processor Type])

        AC_SUBST(OS_NAME)
        AC_SUBST(PROCESSOR_TYPE)
        AC_MSG_RESULT($OS_NAME_VERSION)
  ])


AC_DEFUN([CONFIGURE_HOST],
  [

	#
	#  Set the prfix of the build directory
	#
	HOST_PREFIX="GRIZ4-"
        #
        # Determine the name of the HOST machine
        #
        AC_MSG_CHECKING(for name of host machine)
        HOSTNAME="`(hostname | tr -d '[0-9]' |  tr '[a-z]' '[A-Z]') 2> /dev/null`"

        HOSTDIR="`(echo $SYS_TYPE) 2> /dev/null`"

        # Determine the name of the host directory - this is usually where the 
        # code is installed. We determine this by finding the path for griz
        if test "$HOSTDIR" = ""; then
                 TEMPLIST=`echo ":$TEMPDIR" | tr '[/]' '[ ]'`
                 HOSTDIR=$OS_NAME

                 HOSTDIR="`(type griz) 2> /dev/null`"
                 HOSTDIR="`(ls -l $HOSTDIR | grep mdg) 2> /dev/null`" 

                 for dir in $HOSTDIR; do 
                     TEMPDIR=$dir
                 done

                 TEMPLIST=`echo ":$TEMPDIR" | tr '[/]' '[ ]'`

                 found_platform="False"
                 HOSTDIR=""

                 for dir in $TEMPLIST; do 
                     TEMPDIR=$dir
        
                     if test "$found_platform" = "True"; then
                        HOSTDIR=$dir
                        found_platform="False" 
                     fi
 
                     if test "$dir" = "mdg"; then
                        found_platform="True" 
                     fi
                 done

                 if test "$HOSTDIR" = ""; then
                    TEMPLIST=`echo ":$TEMPDIR" | tr '[/]' '[ ]'`
                    HOSTDIR=$OS_NAME
                 fi
 
                 if test "$HOSTDIR" = "bin"; then
                    HOSTDIR=$HOSTNAME
                 fi
        fi
        
        HOSTDIR="$HOSTDIR-$HOSTNAME"

        AC_DEFINE_UNQUOTED(HOSTDIR, $HOSTDIR, [Host directory name ($SYS_TYPE)])
        AC_SUBST(HOSTDIR)
        AC_DEFINE_UNQUOTED(HOSTNAME, $HOSTNAME, [Host name])
        AC_SUBST(HOSTNAME)
        AC_MSG_RESULT(Hostname=$HOSTNAME)
  ])


AC_DEFUN([CONFIGURE_GPROF],
  [       
        #
        # Set options for running Griz under gprof
        #
        AC_ARG_ENABLE(gprof,
               AC_HELP_STRING(
                      [--enable-gprof],[Compile and load with GPROF]),
		      GPROF_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling: GPROF),
		      GPROF_ENABLE="False")

        if test "$GPROF_ENABLE" = "True"; then
                CFLAGS_TEMP+=" -pg"
                EXE_SUFFIX:="$EXE_SUFFIX_gprof"

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-gprof "  
        fi
  ])

AC_DEFUN([CONFIGURE_BATCH],
  [       
        #
        # Set options for running Griz in Batch
        #
  ])


AC_DEFUN([CONFIGURE_BITMAPS],
  [       
        #
        # Set options for Bitmaps
        #

        BITMAPS_HOME=".."

        AC_ARG_WITH([bitmaps],
               AC_HELP_STRING(
	            [--with-bitmaps=[PATH]],
                    [Use given base PATH for BITMAPS files]),
	            BITMAPS_HOME="${withval}" &&
	            AC_MSG_RESULT("Using BITMAPS Path : $withval")
	           )

        BITMAPS_INCLUDE_PATHS="-I$BITMAPS_HOME"

        # Bitmaps Options
        AC_SUBST(BITMAPS_HOME)
        AC_SUBST(BITMAPS_INCLUDE_PATHS)
  ])


AC_DEFUN([CONFIGURE_EXODUS],
  [       
        #
        # Set options for running with EXODUS support
        #

        AC_ARG_ENABLE([exodus],
               AC_HELP_STRING(
                 [--enable-exodus],[Compile and load with EXODUS]),
		EXODUS_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling: Exodus),
		EXODUS_ENABLE="False")

        EXO_HOME=""
        EXO_DEFINES=""
        EXO_INCLUDE_PATHS=""
        EXO_LIBRARY_PATHS=""

        if test "$EXODUS_ENABLE" = "True"; then
                EXO_HOME="/usr/local/lib"

                AC_ARG_WITH([exodus],
                    AC_HELP_STRING(
	                   [--with-exodus=[PATH]],
                           [Use given base PATH for EXODUS libraries and header files]),
	                   EXO_HOME="${withval}" &&
	                   AC_MSG_RESULT("Using EXODUS Path : $withval")
	                   )

                EXO_DEFINES="-DEXO_SUPPORT"
                EXO_INCLUDE_PATHS="-I$EXODUS_HOME/include"
                EXO_LIBRARY_PATHS="-L$EXODUS_HOME/lib"
                EXO_LIBRARY="-lexoIIc"
                NETCDF_ENABLE="True"
                EXO_TARGET="libgex.so"

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-exodus --with-exodus=$EXODUS_HOME "  
                AC_DEFINE(HAVE_EXODUS, 1, [Exodus library])
        fi

        # Exodus Library Options
        AC_SUBST(EXO_ENABLE)
        AC_SUBST(EXO_HOME)
        AC_SUBST(EXO_DEFINES)
        AC_SUBST(EXO_LIBRARY)
        AC_SUBST(EXO_LIBRARY_PATH)
        AC_SUBST(EXO_INCLUDE_PATH)
        AC_SUBST(EXO_TARGET)
  ])


AC_DEFUN([CONFIGURE_HERSHEY],
  [       
        #
        # Set options for running with HERSHEY support
        #
        HERSHEY_HOME="../HersheyLib"

        AC_ARG_WITH([hershey],
               AC_HELP_STRING(
	            [--with-hershey=[PATH]],
                    [Use given base PATH for HERSHEY libraries and header files]),
	            HERSHEY_HOME="${withval}" &&
		    CONFIG_OPTIONS="$CONFIG_OPTIONS --with-hershey=$HERSHEY_HOME " &&
	            AC_MSG_RESULT("Using HERSHEY Path : $withval")
	           )

	HERSHEY_FONTS_SET="False"
        AC_ARG_WITH([hershey_fonts],
               AC_HELP_STRING(
	            [--with-hershey_fonts=[PATH]],
                    [Use given base PATH for HERSHEY Font files]),
	            HERSHEY_FONTS_PATH="${withval}" &&
	            HERSHEY_FONTS_SET="True" &&
	            CONFIG_OPTIONS="$CONFIG_OPTIONS --with-hershey_fonts=$HERSHEY_FONTS_PATH " &&
		    AC_MSG_RESULT("Using HERSHEY_FONTS Path : $withval")
	           )

        HERSHEY_INCLUDE_PATHS="-I$HERSHEY_HOME/src"
        HERSHEY_LIBRARY_PATHS="-L$HERSHEY_HOME/src"
        HERSHEY_LIBRARY="-lhershey"

        # Hershey Library Options
        AC_SUBST(HERSHEY_FONTS_PATH)
        AC_SUBST(HERSHEY_HOME)
        AC_SUBST(HERSHEY_LIBRARY)
        AC_SUBST(HERSHEY_LIBRARY_PATHS)
        AC_SUBST(HERSHEY_INCLUDE_PATHS)
  ])

AC_DEFUN([CONFIGURE_IMAGE],
  [       
        #
        # Set options for the Image Library
        #
        IMAGE_HOME="../ImageLib"

        AC_ARG_WITH([image],
               AC_HELP_STRING(
	            [--with-image=[PATH]],
                    [Use given base PATH for IMAGE libraries and header files]),
	            IMAGE_HOME="${withval}" &&
                    CONFIG_OPTIONS="$CONFIG_OPTIONS --with-image=$IMAGE_HOME " &&
	            AC_MSG_RESULT("Using IMAGE Path : $withval")
                   )

        IMAGE_INCLUDE_PATHS="-I$IMAGE_HOME"
        IMAGE_LIBRARY_PATHS="-L$IMAGE_HOME"
        IMAGE_LIBRARY="-limage"

        #  IMAGE Library Options
        AC_SUBST(IMAGE_HOME)
        AC_SUBST(IMAGE_LIBRARY)
        AC_SUBST(IMAGE_LIBRARY_PATHS)
        AC_SUBST(IMAGE_INCLUDE_PATHS)
  ])


AC_DEFUN([CONFIGURE_JPEG],
  [    
	JPEG_ENABLE="True"  
	CONFIG_JPEG="True"
        
	#
	# Set options for running with JPEG support
	#
	AC_ARG_ENABLE([nojpeg],       
		AC_HELP_STRING([--enable-jpeg],[Compile and load with JPEG support]),
		JPEG_ENABLE="False" &&  AC_MSG_RESULT(---->Enabling: JPEG),
		JPEG_ENABLE="True")
        
	#
	#  Make JPEG default enabled 
	#  
	
        
	JPEG_HOME="None"
	JPEG_INCLUDE_PATHS=""
	JPEG_LIBRARY_PATHS=""
	JPEG_LIBRARY=""
	JPEG_DEFINES=""

        
	if test "$JPEG_ENABLE" = "True"; then
		AC_ARG_WITH([jpeg],
			AC_HELP_STRING(
				[--with-jpeg=[PATH]],
				[Use given base PATH for JPEG libraries and header files]),
			JPEG_HOME="${withval}" &&
			CONFIG_OPTIONS="$CONFIG_OPTIONS --with-jpeg=$JPEG_HOME " jpeg_set="True" &&
			AC_MSG_RESULT("Using JPEG Path : $withval")
        	   
		)
		if test "$JPEG_HOME" = "None";then
         		JPEG_HOME=".."
			CONFIG_JPEG="true"
		fi

  		JPEG_INCLUDE_PATHS="-I$JPEG_HOME/include"
		JPEG_LIBRARY_PATHS="-L$JPEG_HOME/lib"

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-jpeg "
          
		JPEG_LIBRARY="-ljpeg "

		if test "$OS_NAME" = "OSF1"; then 
			JPEG_LIBRARY_PATHS=" -noso $JPEG_LIBRARY_PATHS"
			JPEG_LIBRARY="$JPEG_LIBRARY -so_archive "
		fi
		
		JPEG_DEFINES="-DJPEG_SUPPORT"
		AC_DEFINE(HAVE_JPEG, 1, [Jpeg library])
        
	fi

        
	#  JPEG Library Options
	AC_SUBST(JPEG_HOME)
	AC_SUBST(CONFIG_JPEG)
	AC_SUBST(JPEG_DEFINES)
	AC_SUBST(JPEG_LIBRARY)
	AC_SUBST(JPEG_LIBRARY_PATHS)
  ])
	
	
AC_DEFUN([CONFIGURE_X11],
  [       
        #
        # Set options for X11 Libraries
        # -lGLw -lXm -lXt -lX11 -lGLU -lGL -lm
	paths="/usr/X11R6 /usr"
	X11_INCLUDE_PATHS=""
	X11_LIBRARY_PATHS=""
	X11_LIBS=" "
	X11_PATH=""
	X11_HOME=""   
	allfound="true"
	   
	   
	required_files="libX11.a libXm.a libXt.a libGL.a libGLU.a libGLw.a"
	for path in $paths; do
		allfound="true"
		for file in $required_files; do
			AC_CHECK_FILE([$path/lib64/$file],
				[
					X11_INCLUDE_PATHS="-I$path/include "
					X11_LIBRARY_PATHS="-L$path/lib64 "
					LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.a||'`
					X11_LIBS="$X11_LIBS -l$LIBNAME "
					X11_PATH="$path/lib64"
					X11_HOME="$path"
				],
				[
					allfound="false"
					X11_LIBS=" "
					break
				]
			)
		done
		if test "$allfound" = "false";then
			allfound="true"
			for file in $required_files; do
				AC_CHECK_FILE([$path/lib/$file],
					[
						X11_INCLUDE_PATHS="-I$path/include "
						X11_LIBRARY_PATHS="-L$path/lib "
						LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.a||'`
						X11_LIBS="$X11_LIBS -l$LIBNAME "
						X11_PATH="$path/lib"
						X11_HOME="$path"
					],
					[
						allfound="false"
						X11_LIBS=" "
						break
					]
				)
			done
			if test "$allfound" = "false";then
				X11_INCLUDE_PATHS=""
	   		X11_LIBRARY_PATHS=""
	   		X11_PATH=""
				X11_HOME=""
			else
				break
			fi
		else
			break
		fi
	done
	
	
	if test "$allfound" = "false";then
		for path in $paths; do
			required_files="libX11.so libXm.so libXt.so libGL.so libGLU.so libGLw.so"
			allfound="true"
			for file in $required_files; do
				AC_CHECK_FILE([$path/lib64/$file],
					[
						X11_INCLUDE_PATHS="-I$path/include "
						X11_LIBRARY_PATHS="-L$path/lib64 "
						LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.so||'`
						X11_LIBS="$X11_LIBS -l$LIBNAME "
						X11_PATH="$path/lib64"
						X11_HOME="$path"
					],
					[
						allfound="false"
						X11_LIBS=" "
						break
					]
				)
			done
			if test "$allfound" = "false";then
				allfound="true"
				for file in $required_files; do
					AC_CHECK_FILE([$path/lib/$file],
						[
							X11_INCLUDE_PATHS="-I$path/include "
							X11_LIBRARY_PATHS="-L$path/lib "
							LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.a||'`
							X11_LIBS="$X11_LIBS -l$LIBNAME "
							X11_PATH="$path/lib"
							X11_HOME="$path"
						],
						[
							allfound="false"
							X11_LIBS=" "
							break
						]
					)
				done
				if test "$allfound" = "false";then
					X11_INCLUDE_PATHS=""
	   			X11_LIBRARY_PATHS=""
	   			X11_PATH=""
					X11_HOME=""
				else
					break
				fi
			else
				break
			fi
		done	
	fi
	   if test "$allfound" = "false"; then
		AC_MSG_ERROR([Could not locate X11 libraries in $paths])
	   else
		AC_MSG_RESULT([X11 libraries validated at $X11_HOME])
	   fi

        AC_SUBST(X11_LIBRARY_PATHS)
        AC_SUBST(X11_INCLUDE_PATHS)
		AC_SUBST(X11_LIBS)
  ])


AC_DEFUN([CONFIGURE_MESA],
  [       
        #
        # Set options for running with the MESA Graphics Library
        #
        AC_ARG_ENABLE([mesa],
               AC_HELP_STRING(
                 [--enable-mesa],[Compile and load with MESA]),
			[
				MESA_ENABLE="True"  
				AC_MSG_RESULT(Enabling : MESA)
			],
			[
				MESA_ENABLE="False"
			]
		)

        	MESA_HOME=""
        	MESA_INCLUDE_PATHS=""
        	MESA_LIBRARY_PATHS=""
        	MESA_LIBRARY=""
        	MESA_DEFINES=""

		AC_ARG_WITH([mesa],
			AC_HELP_STRING([--with-mesa=[PATH]],[Use given base PATH for MESA libraries and header files]),
			[
				AC_CHECK_FILE([${withval}/lib64/libMesaGL.a],
					[
						MESA_LIBRARY_PATHS="-L${withval}/lib64 "
						MESA_HOME="${withval}"
						MESA_INCLUDE_PATHS="-I${withval}/include "
						MESA_ENABLE="True" 
					],
					[
						AC_CHECK_FILE([${withval}/lib/libMesaGL.a],
							[
								MESA_LIBRARY_PATHS="-L${withval}/lib "
								MESA_HOME="${withval}"
								MESA_INCLUDE_PATHS="-I${withval}/include "
								MESA_ENABLE="True" 
							],
							[
								AC_MSG_WARN([No Mesa library at the given location: ${withval}])
								MESA_ENABLE="False"
							]
						)
					]
				)

				AC_MSG_RESULT("Using MESA Path : $withval")
			]
		)
        if test "$MESA_ENABLE" = "True"; then

			MESA_HOME="/usr/X11R6"


                MESA_HOME=""


                MESA_INCLUDE_PATHS="-I$MESA_HOME/include"        
                MESA_LIBRARY_PATHS="-L$MESA_HOME/lib"

                BATCH_ENABLE="False"
                EXE_SUFFIX="$EXE_SUFFIX""_mesa"

                if test "$bits64_ENABLE" = "True"; then
                   MESA_LIBRARY_PATHS="-L$MESA_HOME/lib64"
                fi

#
###############################################################################
#
# Note: On Linux this link works:
#       LIBS = -L../lib -ljpeg -L/usr/local/tools/libpng/lib -lpng -L../HersheyLib/src -lhershey -L../ImageLib \
#              -limage -L/usr/X11R6/lib64 -lXm -lXt -lX11 -lm -lGLU -lOSMesa -lMesaGL -lGLw
#
#
# For SGI Use: -L/home1/corey3/GrizSGI-April0307/GrizDist/GRIZ4-IRIX64-KRAKOV_32bit_osmesa/ext/Mesa-5.0/lib \ 
#               -lMesaGLU -lMesaGL -lGLw -lGLU -lOSMesa -lXm -lXt -lX11 -lm -lPW
#
# For AIX Use:  -L/g/g14/icorey/MDG/Mesa/Mesa-3.5/lib/ -lOSMesa -lGL -lGLU -lGLw \
#               -lOSMesa -lMesaGL
#
###############################################################################
#      
		MESA_LIBRARY="-lMesaGLU -lMesaGL -lGLw -lOSMesa -lXm -lXt -lX11 "

		if test "$OS_NAME" = "Linux" -o \
		        "$OS_NAME" = "HPUX"; then
		        MESA_LIBRARY_PATHS="-L/usr/X11R6/lib64 "
			MESA_LIBRARY="-lXm -lXt -lX11 -lm -lGLU -lOSMesa -lGLw "
		fi

		if test "$OS_NAME" = "AIX"; then
		        MESA_LIBRARY_PATHS="-L/usr/X11R6/lib "
			MESA_LIBRARY=" -lOSMesa -lGL -lGLU -lGLw -lOSMesa -lXm -lXt -lX11 "
		fi

		if test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX"; then
                  	MESA_LIBRARY="-lMesaGLU -lGLw -lMesaGLwM -lOSMesa -lXm -lXt -lX11 "
	        fi

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-mesa --with-mesa=$MESA_HOME "
		AC_DEFINE(HAVE_MESA, 1, [Mesa library])
     fi

        #  MESA Library Options
        AC_SUBST(CONFIG_MESA)
        AC_SUBST(MESA_ENABLE)
        AC_SUBST(MESA_HOME)
        AC_SUBST(MESA_LIBRARY)
        AC_SUBST(MESA_LIBRARY_PATHS)
        AC_SUBST(MESA_INCLUDE_PATHS)


])

AC_DEFUN([CONFIGURE_OSMESA],
  [
	PATH_SET="false"      
	STD_PATHS="/usr /usr/X11R6 /usr/local"
	OSMESA_HOME=""
	OSMESA_INCLUDE_PATHS=""
	OSMESA_LIBRARY_PATHS=""
	OSMESA_LIBRARY=""
	OSMESA_DEFINES=""
	BATCH_DEFINES=""

	AC_ARG_WITH([mesa],        
		[
			AC_HELP_STRING([--with-mesa=[PATH]],[Use given base PATH for MESA libraries and header files])
		],
		[
			OSMESA_HOME="${withval}"
			AC_CHECK_FILE([$OSMESA_HOME/lib64],
				[
					AC_CHECK_FILE([$OSMESA_HOME/lib64/libOSMesa.a],
						[
							OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib64
							PATH_SET="true"
						],
						[
							AC_CHECK_FILE([$OSMESA_HOME/lib64/libOSMesa.so],
								[
									OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib64
									PATH_SET="true"
								]
							)
						]
					)
				],
				[
					AC_CHECK_FILE([$OSMESA_HOME/lib/libOSMesa.a],
						[
							OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib
							PATH_SET="true"
						],
						[
							AC_CHECK_FILE([$OSMESA_HOME/lib/libOSMesa.so],
								[
									OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib
									PATH_SET="true"
								],
								[
									AC_CHECK_FILE([$OSMESA_HOME/libOSMesa.a],
										[
											OSMESA_LIBRARY_PATHS=$OSMESA_HOME
											PATH_SET="true"
										],
										[
											AC_CHECK_FILE([$OSMESA_HOME/libOSMesa.so],
												[
													OSMESA_LIBRARY_PATHS=$OSMESA_HOME
													PATH_SET="true"
												]
											)
										]
									)
								]
							)
						]
					)
				]
			)
			if test "$PATH_SET" = "true"; then
				AC_MSG_RESULT([Using MESA Path : $OSMESA_LIBRARY_PATHS])
			fi
		]
	)
	suffix="a"
	# lets check the usual suspects for Mesa libraries
	if test "$PATH_SET" = "false"; then
		AC_CHECK_FILE([$X11_PATH/libOSMesa.a],
			[
				OSMESA_HOME="$X11_HOME"
				OSMESA_LIBRARY_PATHS="$X11_PATH"
				PATH_SET="true"
			],
			[
				AC_CHECK_FILE([$X11_PATH/libOSMesa.so],
					[
						OSMESA_HOME="$X11_HOME"
						OSMESA_LIBRARY_PATHS="$X11_PATH"
						PATH_SET="true"
						suffix="so"
						
					],
			
					[
						for path in $STD_PATHS; do
						AC_CHECK_FILE([$path/lib64/libOSMesa.a],
							[
								OSMESA_HOME=$path
								OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib64
								PATH_SET="true"
								break
							],
							[
								AC_CHECK_FILE([$path/lib64/libOSMesa.so],
									[
										OSMESA_HOME=$path
										OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib64
										PATH_SET="true"
										suffix="so"
										break
									],
									[
										AC_CHECK_FILE([$path/lib/libOSMesa.a],
											[
												OSMESA_HOME=$path
												OSMESA_LIBRARY_PATHS=$OSMESA_HOME/lib
												PATH_SET="true"
												break
											],
											[
												AC_CHECK_FILE([$path/lib/libOSMesa.so],
													[
														OSMESA_HOME=$path
														OSMESA_LIBRARY_PATHS=$OSMESA_HOME
														PATH_SET="true"
														suffix="so"
														break
													]
												)
											]
										)
									]
								)
							]
						)
						done
					]
				)
			]
		)	
	fi
	NO_MESA_X11_LIBS=""
	if test "$PATH_SET" = "true"; then
		NO_MESA_X11_LIBS=$X11_LIBS
		OSMESA_LIBRARY=" -lOSMesa"
		MESAFILES="libGL.$suffix libGLU.$suffix libGLw.$suffix"
		for file in $MESAFILES; do
			AC_CHECK_FILE([$OSMESA_LIBRARY_PATHS/$file],
				[
					if test "$suffix" = "a";then
						LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.a||'`	
					else
						LIBNAME=`echo $file | sed -e 's|lib||' | sed -e 's|.so||'`	
					fi
					RMNAME="-l$LIBNAME"
					COMMAND="sed -e s|$RMNAME||"
					NO_MESA_X11_LIBS=`echo $NO_MESA_X11_LIBS | $COMMAND` 
					OSMESA_LIBRARY="$OSMESA_LIBRARY -l$LIBNAME "
				]
			)
		done

		OSMESA_INCLUDE_PATHS==$OSMESA_HOME/include
		OSMESA_LIBRARY_PATHS="-L$OSMESA_LIBRARY_PATHS"
		BATCH_DEFINES="-DSERIAL_BATCH"


		CONFIG_OPTIONS="$CONFIG_OPTIONS --with-osmesa=$OSMESA_HOME "
		AC_DEFINE(HAVE_OSMESA, 1, [OSMesa library])
		AC_DEFINE(HAVE_BATCH, 1, [OSMesa library])
		OSMESA_ENABLE="True"
		BATCH_ENABLE="True" 

	else
		OSMESA_ENABLE="False"
		BATCH_ENABLE="False" 
	fi


        
	#  MESA Library Options
         
	AC_SUBST(OSMESA_ENABLE)
	AC_SUBST(OSMESA_HOME)
	AC_SUBST(OSMESA_LIBRARY)
	AC_SUBST(OSMESA_LIBRARY_PATHS)
	AC_SUBST(OSMESA_INCLUDE_PATHS)
	AC_SUBST(BATCH_ENABLE)
	AC_SUBST(BATCH_DEFINES)
	AC_SUBST(NO_MESA_X11_LIBS)
  ])

AC_DEFUN([CONFIGURE_MILI],
  [       
        #
        # Set options for the MILI Library
        #
        MILI_HOME="/usr/apps/mdg"
 
        	AC_ARG_WITH([mili],
               AC_HELP_STRING(
				[--with-mili=[PATH]],
				[Use given base PATH for MILI libraries and header files]
			),
			[
				AC_CHECK_FILE([${withval}/lib/libmili.a],
					[
	            			MILI_HOME="${withval}"
	            			AC_MSG_RESULT([Using MILI Path : $withval])
					],
					[
						AC_MSG_ERROR([No Mili library at $withval/lib])
					]
				)
			],
			[
				AC_CHECK_FILE([$MILI_HOME/lib/libmili.a],
					[
	            			AC_MSG_RESULT([Using MILI Path :$MILI_HOME ])
					],
					[
						AC_MSG_ERROR([No Mili library found. Please set --with-mili=PATH to correct Mili library])
					]
				)
			]
		)

		CONFIG_OPTIONS="$CONFIG_OPTIONS --with-mili=$MILI_HOME "

        	MILI_INCLUDE_PATHS="-I$MILI_HOME/include"
        	MILI_LIBRARY_PATHS="-L$MILI_HOME/lib"

        	MILI_LIBRARY="-leprtf -lmili -ltaurus"

        AC_DEFINE(HAVE_MILI, 1, [Mili library])

        #  MILI Library Options
        AC_SUBST(MILI_HOME)
        AC_SUBST(MILI_LIBRARY)
        AC_SUBST(MILI_LIBRARY_PATHS)
        AC_SUBST(MILI_INCLUDE_PATHS)
  ])


AC_DEFUN([CONFIGURE_NETCDF],
  [       
        #
        # Set options for running with the NETCDF Library
        #
        AC_ARG_ENABLE([netcdf],
               AC_HELP_STRING(
                      [--enable-netcdf],[Compile and load with NETCDF]),
		      NETCDF_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling : NETCDF),
		      NETCDF_ENABLE="False")

        NETCDF_HOME=""
        NETCDF_INCLUDE_PATHS=""
        NETCDF_LIBRARY_PATHS=""
        NETCDF_LIBRARY=""

        if test "$NETCDF_ENABLE" = "True"; then
                NETCDF_HOME="/usr/local/lib"

                AC_ARG_WITH([netcdf],
               AC_HELP_STRING(
	            [--with-netcdf=[PATH]],
                    [Use given base PATH for NETCDF libraries and header files]),
	            NETCDF_HOME="${withval}" &&
	            AC_MSG_RESULT("Using NETCDF Path : $withval")
	        )

                NETCDF_INCLUDE_PATHS="-I$NETCDF_HOME/include"
                NETCDF_LIBRARY_PATHS="-L$NETCDF_HOME/lib"
                NETCDF_LIBRARY="-lnetcdf"

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-netcdf --with-netcdf=$NETCDF_HOME "
                AC_DEFINE(HAVE_NETCDF, 1, [NetCDF library])
        fi

        #  NETCDF Library Options
        AC_SUBST(NETCDF_ENABLE)
        AC_SUBST(NETCDF_HOME)
        AC_SUBST(NETCDF_LIBRARY)
        AC_SUBST(NETCDF_LIBRARY_PATHS)
        AC_SUBST(NETCDF_INCLUDE_PATHS)
  ])


AC_DEFUN([CONFIGURE_PAPI],
  [       
        #
        # Set options for running with PAPI support
        #
        AC_ARG_ENABLE([papi],
               AC_HELP_STRING(
                      [--enable-papi],[Compile and load with PAPI Library]),
		      PAPI_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling : PAPI),
		      PAPI_ENABLE="False")

        PAPI_HOME=""
        PAPI_INCLUDE_PATHS=""
        PAPI_LIBRARY_PATHS=""
        PAPI_LIBRARY=""
        PAPI_DEFINES=""

        if test "$PAPI_ENABLE" = "True"; then
                PAPI_HOME="/usr/local/papi"

                AC_ARG_WITH([papi],
               AC_HELP_STRING(
	            [--with-papi=[PATH]],
                    [Use given base PATH for PAPI libraries and header files]),
	            PAPI_HOME="${withval}" &&
	            AC_MSG_RESULT("Using PAPI Path : $withval")
	          )

                PAPI_INCLUDE_PATHS="-I$PAPI_HOME/include"
                PAPI_LIBRARY_PATHS="-L$PAPI_HOME/lib"
                PAPI_LIBRARY="-lpapi"
                PAPI_DEFINES="-DPAPI_SUPPORT"

		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-papi --with-papi=$PAPI_HOME "
                AC_DEFINE(HAVE_PAPI, 1, [Papi library])
        fi

        # PAPI Library Options
        AC_SUBST(PAPI_ENABLE)
        AC_SUBST(PAPI_HOME)
        AC_SUBST(PAPI_DEFINES)
        AC_SUBST(PAPI_LIBRARY)
        AC_SUBST(PAPI_LIBRARY_PATHS)
  ])


AC_DEFUN([CONFIGURE_PNG],
  [       
        PNG_ENABLE="False"
        #
        # Set options for running with PNG support
        #
        AC_ARG_ENABLE([nopng],
               AC_HELP_STRING(
                      [--enable-png],[Compile and load with PNG Library]),
		      PNG_ENABLE="False" &&  AC_MSG_RESULT(---->Enabling: PNG),
		      PNG_ENABLE="True")
	
	PNG_HOME=""
	PNG_INCLUDE_PATHS=""
	PNG_LIBRARY_PATHS=""
	ZLIB_HOME="None"
	ZLIB_INCLUDE_PATHS=""
	ZLIB_LIBRARY_PATHS=""
	PNG_LIBRARY=""
	PNG_DEFINES=""
	CONFIG_PNG="false"
	CONFIG_ZLIB="false"
        

     
	if test "$PNG_ENABLE" = "True"; then
               
		AC_ARG_WITH([png],
			AC_HELP_STRING(
				[--with-png=[PATH]],
				[Use given base PATH for PNG libraries and header files, "local"
				the png libs included with Griz]),
			PNG_HOME="${withval}" &&
			AC_MSG_RESULT("Using PNG Path : $withval")
		) 
		if test "$PNG_HOME" = "local"; then
			PNG_HOME=".."
			ZLIB_HOME=".."
			CONFIG_ZLIB="true"
			CONFIG_PNG="true"
		elif test -e "/usr/local/tools/libpng/lib/libpng.a" -o -e "/usr/local/tools/libpng/lib/libpng.so"; then
			PNG_HOME="/usr/local/tools/libpng"
			if test "$HOSTNAME" = "GPS";then
				ZLIB_HOME="$PNG_HOME"
				ZLIB_LIBRARY_PATHS="-L$ZLIB_HOME/lib"
				ZLIB_LIBRARY=" -lz "
			fi
		else
			PNG_HOME=".."
			ZLIB_HOME=".."
			CONFIG_ZLIB="true"
			CONFIG_PNG="true"
			
			ZLIB_LIBRARY_PATHS="-L$ZLIB_HOME/lib"
			ZLIB_INCLUDE_PATHS="-I$ZLIB_HOME/include"
			ZLIB_LIBRARY=" -lz "		
		fi
		
		PNG_LIBRARY_PATHS="-L$PNG_HOME/lib"
		PNG_INCLUDE_PATHS="-I$PNG_HOME/include"
		
		PNG_LIBRARY=" -lpng "
		PNG_DEFINES="-DPNG_SUPPORT"
	
		CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-png --with-png=$PNG_HOME "
                AC_DEFINE(HAVE_PNG, 1, [Png library])
          AC_MSG_RESULT([PNG path set to $PNG_HOME ])
	fi
	if test "$OS_NAME" = "OSF1"; then 
		PNG_LIBRARY_PATHS=" -noso $PNG_LIBRARY_PATHS"
		PNG_LIBRARY="$PNG_LIBRARY -so_archive "
		ZLIB_LIBRARY_PATHS=" -noso $ZLIB_LIBRARY_PATHS"
		ZLIB_LIBRARY="$ZLIB_LIBRARY -so_archive "
		
	fi
        
	# PNG Library Options

	AC_SUBST(CONFIG_PNG)
	AC_SUBST(CONFIG_ZLIB)
	AC_SUBST(PNG_ENABLE)
	AC_SUBST(PNG_HOME)
	AC_SUBST(ZLIB_HOME)
	AC_SUBST(PNG_DEFINES)
	AC_SUBST(ZLIB_LIBRARY)
	AC_SUBST(ZLIB_INCLUDE_PATHS)
	AC_SUBST(ZLIB_LIBRARY_PATHS)
	AC_SUBST(PNG_LIBRARY)
	AC_SUBST(PNG_INCLUDE_PATHS)
	AC_SUBST(PNG_LIBRARY_PATHS)
       
  ])

AC_DEFUN([CONFIGURE_64bit],
  [       
        #
        # Set options for running with 64bit support
        #
        AC_ARG_ENABLE([bits64],
               AC_HELP_STRING(
                      [--enable-bits64],[Compile and load with 64bit option]),
		      bits64_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling: 64bit),
		      bits64_ENABLE="False")

        if test "$bits64_ENABLE" = "True"; then
              WORD_SIZE="-64"
              if test "$PROCESSOR_TYPE" = "i386linux"; then
                 WORD_SIZE="-m32"
              fi


              EXE_SUFFIX="$EXE_SUFFIX""_64bit"

	      CONFIG_OPTIONS="$CONFIG_OPTIONS --enable--bits64 "

              AC_MSG_RESULT(** 64-bit Option Enabled **)

        elif test "$bits64_ENABLE" = "False"; then
              if test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX" -o  "$PROCESSOR_TYPE" = "i386linux"; then
                 WORD_SIZE="-n32"
                 if test "$PROCESSOR_TYPE" = "i386linux"; then
                    WORD_SIZE="-m32"
                 fi
                 EXE_SUFFIX="$EXE_SUFFIX""_32bit"
              fi
        else
        	AC_MSG_WARN("****  Unknown 64bit Enable State: $64bit_ENABLE")
        fi
		# Export 64bit Options
        AC_SUBST(64bit_ENABLE)
        AC_SUBST(WORD_SIZE)
        AC_SUBST(EXE_SUFFIX)
  ])


AC_DEFUN([CONFIGURE_INSTALL],
  [       
        #
        # Set the directory names and paths
        #
        INSTALL_HOME="/usr/gapps/mdg"
        TOP_DIR="`pwd`"

        INSTALL_DIR="src"

        AC_SUBST(INSTALL_HOME)
        AC_SUBST(INSTALL_DIR)
        AC_SUBST(PLATFORM_DIRS)
  ])


AC_DEFUN([CONFIGURE_DIRS],
  [       
        #
        # Set the directory names and paths
        #
	GRIZ_HOME="."
	TOP_DIR="`pwd`"

	BIN_DIRS=""
	LIB_DIRS=""
	EXT_DIRS=""
	SRC_DIRS=""
	MISC_DIRS="Bitmaps HersheyLib/data HersheyLib/fonts HersheyLib/hfonts"
	MAKEFILES=""
        #
        # Setup JPEG
        #
	if test "$CONFIG_JPEG" = "true" ;then
		EXT_DIRS=" ext/JPEG "
		SRC_DIRS=" ext/JPEG "
		MAKEFILES=" $HOST_PREFIX$HOSTDIR$EXE_SUFFIX""/ext/Makefile "
		EXT_DIRS="$EXT_DIRS ext/JPEG "
		EXT_BUILDS=" JPEG "
  	 fi

        #
        # Setup PNG
        #
	 if test "$CONFIG_PNG" = "true"; then
		SRC_DIRS=" $SRC_DIRS ext/PNG ext/ZLIB "
		MAKEFILES="$MAKEFILES $HOST_PREFIX$HOSTDIR$EXE_SUFFIX/ext/PNG/Makefile"
		EXT_DIRS="$EXT_DIRS ext/PNG ext/ZLIB "
		EXT_BUILDS=" ZLIB PNG $EXT_BUILDS " 
	fi


        SRC_DIRS="$SRC_DIRS HersheyLib/src ImageLib"
        OBJS_DIRS=""
        INC_DIRS="include"
	   SRC="src"
        AC_ARG_WITH([install],
               AC_HELP_STRING(
	            [--install-path=[PATH]],
                    [Use given base PATH installing Griz4]),
	            INSTALL_DIRS="-I${withval}/include" &&
	            AC_MSG_RESULT(---->Installing Griz st: $withval)
	           )
 	AC_ARG_WITH([usrbuild],
               AC_HELP_STRING(
	            [--with-usrbuild=[PATH]],
                    [Use given PATH for BUILD DIRECTORY]),
	            BASE_DIR="${withval}" &&
	            AC_MSG_RESULT("Setting Build directory to : $withval")
	           )

        AC_SUBST(GRIZ_HOME)
        AC_SUBST(TOP_DIR)
        AC_SUBST(BIN_DIRS)
        AC_SUBST(LIB_DIRS)
        AC_SUBST(SRC_DIRS)
	AC_SUBST(SRC)
        AC_SUBST(INC_DIRS)
	AC_SUBST(INSTALL_DIRS)
	AC_SUBST(EXT_BUILDS)
	AC_SUBST(EXT_DIRS)
	AC_SUBST(BASE_DIR)
  ])

AC_DEFUN([EXTERNAL_CONFIGURATIONS],
  [
	# We run any other external configuration files here. We need this to run after everything else
        # so the directories would be created and the flags set.

	# PNG

	if test "$CONFIG_ZLIB" = "true"; then
		if test ! -e "$HOSTDIR/ext_lib";then 
			mkdir $TOP_DIR/$HOSTDIR/lib
		fi
		if test ! -e "$HOSTDIR/include";then 
			mkdir $TOP_DIR/$HOSTDIR/include
		fi
		cd $TOP_DIR/$HOSTDIR/ext/ZLIB;
		
		AC_MSG_RESULT("Configuring ZLIB")

		configure --prefix=$TOP_DIR/$HOSTDIR

		cd $TOP_DIR;

		cd $TOP_DIR/$HOSTDIR/ext/PNG;
		
		
		AC_MSG_RESULT("Configuring PNG")
		if test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX"; then
			AC_MSG_RESULT("copying scripts/makefile.sgi")
			cp -f scripts/makefile.sgi Makefile.in
		elif test "$OS_NAME" = "AIX";then
			AC_MSG_RESULT("copying scripts/makefile.aix")
			cp -f scripts/makefile.aix Makefile.in
		elif test "$OS_NAME" = "Linux";then
			AC_MSG_RESULT("copying scripts/makefile.linux")
			cp -f scripts/makefile.linux Makefile.in
		elif test "$OS_NAME" = "OSF1";then
			AC_MSG_RESULT("copying scripts/makefile.dec")
			cp -f scripts/makefile.dec Makefile.in
		elif test "$OS_NAME" = "HPUX"   -o "$OS_NAME" = "HP-UX"; then
			AC_MSG_RESULT("copying scripts/makefile.hpux")
			cp -f scripts/makefile.hpux Makefile.in
		else
			AC_MSG_WARN("Trying standard makefile setup makefile.std")
			cp -f scripts/makefile.std Makefile.in
		fi	
				
		cd $TOP_DIR	
	fi
	   #/ PNG
        # JPEG

	if test "$CONFIG_JPEG" = "true";then
		AC_MSG_RESULT("Configuring JPEG")
		cd $TOP_DIR
		if test ! -e "$HOSTDIR/lib";then 
			mkdir $TOP_DIR/$HOSTDIR/lib
		fi
		if test ! -e "$HOSTDIR/include";then 
			mkdir $TOP_DIR/$HOSTDIR/include
		fi
		cd $HOSTDIR/ext/JPEG;
	
		if test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX"; then
			configure MY_CC=cc --prefix=$TOP_DIR/$HOSTDIR 
		else
			configure --prefix=$TOP_DIR/$HOSTDIR 
	     fi
		cd $TOP_DIR	
	fi

        #/ JPEG

	AC_SUBST(PNG_SUBDIRS)
  ])


AC_DEFUN([CONFIGURE_BUILDDIRS],
  [
        #
        # Create the operating-system dependent build directories.
        #
	if test "$BASE_DIR" != ""; then
	    HOSTDIR="GRIZ4-$BASE_DIR-$HOSTNAME"
	else
	    HOSTDIR="$HOST_PREFIX""$HOSTDIR""$EXE_SUFFIX"
        fi

        if test -d "$HOSTDIR"; then
           rm -rf "$HOSTDIR-old"
           mv -f "$HOSTDIR" "$HOSTDIR-old"
        fi     

        if test ! -d "$HOSTDIR"; then
	   	AC_MSG_RESULT("Creating directory: $HOSTDIR")
	   	mkdir "$HOSTDIR"
	   	if test ! -d "$HOSTDIR"; then
	       AC_MSG_WARN("****  Error creating directory:  $HOSTDIR")
	   	fi
        fi
        for dir in $BIN_DIRS $LIB_DIRS $INC_DIRS $SRC_DIRS $SRC $EXT_DIRS $OBJS_DIRS $MISC_DIRS; do
                PATH_NAME=$HOSTDIR
	        set PATH_NAME_TEMP `echo ":$dir" | sed -ne 's/^:\//#/;s/^://;s/\// /g;s/^#/\//;p'`
	        shift
	        for dir_name
	        do
		        PATH_NAME="$PATH_NAME/$dir_name"

		        case "$PATH_NAME" in
		        	-* ) PATH_NAME=./$PATH_NAME ;;
		        esac
		        if test ! -d "$PATH_NAME"; then
			        AC_MSG_RESULT("Creating directory:  $PATH_NAME")
			        mkdir "$PATH_NAME"
			        if test ! -d "$PATH_NAME"; then
			        	AC_MSG_WARN("****  Error creating directory:  $PATH_NAME")
			        fi
		        fi
	        done
        done
	   ROOT_DIR=$TOP_DIR/$HOSTDIR

	   if test "$CONFIG_JPEG" = "true" -o "$CONFIG_PNG" = "true"; then
		cd $ROOT_DIR/ext
		ln -sf $TOP_DIR/ext/Makefile.in
		cd $TOP_DIR
	   fi


	cd $HOSTDIR
        # Make a temporary src and include directory at the top level
        cd $TOP_DIR
        if test ! -e "src"; then
                mkdir src
        fi

        if test ! -e "include"; then
                mkdir include
        fi
        cd src
        ln -sf $TOP_DIR/*.c .     > /dev/null 2>&1


        cd $ROOT_DIR

        cp ../Makefile.Driver Makefile.in
        cp ../Makefile.Library src/Makefile.in

        if test ! -e "src"; then
                mkdir src
        fi

        if test ! -e "include"; then
                mkdir include
        fi

        if test ! -e "Doc"; then
                mkdir Doc
        fi


        # Set up links for source files
        cd $ROOT_DIR
        cd src
        ln -sf $TOP_DIR/*.c .     > /dev/null 2>&1
        ln -sf $TOP_DIR/griz.in . > /dev/null 2>&1
        if test ! -e "GL"; then
                mkdir GL
        fi

	
        cd ..
        cd include
        ln -sf $TOP_DIR/*.h .     > /dev/null 2>&1
        if test ! -e "GL"; then
                mkdir GL
        fi
	cd GL
        ln -sf $TOP_DIR/osmesa.h . > /dev/null 2>&1

        # Set up links for Doc files
        cd $ROOT_DIR
        cd Doc
        ln -sf $TOP_DIR/Doc/* .     > /dev/null 2>&1


        # Set up links for include files
        cd $ROOT_DIR
        cd include
        ln -sf $TOP_DIR/*.h . > /dev/null 2>&1

        cd $ROOT_DIR
        for dir in $INC_DIRS; do
        	dir_tmp="$TOP_DIR/$dir"
        	cd "$dir_tmp"

        	for file in *.h; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
		        	ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
		        fi
        	done
        	for file in *.h.in; do
        		if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
        			ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
        		fi
        	done
        	cd $TOP_DIR
        done

        # Setup links for sources files, makefiles, and include files
        for dir in $SRC_DIRS $SRC; do
	        cd $ROOT_DIR/$dir

	        #   Create a link to the matching include directory for debugging purposes.
	        if test ! -e "include"; then
	        	incl_dir=`echo "$dir" | sed 's/src/include/'`
	        	if  test -e "$TOP_DIR/$HOSTDIR/$incl_dir"; then
	        		ln -s $TOP_DIR/$HOSTDIR/$incl_dir include > /dev/null 2>&1
	        	fi
        	fi

        	if test ! -e "Makefile"; then
        		ln -s $TOP_DIR/$dir/Makefile . > /dev/null 2>&1
	        fi
	        if test ! -e "Makefile.in"; then
	        	ln -s $TOP_DIR/$dir/Makefile.in . > /dev/null 2>&1
	        fi

	        dir_tmp="$TOP_DIR/$dir"
	        cd "$dir_tmp"
	        for file in *.in; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
	        for file in *.c; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
	        for file in *.f; do
        		if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
	        for file in *.F; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
	        for file in *.f90; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
	        for file in *.h; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
		        	ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
		        fi
	        done
	        for file in *_List; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
		for file in *.doc; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done
        done

        # Misc directories do not consist of code - this is stuff like the Bitmaps, Fonts, etc.
        cd $TOP_DIR
        PATH_NAME="$HOSTDIR"
        for dir in $MISC_DIRS; do
	        cd $PATH_NAME/$dir

                # Link to all files in these directories
	        dir_tmp="$TOP_DIR/$dir"
	        cd "$dir_tmp"
	        for file in *; do
	        	if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
	        		ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
	        	fi
	        done

	        cd $TOP_DIR
        done

	cd $TOP_DIR

 	for dir in $EXT_DIRS; do
		dir_tmp="$TOP_DIR/$dir"
          cd $dir_tmp
		for file in *; do
			if test ! -e "$TOP_DIR/$HOSTDIR/$dir/$file" && test -e $dir_tmp/$file; then
				ln -s $dir_tmp/$file $TOP_DIR/$HOSTDIR/$dir/$file > /dev/null 2>&1
			fi
		done
	done

     AC_SUBST(ROOT_DIR)
	
	cd $ROOT_DIR
	echo $CONFIG_CMD > CONFIGURE_OPTIONS
  ])


AC_DEFUN([CONFIGURE_TASKS],
  [
        #
        # Define the number of tasks based on the system type
        #
        TASKS=1;
        if test -e "/usr/sbin/psrinfo"; then
        	TASKS=4
					#"`/usr/sbin/psrinfo -n | cut -d '=' -f2`"
        else
					TASKS=8
        fi
        AC_SUBST(TASKS)
  ])


AC_DEFUN([CONFIGURE_GMAKE],
  [
        #
        # Determine if 'gmake' is available for building the code.
        #
       GMAKE_HOME=""

        AC_ARG_WITH([gmake],
               AC_HELP_STRING(
	            [--with-gmake=[PATH]],
                    [Use given base PATH for gmake tool]),
	            GMAKE_HOME="${withval}/" &&
		    CONFIG_OPTIONS="$CONFIG_OPTIONS --with-gmake=$GMAKE_HOME " &&
	            AC_MSG_RESULT("Using GMAKE Path : $withval")
	           )

        # gmake path
        AC_SUBST(GMAKE_HOME)
  ])


AC_DEFUN([CONFIGURE_GRIZ_LIBRARIES],
  [
        #
        # Define the minimal set of system libraries required
        #

        if test "$MESA_ENABLE" = "True" -o "$OSMESA_ENABLE" = "True"; then
        # If using Mesa OpenGL instead of a native OpenGL, try this:
                SYSLIBS=" -lm"
        else
                SYSLIBS=" -lm"
        fi

        if test "$OS_NAME" = "OSF1"; then
                SYSLIBS="$SYSLIBS -lXext -lXmu "
                LDLIBPATH="/usr/local/lib"

        elif test "$OS_NAME" = "Linux"; then
                SYSLIBS="$SYSLIBS "
                LDLIBPATH="/usr/X11R6/lib"

        elif test "$OS_NAME" = "AIX"; then
                SYSLIBS="$SYSLIBS "
                LDLIBPATH="/usr/local/lib"

        elif test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX"; then
                SYSLIBS="$SYSLIBS -lPW"
                LDLIBPATH="/mdg/lib/irix:/usr/local/lib:/usr/gapps/mdg/bin/grizlib/SGI"

        elif test "$OS_NAME" = "HPUX"   -o "$OS_NAME" = "HP-UX"; then
                SYSLIBS="$SYSLIBS -lXext -lXmu"

        elif test "$OS_NAME" = "SunOS"; then
                SYSLIBS="$SYSLIBS -lXext -lXmu -L/usr/local/lib/gcc-lib/sparc-sun-solaris2.8/2.95.2 -lgcc"
                LDLIBPATH="/mdg/lib/solaris:/usr/openwin/lib:/usr/dt/lib/:/usr/lib:/usr/local/lib:/usr/gapps/mdg/bin/grizlib/SUN"
                MACHINE_TYPE="`(uname -i | grep -v Blade) 2> /dev/null`"
                if test "$MACHINE_TYPE" = ""; then
                   LDLIBPATH="/mdg/lib/solarisBlade:$LDLIBPATH"
                fi
        fi

	LDLIBPATH="`echo $LDLIBPATH | sed -e 's|-O2||'`"
        #
        # Construct the list of MDG Libraries that will be loaded
        #

        MDGLIBS="$HERSHEY_LIBRARY_PATHS $HERSHEY_LIBRARY"
        MDGLIBS="$MDGLIBS $IMAGE_LIBRARY_PATHS $IMAGE_LIBRARY"

        EXTRALIBS="$JPEG_LIBRARY_PATHS $JPEG_LIBRARY"
        EXTRALIBS="$EXTRALIBS $PNG_LIBRARY_PATHS $PNG_LIBRARY $ZLIB_LIBRARY_PATHS $ZLIB_LIBRARY "

        GRIZLIBS="$EXTRALIBS $MDGLIBS $SYSLIBS"
	GRIZLIBS="`echo $GRIZLIBS | sed -e 's|-O2||'`"

        AC_SUBST(GRIZLIBS)
        AC_SUBST(LDLIBPATH)

        AC_MSG_RESULT([**************************************])
        AC_MSG_RESULT([---->Loading with Libraries: $GRIZLIBS])
        AC_MSG_RESULT([**************************************])
  ])


AC_DEFUN([CONFIGURE_COMPILER_FLAGS],
  [
        #
        # Define the compile search order
        #

        CC_DEFINES="$CC_DEFINES \
                $EXO_DEFINES \
                $IMAGE_DEFINES \
                $JPEG_DEFINES \
                $MESA_DEFINES \
                $MILI_DEFINES \
                $NETCDF_DEFINES \
                $PAPI_DEFINES \
                $PNG_DEFINES"
	CC_DEFINES="`echo $CC_DEFINES | sed -e 's|-O2||'`"

        CC_INCLUDE_PATHS="-I/usr/local/include \
                -I/usr/include \
                -I../include "
		if test ! "$JPEG_INCLUDE_PATH" = "-I../include"; then
			CC_INCLUDE_PATHS="$CC_INCLUDE_PATHS $JPEG_INCLUDE_PATH "
		fi
		if test ! "$PNG_INCLUDE_PATH" = "-I../include"; then
			CC_INCLUDE_PATHS="$CC_INCLUDE_PATHS $PNG_INCLUDE_PATH "
		fi
        CC_INCLUDE_PATHS="$CC_INCLUDE_PATHS $MILI_INCLUDE_PATHS \
                $MESA_INCLUDE_PATHS \
                $BITMAPS_INCLUDE_PATHS \
                $EXO_INCLUDE_PATHS \
                $HERSHEY_INCLUDE_PATHS \
                $IMAGE_INCLUDE_PATHS \
                $NETCDF_INCLUDE_PATHS \
                $PAPI_INCLUDE_PATHS \
                $X11_INCLUDE_PATHS"
	
	CC_INCLUDE_PATHS="`echo $CC_INCLUDE_PATHS | sed -e 's|-O2||'`"

        AC_MSG_RESULT(**************************************)
        AC_MSG_RESULT(---->Include Paths: $CC_INCLUDE_PATHS )
        AC_MSG_RESULT(**************************************)

        AC_SUBST(CC_DEFINES)
        AC_SUBST(CC_INCLUDE_PATHS)
  ])


AC_DEFUN([CONFIGURE_COMPILER],
  [
        # First Set non-machine specific defaults

# SGI compile flags; add -mips2 if running on R4000, etc.
#
# Flags used on other platforms:
#
# Sun Solaris           -O2 -Xa
# DEC OSF               -O2
# Cray UNICOS           -h nostdc -Gn
# HPUX                  -Ae
#

# Minimal library set (not complete! See below.)
# LIBS = -lXm -lXt -lX11 -lGLw -lGLU -lGL -lm
#
# For these systems, add libraries...
# SGI                   -PW
# Sun Solaris           -Xext -Xmu
# DEC OSF               -Xext -Xmu
# Cray UNICOS           -Xmu
# HPUX                  -Xext -Xmu
#
# Additional notes on libraries:
# Depending on your platform and environment, it may well be necessary 
# to add appropriate "-I" and/or "-L" switches so all the include 
# files and libraries can be found. Also, order of specification is 
# important.
#

        SHELL="/bin/sh"
        SHELL_ARGS="-ec"

        if test "$OS_NAME" = "OSF1"; then
                AC_CHECK_PROGS(MY_CC, cc c89)
                CC_FLAGS_DEBUG="-g"
                CC_FLAGS_OPT="-g4"
                CC_DEPEND="-M" 

	        AR="ar"
                RANLIB="ranlib"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.dec""$EXE_SUFFIX"

        elif test "$OS_NAME" = "Linux"; then
                AC_CHECK_PROGS(MY_CC, gcc cc)

                SHELL_ARGS=" "

                CC_FLAGS_DEBUG="-g"
                CC_FLAGS_OPT="-g -O3"
                CC_DEPEND="-M" 

	        AR="ar"
                RANLIB="ranlib"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.linux""$EXE_SUFFIX"

                if test  "$SYS_TYPE" = "chaos_3_x86_64_ib" -o "$SYS_TYPE" = "chaos_3_x86_64"; then
                   MY_CC="pathcc"
                   MY_FC="pathf90"
                fi

        elif test "$OS_NAME" = "AIX"; then
                AC_CHECK_PROGS(MY_CC, cc xlc xlr_r)
                CC_FLAGS_DEBUG="-g -bnoobjreorder -qfullpath -DAIX -qMAXMEM=8192 -bmaxdata:0x40000000"
                CC_FLAGS_OPT="-g -O3 -bnoobjreorder -qfullpath -DAIX -qMAXMEM=8192 -bmaxdata:0x40000000"
                CC_DEPEND="-M" 

	        AR="ar"
                RANLIB="ranlib"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.aix""$EXE_SUFFIX"

       elif test "$OS_NAME" = "HPUX" -o "$OS_NAME" = "HP-UX"; then
                AC_CHECK_PROGS(MY_CC, cc)
                CC_FLAGS_DEBUG="-g -Ae"
                CC_FLAGS_OPT="-g -O3 -Ae"
                CC_DEPEND="-M" 

	        AR="ar"
                RANLIB="ranlib"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.hpux""$EXE_SUFFIX"

       elif test "$OS_NAME" = "SunOS"; then
                MACHINE_TYPE="`(uname -i | grep -v Blade) 2> /dev/null`"
                if test "$MACHINE_TYPE" = ""; then
                   MACHINE_TYPE="Blade"
                fi

                AC_MSG_RESULT("Machine_type=$MACHINE_TYPE") 
                AC_CHECK_PROGS(MY_CC, cc)
                CC_FLAGS_DEBUG="-g -Xa"
                CC_FLAGS_OPT="-xO3 -Xa"
                CC_DEPEND="-M" 

	        AR="ar"
                RANLIB="ranlib"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.sol""$EXE_SUFFIX"
                if test "$MACHINE_TYPE" = "Blade"; then
                   GRIZ_EXE="griz4s.sol.blade""$EXE_SUFFIX"
                fi

       elif test "$OS_NAME" = "IRIX64" -o "$OS_NAME" = "IRIX"; then
                AC_CHECK_PROGS(MY_CC, cc)
                CC_FLAGS_DEBUG="-g $WORD_SIZE -ansi -common"
                CC_FLAGS_OPT="-g -O3 $WORD_SIZE -ansi -common"
                CC_DEPEND="-M" 

                CC_DEFINES="-DIRIX -DO2000" 

	        AR="ar"
                RANLIB="echo"
	        AR_FLAGS="-r"

                GRIZ_EXE="griz4s.irix""$EXE_SUFFIX"
        fi

        AC_PROG_CC
        AC_PROG_INSTALL
        AC_PROG_AWK
        AC_PROG_CPP
        AC_PROG_LN_S
        AC_PROG_MAKE_SET
        AC_SUBST(AR)
        AC_SUBST(AR_FLAGS)
        AC_SUBST(RANLIB)
        AC_SUBST(MY_CC)
        AC_SUBST(MY_FC)

        AC_SUBST(SHELL)        
        AC_SUBST(SHELL_ARGS)

        AC_SUBST(GRIZ_EXE)
        AC_SUBST(CC_DEPEND)
        AC_SUBST(CC_FLAGS_DEBUG)
        AC_SUBST(CC_FLAGS_OPT)
        AC_SUBST(LDFLAGS_EXTRA)
  ])


AC_DEFUN([CONFIGURE_COMPILER_IMAGELIB],
  [
        #
        # Define the compile search order for the Image Library
        #

        CC_DEFINES_IMAGELIB="$CC_DEFINES "

	CC_DEFINES_IMAGELIB="`echo $CC_DEFINES_IMAGELIB | sed -e 's|-O2||'`"

        CC_INCLUDE_PATHS_IMAGELIB=" $CC_INCLUDE_PATHS "

        CC_INCLUDE_PATHS_IMAGELIB=" -I/usr/local/include -I/usr/include -I../include "

	CC_INCLUDE_PATHS_IMAGELIB="`echo $CC_INCLUDE_PATHS_IMAGELIB | sed -e 's|-O2||'`"


        AC_SUBST(CC_DEFINES_IMAGELIB)
        AC_SUBST(CC_INCLUDE_PATHS_IMAGELIB)
  ])

AC_DEFUN([CONFIGURE_COMPILER_GRIZ],
  [
        #
        # Define the compile search order
        #

	CC_DEFINES_GRIZ=" $CC_DEFINES "
	CC_DEFINES_GRIZ="`echo $CC_DEFINES_GRIZ | sed -e 's|-O2||'`"
        
	CC_INCLUDE_PATHS_GRIZ=" $CC_INCLUDE_PATHS "
	CC_INCLUDE_PATHS_GRIZ="`echo $CC_INCLUDE_PATHS_GRIZ | sed -e 's|-O2||'`"
        
	AC_SUBST(CC_DEFINES_GRIZ)
	AC_SUBST(CC_INCLUDE_PATHS_GRIZ)
  ])


AC_DEFUN([CONFIGURE_COMPILER_HERSHEYLIB],
  [
        #
        # Define the compile search order
        #
        CC_INCLUDE_PATHS_HERSHEYLIB="$MESA_INCLUDE_PATHS -I/usr/local/include -I/usr/include -I../include "

        HERSHEY_FONTLIB=/usr/local/lib/hershey
        HERSHEY_FONTLIB=$TOP_DIR/$HOSTDIR/HersheyLib/hfonts

        if test "$HERSHEY_FONTS_SET" = "True"; then
                HERSHEY_FONTLIB=$HERSHEY_FONTS_PATH
        fi

        CC_DEFINES_HERSHEYLIB=-DFONTLIB=\\\"$HERSHEY_FONTLIB\\\"

        if test "$OS_NAME" = "OSF1"; then
                CC_DEFINES_HERSHEYLIB="$CC_DEFINES_HERSHEYLIB -DSGI"
        fi

        AC_SUBST(HERSHEY_FONTLIB)
        AC_SUBST(CC_DEFINES_HERSHEYLIB)
        AC_SUBST(CC_INCLUDE_PATHS_HERSHEYLIB)
  ])


AC_DEFUN([CONFIGURE_KLOCWORK],
  [       
        #
        # Set options for KW
        #
	KW_ENABLE="False" 
	KW_HOME="/var/tmp/klocwork/project_sources/griz/"

        AC_ARG_ENABLE([kw],
               AC_HELP_STRING(
                 [--enable-kw],[Compile and load with KlocWorks]),
		KW_ENABLE="True" &&  AC_MSG_RESULT(---->Enabling: KlocWorks),
		KW_ENABLE="False")

        AC_ARG_WITH([kw],
               AC_HELP_STRING(
	            [--with-kw=[PATH]],
                    [Use given base PATH for KW trace file]),
	            KW_HOME="${withval}" &&
	            AC_MSG_RESULT("Using KlocWorks Path : $withval")
	           )

        if test "$KW_ENABLE" = "True"; then
           CONFIG_OPTIONS="$CONFIG_OPTIONS --enable-kw ssss --with-kw=$KW_HOME "	
           EXE_SUFFIX="$EXE_SUFFIX""_KlocWorks"
        fi

        # KW Options
        AC_SUBST(KW_HOME)
        AC_SUBST(KW_ENABLE)
  ])

# End of acinclude.m4