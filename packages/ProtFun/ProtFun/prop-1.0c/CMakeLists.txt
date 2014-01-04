
PROJECT( PROP )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

SET( PROP_BASENAME prop )
SET( PROP_VERSION 1.0c )
SET( PROP_FULLNAME ${PROP_BASENAME}-${PROP_VERSION}  CACHE STRING "Full name of program (to be used during install)" )

SET( LIBRARY_DIR /usr/share/ )
SET( PROP_LIBRARY_DIR ${LIBRARY_DIR}/${PROP_FULLNAME}/ CACHE PATH "Path to install the PROP library" )

MESSAGE( STATUS "Will install prop libraries to: ${PROP_LIBRARY_DIR}" )
MESSAGE( STATUS "Will install prop bash-script to: ${CMAKE_INSTALL_PREFIX}" )

SET( WWW_OUT_DIR /usr/opt/www/pub/CBS/services/ProP-1.0/tmp/ )
MESSAGE( STATUS "Using WWWOUT dir as: ${WWW_OUT_DIR}" )

SET( TEMP_DIR /var/tmp/ CACHE PATH "Temporary directory for ProP to use, should be world writeable." )
MESSAGE( STATUS "Using temporary directory as: ${TEMP_DIR}" )

FIND_PROGRAM( GAWK_DIR gawk )
FIND_PROGRAM( PPM2GIF_DIR ppm2gif )
FIND_PROGRAM( GNUPLOT_DIR gnuplot )

if( NOT PPM2GIF_DIR )
    FIND_PROGRAM( PPM2GIF_DIR ppmtogif )
endif()

if( NOT GAWK_DIR )
    MESSAGE( FATAL_ERROR "gawk not found!" )
endif()

if( NOT PPM2GIF_DIR )
    MESSAGE( STATUS "ppm2gif/ppmtogif not found! Graphics mode will not be available" )
endif()

if( NOT GNUPLOT_DIR )
    MESSAGE( STATUS "gnuplot not found! Graphics mode will not be available" )
endif()

MESSAGE( STATUS "Configuring prop template" )
CONFIGURE_FILE( ./templates/prop.template ./prop @ONLY )

SET( MAN_DIR "/usr/share/man/man1/" CACHE PATH "Manual installation directory" )

OPTION( INSTALL_MANUAL "Install the manual for ProP" )

INSTALL( DIRECTORY ./etc ./syn ./test ./tmp DESTINATION ${PROP_LIBRARY_DIR} )

FILE( MAKE_DIRECTORY ${PROP_LIBRARY_DIR}/bin ${PROP_LIBRARY_DIR}/how )

INSTALL( PROGRAMS ./prop DESTINATION ${CMAKE_INSTALL_PREFIX} )

INSTALL( PROGRAMS ./bin/fasta2dir ./bin/gethow ./bin/getseqstr ./bin/in2how ./bin/in2how+fasta.awk ./bin/mkgraph ./bin/predict_S ./bin/predict_T ./bin/predict_Y ./bin/prop_furin ./bin/prop_total DESTINATION ${PROP_LIBRARY_DIR}/bin )
INSTALL( PROGRAMS ./how/how98_${CMAKE_SYSTEM_NAME} DESTINATION ${PROP_LIBRARY_DIR}/how RENAME how98 )

if( INSTALL_MANUAL )
    INSTALL( FILES ./man/prop.1 DESTINATION ${MAN_DIR} )
endif()

