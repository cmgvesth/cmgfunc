#########################################################################
# -- SIGNALP 4.1 CMAKE INSTALL --                                       #
#########################################################################
# Allows custom installation settings for program signalp (version 4.1) #
# SignalP performs finds Signal Peptide Cleavage sites in the input.    #
# Optionally requests gnuplot and ppm2gif/ppmtogif for graphical data.  #
#########################################################################
project( SIGNALP )
cmake_minimum_required( VERSION 2.6 )

# Allow optional install of man pages
option( INSTALL_MANUAL "Install the Man-Pages for SignalP" )

# Basic Settings
SET( SIGNALP_VERSION 4.1 )
SET( SIGNALP_BASENAME signalp )
SET( SIGNALP_FULLNAME "${SIGNALP_BASENAME}-${SIGNALP_VERSION}" CACHE STRING "Full name of program (to be used during install)" )

# Directory Settings
SET( LIBRARY_DIRECTORY /usr/share/ )
SET( SIGNALP_LIBRARY_DIR ${LIBRARY_DIRECTORY}/${SIGNALP_FULLNAME} CACHE PATH "Path of SIGNALP library install" )
SET( SIGNALP_TEMP_DIR /var/tmp CACHE PATH "Path where SignalP will store temp files")
SET( SIGNALP_WWW_OUT_DIR /usr/opt/www/pub/CBS CACHE PATH "Path for SignalP WWW output")
SET( SIGNALP_MAN_DIRECTORY /usr/share/man/man1/)

# Check optional prerequisites
find_program( GNUPLOT_DIR gnuplot )
find_program( PPM2GIF_DIR ppm2gif )

if( NOT GNUPLOT_DIR )
    message(  "Program 'gnuplot' not found. Graphical mode will not be available")
endif()

if( NOT PPM2GIF_DIR )
    find_program( PPMTOGIF_DIR ppmtogif )
    if( NOT PPM2GIF_DIR )
        message(  "Program 'ppm2gif/ppmtogif' not found. Graphs will not be saved as GIF files")
    endif()
endif()

# PRINT SOME STATUS MESSAGES FOR DEBUGGING SAKE 
message( STATUS "Configuration of project: ${CMAKE_PROJECT_NAME}" )
message(STATUS "Full Name: ${SIGNALP_FULLNAME}")
message( STATUS "Will save library to: ${SIGNALP_LIBRARY_DIR}" )

# Configure the source file
if( NOT EXISTS ./src/signalp.in )
    message( FATAL_ERROR "Could not find required file: src/signalp.in" )
endif()

configure_file( ./src/signalp.in ./${SIGNALP_FULLNAME} @ONLY )

# Install Commands
install( DIRECTORY ./bin ./lib ./syn ./test DESTINATION ${SIGNALP_LIBRARY_DIR} )
install( PROGRAMS "./${SIGNALP_FULLNAME}" DESTINATION ${CMAKE_INSTALL_PREFIX} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ RENAME ${SIGNALP_BASENAME})

if( INSTALL_MANUAL )
    install( FILES ./signalp.1 DESTINATION ${SIGNALP_MAN_DIRECTORY} RENAME ${SIGNALP_BASENAME}.1 )
endif()
