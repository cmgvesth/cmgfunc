
PROJECT( TMHMM )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

SET( TMHMM_BASENAME tmhmm )
SET( TMHMM_VERSION 2.0c )
SET( TMHMM_FULLNAME ${TMHMM_BASENAME}-${TMHMM_VERSION}  CACHE STRING "Full name of program (to be used during install)" )

SET( LIBRARY_DIR /usr/share/ )
SET( TMHMM_LIBRARY_DIR ${LIBRARY_DIR}/${TMHMM_FULLNAME}/ CACHE PATH "Path to install the PROP library" )

MESSAGE( STATUS "Will install tmhmm libraries to: ${PROP_LIBRARY_DIR}" )
MESSAGE( STATUS "Will install tmhmm bin to: ${CMAKE_INSTALL_PREFIX}" )

SET( WWW_OUT_DIR /usr/opt/www/pub/CBS/services/ProP-1.0/tmp/ )
MESSAGE( STATUS "Using WWWOUT dir as: ${WWW_OUT_DIR}" )

SET( TEMP_DIR /var/tmp/ CACHE PATH "Temporary directory for TMHMM to use, should be world writeable." )

MESSAGE( STATUS "Using temporary directory as: ${TEMP_DIR}" )
MESSAGE( STATUS "Configuring tmhmm template" )

find_program( GNUPLOT_DIR gnuplot )
find_program( PPM2GIF_DIR ppm2gif )
find_program( GREP_DIR grep )
find_program( PERL_DIR perl )
find_program( GS_DIR gs )

if( NOT GS_DIR )
    message( ERROR "Program 'gs' not found!")
endif()

if( NOT PERL_DIR )
    message( ERROR "Program 'perl' not found!")
endif()

if( NOT GREP_DIR )
    message( ERROR "Program 'grep' not found!")
endif()

if( NOT GNUPLOT_DIR )
    message( ERROR "Program 'gnuplot' not found. Graphical mode will not be available")
endif()

if( NOT PPM2GIF_DIR )
    find_program( PPM2GIF_DIR ppmtogif )
    if( NOT PPM2GIF_DIR )
        message( ERROR  "Program 'ppm2gif/ppmtogif' not found. Graphs will not be saved as GIF files")
    endif()
endif()

CONFIGURE_FILE( ./templates/tmhmm.in ./bin/tmhmm @ONLY )
CONFIGURE_FILE( ./templates/tmhmmformat.in ./bin/tmhmmformat.pl @ONLY )

SET( MAN_DIR "/usr/share/man/man1/" CACHE PATH "Manual installation directory" )

OPTION( INSTALL_MANUAL "Install the manual for tmhmm" )

INSTALL( DIRECTORY ./bin ./lib DESTINATION ${TMHMM_LIBRARY_DIR} )
INSTALL( PROGRAMS ./bin/tmhmmformat.pl ./bin/decodeanhmm.Linux_i386 ./bin/decodeanhmm.Linux_i686 ./bin/decodeanhmm.Linux_ia64 ./bin/decodeanhmm.Linux_x86_64 DESTINATION ${TMHMM_LIBRARY_DIR}/bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ)

INSTALL( PROGRAMS ./bin/tmhmm DESTINATION ${CMAKE_INSTALL_PREFIX} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ )

if( INSTALL_MANUAL )
    INSTALL( FILES ./TMHMM2.0.html DESTINATION ${MAN_DIR} )
endif()

