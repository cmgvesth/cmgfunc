
PROJECT( APE )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

SET( APE_BASENAME ape )
SET( APE_VERSION 1.0 )
SET( APE_FULLNAME ${APE_BASENAME}-${APE_VERSION}  CACHE STRING "Full name of program (to be used during install)" )

SET( LIBRARY_DIR /usr/share/ )
SET( APE_LIBRARY_DIR ${LIBRARY_DIR}/${APE_FULLNAME}/ CACHE PATH "Path to install the APE library" )

MESSAGE( STATUS "Will install ape libraries to: ${APE_LIBRARY_DIR}" )
MESSAGE( STATUS "Will install ape bin to: ${CMAKE_INSTALL_PREFIX}" )

SET( WWW_OUT_DIR /usr/opt/www/pub/CBS/services/ape-1.0/tmp/ )
MESSAGE( STATUS "Using WWWOUT dir as: ${WWW_OUT_DIR}" )

SET( TEMP_DIR /var/tmp/ CACHE PATH "Temporary directory for ProP to use, should be world writeable." )

MESSAGE( STATUS "Using temporary directory as: ${TEMP_DIR}" )
MESSAGE( STATUS "Configuring ape template" )
CONFIGURE_FILE( ./ape.in ./ape @ONLY )

SET( MAN_DIR "/usr/share/man/man1/" CACHE PATH "Manual installation directory" )

OPTION( INSTALL_MANUAL "Install the manual for ape" )

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${APE_LIBRARY_DIR}/bin)
INSTALL( DIRECTORY ./clp ./tmp ./etc ./pred ./disp ./graphics ./test DESTINATION ${APE_LIBRARY_DIR} )
INSTALL( PROGRAMS ./bin/loadenv.awk ./bin/nnhowplayer6_Linux.i686 ./bin/nnhowplayer6_Linux.ia64 ./bin/nnhowplayer6_Linux.x86_64 ./bin/paste.awk ./bin/pp.awk ./bin/seq2seq.awk DESTINATION ${APE_LIBRARY_DIR}/bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ )

INSTALL( PROGRAMS ./ape DESTINATION ${CMAKE_INSTALL_PREFIX} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ )

if( INSTALL_MANUAL )
    INSTALL( FILES ./doc/netphos.1 DESTINATION ${MAN_DIR} )
endif()

