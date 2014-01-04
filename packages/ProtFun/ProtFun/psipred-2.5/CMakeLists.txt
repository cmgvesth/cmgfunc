# CMAKE FOR PSIPRED

PROJECT( PSIPRED )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

# There is no manual for this source build...
ADD_EXECUTABLE( ./bin/psipred ./src/sspred_avpred.c ./src/ssdefs.h ./src/sspred_net.h )
ADD_EXECUTABLE( ./bin/psipass2 ./src/sspred_hmulti.c ./src/ssdefs.h ./src/sspred_net.h )
ADD_EXECUTABLE( ./bin/seq2mtx ./src/seq2mtx.c )
ADD_EXECUTABLE( ./bin/pfilt ./src/pfilt.c )

TARGET_LINK_LIBRARIES( ./bin/psipred m )
TARGET_LINK_LIBRARIES( ./bin/psipass2 m )
TARGET_LINK_LIBRARIES( ./bin/seq2mtx m )
TARGET_LINK_LIBRARIES( ./bin/pfilt m )

MESSAGE( STATUS "Successfully configured make for targets: psipred psipass2 seq2mtx pfilt" )

# Basic Settings
SET( PSIPRED_VERSION 2.5 )
SET( PSIPRED_BASENAME psipred)
SET( PSIPRED_FULLNAME "${PSIPRED_BASENAME}-${PSIPRED_VERSION}" CACHE STRING "Full name of program (to be used during install)" )

# Directory Settings
SET( LIBRARY_DIRECTORY /usr/share/ )
SET( PSIPRED_LIBRARY_DIR ${LIBRARY_DIRECTORY}/${PSIPRED_FULLNAME} CACHE PATH "Path of PSIPRED library install" )
SET( NCBI_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Path to your NCBI-tools install (blastpgp, etc)" )

MESSAGE( STATUS "Will install psipred library to ${PSIPRED_LIBRARY_DIR}" )
SET( EXEC_DIR ${PSIPRED_LIBRARY_DIR}/bin )
SET( DATA_DIR ${PSIPRED_LIBRARY_DIR}/data )

FIND_PROGRAM( NCBI_BLAST_DIR blastpgp PATHS ${NCBI_DIR})

if( NOT NCBI_BLAST_DIR )
	MESSAGE( FATAL_ERROR "Could not locate NCBI-Blast install directory." )
endif()

MESSAGE( STATUS "Configuring templates for psipred." )
configure_file( ./templates/runpsipred.template ./runpsipred @ONLY )
configure_file( ./templates/runpsipred_single.template ./runpsipred_single @ONLY )

INSTALL( DIRECTORY ./bin ./data DESTINATION ${PSIPRED_LIBRARY_DIR} )
INSTALL( PROGRAMS ./runpsipred ./runpsipred_single DESTINATION ${CMAKE_INSTALL_PREFIX} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ )

