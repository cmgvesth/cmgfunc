#########################################################################
# -- SEG/NSEG/PSEG CMAKE INSTALL --                                     #
#########################################################################

project( SEG_NSEG_PSEG )
cmake_minimum_required( VERSION 2.6 )

# Allow optional install of man pages
option( INSTALL_MANUAL "Install the Man-Pages for SEG" )



# Add Targets
add_executable( seg ./src/seg/seg.c ./src/seg/lnfac.h ./src/seg/genwin.h ./src/seg/genwin.c )
add_executable( nseg ./src/nseg/nseg.c ./src/nseg/lnfac.h ./src/nseg/genwin.h ./src/nseg/genwin.c )
add_executable( nmerge ./src/nseg/nmerge.c ./src/nseg/genwin.h ./src/nseg/genwin.c )
add_executable( pseg ./src/pseg/pseg.c ./src/pseg/lnfac.h ./src/pseg/pgenwin.h ./src/pseg/pgenwin.c )
add_executable( pmerge ./src/pseg/pmerge.c ./src/pseg/pgenwin.h ./src/pseg/pgenwin.c )

target_link_libraries( seg m )
target_link_libraries( nseg m )
target_link_libraries( nmerge m )
target_link_libraries( pseg m )
target_link_libraries( pmerge m )

message( STATUS "Will install targets to $(CMAKE_INSTALL_PREFIX)" )
message( STATUS "Targets are: seg nseg nmerge pseg pmerge" )

# Install Executables
install( TARGETS seg nseg nmerge pseg pmerge DESTINATION ${CMAKE_INSTALL_PREFIX} )
