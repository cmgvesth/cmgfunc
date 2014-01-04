# Install script for directory: /home/cmgfunc/ProtFun/ProtFun/ape-1.0

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/cmgfunc/ProtFun/")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/ape-1.0/clp;/home/cmgfunc/ProtFun/ape-1.0/tmp;/home/cmgfunc/ProtFun/ape-1.0/etc;/home/cmgfunc/ProtFun/ape-1.0/pred;/home/cmgfunc/ProtFun/ape-1.0/disp;/home/cmgfunc/ProtFun/ape-1.0/graphics;/home/cmgfunc/ProtFun/ape-1.0/test")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun/ape-1.0" TYPE DIRECTORY FILES
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./clp"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./tmp"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./etc"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./pred"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./disp"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./graphics"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./test"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/ape-1.0/bin/loadenv.awk;/home/cmgfunc/ProtFun/ape-1.0/bin/nnhowplayer6_Linux.i686;/home/cmgfunc/ProtFun/ape-1.0/bin/nnhowplayer6_Linux.ia64;/home/cmgfunc/ProtFun/ape-1.0/bin/nnhowplayer6_Linux.x86_64;/home/cmgfunc/ProtFun/ape-1.0/bin/paste.awk;/home/cmgfunc/ProtFun/ape-1.0/bin/pp.awk;/home/cmgfunc/ProtFun/ape-1.0/bin/seq2seq.awk")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun/ape-1.0/bin" TYPE PROGRAM PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ FILES
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/loadenv.awk"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/nnhowplayer6_Linux.i686"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/nnhowplayer6_Linux.ia64"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/nnhowplayer6_Linux.x86_64"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/paste.awk"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/pp.awk"
    "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./bin/seq2seq.awk"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/ape")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun" TYPE PROGRAM PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_EXECUTE WORLD_READ FILES "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/./ape")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/cmgfunc/ProtFun/ProtFun/ape-1.0/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
