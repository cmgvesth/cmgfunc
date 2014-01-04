# Install script for directory: /home/cmgfunc/ProtFun/ProtFun/prop-1.0c

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/cmgfunc/ProtFun")
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
   "/home/cmgfunc/ProtFun/prop-1.0c/etc;/home/cmgfunc/ProtFun/prop-1.0c/syn;/home/cmgfunc/ProtFun/prop-1.0c/test;/home/cmgfunc/ProtFun/prop-1.0c/tmp")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun/prop-1.0c" TYPE DIRECTORY FILES
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./etc"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./syn"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./test"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./tmp"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/prop")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun" TYPE PROGRAM FILES "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./prop")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/prop-1.0c/bin/fasta2dir;/home/cmgfunc/ProtFun/prop-1.0c/bin/gethow;/home/cmgfunc/ProtFun/prop-1.0c/bin/getseqstr;/home/cmgfunc/ProtFun/prop-1.0c/bin/in2how;/home/cmgfunc/ProtFun/prop-1.0c/bin/in2how+fasta.awk;/home/cmgfunc/ProtFun/prop-1.0c/bin/mkgraph;/home/cmgfunc/ProtFun/prop-1.0c/bin/predict_S;/home/cmgfunc/ProtFun/prop-1.0c/bin/predict_T;/home/cmgfunc/ProtFun/prop-1.0c/bin/predict_Y;/home/cmgfunc/ProtFun/prop-1.0c/bin/prop_furin;/home/cmgfunc/ProtFun/prop-1.0c/bin/prop_total")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun/prop-1.0c/bin" TYPE PROGRAM FILES
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/fasta2dir"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/gethow"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/getseqstr"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/in2how"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/in2how+fasta.awk"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/mkgraph"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/predict_S"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/predict_T"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/predict_Y"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/prop_furin"
    "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./bin/prop_total"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/cmgfunc/ProtFun/prop-1.0c/how/how98")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/home/cmgfunc/ProtFun/prop-1.0c/how" TYPE PROGRAM RENAME "how98" FILES "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/./how/how98_Linux")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/cmgfunc/ProtFun/ProtFun/prop-1.0c/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
