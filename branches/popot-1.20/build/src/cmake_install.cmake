# Install script for directory: /home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
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
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}/usr/local/lib/libpopot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libpopot.so")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/lib/libpopot.so"
         RPATH "")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libpopot.so")
FILE(INSTALL DESTINATION "/usr/local/lib" TYPE SHARED_LIBRARY FILES "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/build/src/libpopot.so")
  IF(EXISTS "$ENV{DESTDIR}/usr/local/lib/libpopot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libpopot.so")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/lib/libpopot.so")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/popot/standard_pso.h;/usr/local/include/popot/initializers.h;/usr/local/include/popot/neighborhood.h;/usr/local/include/popot/topology.h;/usr/local/include/popot/tools.h;/usr/local/include/popot/individuals.h;/usr/local/include/popot/rng_generators.h;/usr/local/include/popot/problems.h;/usr/local/include/popot/exceptions.h;/usr/local/include/popot/maths.h;/usr/local/include/popot/algorithm.h;/usr/local/include/popot/popot.h")
FILE(INSTALL DESTINATION "/usr/local/include/popot" TYPE FILE FILES
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/standard_pso.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/initializers.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/neighborhood.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/topology.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/tools.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/individuals.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/rng_generators.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/problems.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/exceptions.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/maths.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/algorithm.h"
    "/home/fix_jer/Developpement/Google/popot/branches/popot-1.20/src/popot.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

