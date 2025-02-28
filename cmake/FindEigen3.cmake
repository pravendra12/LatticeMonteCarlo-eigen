# - Try to find Eigen3 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen3 3.1.2)
# to require version 3.1.2 or newer of Eigen3.
#
# Once done this will define

#  EIGEN3_FOUND - system has eigen lib with correct version
#  EIGEN3_INCLUDE_DIRS - the eigen include directory
#  EIGEN3_VERSION - eigen version
#
#  and the following imported target:

#  Eigen3::Eigen - The header-only Eigen library
#
# This module reads hints about search locations from
# the following environment variables:
#
# EIGEN3_ROOT
# EIGEN3_ROOT_DIR

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.
if (NOT Eigen3_FIND_VERSION)
    if (NOT Eigen3_FIND_VERSION_MAJOR)
        set(Eigen3_FIND_VERSION_MAJOR 2)
    endif ()
    if (NOT Eigen3_FIND_VERSION_MINOR)
        set(Eigen3_FIND_VERSION_MINOR 91)
    endif ()
    if (NOT Eigen3_FIND_VERSION_PATCH)
        set(Eigen3_FIND_VERSION_PATCH 0)
    endif ()
    set(Eigen3_FIND_VERSION "${Eigen3_FIND_VERSION_MAJOR}.${Eigen3_FIND_VERSION_MINOR}.${Eigen3_FIND_VERSION_PATCH}")
endif ()

macro(_eigen3_check_version)
    file(READ "${EIGEN3_INCLUDE_DIRS}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)
    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
    set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
    set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
    set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")
    set(EIGEN3_VERSION ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION})
    if (${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})
        set(EIGEN3_VERSION_OK FALSE)
    else ()
        set(EIGEN3_VERSION_OK TRUE)
    endif ()
    if (NOT EIGEN3_VERSION_OK)
        message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIRS}, "
                "but at least version ${Eigen3_FIND_VERSION} is required")
    endif ()
endmacro()

if (EIGEN3_INCLUDE_DIRS)
    # in cache already
    _eigen3_check_version()
    set(EIGEN3_FOUND ${EIGEN3_VERSION_OK})
else ()
    # search first if an Eigen3Config.cmake is available in the system,
    # if successful this would set EIGEN3_INCLUDE_DIRS and the rest of
    # the script will work as usual
    find_package(Eigen3 ${Eigen3_FIND_VERSION} NO_MODULE QUIET)
    if (NOT EIGEN3_INCLUDE_DIRS)
        find_path(EIGEN3_INCLUDE_DIRS NAMES signature_of_eigen3_matrix_library
                  HINTS
                  ENV EIGEN3_ROOT
                  ENV EIGEN3_ROOT_DIR
                  PATHS
                  ${KDE4_INCLUDE_DIR}
                  ${CMAKE_INSTALL_PREFIX}/include/eigen3
                  ${CMAKE_INSTALL_PREFIX}/include
                  /usr/local/include/eigen3
                  /usr/local/include
                  /usr/include/eigen3
                  /usr/include
                  /sw/pkgs/arc/eigen/eigen-3.4.0
                  ~/.local/include/eigen3
                  ~/.local/include/
                  PATH_SUFFIXES eigen3 eigen
                  )
    endif ()
    if (EIGEN3_INCLUDE_DIRS)
        _eigen3_check_version()
    endif ()
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIRS EIGEN3_VERSION_OK)
    mark_as_advanced(EIGEN3_INCLUDE_DIRS)
endif ()

if (EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIRS}")
endif ()
