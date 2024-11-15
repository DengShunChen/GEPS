# Find the DMS includes and library
#
# Usage: find_package(DMS [REQUIRED] [QUIET])
#
# It sets the following variables: DMS_FOUND               ... true if DMS is
# found on the system; if false, you cannot build anything that requires DMS
# DMS_INCLUDE_DIR         ... DMS include directory where to locate DMS.h file
# DMS_LIBRARY             ... full path to DMS library to link against to use
# DMS
#
# The following variables will be checked by the function DMS_USE_STATIC_LIBS
# ... if true, only static libraries are found

# If environment variable DMSDIR is specified, it has same effect as DMS_ROOT
if(NOT DMS_ROOT)
  set(DMS_ROOT $ENV{DMS})
endif()
message(STATUS "Found DMS root directory at ${DMS_ROOT}")

# Check whether to search static or dynamic libs
set(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})

set(DMS_USE_STATIC_LIBS NO)
if(${DMS_USE_STATIC_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
  # set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
endif()

set(DMS_LIBRARIES)
foreach(_lib IN ITEMS rdms gdbm)

  # find DMS library
  find_library(
    DMS_LIBRARY_${_lib}
    NAMES ${_lib}
    PATHS /users/xa09/pkg/x86_64/dms38key/lib
          ${DMS_ROOT}/lib
          $ENV{DMS_ROOT}/lib
          ${DMS_ROOT}/lib64
          $ENV{DMS}
          $ENV{DMS}/lib
          $ENV{DMS}/.libs
          /usr/local/lib
          /usr/lib
          /opt/DMS/lib
    DOC "Specify the dms library here."
    NO_DEFAULT_PATH)
  message(STATUS "  - DMS_LIBRARY: ${DMS_LIBRARY_${_lib}}")
  list(APPEND DMS_LIBRARIES ${DMS_LIBRARY_${_lib}})

endforeach()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})

if(DMS_LIBRARIES)
  set(DMS_FOUND 1)
  if(NOT DMS_FIND_QUIETLY)
    message(STATUS "Looking for DMS... found ${DMS_LIBRARIES}")
  endif(NOT DMS_FIND_QUIETLY)
else()
  set(DMS_FOUND 0)
  if(DMS_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for DMS... not found")
  else(DMS_FIND_REQUIRED)
    if(NOT DMS_FIND_QUIETLY)
      message(STATUS "Looking for DMS... not found")
    endif(NOT DMS_FIND_QUIETLY)
  endif(DMS_FIND_REQUIRED)
endif()

# set(DMS_LIBRARIES ${DMS_LIBRARY})

mark_as_advanced(DMS_FOUND DMS_INCLUDE_DIR DMS_LIBRARY DMS_LIBRARIES)
