cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
# This will define the following variables
#
#    NEUT_FOUND
#
# and the following imported targets
#
#     NEUT::Core
#

# Colorful messaging facilities
if(NOT COMMAND cmessage)
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(CM_ColourReset "${Esc}[m")
    set(CM_ColourBold "${Esc}[1m")
    set(CM_Red "${Esc}[31m")
    set(CM_Green "${Esc}[32m")
    set(CM_Yellow "${Esc}[33m")
    set(CM_Blue "${Esc}[34m")
    set(CM_Magenta "${Esc}[35m")
    set(CM_Cyan "${Esc}[36m")
    set(CM_White "${Esc}[37m")
    set(CM_BoldRed "${Esc}[1;31m")
    set(CM_BoldGreen "${Esc}[1;32m")
    set(CM_BoldYellow "${Esc}[1;33m")
    set(CM_BoldBlue "${Esc}[1;34m")
    set(CM_BoldMagenta "${Esc}[1;35m")
    set(CM_BoldCyan "${Esc}[1;36m")
    set(CM_BoldWhite "${Esc}[1;37m")
  endif()

  message(STATUS "Setting up colored messages...")

  function(cmessage)
    list(GET ARGV 0 MessageType)
    if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldRed}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldYellow}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL AUTHOR_WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldCyan}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL STATUS)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL CACHE)        
      list(REMOVE_AT ARGV 0)
      message(-- "${CM_Blue}[CACHE]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL DEBUG)
      list(REMOVE_AT ARGV 0)
      if(BUILD_DEBUG_MSGS)
        message("${CM_Magenta}[DEBUG]:${CM_ColourReset} ${ARGV}")
      endif()
    else()
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    endif()
  endfunction()
endif()


find_package(PkgConfig REQUIRED)
pkg_check_modules(NEUT QUIET NEUT)

find_path(NEUT_INCLUDE_DIR
  NAMES necardC.h
  PATHS ${NEUT_INCLUDE_DIRS}
)

pkg_get_variable(NEUT_PREFIX NEUT prefix)
pkg_get_variable(NEUT_EXE_LINKFLAGS NEUT EXE_LINKFLAGS)

find_path(MANUAL_NEUT_BIN_DIR
  NAMES neutroot2
  PATHS ${NEUT_PREFIX}
)

find_library(MANUAL_NEUT_LIBRARY_DIR
  NAMES libNEUT.a
  PATHS ${NEUT_LIBRARY_DIRS}
)

set(NEUT_VERSION ${NEUT_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NEUT
    REQUIRED_VARS NEUT_INCLUDE_DIR MANUAL_NEUT_BIN_DIR MANUAL_NEUT_LIBRARY_DIR
    VERSION_VAR NEUT_VERSION
)

if(NEUT_FOUND)

  include(CheckPIESupported)
  check_pie_supported(LANGUAGES C CXX)

  #We should let CMake set this
  list(REMOVE_ITEM NEUT_CFLAGS_OTHER "-fPIC")

  cmessage(STATUS "NEUT Found: ${NEUT_PREFIX} ")
  cmessage(STATUS "    NEUT_INCLUDE_DIRS: ${NEUT_INCLUDE_DIRS}")
  cmessage(STATUS "    NEUT_CFLAGS_OTHER: ${NEUT_CFLAGS_OTHER}")
  cmessage(STATUS "    NEUT_LIBRARY_DIRS: ${NEUT_LIBRARY_DIRS}")
  cmessage(STATUS "    NEUT_LIBRARIES: ${NEUT_LIBRARIES}")
  cmessage(STATUS "    NEUT_LINK_OPTIONS: ${NEUT_LINK_OPTIONS}")

  #This horrow show lets us add the exe-flag only for executable targets.
  SET(NEUT_LINK_OPTIONS ${NEUT_LDFLAGS_OTHER})
  foreach(opt ${NEUT_EXE_LINKFLAGS})
    LIST(APPEND NEUT_LINK_OPTIONS $<IF:$<STREQUAL:"$<TARGET_PROPERTY:TYPE>","EXECUTABLE">,${opt},>)
  endforeach()

  if(NOT TARGET NEUT::Core)
      add_library(NEUT::Core INTERFACE IMPORTED)
      set_target_properties(NEUT::Core PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${NEUT_INCLUDE_DIRS}"
          INTERFACE_COMPILE_OPTIONS "${NEUT_CFLAGS_OTHER}"
          INTERFACE_LINK_DIRECTORIES "${NEUT_LIBRARY_DIRS}"
          INTERFACE_LINK_LIBRARIES "${NEUT_LIBRARIES}"
          INTERFACE_LINK_OPTIONS "${NEUT_LINK_OPTIONS}"
      )
  endif()

endif()