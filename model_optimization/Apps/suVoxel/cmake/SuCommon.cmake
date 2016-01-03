
# set directory structures for the different platforms
if (WIN32)  
  set (SU_PROJECT_LIBDIR "lib")
  set (SU_PROJECT_BINDIR ".")
  set (SU_PROJECT_PLUGINDIR "Plugins")
  if (NOT EXISTS ${CMAKE_BINARY_DIR}/Build/${SU_PROJECT_LIBDIR})
    file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Build/${SU_PROJECT_LIBDIR})
  endif ()
endif (WIN32) 
# ========================================================================
# Macros & Functions
# ========================================================================
## extended version of add_executable that also copies output to out Build directory
FUNCTION (SU_ADD_EXECUTABLE _target)
  add_executable (${_target} ${ARGN})

  # set common target properties defined in common.cmake
  SU_SET_TARGET_PROPS (${_target})

  if (WIN32)
    add_custom_command (TARGET ${_target} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E
                        copy_if_different
                          $<TARGET_FILE:${_target}>
                          ${CMAKE_BINARY_DIR}/Build/${SU_PROJECT_BINDIR}/$<TARGET_FILE_NAME:${_target}>)
  endif (WIN32)
  
  ##install (TARGETS ${_target} DESTINATION ${SU_PROJECT_BINDIR})
endfunction ()


# sets default build properties
MACRO (SU_SET_TARGET_PROPS target)
  if (WIN32)
    set_target_properties (
      ${target} PROPERTIES
      BUILD_WITH_INSTALL_RPATH 1
      SKIP_BUILD_RPATH 0
    )  
  endif ()
  
endmacro ()