install(
  TARGETS tuvxlib
  EXPORT tuvxlibTargets
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
)


install(
  CODE "file(GLOB_RECURSE MODULE_FILES \"${TUVX_MOD_DIR}/*.mod\")"
  CODE "file(GLOB_RECURSE SUBMOD_FILES \"${TUVX_MOD_DIR}/*.smod\")"
  CODE "file(INSTALL \${MODULE_FILES} DESTINATION ${INSTALL_MOD_DIR})"
  CODE "file(INSTALL \${SUBMOD_FILES} DESTINATION ${INSTALL_MOD_DIR})"
  )

install(EXPORT tuvxlibTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${CMAKE_Fortran_COMPILER_ID}
  NAMESPACE musica::
  FILE tuvx-${CMAKE_Fortran_COMPILER_ID}-config.cmake
)

install(EXPORT tuvxlibTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${CMAKE_Fortran_COMPILER_ID}
  NAMESPACE musica::
  FILE tuvx-${CMAKE_Fortran_COMPILER_ID}Config.cmake 
)

string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} compiler_id)
install(EXPORT tuvxlibTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${CMAKE_Fortran_COMPILER_ID}
  NAMESPACE musica::
  FILE tuvx-${compiler_id}-config.cmake
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()