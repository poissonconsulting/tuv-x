include(CMakePackageConfigHelpers)

install(
  TARGETS 
    tuvx tuv-x
  EXPORT 
    tuvx_Exports
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
)

install(
  DIRECTORY 
    ${TUVX_MOD_DIR}/
  DESTINATION 
    ${INSTALL_MOD_DIR}
  FILES_MATCHING PATTERN "*.mod"
)

install(
  DIRECTORY 
    ${PROJECT_SOURCE_DIR}/data/
  DESTINATION 
    ${INSTALL_PREFIX}/data
)

install(
  DIRECTORY 
    ${PROJECT_SOURCE_DIR}/examples/
  DESTINATION 
    ${INSTALL_PREFIX}/examples
)

install(
  EXPORT 
    tuvx_Exports
  DESTINATION 
    ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${PROJECT_VERSION}
  NAMESPACE musica::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/tuvxConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/tuvxConfig.cmake"
  INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${PROJECT_VERSION}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/tuvxConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/tuvxConfig.cmake
    ${PROJECT_BINARY_DIR}/tuvxConfigVersion.cmake
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/tuvx-${PROJECT_VERSION}
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()