# Add library subdirectories
add_subdirectory(input)
add_subdirectory(vtkio)
add_subdirectory(meshgen)
add_subdirectory(polymesh)

# Add applications subdirectory (executables)
add_subdirectory(apps)

# =============================================================================
# Installation Rules
# =============================================================================

# Install core libraries with export configuration
install(TARGETS meshgen polymesh vtkio
    EXPORT FVM_meshTargets
    ARCHIVE DESTINATION ${FVM_INSTALL_LIBDIR} COMPONENT Development
    LIBRARY DESTINATION ${FVM_INSTALL_LIBDIR} COMPONENT Runtime
    RUNTIME DESTINATION ${FVM_INSTALL_BINDIR} COMPONENT Runtime  # For Windows DLLs
    INCLUDES DESTINATION ${FVM_INSTALL_INCLUDEDIR}
)

# Install yaml_input separately (not exported due to yaml-cpp dependency)
install(TARGETS yaml_input
    ARCHIVE DESTINATION ${FVM_INSTALL_LIBDIR} COMPONENT Development
    LIBRARY DESTINATION ${FVM_INSTALL_LIBDIR} COMPONENT Runtime
    RUNTIME DESTINATION ${FVM_INSTALL_BINDIR} COMPONENT Runtime
)

# Install header files
install(DIRECTORY common/
    DESTINATION ${FVM_INSTALL_INCLUDEDIR}/common
    COMPONENT Development
    FILES_MATCHING PATTERN "*.hpp"
)
install(DIRECTORY meshgen/
    DESTINATION ${FVM_INSTALL_INCLUDEDIR}/meshgen
    COMPONENT Development
    FILES_MATCHING PATTERN "*.hpp"
)
install(DIRECTORY polymesh/
    DESTINATION ${FVM_INSTALL_INCLUDEDIR}/polymesh
    COMPONENT Development
    FILES_MATCHING PATTERN "*.hpp"
)
install(DIRECTORY vtkio/
    DESTINATION ${FVM_INSTALL_INCLUDEDIR}/vtkio
    COMPONENT Development
    FILES_MATCHING PATTERN "*.hpp"
)

# Export targets for use by other projects
install(EXPORT FVM_meshTargets
    FILE FVM_meshTargets.cmake
    NAMESPACE FVM::
    DESTINATION ${FVM_INSTALL_CMAKEDIR}
    COMPONENT Development
)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/FVM_meshConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Generate package config file from template
configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/FVM_meshConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/FVM_meshConfig.cmake"
    INSTALL_DESTINATION ${FVM_INSTALL_CMAKEDIR}
)

# Install CMake configuration files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FVM_meshConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FVM_meshConfigVersion.cmake"
    DESTINATION ${FVM_INSTALL_CMAKEDIR}
    COMPONENT Development
)

# Install FindGmsh and FindMETIS for consumers who might need them
install(FILES
    "${CMAKE_SOURCE_DIR}/cmake/FindGmsh.cmake"
    "${CMAKE_SOURCE_DIR}/cmake/FindMETIS.cmake"
    DESTINATION ${FVM_INSTALL_CMAKEDIR}/modules
    COMPONENT Development
)

# Copy Gmsh DLL to install bin directory (Windows only)
if(WIN32 AND GMSH_DLL)
    install(FILES "${GMSH_DLL}"
        DESTINATION ${FVM_INSTALL_BINDIR}
        COMPONENT Runtime
    )
endif()

# Copy METIS DLL if available (Windows only)
if(WIN32 AND HAVE_METIS AND METIS_DLL)
    install(FILES "${METIS_DLL}"
        DESTINATION ${FVM_INSTALL_BINDIR}
        COMPONENT Runtime
    )
endif()
