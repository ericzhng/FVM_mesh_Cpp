# CMakeLists.txt for executable applications

# =============================================================================
# Mesh Generator Application
# =============================================================================
add_executable(fvm_meshgen meshgen_main.cpp)

target_include_directories(fvm_meshgen PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${Gmsh_INCLUDE_DIR}
)

target_link_libraries(fvm_meshgen
    PRIVATE
        FVM::meshgen
        FVM::vtkio
)

# =============================================================================
# Mesh Partitioner Application
# =============================================================================
add_executable(fvm_partitioner partitioner_main.cpp)

target_include_directories(fvm_partitioner PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(fvm_partitioner
    PRIVATE
        FVM::polymesh
        FVM::vtkio
)

# =============================================================================
# Unified Mesh Generator Application (with YAML input)
# =============================================================================
add_executable(whole_mesh_main whole_mesh_main.cpp)

target_include_directories(whole_mesh_main PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${Gmsh_INCLUDE_DIR}
    ${CLI11_INCLUDE_DIR}
)

target_link_libraries(whole_mesh_main
    PRIVATE
        FVM::yaml_input
        FVM::meshgen
        FVM::polymesh
        FVM::vtkio
        CLI11::CLI11
)

# =============================================================================
# DLL Copy Commands (Windows only)
# =============================================================================

# Common DLLs for all executables
set(FVM_EXECUTABLES fvm_meshgen fvm_partitioner whole_mesh_main)

foreach(TARGET_EXE ${FVM_EXECUTABLES})
    # Copy Gmsh DLL
    if(WIN32 AND GMSH_DLL)
        add_custom_command(TARGET ${TARGET_EXE} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${GMSH_DLL}"
                "$<TARGET_FILE_DIR:${TARGET_EXE}>"
            COMMENT "Copying Gmsh DLL to ${TARGET_EXE} output directory"
        )
    endif()

    # Copy METIS DLL
    if(WIN32 AND HAVE_METIS AND METIS_DLL)
        add_custom_command(TARGET ${TARGET_EXE} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${METIS_DLL}"
                "$<TARGET_FILE_DIR:${TARGET_EXE}>"
            COMMENT "Copying METIS DLL to ${TARGET_EXE} output directory"
        )
    endif()

    # Copy vtkio DLL (shared builds only - needed by all apps)
    if(WIN32 AND BUILD_SHARED)
        add_custom_command(TARGET ${TARGET_EXE} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:vtkio>"
                "$<TARGET_FILE_DIR:${TARGET_EXE}>"
            COMMENT "Copying vtkio DLL to ${TARGET_EXE} output directory"
        )
    endif()
endforeach()

# fvm_meshgen needs meshgen DLL
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET fvm_meshgen POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:meshgen>"
            "$<TARGET_FILE_DIR:fvm_meshgen>"
        COMMENT "Copying meshgen DLL to fvm_meshgen output directory"
    )
endif()

# fvm_partitioner needs polymesh DLL
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET fvm_partitioner POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:polymesh>"
            "$<TARGET_FILE_DIR:fvm_partitioner>"
        COMMENT "Copying polymesh DLL to fvm_partitioner output directory"
    )
endif()

# whole_mesh_main needs meshgen, polymesh, and yaml_input DLLs
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET whole_mesh_main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:meshgen>"
            "$<TARGET_FILE_DIR:whole_mesh_main>"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:polymesh>"
            "$<TARGET_FILE_DIR:whole_mesh_main>"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:yaml_input>"
            "$<TARGET_FILE_DIR:whole_mesh_main>"
        COMMENT "Copying DLLs to whole_mesh_main output directory"
    )
endif()

# =============================================================================
# Installation Rules
# =============================================================================
install(TARGETS fvm_meshgen fvm_partitioner whole_mesh_main
    RUNTIME DESTINATION ${FVM_INSTALL_BINDIR}
    COMPONENT Runtime
)
