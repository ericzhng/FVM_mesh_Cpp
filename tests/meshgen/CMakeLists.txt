# Meshgen tests

# =============================================================================
# Geometry Tests
# =============================================================================

add_executable(geometry_tests
    geometry_test.cpp
)

set_target_properties(geometry_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

target_include_directories(geometry_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${Gmsh_INCLUDE_DIR}
)

target_link_libraries(geometry_tests
    PRIVATE
        meshgen
        Gmsh::Gmsh
        GTest::gtest_main
)

if(WIN32 AND GMSH_DLL)
    add_custom_command(TARGET geometry_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GMSH_DLL}"
            "$<TARGET_FILE_DIR:geometry_tests>"
        COMMENT "Copying Gmsh DLL to test directory"
    )
endif()

add_test(NAME geometry_tests
    COMMAND geometry_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

# =============================================================================
# MeshGenerator Tests
# =============================================================================

add_executable(mesh_generator_tests
    mesh_generator_test.cpp
)

set_target_properties(mesh_generator_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

target_include_directories(mesh_generator_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${Gmsh_INCLUDE_DIR}
)

target_link_libraries(mesh_generator_tests
    PRIVATE
        meshgen
        Gmsh::Gmsh
        GTest::gtest_main
)

if(WIN32 AND GMSH_DLL)
    add_custom_command(TARGET mesh_generator_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GMSH_DLL}"
            "$<TARGET_FILE_DIR:mesh_generator_tests>"
        COMMENT "Copying Gmsh DLL to test directory"
    )
endif()

add_test(NAME mesh_generator_tests
    COMMAND mesh_generator_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)
