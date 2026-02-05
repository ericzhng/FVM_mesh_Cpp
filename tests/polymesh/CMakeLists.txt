# Polymesh module tests

# =============================================================================
# PolyMesh Tests
# =============================================================================

add_executable(poly_mesh_tests
    poly_mesh_test.cpp
)

set_target_properties(poly_mesh_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

target_include_directories(poly_mesh_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${Gmsh_INCLUDE_DIR}
)

target_link_libraries(poly_mesh_tests
    PRIVATE
        FVM::polymesh
        Gmsh::Gmsh
        GTest::gtest_main
)

# Copy DLLs to test directory (Windows shared builds)
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET poly_mesh_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:polymesh>"
            "$<TARGET_FILE_DIR:poly_mesh_tests>"
        COMMENT "Copying polymesh DLL to poly_mesh_tests directory"
    )
endif()

if(WIN32 AND GMSH_DLL)
    add_custom_command(TARGET poly_mesh_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GMSH_DLL}"
            "$<TARGET_FILE_DIR:poly_mesh_tests>"
        COMMENT "Copying Gmsh DLL to poly_mesh_tests directory"
    )
endif()

add_test(NAME poly_mesh_tests
    COMMAND poly_mesh_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

# =============================================================================
# Reorder Tests
# =============================================================================

add_executable(reorder_tests
    reorder_test.cpp
)

set_target_properties(reorder_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)

target_include_directories(reorder_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${Gmsh_INCLUDE_DIR}
)

target_link_libraries(reorder_tests
    PRIVATE
        FVM::polymesh
        Gmsh::Gmsh
        GTest::gtest_main
)

# Copy DLLs to test directory (Windows shared builds)
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET reorder_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:polymesh>"
            "$<TARGET_FILE_DIR:reorder_tests>"
        COMMENT "Copying polymesh DLL to reorder_tests directory"
    )
endif()

if(WIN32 AND GMSH_DLL)
    add_custom_command(TARGET reorder_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GMSH_DLL}"
            "$<TARGET_FILE_DIR:reorder_tests>"
        COMMENT "Copying Gmsh DLL to reorder_tests directory"
    )
endif()

add_test(NAME reorder_tests
    COMMAND reorder_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>
)
