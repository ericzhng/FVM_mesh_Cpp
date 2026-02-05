# VTK I/O Tests

add_executable(vtkio_tests
    vtk_reader_test.cpp
)

target_include_directories(vtkio_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(vtkio_tests
    PRIVATE
        FVM::vtkio
        GTest::gtest_main
)

# Inject the path as a macro named TEST_DATA_DIR
target_compile_definitions(vtkio_tests PRIVATE 
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

# Copy vtkio DLL to test directory (Windows shared builds)
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET vtkio_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:vtkio>"
            "$<TARGET_FILE_DIR:vtkio_tests>"
        COMMENT "Copying vtkio DLL to vtkio_tests directory"
    )
endif()

# Use add_test instead of gtest_discover_tests to avoid DLL issues during build
add_test(NAME vtkio_tests COMMAND vtkio_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
