# Input module tests

# =============================================================================
# Expression Evaluator Tests
# =============================================================================

add_executable(expr_evaluator_tests
    expr_evaluator_test.cpp
)

target_include_directories(expr_evaluator_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(expr_evaluator_tests
    PRIVATE
        FVM::yaml_input
        GTest::gtest_main
)

# Copy DLL to test directory (Windows shared builds)
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET expr_evaluator_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:yaml_input>"
            "$<TARGET_FILE_DIR:expr_evaluator_tests>"
        COMMENT "Copying yaml_input DLL to expr_evaluator_tests directory"
    )
endif()

add_test(NAME expr_evaluator_tests COMMAND expr_evaluator_tests)

# =============================================================================
# Input Parser Tests
# =============================================================================

add_executable(input_parser_tests
    input_parser_test.cpp
)

target_include_directories(input_parser_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(input_parser_tests
    PRIVATE
        FVM::yaml_input
        GTest::gtest_main
)

# Copy DLL to test directory (Windows shared builds)
if(WIN32 AND BUILD_SHARED)
    add_custom_command(TARGET input_parser_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:yaml_input>"
            "$<TARGET_FILE_DIR:input_parser_tests>"
        COMMENT "Copying yaml_input DLL to input_parser_tests directory"
    )
endif()

add_test(NAME input_parser_tests COMMAND input_parser_tests)
