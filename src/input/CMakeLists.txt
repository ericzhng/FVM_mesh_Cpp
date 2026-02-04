# Common Library - Input parsing and expression evaluation

set(INPUT_SOURCES
    input_config.cpp
    input_parser.cpp
    expr_evaluator.cpp
)

set(INPUT_HEADERS
    input_config.hpp
    input_parser.hpp
    expr_evaluator.hpp
)

# Create library (STATIC or SHARED based on FVM_LIBRARY_TYPE)
add_library(yaml_input ${FVM_LIBRARY_TYPE} ${INPUT_SOURCES} ${INPUT_HEADERS})

# Set library properties
target_include_directories(yaml_input
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

# Link yaml-cpp (required for YAML parsing)
target_link_libraries(yaml_input PUBLIC yaml-cpp)

# Handle export definitions for shared/static builds
if(BUILD_SHARED)
    # Building shared library - define EXPORTS macro
    target_compile_definitions(yaml_input PRIVATE FVM_EXPORTS)

    # Set shared library properties
    set_target_properties(yaml_input PROPERTIES
        # Version information
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        # Windows: don't add "lib" prefix
        PREFIX ""
        # Enable visibility on Unix
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
    )

    # On Windows, ensure runtime library is copied to output directory
    if(WIN32)
        set_target_properties(yaml_input PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        )
    endif()
else()
    # Building static library - define STATIC macro
    target_compile_definitions(yaml_input PUBLIC FVM_STATIC)
endif()

# Set C++ standard
target_compile_features(yaml_input PUBLIC cxx_std_17)

# Export alias for easier usage
add_library(FVM::yaml_input ALIAS yaml_input)
