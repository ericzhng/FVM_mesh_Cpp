#pragma once

/**
 * @file input_parser.hpp
 * @brief YAML parser for mesh generation input files.
 *
 * This module provides functions to parse YAML configuration files
 * into InputConfig structures for mesh generation.
 */

#include "common/fvm_export.hpp"
#include "input_config.hpp"

#include <string>

namespace fvm
{

    /**
     * @brief Parse a YAML configuration file.
     * @param filepath Path to the YAML file
     * @return Parsed InputConfig structure
     * @throws std::runtime_error if file cannot be read or parsed
     */
    FVM_API InputConfig parseYamlFile(const std::string &filepath);

    /**
     * @brief Parse YAML content from a string.
     * @param yamlContent YAML content as string
     * @return Parsed InputConfig structure
     * @throws std::runtime_error if content cannot be parsed
     */
    FVM_API InputConfig parseYamlString(const std::string &yamlContent);

} // namespace fvm
