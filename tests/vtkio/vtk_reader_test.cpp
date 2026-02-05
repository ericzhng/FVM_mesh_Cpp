#include <gtest/gtest.h>
#include "vtkio/vtk_reader.hpp"
#include "vtkio/vtk_writer.hpp"
#include <filesystem>
#include <fstream>
#include <cstdio>

namespace fs = std::filesystem;

// Simple test that doesn't use std::filesystem
TEST(VTKReaderBasic, CanCreateReader)
{
    // Just verify the header compiles and links correctly
    EXPECT_TRUE(true);
}

TEST(VTKReaderBasic, ReadNonExistentFileThrows)
{
    EXPECT_THROW(fvm::VTKReader::readVTK("nonexistent_file_xyz.vtk"), std::runtime_error);
}

TEST(VTKReaderBasic, UnsupportedFormatThrows)
{
    EXPECT_THROW(fvm::VTKReader::read("test.xyz"), std::runtime_error);
}

// Test reading actual files - only if they exist
TEST(VTKReaderFile, ReadVTKFile)
{
    fs::path testFile = fs::path(TEST_DATA_DIR) / "sample_rect_mesh.vtk";

    std::ifstream f(testFile);
    if (!f.good())
    {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    auto mesh = fvm::VTKReader::readVTK(testFile.string());
    EXPECT_GT(mesh.nodes.size(), 0u);
    EXPECT_GT(mesh.elements.size(), 0u);
    EXPECT_EQ(mesh.elements.size(), mesh.elementTypes.size());
}

TEST(VTKReaderFile, ReadVTUFile)
{
    fs::path testFile = fs::path(TEST_DATA_DIR) / "sample_rect_mesh.vtu";

    std::ifstream f(testFile);
    if (!f.good())
    {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    auto mesh = fvm::VTKReader::readVTU(testFile.string());
    EXPECT_GT(mesh.nodes.size(), 0u);
    EXPECT_GT(mesh.elements.size(), 0u);
    EXPECT_EQ(mesh.elements.size(), mesh.elementTypes.size());
}

TEST(VTKReaderFile, RoundTripVTK)
{
    fs::path testFile = fs::path(TEST_DATA_DIR) / "sample_rect_mesh.vtk";

    std::ifstream f(testFile);
    if (!f.good())
    {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    // Read original
    auto mesh1 = fvm::VTKReader::readVTK(testFile.string());

    // Write to temp file
    fs::path tempFile = fs::path(TEST_DATA_DIR) / "test_roundtrip_temp.vtk";
    fvm::VTKWriter::writeVTK(mesh1, tempFile.string());

    // Read back
    auto mesh2 = fvm::VTKReader::readVTK(tempFile.string());

    // Compare
    EXPECT_EQ(mesh1.nodes.size(), mesh2.nodes.size());
    EXPECT_EQ(mesh1.elements.size(), mesh2.elements.size());
    EXPECT_EQ(mesh1.elementTypes.size(), mesh2.elementTypes.size());

    // Compare first few nodes
    for (size_t i = 0; i < std::min(mesh1.nodes.size(), size_t(10)); ++i)
    {
        EXPECT_NEAR(mesh1.nodes[i][0], mesh2.nodes[i][0], 1e-6);
        EXPECT_NEAR(mesh1.nodes[i][1], mesh2.nodes[i][1], 1e-6);
        EXPECT_NEAR(mesh1.nodes[i][2], mesh2.nodes[i][2], 1e-6);
    }
}

TEST(VTKReaderFile, RoundTripVTU)
{
    fs::path testFile = fs::path(TEST_DATA_DIR) / "sample_rect_mesh.vtu";

    std::ifstream f(testFile);
    if (!f.good())
    {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    // Read original
    auto mesh1 = fvm::VTKReader::readVTU(testFile.string());

    // Write to temp file
    fs::path tempFile = fs::path(TEST_DATA_DIR) / "test_roundtrip_temp.vtu";
    fvm::VTKWriter::writeVTU(mesh1, tempFile.string());

    // Read back
    auto mesh2 = fvm::VTKReader::readVTU(tempFile.string());

    // Compare
    EXPECT_EQ(mesh1.nodes.size(), mesh2.nodes.size());
    EXPECT_EQ(mesh1.elements.size(), mesh2.elements.size());
}
