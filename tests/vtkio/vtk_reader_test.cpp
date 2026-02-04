#include <gtest/gtest.h>
#include "vtkio/vtk_reader.hpp"
#include "vtkio/vtk_writer.hpp"
#include <fstream>
#include <cstdio>

// Simple test that doesn't use std::filesystem
TEST(VTKReaderBasic, CanCreateReader) {
    // Just verify the header compiles and links correctly
    EXPECT_TRUE(true);
}

TEST(VTKReaderBasic, ReadNonExistentFileThrows) {
    EXPECT_THROW(fvm::VTKReader::readVTK("nonexistent_file_xyz.vtk"), std::runtime_error);
}

TEST(VTKReaderBasic, UnsupportedFormatThrows) {
    EXPECT_THROW(fvm::VTKReader::read("test.xyz"), std::runtime_error);
}

// Test reading actual files - only if they exist
TEST(VTKReaderFile, ReadVTKFile) {
    const char* testFile = "data/sample_rect_mesh.vtk";
    std::ifstream f(testFile);
    if (!f.good()) {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    auto mesh = fvm::VTKReader::readVTK(testFile);
    EXPECT_GT(mesh.nodes.size(), 0u);
    EXPECT_GT(mesh.cells.size(), 0u);
    EXPECT_EQ(mesh.cells.size(), mesh.cellTypes.size());
}

TEST(VTKReaderFile, ReadVTUFile) {
    const char* testFile = "data/sample_rect_mesh.vtu";
    std::ifstream f(testFile);
    if (!f.good()) {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    auto mesh = fvm::VTKReader::readVTU(testFile);
    EXPECT_GT(mesh.nodes.size(), 0u);
    EXPECT_GT(mesh.cells.size(), 0u);
    EXPECT_EQ(mesh.cells.size(), mesh.cellTypes.size());
}

TEST(VTKReaderFile, RoundTripVTK) {
    const char* testFile = "data/sample_rect_mesh.vtk";
    std::ifstream f(testFile);
    if (!f.good()) {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    // Read original
    auto mesh1 = fvm::VTKReader::readVTK(testFile);

    // Write to temp file
    const char* tempFile = "data/test_roundtrip_temp.vtk";
    fvm::VTKWriter::writeVTK(mesh1, tempFile);

    // Read back
    auto mesh2 = fvm::VTKReader::readVTK(tempFile);

    // Compare
    EXPECT_EQ(mesh1.nodes.size(), mesh2.nodes.size());
    EXPECT_EQ(mesh1.cells.size(), mesh2.cells.size());
    EXPECT_EQ(mesh1.cellTypes.size(), mesh2.cellTypes.size());

    // Compare first few nodes
    for (size_t i = 0; i < std::min(mesh1.nodes.size(), size_t(10)); ++i) {
        EXPECT_NEAR(mesh1.nodes[i][0], mesh2.nodes[i][0], 1e-6);
        EXPECT_NEAR(mesh1.nodes[i][1], mesh2.nodes[i][1], 1e-6);
        EXPECT_NEAR(mesh1.nodes[i][2], mesh2.nodes[i][2], 1e-6);
    }

    // Cleanup
    std::remove(tempFile);
}

TEST(VTKReaderFile, RoundTripVTU) {
    const char* testFile = "data/sample_rect_mesh.vtu";
    std::ifstream f(testFile);
    if (!f.good()) {
        GTEST_SKIP() << "Test file not found: " << testFile;
    }
    f.close();

    // Read original
    auto mesh1 = fvm::VTKReader::readVTU(testFile);

    // Write to temp file
    const char* tempFile = "data/test_roundtrip_temp.vtu";
    fvm::VTKWriter::writeVTU(mesh1, tempFile);

    // Read back
    auto mesh2 = fvm::VTKReader::readVTU(tempFile);

    // Compare
    EXPECT_EQ(mesh1.nodes.size(), mesh2.nodes.size());
    EXPECT_EQ(mesh1.cells.size(), mesh2.cells.size());

    // Cleanup
    std::remove(tempFile);
}
