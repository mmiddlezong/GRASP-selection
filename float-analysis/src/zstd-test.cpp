#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <filesystem>
#include <vector>
#include <zstd.h>

namespace fs = std::filesystem;

void ZSTD_compressFile(const fs::path& filePath, int level) {
    // Open file and read content
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        std::cerr << "Failed to read file: " << filePath << std::endl;
        return;
    }

    // Prepare for compression
    size_t compressedSize = ZSTD_compressBound(size); // Maximum compressed size
    std::vector<char> compressed(compressedSize);

    // Measure compression time
    auto start = std::chrono::high_resolution_clock::now();
    size_t const cSize = ZSTD_compress(compressed.data(), compressedSize, buffer.data(), size, level); // Level 1 compression
    auto end = std::chrono::high_resolution_clock::now();

    if (ZSTD_isError(cSize)) {
        std::cerr << "Compression error: " << ZSTD_getErrorName(cSize) << std::endl;
        return;
    }

    std::chrono::duration<double, std::milli> elapsed = end - start;
    double throughput = size / (elapsed.count() / 1000.0) / 1024.0 / 1024.0; // Throughput in MB/s
    double ratio = static_cast<double>(size) / cSize;

    // Print in CSV format: FileName, Size, CompressedSize, Throughput(MB/s), CompressionRatio
    std::cout << filePath.filename() << "," << size << "," << cSize << "," << throughput << "," << ratio << std::endl;
}

int main() {
    fs::path realDatasetsDir("/home/appllo/projects/sdr-projects/huffman/real-datasets/");
    std::vector<fs::path> directoryPaths = {realDatasetsDir / "NYX", realDatasetsDir / "SCALE" };

    std::cout << "Input compression level: ";
    int level;
    std::cin >> level;
    std::cout << level << std::endl;

    for (const std::string& directoryPath : directoryPaths) {
        std::cout << "\n" << directoryPath << " ---------------------------------\n\n";
        // Print CSV header
        std::cout << "FileName,Size,CompressedSize,Throughput(MB/s),CompressionRatio\n";

        // Loop through all files in the directory
        for (const auto& entry : fs::directory_iterator(directoryPath)) {
            if (entry.is_regular_file()) { // Check if it is a normal file
                ZSTD_compressFile(entry.path(), level);
            }
        }
    }

    return 0;
}
