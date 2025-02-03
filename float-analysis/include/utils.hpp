#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <cassert>
#include <chrono>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <zstd.h>

// ---------------------------------------
// UTIL FUNCTIONS FOR FILE I/O WITH FLOATS
// ---------------------------------------

namespace fs = std::filesystem;

size_t getFileSize(const fs::path &filePath) {
    std::ifstream in(filePath, std::ios::ate | std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "File could not be opened.\n";
        return 0;
    }
    return in.tellg();
}
char *readBigEndianBytes(fs::path filename) {
    std::ifstream inputBinaryFile(filename);
    if (!inputBinaryFile) {
        std::cerr << "Unable to open file input.bin";
        exit(1);
    }
    size_t size = getFileSize(filename);
    char *bytes = new char[size];
    inputBinaryFile.read(bytes, size);
    inputBinaryFile.close();
    assert(size % 4 == 0);
    for (long long i = 0; i < size / 4; i++) {
        // Reverse endianness
        char byte0 = bytes[4 * i];
        char byte1 = bytes[4 * i + 1];
        char byte2 = bytes[4 * i + 2];
        char byte3 = bytes[4 * i + 3];
        bytes[4 * i] = byte3;
        bytes[4 * i + 1] = byte2;
        bytes[4 * i + 2] = byte1;
        bytes[4 * i + 3] = byte0;
    }
    return bytes;
}
char *readBytes(fs::path filename) {
    std::ifstream inputBinaryFile(filename);
    if (!inputBinaryFile) {
        std::cerr << "Unable to open file input.bin";
        exit(1);
    }
    size_t size = getFileSize(filename);
    char *bytes = new char[size];
    inputBinaryFile.read(bytes, size);
    return bytes;
}
std::vector<float> readFloatsFromFile(const fs::path &filename) {
    // Open the file for binary input
    std::ifstream inputBinaryFile(filename, std::ios::binary);
    if (!inputBinaryFile) {
        std::cerr << "Unable to open file " << filename << std::endl;
        exit(1);
    }

    // Determine the size of the file
    inputBinaryFile.seekg(0, std::ios::end);
    std::size_t fileSize = inputBinaryFile.tellg();
    inputBinaryFile.seekg(0, std::ios::beg);

    // Calculate the number of float elements
    std::size_t numFloats = fileSize / sizeof(float);

    // Read the floats from the file
    std::vector<float> floats(numFloats);
    if (numFloats > 0) {
        inputBinaryFile.read(reinterpret_cast<char *>(floats.data()), numFloats * sizeof(float));
    }

    // Check for reading errors
    if (inputBinaryFile.fail()) {
        std::cerr << "Error occurred while reading from file " << filename << std::endl;
        exit(1);
    }

    return floats;
}
void writeFloatsToFile(const fs::path &filename, const std::vector<float> &floats) {
    // Open the file for binary output
    std::ofstream outputBinaryFile(filename, std::ios::binary);
    if (!outputBinaryFile) {
        std::cerr << "Unable to open file " << filename << std::endl;
        exit(1);
    }

    // Write the floats to the file
    if (!floats.empty()) {
        outputBinaryFile.write(reinterpret_cast<const char *>(floats.data()), floats.size() * sizeof(float));
    }

    // Check for writing errors
    if (outputBinaryFile.fail()) {
        std::cerr << "Error occurred while writing to file " << filename << std::endl;
        exit(1);
    }
}

// ---------------------------------------------------
// UTIL FUNCTIONS FOR WORKING WITH IEEE REPRESENTATION
// ---------------------------------------------------

uint8_t getExponentFromFloat(float value) {
    unsigned int binaryRepresentation = *reinterpret_cast<unsigned int *>(&value);
    return (binaryRepresentation >> 23) & 0xFF;
}
// Function to print the binary representation of a char
void debug_printBinary(char c) {
    // Iterate over each bit in the char (from most significant bit to least significant bit)
    for (int i = 7; i >= 0; --i) {
        // Extract the bit using bitwise AND and shift
        std::cout << ((c >> i) & 1);
    }
    std::cout << std::endl;
}

// ---------------------------------------------------------
// UTIL FUNCTIONS FOR FAST READING/WRITING BITS FROM BUFFERS
// ---------------------------------------------------------

inline bool readBit(uint8_t *data, long long bitIndex) {
    // Warning: this code can be unsafe with a bad bitIndex
    return (data[bitIndex / 8] >> (7 - (bitIndex % 8))) & 1;
}
inline uint8_t readByteFromBuffer(uint8_t *data, long long startingBitIndex) {
    uint8_t result = 0;
    int mod = startingBitIndex % 8;
    long long div = startingBitIndex / 8;
    // Read until end of byte
    result |= data[div] << mod;
    // Read beginning of next byte
    if (mod != 0) {
        result |= data[div + 1] >> (8 - mod);
    }
    return result;
}
inline uint8_t readByteFromBuffer(uint8_t *data, long long startingBitIndex, int numBits, bool enforceMask = false) {
    uint8_t result = 0;
    int mod = startingBitIndex % 8;
    long long div = startingBitIndex / 8;
    // Read until end of byte
    result |= data[div] << mod;
    // Read beginning of next byte
    if (numBits > 8 - mod) {
        result |= data[div + 1] >> (8 - mod);
    }
    if (enforceMask && numBits < 8) {
        // Mask from the left according to numBits
        // For example, if result is 11010101 and numBits is 4, the result should be 11010000
        result &= ((1 << numBits) - 1) << (8 - numBits);
    }
    return result;
}
inline void writeBit(char *data, long long bitIndex, bool bit) {
    // Warning: this code can be unsafe with a bad bitIndex
    // This code also assumes data is initialized at 0.
    data[bitIndex / 8] |= (bit << (7 - (bitIndex % 8)));
}
inline void writeByteToBuffer(uint8_t byte, char *buffer, long long startingBitIndex) {
    int mod = startingBitIndex % 8;
    long long div = startingBitIndex / 8;
    // Write until end of byte
    buffer[div] |= byte >> mod;
    // Write beginning of next byte
    if (mod != 0) {
        buffer[div + 1] |= byte << (8 - mod);
    }
}
inline void writeByteToBuffer(uint8_t byte, char *buffer, long long startingBitIndex, int numBits) {
    int mod = startingBitIndex % 8;
    long long div = startingBitIndex / 8;
    // Write until end of byte
    buffer[div] |= byte >> mod;
    // Write beginning of next byte
    if (numBits > 8 - mod) {
        buffer[div + 1] |= byte << (8 - mod);
    }
}

// -----------------------------------
// OTHER STUFF
// -----------------------------------

std::string getCurrentTimeFormatted() {
    // Get current time as a time_t object
    time_t t = time(nullptr);
    // Convert to local time
    tm *localTime = localtime(&t);

    // Get current time including milliseconds
    auto now = std::chrono::system_clock::now();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Format the time
    char buffer[20];
    strftime(buffer, sizeof(buffer), "%Y-%m-%d-%H-%M-%S", localTime);

    // Combine the formatted time with milliseconds
    std::ostringstream oss;
    oss << buffer << "-" << std::setw(3) << std::setfill('0') << milliseconds.count();

    return oss.str();
}
size_t ZSTD_compressFile(const fs::path &filePath, int level) {
    // Open file and read content
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        std::cerr << "Failed to read file: " << filePath << std::endl;
        return 0;
    }

    // Prepare for compression
    size_t compressedSize = ZSTD_compressBound(size); // Maximum compressed size
    std::vector<char> compressed(compressedSize);

    size_t const cSize = ZSTD_compress(compressed.data(), compressedSize, buffer.data(), size, level); // Level 1 compression

    if (ZSTD_isError(cSize)) {
        std::cerr << "Compression error: " << ZSTD_getErrorName(cSize) << std::endl;
        return 0;
    }

    return cSize;
}

#endif // UTILS_HPP