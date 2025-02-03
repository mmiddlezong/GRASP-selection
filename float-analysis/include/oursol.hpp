#ifndef OURSOL_HPP
#define OURSOL_HPP

#include "utils.hpp"
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

// Returns number of bytes in encodedData
long long groupBitsWithLoss(char *bytes, char *encodedData, long long n, float absoluteErrorBound) {
    // Compute the exponent part of a float representing absoluteErrorBound
    // BELOW CODE PROBABLY DEPENDS ON LITTLE ENDIAN
    uint8_t errorBoundExp = getExponentFromFloat(absoluteErrorBound);
    long long curMantissaIndex = 9 * n;
    for (long long i = 0; i < n; i++) {
        // seeeeeee
        // emmmmmmm
        // mmmmmmmm
        // mmmmmmmm
        long long startingByte = 4 * i;

        // Write sign bit
        writeBit(encodedData, i, (bytes[startingByte + 3] & 0x80) != 0);

        // Write exponent bits
        uint8_t exponentByte = bytes[startingByte + 3] << 1;
        // exponentByte: eeeeeee0
        exponentByte |= (bytes[startingByte + 2] >> 7 & 0x01);
        // exponentByte: eeeeeeee
        writeByteToBuffer(exponentByte, encodedData, n + i * 8);

        // Compute how many mantissa bits to store
        uint8_t valueExp = (bytes[startingByte + 3] << 1) & 0xFE;
        valueExp |= ((bytes[startingByte + 2] >> 7 & 0x01));
        int numMantissaBits = std::max(0, valueExp - errorBoundExp); // Clamp to [0, 23]
        numMantissaBits = std::min(23, numMantissaBits);

        // Directly write the mantissa bits, considering the offset
        if (numMantissaBits > 0) {
            uint8_t modifiedByte = bytes[startingByte + 2] << 1; // Remove the e bit from this byte
            // modifiedByte: mmmmmmm0
            int m = std::min(numMantissaBits, 7);
            modifiedByte &= (1 << m) - 1 << (8 - m); // Keep the first m bits
            // modifiedByte: mmmm0000 (if numMantissaBits is 4)
            writeByteToBuffer(modifiedByte, encodedData, curMantissaIndex, m);

            if (numMantissaBits >= 8) {
                uint8_t mantissaByte = bytes[startingByte + 1];
                m = std::min(numMantissaBits - 7, 8);
                mantissaByte &= (1 << m) - 1 << (8 - m);

                writeByteToBuffer(mantissaByte, encodedData, curMantissaIndex + 7, m);
            }
            if (numMantissaBits >= 16) {
                uint8_t mantissaByte = bytes[startingByte];
                m = std::min(numMantissaBits - 15, 8);
                mantissaByte &= (1 << m) - 1 << (8 - m);

                writeByteToBuffer(mantissaByte, encodedData, curMantissaIndex + 15, m);
            }
            curMantissaIndex += numMantissaBits;
        }
    }
    return (curMantissaIndex + 7) / 8;
}
void saveToFile(long long originalFileSize, float absoluteErrorBound, char *bytes, fs::path &filePath, long long numBytes) {
    std::ofstream file(filePath, std::ios::binary);
    if (file.is_open()) {
        // Write the buffer to the file
        file.write(reinterpret_cast<const char *>(&originalFileSize), sizeof(originalFileSize));
        file.write(reinterpret_cast<const char *>(&absoluteErrorBound), sizeof(absoluteErrorBound));
        file.write(bytes, numBytes);
        file.close();
    } else {
        std::cerr << "Failed to open file for writing: " << filePath << std::endl;
    }
}
void decodeFile(std::vector<float> &floats, std::string filename) {
    std::ifstream file(filename, std::ios::binary);

    // Get the size of the file
    file.seekg(0, std::ios::end);
    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // Read in the original file size
    long long originalFileSize;
    file.read(reinterpret_cast<char *>(&originalFileSize), sizeof(originalFileSize));
    // Read in the absolute error bound
    float absoluteErrorBound;
    file.read(reinterpret_cast<char *>(&absoluteErrorBound), sizeof(absoluteErrorBound));

    int headerSize = 12;
    
    // Read the remainder of file into this
    uint8_t *encodedBuffer = new uint8_t[fileSize - headerSize];
    file.read(reinterpret_cast<char *>(encodedBuffer), fileSize - headerSize);
    file.close();

    long long n = originalFileSize / 4;

    // BELOW CODE DEPENDS ON LITTLE ENDIAN
    uint8_t errorBoundExp = getExponentFromFloat(absoluteErrorBound);
    long long curMantissaIndex = 9 * n;
    for (long long i = 0; i < n; i++) {
        uint8_t byte1 = (readBit(encodedBuffer, i) << 7) & 0x80;
        byte1 |= (readByteFromBuffer(encodedBuffer, n + i * 8, 7) >> 1) & 0x7F;

        uint8_t byte2 = (readBit(encodedBuffer, n + i * 8 + 7) << 7) & 0x80;

        // Compute number of mantissa bits to read
        uint8_t valueExp = (byte1 << 1) & 0xFE;
        valueExp |= ((byte2 >> 7) & 0x01);
        int numMantissaBits = std::max(0, valueExp - errorBoundExp); // Clamp to [0, 23]
        numMantissaBits = std::min(23, numMantissaBits);

        uint8_t byte3 = 0;
        uint8_t byte4 = 0;
        // Directly read the mantissa bits, considering the offset
        if (numMantissaBits > 0) {
            byte2 |= (readByteFromBuffer(encodedBuffer, curMantissaIndex, std::min(numMantissaBits, 7), true) >> 1) & 0x7F;
            if (numMantissaBits >= 8) {
                byte3 |= readByteFromBuffer(encodedBuffer, curMantissaIndex + 7, std::min(numMantissaBits - 7, 8), true);
            }
            if (numMantissaBits >= 16) {
                byte4 |= readByteFromBuffer(encodedBuffer, curMantissaIndex + 15, std::min(numMantissaBits - 15, 8), true);
            }
            curMantissaIndex += numMantissaBits;
        }

        // Construct the float from the bytes
        float result;
        unsigned char *bytePtr = reinterpret_cast<unsigned char *>(&result);

        // Little endian order
        bytePtr[0] = byte4;
        bytePtr[1] = byte3;
        bytePtr[2] = byte2;
        bytePtr[3] = byte1;

        // Remove the bias towards 0s in the compressed mantissa
        if (result > absoluteErrorBound / 2) {
            result += absoluteErrorBound / 2;
        } else if (result < -absoluteErrorBound / 2) {
            result -= absoluteErrorBound / 2;
        }
        floats.push_back(result);
    }
    delete[] encodedBuffer;
}

#endif // OURSOL_HPP