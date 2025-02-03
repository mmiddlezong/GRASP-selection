#include "oursol.hpp"
#include "utils.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <libpressio_ext/cpp/libpressio.h>
#include <libpressio_ext/io/posix.h>
#include <math.h>
#include <numeric>
#include <stdio.h>
namespace fs = std::filesystem;

void testCompressor(pressio_data *inputData, std::string compressorName, double absoluteErrorBound,
                    double &compressionRatio, double &totalTime) {
    pressio library;
    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});
    pressio_data output = pressio_data::owning(inputData->dtype(), inputData->dimensions());
    pressio_compressor compressor = library.get_compressor(compressorName);
    if (!compressor) {
        std::cerr << "Compressor not found.\n";
        return;
    }
    const char *metrics_ids[] = {"precise_time", "size", "error_stat"};
    pressio_metrics metrics = library.get_metrics(std::begin(metrics_ids), std::end(metrics_ids));
    compressor->set_options({
        {"pressio:abs", absoluteErrorBound},
    });
    compressor->set_metrics(metrics);
    if (compressor->compress(inputData, &compressed)) {
        std::cerr << compressor->error_msg() << std::endl;
        return;
    }
    if (compressor->decompress(&compressed, &output)) {
        std::cerr << compressor->error_msg() << std::endl;
        return;
    }

    pressio_options results = compressor->get_metrics_results();
    results.get("size:compression_ratio", &compressionRatio);
    uint32_t compressionTimeMicroseconds = 0;
    results.get("precise_time:compress", &compressionTimeMicroseconds);
    uint32_t decompressionTimeMicroseconds = 0;
    results.get("precise_time:decompress", &decompressionTimeMicroseconds);
    totalTime = (double)(compressionTimeMicroseconds + decompressionTimeMicroseconds) / 1000000.0;
}
enum SamplingMethod { MULTIBLOCK, CHUNK, SZ3 };
void predictOptimalCompressor(const fs::path &datasetDir, const std::vector<size_t> &dims, float relativeErrorBound,
                              double sampleFraction, double f, SamplingMethod samplingMethod, bool outputResults = true,
                              bool outputLog = true, bool outputSampleData = false) {
    fs::path outDir = "out" / datasetDir.filename();
    fs::create_directories(outDir);
    // Output optimal compressor results
    std::ofstream resultsFile;
    if (outputResults) {
        resultsFile.open(outDir / ("algo-select-" + getCurrentTimeFormatted() + ".csv"));
        resultsFile << "dataset,filename,filesize,samplingMethod,sampleSize,sampleFraction,relativeErrorBound,f,"
                       "SZ3CompressionRatio,SZ3Throughput,SZ3Efficiency,ZFPCompressionRatio,"
                       "ZFPThroughput,ZFPEfficiency,OurSolCompressionRatio,OurSolThroughput,"
                       "OurSolEfficiency,optimalCompressor,algoSelectTime\n";
    }

    std::string datasetName = datasetDir.filename().string();
    std::cout << "------------------- Predicting optimal compressor for " << datasetName << " -------------------"
              << std::endl;
    std::cout << "Relative error bound: " << relativeErrorBound << std::endl;
    std::cout << "Sample fraction: " << sampleFraction << std::endl;
    std::cout << "f: " << f << std::endl;
    std::cout << "Dimensions: " << dims[0] << "x" << dims[1] << std::endl;

    std::ofstream log;
    if (outputLog) {
        log.open(outDir / ("algo-select-" + getCurrentTimeFormatted() + ".log"));
        log << "------------------- Predicting optimal compressor for " << datasetName << " -------------------"
            << std::endl;
        log << "Relative error bound: " << relativeErrorBound << std::endl;
        log << "Sample fraction: " << sampleFraction << std::endl;
        log << "f: " << f << std::endl;
        log << "Dimensions: " << dims[0] << "x" << dims[1] << std::endl;
    }

    // Measure time to compute sampling indices
    auto startSample = std::chrono::high_resolution_clock::now();

    std::vector<size_t> sampleIndices;
    size_t sampleNumFloats;
    size_t sampleNumBytes;
    std::vector<size_t> sampleDims(dims.size());

    if (sampleFraction > 0.999) {
        for (size_t i = 0; i < dims[0] * dims[1]; i++) {
            sampleIndices.push_back(i);
        }
        sampleNumFloats = dims[0] * dims[1];
        sampleNumBytes = 4 * sampleNumFloats;
        sampleDims = dims;
    } else if (samplingMethod == MULTIBLOCK) {
        std::cout << "Using multi-block sampling" << std::endl;
        // Find aspect ratio
        assert(dims.size() == 2);        // We only support 2D datasets
        int rows = 0;                    // how many rows of sample blocks
        int cols = 0;                    // how many columns of sample blocks
        int preferredMinBlockSize = 100; // a rough estimate of block size to start with, the actual
                                         // block size is adaptive to the size of the dataset
        while (true) {
            rows = dims[1] / preferredMinBlockSize; // Corresponds to dims[1]
            cols = dims[0] / preferredMinBlockSize; // Corresponds to dims[0]

            rows = std::max(rows, 1); // clamp to at least 1 row and col
            cols = std::max(cols, 1);

            if (rows * cols > 20 && rows > 1 && cols > 1) {
                // Too many blocks
                preferredMinBlockSize *= 1.01; // adaptively increase the block size
            } else {
                break;
            }
        }
        // Pre-compute indices to sample
        size_t actualBlockSize = std::min(dims[0] / cols, dims[1] / rows);
        size_t sampleBlockSize = sqrt(sampleFraction * dims[0] * dims[1] / (rows * cols));
        // Clamp sample block size to be at least 1 and at most the actual block size
        sampleBlockSize = std::max((int)sampleBlockSize, 1);
        sampleBlockSize = std::min((int)sampleBlockSize, (int)actualBlockSize);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                size_t topLeftX = j * actualBlockSize + (actualBlockSize - sampleBlockSize) / 2;
                size_t topLeftY = i * actualBlockSize + (actualBlockSize - sampleBlockSize) / 2;
                for (size_t x = topLeftX; x < topLeftX + sampleBlockSize; x++) {
                    for (size_t y = topLeftY; y < topLeftY + sampleBlockSize; y++) {
                        sampleIndices.push_back(y * dims[0] + x);
                    }
                }
            }
        }
        std::sort(sampleIndices.begin(), sampleIndices.end());

        sampleNumFloats = sampleBlockSize * sampleBlockSize * rows * cols;
        sampleNumBytes = 4 * sampleNumFloats;
        // std::cout << "Sample size: " << sampleSize << std::endl;
        sampleDims = {sampleBlockSize * cols, sampleBlockSize * rows};

        std::cout << "Sample aspect ratio: " << rows << "x" << cols << std::endl;
        std::cout << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        if (outputLog) {
            log << "Sample aspect ratio: " << rows << "x" << cols << std::endl;
            log << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        }
    } else if (samplingMethod == CHUNK) {
        std::cout << "Using single block sampling" << std::endl;
        // Find aspect ratio
        assert(dims.size() == 2); // We only support 2D datasets
        size_t sampleHeight = dims[1] * sqrt(sampleFraction);
        size_t sampleWidth = dims[0] * sqrt(sampleFraction);
        sampleHeight = std::max(sampleHeight, (size_t)1);
        sampleWidth = std::max(sampleWidth, (size_t)1);
        size_t topLeftX = (dims[0] - sampleWidth) / 2;
        size_t topLeftY = (dims[1] - sampleHeight) / 2;
        for (size_t x = topLeftX; x < topLeftX + sampleWidth; x++) {
            for (size_t y = topLeftY; y < topLeftY + sampleHeight; y++) {
                sampleIndices.push_back(y * dims[0] + x);
            }
        }
        std::sort(sampleIndices.begin(), sampleIndices.end());
        sampleNumFloats = sampleHeight * sampleWidth;
        sampleNumBytes = 4 * sampleNumFloats;
        sampleDims = {sampleWidth, sampleHeight};
        std::cout << "Sample dimensions: " << sampleHeight << "x" << sampleWidth << std::endl;
        std::cout << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        if (outputLog) {
            log << "Sample dimensions: " << sampleHeight << "x" << sampleWidth << std::endl;
            log << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        }
    } else if (samplingMethod == SZ3) {
        std::cout << "Using SZ3 sampling" << std::endl;
        size_t num = std::accumulate(dims.begin(), dims.end(), (size_t)1, std::multiplies<size_t>());

        size_t dmin = *std::min_element(dims.begin(), dims.end());
        size_t sampling_block = dmin;

        while (true) {
            size_t sample_n = 1;
            for (auto dim : dims) {
                sample_n *= dim / dmin * 2 * sampling_block;
            }
            float cur_sampling_ratio = sample_n * 1.0 / num;
            if (cur_sampling_ratio > sampleFraction) {
                sampling_block--;
            } else {
                break;
            }
        }

        if (sampling_block * 2 > dmin) {
            sampling_block = dmin / 2;
        }
        size_t b0 = dims[1] / dmin;
        size_t b1 = dims[0] / dmin;
        sampleDims[1] = b0 * 2 * sampling_block;
        sampleDims[0] = b1 * 2 * sampling_block;
        size_t di, dj;
        sampleNumFloats = sampleDims[0] * sampleDims[1];
        sampleNumBytes = 4 * sampleNumFloats;
        for (size_t bi = 0; bi < b0; bi++) {
            for (size_t bj = 0; bj < b1; bj++) {
                for (size_t i = 0; i < 2 * sampling_block; i++) {
                    for (size_t j = 0; j < 2 * sampling_block; j++) {
                        di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                        dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                        sampleIndices.push_back((bi * dmin + di) * dims[0] + bj * dmin + dj);
                    }
                }
            }
        }
        std::sort(sampleIndices.begin(), sampleIndices.end());
        std::cout << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        if (outputLog) {
            log << "Sample size: " << sampleNumBytes << " bytes" << std::endl;
        }
    }

    auto endSample = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSample = endSample - startSample;
    std::cout << "Sampling time: " << elapsedSample.count() << "s" << std::endl;
    if (outputLog) {
        log << "Sampling time: " << elapsedSample.count() << "s" << std::endl;
    }

    // Output sample indices to a file
    if (outputSampleData) {
        std::ofstream sampleFile(outDir / ("sample_indices-" + std::to_string(static_cast<int>(samplingMethod)) + "-" +
                                           getCurrentTimeFormatted() + ".txt"));
        for (size_t i : sampleIndices) {
            sampleFile << i << " ";
        }
        sampleFile << std::endl;
    }

    std::vector<std::string> files;

    for (const auto &entry : fs::directory_iterator(datasetDir)) {
        if (entry.is_regular_file()) {
            files.push_back(entry.path().filename().string());
        }
    }

    pressio library;
    pressio_data metadata = pressio_data::empty(pressio_float_dtype, dims);

    // Store data for aggregates
    std::vector<float> SZ3CompressionRatios;
    std::vector<float> SZ3Throughputs;
    std::vector<float> ZFPCompressionRatios;
    std::vector<float> ZFPThroughputs;
    std::vector<float> OurSolCompressionRatios;
    std::vector<float> OurSolThroughputs;
    int sz3Count = 0;
    int zfpCount = 0;
    int ourSolCount = 0;

    for (const std::string &filename : files) {
        fs::path inputPath = datasetDir / filename;
        std::cout << "------------------- Processing " << inputPath << std::endl;
        if (outputLog) {
            log << "------------------- Processing " << inputPath << std::endl;
        }

        std::ifstream inputFile(inputPath, std::ios::binary);
        if (!inputFile.is_open()) {
            std::cerr << "File at " << inputPath << " could not be opened. Error: " << std::strerror(errno) << "\n";
            continue;
        }
        size_t filesize = getFileSize(inputPath);
        char *data = new char[filesize];
        inputFile.read(data, filesize);
        inputFile.close();

        // Calculate value range
        std::vector<float> originalFloats = readFloatsFromFile(inputPath);
        float minFloat = *std::min_element(originalFloats.begin(), originalFloats.end());
        float maxFloat = *std::max_element(originalFloats.begin(), originalFloats.end());
        float range = maxFloat - minFloat;

        // Track total time spent (does not include finding the absolute error bound through finding
        // range)
        auto start = std::chrono::high_resolution_clock::now();

        // Sample data
        std::vector<float> sampleDataPoints;
        for (size_t i : sampleIndices) {
            sampleDataPoints.push_back((reinterpret_cast<float *>(data))[i]);
        }
        if (outputSampleData) {
            std::cout << "Writing sample data to file" << std::endl;
            std::cout << "Sample dims: " << sampleDims[0] << "x" << sampleDims[1] << std::endl;
            if (outputLog) {
                log << "Writing sample data to file" << std::endl;
                log << "Sample dims: " << sampleDims[0] << "x" << sampleDims[1] << std::endl;
            }
            writeFloatsToFile(outDir /
                                  (filename + "-" + std::to_string(static_cast<int>(samplingMethod)) + "-" +
                                   std::to_string(sampleDims[0]) + "x" + std::to_string(sampleDims[1]) + ".sample"),
                              sampleDataPoints);
        }
        char *sampleData = reinterpret_cast<char *>(sampleDataPoints.data());
        delete[] data;

        // Run each compressor on sampled data
        pressio_data inputData = pressio_data::copy(pressio_float_dtype, sampleData, sampleDims);

        double absoluteErrorBound = (double)relativeErrorBound * (double)range;

        // SZ3
        double SZ3CompressionRatio;
        double SZ3TotalTime;
        testCompressor(&inputData, "sz3", absoluteErrorBound, SZ3CompressionRatio, SZ3TotalTime);
        double SZ3Throughput = (double)sampleNumBytes / SZ3TotalTime;

        // ZFP
        double ZFPCompressionRatio;
        double ZFPTotalTime;
        testCompressor(&inputData, "zfp", absoluteErrorBound, ZFPCompressionRatio, ZFPTotalTime);
        double ZFPThroughput = (double)sampleNumBytes / ZFPTotalTime;

        // OurSol
        char *encodedData = new char[sampleNumBytes]();

        // Group bits
        auto c0 = std::chrono::high_resolution_clock::now();
        long long numBytesToWrite = groupBitsWithLoss(sampleData, encodedData, sampleNumFloats, absoluteErrorBound);
        auto c1 = std::chrono::high_resolution_clock::now();

        // Save to file
        fs::path outputPath = outDir / (filename + ".oursol");
        saveToFile(sampleNumBytes, absoluteErrorBound, encodedData, outputPath, numBytesToWrite);

        auto c4 = std::chrono::high_resolution_clock::now();
        size_t compressedSize = ZSTD_compressFile(outputPath, 1);
        auto c5 = std::chrono::high_resolution_clock::now();

        delete[] encodedData;

        // Decoding
        std::vector<float> decodedFloats;

        auto c2 = std::chrono::high_resolution_clock::now();
        decodeFile(decodedFloats, outputPath);
        auto c3 = std::chrono::high_resolution_clock::now();
        fs::remove(outputPath);

        // Calculate time in seconds
        double OurSolCompressionTime =
            std::chrono::duration<double>(c1 - c0).count() + std::chrono::duration<double>(c5 - c4).count();
        double OurSolDecompressionTime = std::chrono::duration<double>(c3 - c2).count();
        double OurSolTotalTime = OurSolCompressionTime + OurSolDecompressionTime;

        // Print each time for debug
        // std::cout << "group time: " << std::chrono::duration<double>(c1 - c0).count() << "s" << std::endl;
        // std::cout << "zstd time: " << std::chrono::duration<double>(c5 - c4).count() << "s" << std::endl;
        // std::cout << "decompression time: " << OurSolDecompressionTime << "s" << std::endl;

        // std::cout << sampleNumBytes << "\n";
        // std::cout << (int)(sampleFraction * filesize) << "\n";
        double OurSolCompressionRatio = (double)sampleNumBytes / compressedSize;
        double OurSolThroughput = (double)sampleNumBytes / OurSolTotalTime;

        // Store in list
        SZ3CompressionRatios.push_back(SZ3CompressionRatio);
        SZ3Throughputs.push_back(SZ3Throughput);
        ZFPCompressionRatios.push_back(ZFPCompressionRatio);
        ZFPThroughputs.push_back(ZFPThroughput);
        OurSolCompressionRatios.push_back(OurSolCompressionRatio);
        OurSolThroughputs.push_back(OurSolThroughput);

        // Calculate metric
        double SZ3Efficiency = 1.0 - ((f * 1000000) / SZ3Throughput) - (1.0 / SZ3CompressionRatio);
        double ZFPEfficiency = 1.0 - ((f * 1000000) / ZFPThroughput) - (1.0 / ZFPCompressionRatio);
        double OurSolEfficiency = 1.0 - ((f * 1000000) / OurSolThroughput) - (1.0 / OurSolCompressionRatio);

        // Output optimal compressor using max
        std::string optimalCompressor;
        if (SZ3Efficiency > ZFPEfficiency && SZ3Efficiency > OurSolEfficiency) {
            optimalCompressor = "SZ3";
            sz3Count++;
        } else if (ZFPEfficiency > SZ3Efficiency && ZFPEfficiency > OurSolEfficiency) {
            optimalCompressor = "ZFP";
            zfpCount++;
        } else {
            optimalCompressor = "OurSol";
            ourSolCount++;
        }

        // End timer
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        // Calculate algo select throughput
        double algoSelectThroughput = (double)filesize / elapsed.count();

        // Print/log stuff
        std::cout << "SZ3: " << SZ3TotalTime << "s, " << SZ3Throughput / 1000000 << " MB/s, " << SZ3CompressionRatio
                  << " compression ratio" << std::endl;
        std::cout << "ZFP: " << ZFPTotalTime << "s, " << ZFPThroughput / 1000000 << " MB/s, " << ZFPCompressionRatio
                  << " compression ratio" << std::endl;
        std::cout << "OurSol: " << OurSolTotalTime << "s, " << OurSolThroughput / 1000000 << " MB/s, "
                  << OurSolCompressionRatio << " compression ratio" << std::endl;
        std::cout << "SZ3 efficiency: " << SZ3Efficiency << (optimalCompressor == "SZ3" ? " <--" : "") << std::endl;
        std::cout << "ZFP efficiency: " << ZFPEfficiency << (optimalCompressor == "ZFP" ? " <--" : "") << std::endl;
        std::cout << "OurSol efficiency: " << OurSolEfficiency << (optimalCompressor == "OurSol" ? " <--" : "")
                  << std::endl;
        std::cout << "Optimal compressor for " << filename << " is " << optimalCompressor << "." << std::endl;
        std::cout << "Total time: " << elapsed.count() << "s" << std::endl;
        std::cout << "Algorithm selection throughput: " << algoSelectThroughput / 1000000 << " MB/s" << std::endl;
        if (outputLog) {
            log << "SZ3: " << SZ3TotalTime << "s, " << SZ3Throughput / 1000000 << " MB/s, " << SZ3CompressionRatio
                << " compression ratio" << std::endl;
            log << "ZFP: " << ZFPTotalTime << "s, " << ZFPThroughput / 1000000 << " MB/s, " << ZFPCompressionRatio
                << " compression ratio" << std::endl;
            log << "OurSol: " << OurSolTotalTime << "s, " << OurSolThroughput / 1000000 << " MB/s, "
                << OurSolCompressionRatio << " compression ratio" << std::endl;
            log << "SZ3 efficiency: " << SZ3Efficiency << (optimalCompressor == "SZ3" ? " <--" : "") << std::endl;
            log << "ZFP efficiency: " << ZFPEfficiency << (optimalCompressor == "ZFP" ? " <--" : "") << std::endl;
            log << "OurSol efficiency: " << OurSolEfficiency << (optimalCompressor == "OurSol" ? " <--" : "")
                << std::endl;
            log << "Optimal compressor for " << filename << " is " << optimalCompressor << "." << std::endl;
            log << "Total time: " << elapsed.count() << "s" << std::endl;
            log << "Algorithm selection throughput: " << algoSelectThroughput / 1000000 << " MB/s" << std::endl;
        }

        // Output results
        if (outputResults) {
            resultsFile << datasetName << "," << filename << "," << filesize << "," << samplingMethod << ","
                        << sampleNumBytes << "," << sampleFraction << "," << relativeErrorBound << "," << f << ","
                        << SZ3CompressionRatio << "," << SZ3Throughput / 1000000 << "," << SZ3Efficiency << ","
                        << ZFPCompressionRatio << "," << ZFPThroughput / 1000000 << "," << ZFPEfficiency << ","
                        << OurSolCompressionRatio << "," << OurSolThroughput / 1000000 << "," << OurSolEfficiency << ","
                        << optimalCompressor << "," << elapsed.count() << "\n";
        }
    }
    // Find harmonic mean of compression ratios and throughputs
    double SZ3CompressionRatioHarmonicMean = 0;
    for (float ratio : SZ3CompressionRatios) {
        SZ3CompressionRatioHarmonicMean += 1.0 / ratio;
    }
    SZ3CompressionRatioHarmonicMean = (double)SZ3CompressionRatios.size() / SZ3CompressionRatioHarmonicMean;
    double SZ3ThroughputHarmonicMean = 0;
    for (float throughput : SZ3Throughputs) {
        SZ3ThroughputHarmonicMean += 1.0 / throughput;
    }
    SZ3ThroughputHarmonicMean = (double)SZ3Throughputs.size() / SZ3ThroughputHarmonicMean;

    double ZFPCompressionRatioHarmonicMean = 0;
    for (float ratio : ZFPCompressionRatios) {
        ZFPCompressionRatioHarmonicMean += 1.0 / ratio;
    }
    ZFPCompressionRatioHarmonicMean = (double)ZFPCompressionRatios.size() / ZFPCompressionRatioHarmonicMean;
    double ZFPThroughputHarmonicMean = 0;
    for (float throughput : ZFPThroughputs) {
        ZFPThroughputHarmonicMean += 1.0 / throughput;
    }
    ZFPThroughputHarmonicMean = (double)ZFPThroughputs.size() / ZFPThroughputHarmonicMean;

    double OurSolCompressionRatioHarmonicMean = 0;
    for (float ratio : OurSolCompressionRatios) {
        OurSolCompressionRatioHarmonicMean += 1.0 / ratio;
    }
    OurSolCompressionRatioHarmonicMean = (double)OurSolCompressionRatios.size() / OurSolCompressionRatioHarmonicMean;
    double OurSolThroughputHarmonicMean = 0;
    for (float throughput : OurSolThroughputs) {
        OurSolThroughputHarmonicMean += 1.0 / throughput;
    }
    OurSolThroughputHarmonicMean = (double)OurSolThroughputs.size() / OurSolThroughputHarmonicMean;

    // Print aggregate results
    std::string mostSelectedCompressor = (sz3Count > zfpCount ? (sz3Count > ourSolCount ? "SZ3" : "OurSol")
                                                              : (zfpCount > ourSolCount ? "ZFP" : "OurSol"));
    std::cout << "------------------- Aggregate Results -------------------" << std::endl;
    std::cout << "SZ3: " << SZ3CompressionRatioHarmonicMean << " compression ratio, "
              << SZ3ThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
    std::cout << "ZFP: " << ZFPCompressionRatioHarmonicMean << " compression ratio, "
              << ZFPThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
    std::cout << "OurSol: " << OurSolCompressionRatioHarmonicMean << " compression ratio, "
              << OurSolThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
    std::cout << "Most selected compressor: " << mostSelectedCompressor << std::endl;
    if (outputLog) {
        log << "------------------- Aggregate Results -------------------" << std::endl;
        log << "SZ3: " << SZ3CompressionRatioHarmonicMean << " compression ratio, "
            << SZ3ThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
        log << "ZFP: " << ZFPCompressionRatioHarmonicMean << " compression ratio, "
            << ZFPThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
        log << "OurSol: " << OurSolCompressionRatioHarmonicMean << " compression ratio, "
            << OurSolThroughputHarmonicMean / 1000000 << " MB/s throughput" << std::endl;
        log << "Most selected compressor: " << mostSelectedCompressor << std::endl;
    }
    if (outputResults) {
        // Calculate total file size
        size_t totalFileSize = dims[0] * dims[1] * 4 * files.size();
        size_t totalSampleSize = sampleNumBytes * files.size();
        resultsFile << datasetName << "," << datasetName << "," << totalFileSize << "," << samplingMethod << ","
                    << totalSampleSize << "," << sampleFraction << "," << relativeErrorBound << "," << f << ","
                    << SZ3CompressionRatioHarmonicMean << "," << SZ3ThroughputHarmonicMean / 1000000 << "," << 0 << ","
                    << ZFPCompressionRatioHarmonicMean << "," << ZFPThroughputHarmonicMean / 1000000 << "," << 0 << ","
                    << OurSolCompressionRatioHarmonicMean << "," << OurSolThroughputHarmonicMean / 1000000 << "," << 0
                    << "," << mostSelectedCompressor << ",0\n";
    }
}
struct Dataset {
    fs::path path;
    std::vector<size_t> dimensions;
};
int main(int, char **) {
    std::vector<float> relativeErrorBounds = {1E-6};
    std::vector<float> fSpeeds = {40.0};
    std::vector<float> sampleFractions = {0.1, 0.01, 0.001};
    fs::path dataPath(DATADIR);
    std::vector<SamplingMethod> samplingMethods = {MULTIBLOCK, CHUNK, SZ3};
    // std::vector<Dataset> datasets = {{dataPath / "CESM-ATM", {3600, 1800}}, {dataPath / "EXAALT", {3137, 5423}},
    // {dataPath / "ISABEL", {500, 50000}}}; std::vector<Dataset> datasets = {{dataPath / "ISABEL", {500, 50000}}};
    // std::vector<Dataset> datasets = {{dataPath / "NYX", {512, 512 * 512}}, {dataPath / "SCALE", {1200, 1200 * 98}}};
    // std::vector<Dataset> datasets = {{dataPath / "CESM-ATM", {3600, 1800}}};
    std::vector<Dataset> datasets = {{dataPath / "EXAALT", {3137, 5423}}};
    for (const Dataset &dataset : datasets) {
        for (SamplingMethod samplingMethod : samplingMethods) {
            for (float relativeErrorBound : relativeErrorBounds) {
                for (float f : fSpeeds) {
                    for (float sampleFraction : sampleFractions) {
                        predictOptimalCompressor(dataset.path, dataset.dimensions, relativeErrorBound, sampleFraction,
                                                 f, samplingMethod, true, true, true);
                    }
                }
            }
        }
    }
    return 0;
}
