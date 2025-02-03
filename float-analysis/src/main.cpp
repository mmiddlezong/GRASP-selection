#include "oursol.hpp"

void compressDataset(const fs::path &datasetDir, float relativeErrorBound, bool outputErrors = false) {
    std::vector<std::string> dataFiles;

    fs::path outDir = "out" / datasetDir.filename();
    fs::create_directories(outDir);
    for (const auto &entry : fs::directory_iterator(datasetDir)) {
        if (entry.is_regular_file()) {
            dataFiles.push_back(entry.path().filename().string());
        }
    }
    // Declare totals
    long long datasetSize = 0;
    long long compressedDatasetSize = 0;
    long long justZSTDDatasetSize = 0;
    float totalCompressionTime = 0.0f;
    float totalDecompressionTime = 0.0f;
    std::vector<double> normedRMSEs;
    std::vector<float> fileCompressionRatios;

    // Dataset compression log
    std::string compressionLogName = "compression-" +
                                     getCurrentTimeFormatted() +
                                     ".log"; // Based on current time
    std::ofstream compressionLog(outDir / compressionLogName);
    if (!compressionLog.is_open()) {
        throw std::runtime_error("File could not be opened");
    }

    for (const std::string &filename : dataFiles) {
        fs::path inputPath = datasetDir / filename;
        fs::path outputPath = outDir / (filename + "-grouped.bin");

        auto c0 = std::chrono::high_resolution_clock::now();

        char *bytes = readBytes(inputPath);
        std::vector<float> originalFloats = readFloatsFromFile(inputPath);
        float minFloat = *std::min_element(originalFloats.begin(), originalFloats.end());
        float maxFloat = *std::max_element(originalFloats.begin(), originalFloats.end());
        float range = maxFloat - minFloat;
        float absoluteErrorBound = relativeErrorBound * range;

        auto c1 = std::chrono::high_resolution_clock::now();

        long long len = getFileSize(inputPath);
        long long n = len / 4;
        char *encodedData = new char[4 * n]();
        long long numBytesToWrite = groupBitsWithLoss(bytes, encodedData, n, absoluteErrorBound);

        auto c2 = std::chrono::high_resolution_clock::now();

        saveToFile(len, absoluteErrorBound, encodedData, outputPath, numBytesToWrite);

        auto c3 = std::chrono::high_resolution_clock::now();

        size_t originalSize = getFileSize(inputPath);
        size_t compressedSize = ZSTD_compressFile(outputPath, 1);

        auto c4 = std::chrono::high_resolution_clock::now();

        size_t justZSTDSize = ZSTD_compressFile(inputPath, 1);

        auto c5 = std::chrono::high_resolution_clock::now();

        // Decoding
        std::vector<float> decodedFloats;
        decodeFile(decodedFloats, outputPath);

        auto c6 = std::chrono::high_resolution_clock::now();

        // Verify the compression was lossless
        assert(decodedFloats.size() == originalFloats.size());

        if (outputErrors) {
            fs::path errorsOutputPath = outDir / (filename + "-1E-" + std::to_string(-static_cast<int>(round(log10(relativeErrorBound)))) + "-errors.dat");
            std::vector<float> errors(decodedFloats.size());
            for (long long i = 0; i < decodedFloats.size(); i++) {
                if (std::abs(decodedFloats[i] - originalFloats[i]) >= absoluteErrorBound) {
                    std::cout << "Decoded float: " << decodedFloats[i] << std::endl;
                    std::cout << "Original float: " << originalFloats[i] << std::endl;
                    std::cout << "Error: " << std::abs(decodedFloats[i] - originalFloats[i]) << " >= " << absoluteErrorBound << std::endl;
                }
                assert(std::abs(decodedFloats[i] - originalFloats[i]) < absoluteErrorBound);
                errors[i] = decodedFloats[i] - originalFloats[i];
            }
            writeFloatsToFile(errorsOutputPath, errors);
        } else {
            for (long long i = 0; i < decodedFloats.size(); i++) {
                if (std::abs(decodedFloats[i] - originalFloats[i]) >= absoluteErrorBound) {
                    std::cout << "Decoded float: " << decodedFloats[i] << std::endl;
                    std::cout << "Original float: " << originalFloats[i] << std::endl;
                    std::cout << "Error: " << std::abs(decodedFloats[i] - originalFloats[i]) << " >= " << absoluteErrorBound << std::endl;
                }
                // assert(std::abs(decodedFloats[i] - originalFloats[i]) < absoluteErrorBound);
            }
        }

        auto c7 = std::chrono::high_resolution_clock::now();

        // Calculate PSNR
        double mse = 0.0f;
        for (long long i = 0; i < decodedFloats.size(); i++) {
            mse += (decodedFloats[i] - originalFloats[i]) * (decodedFloats[i] - originalFloats[i]);
        }
        mse /= decodedFloats.size();

        double normedRMSE = sqrt(mse) / range;
        double psnr = 20 * log10(range) - 10 * log10(mse);

        auto c8 = std::chrono::high_resolution_clock::now();

        // Get time in ms
        auto readTime = std::chrono::duration_cast<std::chrono::microseconds>(c1 - c0).count() / 1000.0f;
        auto groupTime = std::chrono::duration_cast<std::chrono::microseconds>(c2 - c1).count() / 1000.0f;
        auto writeTime = std::chrono::duration_cast<std::chrono::microseconds>(c3 - c2).count() / 1000.0f;
        auto ZSTDTime = std::chrono::duration_cast<std::chrono::microseconds>(c4 - c3).count() / 1000.0f;
        auto calculatingJustZSTDTime = std::chrono::duration_cast<std::chrono::microseconds>(c5 - c4).count() / 1000.0f;
        auto decodingTime = std::chrono::duration_cast<std::chrono::microseconds>(c6 - c5).count() / 1000.0f;
        auto errorProcessingTime = std::chrono::duration_cast<std::chrono::microseconds>(c7 - c6).count() / 1000.0f;
        auto psnrTime = std::chrono::duration_cast<std::chrono::microseconds>(c8 - c7).count() / 1000.0f;
        auto compressionTime = groupTime + ZSTDTime;

        // Update running totals for the entire dataset
        datasetSize += originalSize;
        compressedDatasetSize += compressedSize;
        totalCompressionTime += compressionTime;
        totalDecompressionTime += decodingTime;
        justZSTDDatasetSize += justZSTDSize;
        normedRMSEs.push_back(normedRMSE);
        fileCompressionRatios.push_back((float)originalSize / compressedSize);

        compressionLog << "Compressed " << inputPath << ": size " << originalSize << " -> " << compressedSize << " (" << justZSTDSize << " with just ZSTD)"
                       << ", time " << compressionTime << " ms"
                       << ", PSNR " << psnr << std::endl;
        compressionLog << "Profiling: " << readTime << " ms (read) + "
                       << groupTime << " ms (group) + "
                       << writeTime << " ms (write) + "
                       << ZSTDTime << " ms (ZSTD) + "
                       << calculatingJustZSTDTime << " ms (just ZSTD) + "
                       << decodingTime << " ms (decoding) + "
                       << errorProcessingTime << " ms (error processing) + "
                       << psnrTime << " ms (PSNR)" << std::endl;
        compressionLog << "Other metrics: normed RMSE " << normedRMSE << ", bitrate " << 32.0f * (float)compressedSize / originalSize << " bits/float, ratio " << (float)originalSize / compressedSize << std::endl;

        delete[] bytes;
        delete[] encodedData;
        fs::remove(outputPath);
    }

    // Calculate the harmonic mean of PSNRs
    double overallNormedRMSESquared = 0.0f;
    for (double normedRMSE : normedRMSEs) {
        overallNormedRMSESquared += normedRMSE * normedRMSE;
    }
    overallNormedRMSESquared /= normedRMSEs.size();
    double overallNormedRMSE = sqrt(overallNormedRMSESquared);
    double overallPSNR = -20 * log10(overallNormedRMSE);

    // Find the file with the file compression ratio most similar to overall
    float overallCompressionRatio = (float)datasetSize / compressedDatasetSize;
    float minDiff = 1E9;
    int minDiffIndex = -1;
    for (int i = 0; i < fileCompressionRatios.size(); i++) {
        float diff = std::abs(fileCompressionRatios[i] - overallCompressionRatio);
        if (diff < minDiff) {
            minDiff = diff;
            minDiffIndex = i;
        }
    }
    // Fetch the file name using the dataFiles
    std::string mostSimilarFile = dataFiles[minDiffIndex];

    compressionLog << "FINISHED! Summary:\n";
    compressionLog << "- Dataset size: " << datasetSize << " B\n";
    compressionLog << "- Compressed dataset size: " << compressedDatasetSize
                   << " B\n";
    compressionLog << "- Compression time: " << totalCompressionTime << " ms\n";
    compressionLog << "- Decompression time: " << totalDecompressionTime << " ms\n";
    compressionLog << "- Relative error bound: " << relativeErrorBound << "\n";
    compressionLog << "- Most representative file: " << mostSimilarFile << "\n";
    compressionLog << "METRICS:\n";
    compressionLog << "- Overall compression ratio: "
                   << overallCompressionRatio << " (" << (float)datasetSize / justZSTDDatasetSize << " with just ZSTD)\n";
    compressionLog << "- Overall throughput (group + ZSTD): "
                   << (float)datasetSize / totalCompressionTime << " KB/s\n";
    compressionLog << "- Overall decompression speed: "
                   << (float)datasetSize / totalDecompressionTime << " KB/s\n";
    compressionLog << "- Overall PSNR: " << overallPSNR << "\n";
    compressionLog << "- Overall bitrate: " << 32.0f * (float)compressedDatasetSize / datasetSize << " bits/float\n";
    compressionLog.close();
}

int main() {
    fs::path realDatasetsDir("/home/appllo/projects/sdr-projects/huffman/real-datasets/");
    std::vector<float> relativeErrorBounds = {1E-5, 1E-6, 1E-7, 1E-8};
    std::vector<fs::path> datasets = { realDatasetsDir / "CESM-ATM" };
    // std::vector<fs::path> datasets = {realDatasetsDir / "ISABEL", realDatasetsDir / "NYX", realDatasetsDir / "SCALE", realDatasetsDir / "QMCPACK"};
    for (float relativeErrorBound : relativeErrorBounds) {
        for (const fs::path &dataset : datasets) {
            compressDataset(dataset, relativeErrorBound, false);
        }
    }
    return 0;
}