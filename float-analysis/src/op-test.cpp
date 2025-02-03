#include <iostream>
#include <chrono>

void writeBit(char *data, int bitIndex, bool bit) {
    // Warning: this code can be unsafe with a bad bitIndex
    // This code also assumes data is initialized at 0.
    data[bitIndex / 8] |= (bit << (7 - (bitIndex % 8)));
}

int main() {
    char *buffer = new char[8000000];
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 8000000; i++) {
        writeBit(buffer, i, true);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    auto ms = elapsed.count();
    std::cout << "Time taken for 8,000,000 bit operations: " << ms << " ms" << std::endl;
    std::cout << "Throughput: " << (1000 / ms) << " MB/s" << std::endl;
    delete[] buffer;
    return 0;
}
