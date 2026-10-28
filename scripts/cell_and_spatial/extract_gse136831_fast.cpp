#include <zlib.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Fast streaming extractor for the column-sorted GSE136831 Matrix Market file.
// Usage: extract_gse136831_fast matrix.mtx.gz output_uint32.bin

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " matrix.mtx.gz output.bin\n";
        return 2;
    }

    constexpr int n_rows = 45947;
    constexpr int n_cols = 312928;
    constexpr int n_targets = 24;
    constexpr std::uint64_t expected_entries = 692789348ULL;
    const int target_rows[n_targets] = {
        8484, 1002, 8554, 9889, 11493, 5011, 11879, 5645,
        12062, 4521, 4710, 435, 7892, 10219, 6881, 4592,
        438, 5640, 3982, 2063, 8160, 3469, 8291, 15103
    };

    std::vector<int> row_to_target(n_rows + 1, -1);
    for (int i = 0; i < n_targets; ++i) row_to_target[target_rows[i]] = i;
    std::vector<std::uint32_t> output(static_cast<std::size_t>(n_targets) * n_cols, 0);

    gzFile in = gzopen(argv[1], "rb");
    if (!in) {
        std::cerr << "cannot open input: " << argv[1] << "\n";
        return 3;
    }
    gzbuffer(in, 8 * 1024 * 1024);

    constexpr std::size_t buf_size = 16 * 1024 * 1024;
    std::vector<unsigned char> buffer(buf_size);
    int header_newlines = 0;
    std::uint64_t number = 0;
    bool in_number = false;
    int field = 0;
    std::uint64_t triple[3] = {0, 0, 0};
    std::uint64_t parsed = 0;
    std::uint64_t retained = 0;

    while (true) {
        int got = gzread(in, buffer.data(), static_cast<unsigned int>(buffer.size()));
        if (got < 0) {
            int errnum = 0;
            std::cerr << "gzread failed: " << gzerror(in, &errnum) << "\n";
            gzclose(in);
            return 4;
        }
        if (got == 0) break;

        for (int i = 0; i < got; ++i) {
            unsigned char ch = buffer[i];
            if (header_newlines < 2) {
                if (ch == '\n') ++header_newlines;
                continue;
            }
            if (ch >= '0' && ch <= '9') {
                number = number * 10 + (ch - '0');
                in_number = true;
            } else if (in_number) {
                triple[field++] = number;
                number = 0;
                in_number = false;
                if (field == 3) {
                    const int row = static_cast<int>(triple[0]);
                    const int col = static_cast<int>(triple[1]);
                    const std::uint32_t count = static_cast<std::uint32_t>(triple[2]);
                    if (row >= 1 && row <= n_rows) {
                        const int target = row_to_target[row];
                        if (target >= 0 && col >= 1 && col <= n_cols) {
                            output[static_cast<std::size_t>(target) * n_cols + (col - 1)] += count;
                            ++retained;
                        }
                    }
                    field = 0;
                    ++parsed;
                    if (parsed % 50000000ULL == 0) {
                        std::cerr << "parsed=" << parsed << "/" << expected_entries
                                  << " retained=" << retained << "\n";
                    }
                }
            }
        }
    }
    gzclose(in);

    // Process a final number if the file does not end in whitespace.
    if (in_number) {
        triple[field++] = number;
        if (field == 3) {
            const int row = static_cast<int>(triple[0]);
            const int col = static_cast<int>(triple[1]);
            if (row >= 1 && row <= n_rows && row_to_target[row] >= 0 && col >= 1 && col <= n_cols) {
                output[static_cast<std::size_t>(row_to_target[row]) * n_cols + (col - 1)] +=
                    static_cast<std::uint32_t>(triple[2]);
                ++retained;
            }
            ++parsed;
        }
    }

    if (parsed != expected_entries) {
        std::cerr << "entry count mismatch: parsed=" << parsed
                  << " expected=" << expected_entries << "\n";
        return 5;
    }

    std::ofstream out(argv[2], std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "cannot open output: " << argv[2] << "\n";
        return 6;
    }
    out.write(reinterpret_cast<const char*>(output.data()),
              static_cast<std::streamsize>(output.size() * sizeof(std::uint32_t)));
    out.close();
    if (!out) {
        std::cerr << "failed writing output\n";
        return 7;
    }
    std::cerr << "complete parsed=" << parsed << " retained=" << retained
              << " bytes=" << output.size() * sizeof(std::uint32_t) << "\n";
    return 0;
}
