#include <zlib.h>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Streaming extractor for 33 fixed myeloid-state genes from GSE136831.
// Usage: extract_gse136831_myeloid_modules matrix.mtx.gz output_uint32.bin

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " matrix.mtx.gz output.bin\n";
        return 2;
    }
    constexpr int n_rows = 45947;
    constexpr int n_cols = 312928;
    constexpr int n_targets = 33;
    constexpr std::uint64_t expected_entries = 692789348ULL;
    const int target_rows[n_targets] = {
        11858, 4986, 373, 13895, 12781, 26020, 12449, 12448, 9640,
        7892, 10219, 1589, 529, 10200, 2530, 1804,
        6777, 3284, 1656, 30184, 2823, 642, 14463, 6716,
        4592, 2319, 2258, 8231, 6274, 5830, 1967, 240, 6271
    };
    std::vector<int> row_to_target(n_rows + 1, -1);
    for (int i = 0; i < n_targets; ++i) row_to_target[target_rows[i]] = i;
    std::vector<std::uint32_t> output(static_cast<std::size_t>(n_targets) * n_cols, 0);

    gzFile in = gzopen(argv[1], "rb");
    if (!in) { std::cerr << "cannot open input\n"; return 3; }
    gzbuffer(in, 8 * 1024 * 1024);
    constexpr std::size_t buf_size = 16 * 1024 * 1024;
    std::vector<unsigned char> buffer(buf_size);
    int header_newlines = 0, field = 0;
    std::uint64_t number = 0, triple[3] = {0, 0, 0}, parsed = 0, retained = 0;
    bool in_number = false;
    while (true) {
        int got = gzread(in, buffer.data(), static_cast<unsigned int>(buffer.size()));
        if (got < 0) { int e = 0; std::cerr << gzerror(in, &e) << "\n"; gzclose(in); return 4; }
        if (got == 0) break;
        for (int i = 0; i < got; ++i) {
            unsigned char ch = buffer[i];
            if (header_newlines < 2) { if (ch == '\n') ++header_newlines; continue; }
            if (ch >= '0' && ch <= '9') { number = number * 10 + (ch - '0'); in_number = true; }
            else if (in_number) {
                triple[field++] = number; number = 0; in_number = false;
                if (field == 3) {
                    int row = static_cast<int>(triple[0]);
                    int col = static_cast<int>(triple[1]);
                    if (row >= 1 && row <= n_rows && col >= 1 && col <= n_cols) {
                        int target = row_to_target[row];
                        if (target >= 0) {
                            output[static_cast<std::size_t>(target) * n_cols + col - 1] += static_cast<std::uint32_t>(triple[2]);
                            ++retained;
                        }
                    }
                    field = 0; ++parsed;
                    if (parsed % 50000000ULL == 0) std::cerr << "parsed=" << parsed << "/" << expected_entries << " retained=" << retained << "\n";
                }
            }
        }
    }
    gzclose(in);
    if (in_number) {
        triple[field++] = number;
        if (field == 3) {
            int row = static_cast<int>(triple[0]); int col = static_cast<int>(triple[1]);
            if (row >= 1 && row <= n_rows && col >= 1 && col <= n_cols && row_to_target[row] >= 0) {
                output[static_cast<std::size_t>(row_to_target[row]) * n_cols + col - 1] += static_cast<std::uint32_t>(triple[2]);
                ++retained;
            }
            ++parsed;
        }
    }
    if (parsed != expected_entries) { std::cerr << "entry mismatch " << parsed << "\n"; return 5; }
    std::ofstream out(argv[2], std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(output.data()), static_cast<std::streamsize>(output.size() * sizeof(std::uint32_t)));
    out.close();
    if (!out) return 6;
    std::cerr << "complete parsed=" << parsed << " retained=" << retained << " bytes=" << output.size() * sizeof(std::uint32_t) << "\n";
    return 0;
}
