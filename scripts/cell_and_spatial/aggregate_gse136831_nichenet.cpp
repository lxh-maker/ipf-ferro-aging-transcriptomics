#include <zlib.h>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Full-transcriptome pseudobulk aggregator for GSE136831.
// Categories: 0 ATII, 1 Club/Goblet secretory, 2 Myeloid; IPF donors only.
// Usage: aggregate_gse136831_nichenet matrix.gz donor.bin bits.bin counts.bin detected.bin

template <typename T>
bool read_binary(const char* path, std::vector<T>& values) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    in.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    return static_cast<bool>(in);
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: " << argv[0]
                  << " matrix.gz donor.int16.bin bits.uint8.bin counts.bin detected.bin\n";
        return 2;
    }
    constexpr int n_rows = 45947;
    constexpr int n_cols = 312928;
    constexpr int n_donors = 32;
    constexpr int n_categories = 3;
    constexpr std::uint64_t expected_entries = 692789348ULL;

    std::vector<std::int16_t> donor(n_cols);
    std::vector<std::uint8_t> bits(n_cols);
    if (!read_binary(argv[2], donor) || !read_binary(argv[3], bits)) {
        std::cerr << "failed to read cell map files\n";
        return 3;
    }

    const std::size_t block = static_cast<std::size_t>(n_rows) * n_donors;
    std::vector<std::uint64_t> counts(static_cast<std::size_t>(n_categories) * block, 0);
    std::vector<std::uint32_t> detected(static_cast<std::size_t>(n_categories) * block, 0);

    gzFile in = gzopen(argv[1], "rb");
    if (!in) {
        std::cerr << "cannot open matrix\n";
        return 4;
    }
    gzbuffer(in, 8 * 1024 * 1024);
    std::vector<unsigned char> buffer(16 * 1024 * 1024);
    int header_newlines = 0;
    std::uint64_t number = 0, triple[3] = {0, 0, 0}, parsed = 0, retained = 0;
    bool in_number = false;
    int field = 0;

    auto process_triple = [&]() {
        const int row = static_cast<int>(triple[0]) - 1;
        const int col = static_cast<int>(triple[1]) - 1;
        const std::uint64_t value = triple[2];
        if (row < 0 || row >= n_rows || col < 0 || col >= n_cols) return;
        const int d = donor[col];
        if (d < 0 || d >= n_donors) return;
        const std::uint8_t cell_bits = bits[col];
        for (int cat = 0; cat < n_categories; ++cat) {
            if ((cell_bits & (1u << cat)) == 0) continue;
            const std::size_t idx = static_cast<std::size_t>(cat) * block
                                  + static_cast<std::size_t>(row) * n_donors + d;
            counts[idx] += value;
            detected[idx] += 1;
            ++retained;
        }
    };

    while (true) {
        int got = gzread(in, buffer.data(), static_cast<unsigned int>(buffer.size()));
        if (got < 0) {
            int errnum = 0;
            std::cerr << "gzread failed: " << gzerror(in, &errnum) << "\n";
            gzclose(in);
            return 5;
        }
        if (got == 0) break;
        for (int i = 0; i < got; ++i) {
            const unsigned char ch = buffer[i];
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
                    process_triple();
                    field = 0;
                    ++parsed;
                    if (parsed % 100000000ULL == 0)
                        std::cerr << "parsed=" << parsed << "/" << expected_entries
                                  << " retained_category_entries=" << retained << "\n";
                }
            }
        }
    }
    if (in_number) {
        triple[field++] = number;
        if (field == 3) {
            process_triple();
            ++parsed;
        }
    }
    gzclose(in);
    if (parsed != expected_entries) {
        std::cerr << "entry count mismatch: " << parsed << " != " << expected_entries << "\n";
        return 6;
    }

    std::ofstream count_out(argv[4], std::ios::binary | std::ios::trunc);
    count_out.write(reinterpret_cast<const char*>(counts.data()),
                    static_cast<std::streamsize>(counts.size() * sizeof(std::uint64_t)));
    count_out.close();
    std::ofstream detect_out(argv[5], std::ios::binary | std::ios::trunc);
    detect_out.write(reinterpret_cast<const char*>(detected.data()),
                     static_cast<std::streamsize>(detected.size() * sizeof(std::uint32_t)));
    detect_out.close();
    if (!count_out || !detect_out) {
        std::cerr << "failed to write output\n";
        return 7;
    }
    std::cerr << "complete parsed=" << parsed << " retained_category_entries=" << retained
              << " counts_bytes=" << counts.size() * sizeof(std::uint64_t)
              << " detected_bytes=" << detected.size() * sizeof(std::uint32_t) << "\n";
    return 0;
}
