#include "normalize_reuse_fixture.h"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>

static std::uint64_t bits(double value) {
    std::uint64_t raw;
    static_assert(sizeof(raw) == sizeof(value));
    std::memcpy(&raw, &value, sizeof(raw));
    return raw;
}

int main(int argc, char** argv) {
    const bool parity = argc == 3 && std::string_view(argv[1]) == "--parity";
    if (argc != 2 && !parity) return 2;
    const std::filesystem::path root(argv[parity ? 2 : 1]);
    reuse_fixture::make_inputs(root, false);
    const auto rows = reuse_fixture::mixed_rows();
    constexpr int iterations = 256;
    using Clock = std::chrono::steady_clock;
    if (!parity) std::cout << std::setprecision(17) << "format\tcase\trows\tseconds\tstatus\n";
    for (int format : {2, 3}) {
        const auto index = reuse_fixture::prepare(root, format);
        if (parity) {
            mss::Reference ref(index, (root / "tiny.fa").string());
            for (const auto& item : rows) {
                const auto r = ref.harmonise(item.second);
                std::cout << format << ' ' << item.first << ' ' << std::quoted(r.status) << ' '
                          << std::quoted(r.chr) << ' ' << r.pos << ' ' << std::quoted(r.ref) << ' '
                          << std::quoted(r.alt) << ' ' << std::quoted(r.id) << ' '
                          << bits(r.effect) << ' ' << bits(r.eaf) << ' ' << r.candidate_count
                          << ' ' << r.swapped << ' ' << r.complemented << '\n';
            }
            continue;
        }
        const auto start_open = Clock::now();
        mss::Reference ref(index, (root / "tiny.fa").string());
        const std::chrono::duration<double> open_time = Clock::now() - start_open;
        std::cout << format << "\topen\t0\t" << open_time.count() << "\t-\n";
        for (const auto& item : rows) {
            std::string status;
            const auto start = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                auto result = ref.harmonise(item.second);
                status = result.status;
            }
            const std::chrono::duration<double> elapsed = Clock::now() - start;
            std::cout << format << '\t' << item.first << '\t' << iterations << '\t'
                      << elapsed.count() << '\t' << status << '\n';
        }
    }
    return std::cout.flush() ? 0 : 1;
}
