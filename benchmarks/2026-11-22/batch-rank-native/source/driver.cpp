#include <preservation_core.h>
#include <array>
#include <cstring>
#include <chrono>
#include <numeric>
#include <random>
#include <iostream>

static void write_file(const std::string& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << text;
    out.close();
}

#include "mixed_fixture.inc"

static std::uint64_t bits(double value) {
    static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == 8);
    std::uint64_t result;
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

// The pinned panel and fixture have tab/newline-free text. Empty strings remain empty.
// Decimal uint64 encodes the exact binary64 bit pattern, including signed zero and NaN.
static void emit(std::size_t ordinal, const std::string& label, const mss::Result& r) {
    std::cout << ordinal << '\t' << label << '\t' << r.status << '\t' << r.chr
              << '\t' << r.pos << '\t' << r.ref << '\t' << r.alt << '\t' << r.id
              << '\t' << bits(r.effect) << '\t' << bits(r.eaf) << '\t' << r.swapped
              << '\t' << r.complemented << '\t' << r.candidate_count << '\n';
}

constexpr std::size_t panel_batch = 100000;
constexpr unsigned shuffle_seed = 257;

struct Run {
    bool shuffled;
    std::minstd_rand generator{shuffle_seed};
    std::size_t rows = 0, batches = 0, last_batch = 0;
    std::uint64_t order_hash = UINT64_C(14695981039346656037);
    std::chrono::nanoseconds elapsed{};
#ifdef MSS_BATCH_RANK_BENCH
    std::size_t prefetch_attempts=0, prefetch_admitted=0, prefetch_peak=0;
    std::size_t gather_peak_estimate=0, gather_max_keys=0, gather_max_key_bytes=0;
#endif
};

// Fisher–Yates with a fixed seed and explicit bounded draw for library-independent order.
static void shuffle_order(std::vector<std::size_t>& order, std::minstd_rand& random) {
    constexpr std::uint64_t range = std::uint64_t(std::minstd_rand::max()) - std::minstd_rand::min() + 1;
    for (auto bound = order.size(); bound > 1; --bound) {
        const auto limit = range - range % bound;
        std::uint64_t draw;
        do { draw = random() - random.min(); } while (draw >= limit);
        std::swap(order[bound - 1], order[draw % bound]);
    }
}

static void process_batch(mss::Reference& reference, const std::vector<mss::Input>& inputs,
                          Run& run, const std::vector<std::string>& labels = {}) {
    const auto count = inputs.size();
    if (!count || (!labels.empty() && labels.size() != count))
        throw std::runtime_error("Invalid batch extent");
    std::vector<std::size_t> order(count), inverse(count, count);
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (run.shuffled) shuffle_order(order, run.generator);
    for (std::size_t step = 0; step < count; ++step) {
        const auto original = order[step];
        if (original >= count || inverse[original] != count)
            throw std::runtime_error("Invalid query permutation");
        inverse[original] = step;
        run.order_hash ^= run.rows + original + 1;
        run.order_hash *= UINT64_C(1099511628211);
    }

    std::vector<mss::Result> results(count);
    const auto start = std::chrono::steady_clock::now();
#ifdef MSS_BATCH_RANK_BENCH
    {
        // Borrowed caller keys are separate from the core's incremental arena.
        std::vector<std::string> keys;
        keys.reserve(2 * count);
        auto add = [&](std::string key) {
            if (key.size() <= 224) {
                run.gather_max_key_bytes = std::max(run.gather_max_key_bytes, key.size());
                keys.push_back(std::move(key));
            }
        };
        for (const auto& input : inputs) {
            if (!mss::missing(input.chr) && input.pos > 0)
                add("C:" + input.chr + ":" + std::to_string(input.pos));
            if (!mss::missing(input.id)) add("I:" + input.id);
        }
        std::size_t estimate = keys.capacity() * sizeof(std::string);
        for (const auto& key : keys) estimate += key.capacity() + 1;
        run.gather_peak_estimate = std::max(run.gather_peak_estimate, estimate);
        run.gather_max_keys = std::max(run.gather_max_keys, keys.size());
        ++run.prefetch_attempts;
        if (reference.prefetch_rank_keys(keys, 32U * 1024 * 1024)) ++run.prefetch_admitted;
        run.prefetch_peak = std::max(run.prefetch_peak, reference.prefetch_peak_bytes());
#endif
    for (std::size_t step = 0; step < count; ++step)
        results[step] = reference.harmonise(inputs[order[step]]);
#ifdef MSS_BATCH_RANK_BENCH
        reference.clear_prefetched_ranks();
        if (reference.prefetch_owned_bytes() || reference.prefetched_rank_count())
            throw std::runtime_error("Batch rank directory was not cleared");
    } // Includes key-vector/string destruction in the candidate clock.
#endif
    const auto end = std::chrono::steady_clock::now();
    run.elapsed += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    for (std::size_t i = 0; i < count; ++i)
        emit(run.rows + i + 1, labels.empty() ? inputs[i].id : labels[i], results[inverse[i]]);
    run.rows += count;
    ++run.batches;
    run.last_batch = count;
}

static void process_panel(mss::Reference& reference, Run& run, std::size_t capacity,
                          std::size_t expected_rows) {
    std::array<char, 65537> buffer;
    auto read_line = [&] {
        if (!std::cin.getline(buffer.data(), buffer.size()))
            throw std::runtime_error("Missing or oversized matched panel line");
        return std::string(buffer.data());
    };
    if (read_line() != "SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\tSE\tP")
        throw std::runtime_error("Unexpected matched input header");

    std::vector<mss::Input> inputs;
    inputs.reserve(capacity);
    while (std::cin.peek() != std::char_traits<char>::eof()) {
        if (run.rows + inputs.size() >= expected_rows)
            throw std::runtime_error("Unexpected matched panel tail");
        const auto fields = mss::split(read_line(), '\t');
        if (fields.size() != 9) throw std::runtime_error("Unexpected matched input columns");
        inputs.emplace_back(fields[1], std::stoll(fields[2]), fields[0], fields[3], fields[4],
                            std::stod(fields[5]), std::stod(fields[6]), "EAF");
        if (inputs.size() == capacity) {
            process_batch(reference, inputs, run);
            inputs.clear();
        }
    }
    if (std::cin.bad() || !std::cin.eof() || run.rows + inputs.size() != expected_rows)
        throw std::runtime_error("Incomplete matched input");
    if (!inputs.empty()) process_batch(reference, inputs, run);
}

int main(int argc, char** argv) {
    try {
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        const std::string mode = argc > 1 ? argv[1] : "";
        if (argc == 3 && mode == "prepare-fixture") {
            const std::string root = argv[2];
            if (!std::filesystem::create_directory(root)) throw std::runtime_error("Fixture directory exists");
            const auto fasta = root + "/tiny.fa", sites = root + "/sites.tsv";
            write_fixture(fasta, sites);
            for (int format : {2, 3})
                mss::prepare(fasta, sites, root + "/format" + std::to_string(format), 100, format);
            return 0;
        }
        if (argc != 5 || (mode != "panel" && mode != "fixture" && mode != "short" && mode != "empty"))
            throw std::runtime_error("Use panel|fixture|short|empty INDEX FASTA matched|shuffled");
        const std::string order = argv[4];
        if (order != "matched" && order != "shuffled") throw std::runtime_error("Unknown query order");
        Run run{order == "shuffled"};
        mss::Reference reference(argv[2], argv[3]);
        const std::size_t capacity = mode == "panel" ? panel_batch : 5;
        if (mode == "fixture") {
            std::vector<mss::Input> inputs;
            std::vector<std::string> labels;
            for (const auto& item : fixture_cases()) {
                labels.push_back(item.first);
                inputs.push_back(item.second);
                if (inputs.size() == capacity) {
                    process_batch(reference, inputs, run, labels);
                    inputs.clear();
                    labels.clear();
                }
            }
            if (!inputs.empty()) process_batch(reference, inputs, run, labels);
        } else {
            const std::size_t rows = mode == "panel" ? 1000000 : mode == "short" ? 7 : 0;
            process_panel(reference, run, capacity, rows);
        }
        std::cout.flush();
        if (!std::cout) throw std::runtime_error("Scientific result output failed");
        std::cerr << "{\"mode\":\"" << mode << "\",\"order\":\"" << order << "\",\"seed\":" << shuffle_seed
                  << ",\"capacity\":" << capacity << ",\"rows\":" << run.rows << ",\"batches\":" << run.batches
                  << ",\"last_batch\":" << run.last_batch << ",\"order_hash\":" << run.order_hash
                  << ",\"harmonise_ns\":" << run.elapsed.count()
#ifdef MSS_BATCH_RANK_BENCH
                  << ",\"batch_rank\":{\"attempts\":" << run.prefetch_attempts
                  << ",\"admitted\":" << run.prefetch_admitted
                  << ",\"arena_peak_bytes\":" << run.prefetch_peak
                  << ",\"gather_estimate_bytes\":" << run.gather_peak_estimate
                  << ",\"gather_max_keys\":" << run.gather_max_keys
                  << ",\"gather_max_key_bytes\":" << run.gather_max_key_bytes << "}"
#endif
                  << "}\n";
        return std::cerr ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
