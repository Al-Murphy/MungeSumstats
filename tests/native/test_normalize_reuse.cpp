#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>

// The ephemeral header has one hook inside normalize itself. The phase markers
// only label entries. An attempted comparison that reuses a value records zero.
namespace reuse_diag {
enum Phase { seed, hit, proof, other };
Phase phase = other;
std::array<unsigned, 4> calls{};
bool measuring = false, selected = false, bad_operand = false;
void normal(const std::string& chr, std::int64_t pos, const std::string& a, const std::string& b) {
    if (!measuring) return;
    const auto ordinal = calls.at(phase)++;
    if (!selected) return;
    constexpr char refs[] = "TCAG", alts[] = "CTGA";
    const bool extra = phase == other;
    if (chr != "1" || pos != 7 || ordinal >= (extra ? 1U : 4U) ||
        a != std::string(1, extra ? 'C' : refs[ordinal]) ||
        b != std::string(1, extra ? 'T' : alts[ordinal])) bad_operand = true;
}
void remaining(const std::size_t* p) {
    if (p) std::cout << *p; else std::cout << '-';
}
void lookup(const std::string& key, std::size_t cap, const std::size_t* p) {
    std::cout << "lookup " << std::quoted(key) << ' ' << cap << ' ';
    remaining(p);
    std::cout << '\n';
}
void occurrence(std::uint64_t rank, std::size_t earlier) {
    std::cout << "occurrence " << rank << ' ' << earlier << '\n';
}
void decode(std::uint64_t offset, std::string_view row) {
    std::cout << "decode " << offset << ' ' << std::quoted(std::string(row)) << '\n';
}
void debit(const std::size_t* p) {
    std::cout << "debit ";
    remaining(p);
    std::cout << '\n';
}
template<class Candidates> void returned(const Candidates& found) {
    std::cout << "candidates " << found.size() << '\n';
    for (const auto& c : found) {
        std::cout << std::quoted(c.chr) << ' ' << c.pos << ' ' << std::quoted(c.ref) << ' '
                  << std::quoted(c.alt) << ' ' << std::quoted(c.id) << ' ' << c.original_pos
                  << ' ' << std::quoted(c.original_ref) << ' ' << std::quoted(c.alleles);
        std::cout << " identifiers";
        for (const auto& id : c.identifiers) std::cout << ' ' << std::quoted(id);
        std::cout << " positions";
        for (const auto p : c.original_positions) std::cout << ' ' << p;
        std::cout << '\n';
    }
}
} // namespace reuse_diag
#include "normalize_reuse_fixture.h"

namespace {
using namespace reuse_fixture;
constexpr auto candidate_limit = "Reference candidate limit exceeded; no arbitrary candidate selected";
constexpr auto closure_limit = "Reference context closure limit exceeded; no arbitrary context selected";

void require(bool ok, const std::string& detail) {
    if (!ok) throw std::logic_error(detail);
}
std::uint64_t bits(double x) {
    std::uint64_t n;
    static_assert(sizeof(n) == sizeof(x));
    std::memcpy(&n, &x, sizeof(n));
    return n;
}
void result(const mss::Result& r) {
    std::cout << "result " << std::quoted(r.status) << ' ' << std::quoted(r.chr) << ' '
              << r.pos << ' ' << std::quoted(r.ref) << ' ' << std::quoted(r.alt) << ' '
              << std::quoted(r.id) << ' ' << bits(r.effect) << ' ' << bits(r.eaf) << ' '
              << r.candidate_count << ' ' << r.swapped << ' ' << r.complemented << '\n';
}
void error(const std::exception& e) {
    std::cout << "error " << typeid(e).name() << ' ' << std::quoted(e.what()) << '\n';
}
void row(mss::Reference& ref, const Row& item, std::size_t cap, int format,
         std::ostream& counts, const char* expected = nullptr, bool expect_error = false) {
    using namespace reuse_diag;
    std::cout << "row " << std::quoted(item.first) << ' ' << cap << '\n';
    phase = other;
    calls.fill(0);
    measuring = true;
    selected = item.first == "counted_swap";
    bad_operand = false;
    try {
        const auto r = ref.harmonise(item.second, cap);
        result(r);
        require(!expect_error && (!expected || r.status == expected), item.first + ": status");
        if (selected) {
            require(r.chr == "1" && r.pos == 7 && r.ref == "C" && r.alt == "T" &&
                    r.id == "unique" && r.swapped && !r.complemented && r.candidate_count == 1 &&
                    bits(r.effect) == bits(-0.0) && bits(r.eaf) == bits(.75), "counted swap fields");
        }
    } catch (const std::runtime_error& e) {
        error(e);
        require(expect_error && expected && e.what() == std::string(expected),
                item.first + ": exception");
    }
    measuring = false;
    require(!bad_operand, "unexpected normalization operands at counted SNP");
    if (selected) {
        counts << format;
        for (auto n : calls) counts << ' ' << n;
        counts << '\n';
    }
    // These branches have no eligible stored result. Compare the actual work
    // for them too, including every >64 KiB normalized allele that must bypass.
    if (item.second.source_contract || item.second.chr.empty() || item.second.pos <= 0 ||
        item.first == "oversized_allele") {
        std::cout << "fallback_calls";
        for (auto n : calls) std::cout << ' ' << n;
        std::cout << '\n';
    }
}

using Count = std::array<unsigned, 4>;
std::array<Count, 2> read_counts(const char* path) {
    require(std::filesystem::file_size(path) <= 128, "unexpected work count file size");
    std::ifstream in(path);
    std::array<Count, 2> counts{};
    for (int format : {2, 3}) {
        int observed = 0;
        require(bool(in >> observed) && observed == format, "missing work counts for format");
        for (auto& n : counts[format - 2]) require(bool(in >> n), "missing phase count");
    }
    in >> std::ws;
    require(in.eof(), "unexpected trailing work counts");
    return counts;
}

int compare_work(const char* control_path, const char* treatment_path) {
    const auto control = read_counts(control_path), treatment = read_counts(treatment_path);
    bool reduced = true;
    for (int format : {2, 3}) {
        const auto& a = control[format - 2];
        const auto& b = treatment[format - 2];
        require(a == Count{4, 4, 4, 1}, "control did not traverse all three stages and swapped EAF");
        const auto saved = int(a[1] + a[2]) - int(b[1] + b[2]);
        std::cout << "format " << format << " saved " << saved << " expected 8\n";
        reduced = reduced && saved == 8 && b == Count{4, 0, 0, 1};
    }
    return reduced ? 0 : 1;
}

// One readable record followed by a corrupt record under the same key, as in
// tests/standalone/lookup_parity.cpp at the qualified source revision.
void corruption(const fs::path& root, int format) {
    const auto dir = root / ("corrupt" + std::to_string(format));
    fs::create_directory(dir);
    const auto fasta = (root / "tiny.fa").string();
    mss::Fasta(fasta).save((dir / "fasta.fai").string());
    write_file(dir / "core-version", "mss-preservation-" + std::to_string(format) + "\n");
    const std::string good = format == 3 ? "A\t1\t7\tC\tT\tr\t7\n" : "1\t7\tC\tT\tr\t7\tC\tT\n";
    write_file(dir / "records", good + "broken record\n");
    const std::string first = "I:r\t0\n";
    write_file(dir / "keys", first + "I:r\t" + std::to_string(good.size()) + "\n");
    {
        std::ofstream out(dir / "offsets", std::ios::binary);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        mss::write64(out, 0);
        mss::write64(out, first.size());
        out.close();
    }
    mss::Reference ref(dir.string(), fasta);
    for (unsigned cap : {1U, 2U}) for (unsigned initial : {0U, 1U, 2U}) {
        std::size_t work = initial;
        std::cout << "corruption " << cap << ' ' << work << '\n';
        try {
            ref.lookup("I:r", cap, &work);
            throw std::logic_error("corrupt query unexpectedly succeeded");
        } catch (const std::runtime_error& e) {
            const char* wanted = initial == 0 ? closure_limit : cap == 1 ? candidate_limit :
                initial == 1 ? closure_limit : format == 3 ?
                "Typed reference record has an unknown tag or field count" : "Corrupt reference record";
            const auto left = initial == 2 && cap == 1 ? 1U : 0U;
            error(e);
            std::cout << "remaining " << work << '\n';
            require(typeid(e) == typeid(std::runtime_error) && std::string(e.what()) == wanted &&
                    work == left, "cap/debit/decode precedence");
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--compare-work")
        return compare_work(argv[2], argv[3]);
    if (argc != 2) return 2;
    const fs::path root(argv[1]);
    make_inputs(root, true);
    std::ofstream counts(root / "work.txt");
    counts.exceptions(std::ios::failbit | std::ios::badbit);
    for (int format : {2, 3}) {
        std::cout << "format " << format << '\n';
        auto path = reuse_fixture::prepare(root, format);
        mss::Reference ref(path, (root / "tiny.fa").string());
        std::cout << "case candidate_context\n";
        for (const auto& key : {"I:shared", "I:duplicate", "C:3:1", "S:3:30"})
            ref.lookup(key);
        std::cout << "case mixed_science\n";
        for (const auto& item : mixed_rows()) row(ref, item, 100, format, counts);

        auto check = [&](const char* name, mss::Input q, const char* expected = nullptr,
                         std::size_t cap = 100, bool fails = false) {
            row(ref, {name, std::move(q)}, cap, format, counts, expected, fails);
        };
        std::cout << "case orientation_frequency\n";
        check("counted_swap", {"1",7,"unique","T","C",0.0,.25,"EAF"}, "swap");
        check("signed_zero_forward", {"1",7,"unique","C","T",-0.0,-0.0,"EAF"}, "match");
        check("swapped_complement", {"1",7,"unique","A","G",1,1,"EAF"}, "swapped_complement");
        check("masked_frequency", {"1",8,"old8","T","G",.5,.25,"EAF"}, "swap");
        check("maf", {"1",7,"unique","C","T",1,.25,"MAF"}, "match");
        check("bad_eaf", {"1",7,"unique","C","T",1,1.5,"EAF"}, "match");
        std::cout << "case anchor_alias_closure\n";
        check("alias_bridge", {"",0,"aliasBridge","C","CA",1,.25,"EAF"}, "source_coordinates_missing");
        check("canonical_bridge", {"3",200,"aliasBridge","AA","A",1,.25,"EAF"});
        check("id_alternate_contigs", {"",1,"duplicate","A","C",1,.25,"EAF"});
        check("id_only_distinct_anchors", {"",0,"shared","G","A",1,.25,"EAF"});
        check("conflict", {"1",7,"rs8","C","T",1,.25,"EAF"}, "coordinate_id_conflict");
        check("absent", {"1",13,"","T","A",1,.25,"EAF"}, "unmapped");
        std::cout << "case original_symbolic_site\n";
        check("symbolic_context", {"4",3,"symbolic","G","A",1,.25,"EAF"}, "source_allele_context_required");
        std::cout << "case large_operands\n";
        check("oversized_allele", {"1",7,"unique","C",std::string(65537,'A'),1,.25,"EAF"}, "allele_mismatch");
        check("oversized_id", {"1",7,std::string(65537,'x'),"C","T",1,.25,"EAF"}, "match");
        std::cout << "case proof_limits\n";
        check("occurrence_cap", {"1",8,"old8","G","T",1,.25,"EAF"}, candidate_limit, 1, true);
        check("proof_work_exhausted", {"1",9,"budget","T","C",1,.25,"EAF"}, closure_limit, 1, true);
        check("proof_work_available", {"1",9,"budget","T","C",1,.25,"EAF"}, "match", 2);
        check("closure_overflow", {"1",7,"unique","C","T",1,.25,"EAF"},
              "Reference context closure limit is not representable",
              std::numeric_limits<std::size_t>::max(), true);
        std::cout << "case corrupt_equal_run\n";
        corruption(root, format);
    }
    counts.close();
    return std::cout.flush() ? 0 : 1;
}
