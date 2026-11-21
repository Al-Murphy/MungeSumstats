#ifndef MSS_TEST_NORMALIZE_REUSE_FIXTURE_H
#define MSS_TEST_NORMALIZE_REUSE_FIXTURE_H
#include <preservation_core.h>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace reuse_fixture {
namespace fs = std::filesystem;
using Row = std::pair<std::string, mss::Input>;

inline void write_file(const fs::path& path, const std::string& data) {
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << data;
    out.close();
}

inline void make_inputs(const fs::path& root, bool diagnostic) {
    if (!fs::create_directory(root))
        throw std::runtime_error("Fixture directory must be fresh");
    // The timing inputs and case order come from the November 18 mixed driver.
    write_file(root / "tiny.fa", ">1\nCAAAAACGTACGT\n>2\nACGT\n>3\nC" +
               std::string(256, 'A') + "G\n>4\nACGTACGT\n");
    std::string sites = "chr\tpos\tref\talt\tid\taliases\n"
        "1\t8\tG\tA,T\trs8\told8;merged8\n"
        "1\t4\tAA\tA\trsDel\t.\n"
        "1\t10\tA\tC\trsDup\tshared\n"
        "1\t12\tG\tA\trsDup\tshared\n"
        "2\t2\tC\tG\tpalindrome\t.\n"
        "2\t1\tA\tC\tduplicate\t.\n"
        "1\t10\tA\tC\tduplicate\t.\n"
        "2\t3\tG\tA\tsplitThird\t.\n"
        "2\t3\tG\tC\tsplitThird2\t.\n"
        "3\t128\tAA\tA\trsRepeat\t.\n"
        "4\t1\tA\tC,G,T,AC,AG,AT,AAC,AAG,AAT,AAAC,AAAG,AAAT,<DEL>\tdense\toldDense\n";
    if (diagnostic)
        sites += "1\t7\tC\tT\tunique\t.\n"
                 "1\t9\tT\tC,T,T,T,T,T,T,T,T,T,T,T\tbudget\t.\n"
                 "3\t30\tAA\tA\tbridge\taliasBridge\n"
                 "3\t30\tAA\t<DEL>\tsymbolicOnly\t.\n"
                 "4\t3\tG\tA,<DEL>\tsymbolic\t.\n";
    write_file(root / "sites.tsv", sites);
}

inline std::vector<Row> mixed_rows() {
    std::vector<Row> rows;
    auto add = [&](const char* name, const char* chr, int pos, const char* id,
                   const char* other, const char* effect) {
        rows.emplace_back(name, mss::Input{chr, pos, id, other, effect, .125, .25, "EAF"});
    };
    add("simple_snp", "2", 1, "", "A", "C");
    add("simple_swap", "2", 1, "", "C", "A");
    add("synonym_multiallelic", "1", 8, "old8", "G", "T");
    add("complement_multiallelic", "1", 8, "rs8", "C", "A");
    add("alt_alt_multiallelic", "1", 8, "rs8", "A", "T");
    add("palindromic", "2", 2, "palindrome", "C", "G");
    add("alias_multiple_events", "1", 10, "shared", "A", "C");
    add("missing_coordinate", "", 0, "duplicate", "A", "C");
    add("generic_indel", "1", 4, "", "AA", "A");
    add("repeat_indel", "3", 200, "rsRepeat", "AA", "A");
    auto declared = rows.back();
    declared.first = "declared_repeat_indel";
    declared.second.source_contract = true;
    declared.second.source_ref = "AA";
    declared.second.source_alt = "A";
    declared.second.source_strand = "+";
    rows.push_back(declared);
    add("dense_multiallelic", "4", 1, "oldDense", "A", "C");
    return rows;
}

inline std::string prepare(const fs::path& root, int format) {
    const auto path = (root / ("format" + std::to_string(format))).string();
    mss::prepare((root / "tiny.fa").string(), (root / "sites.tsv").string(), path, 100, format);
    return path;
}
} // namespace reuse_fixture
#endif
