// Build against either header using -I <worktree>/src. Compare parity stdout
// byte-for-byte between builds. The pressure mode exits 1 on allocation churn.
#ifdef NDEBUG
#error "lookup_parity requires assertions enabled"
#endif

#include <preservation_core.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <new>
#include <random>

static std::size_t allocation_calls = 0;
static bool counting = false;
void* operator new(std::size_t n) {
  if (counting) ++allocation_calls;
  if (void* p = std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

namespace {
using namespace mss;
namespace fs = std::filesystem;

void write(const fs::path& path, const std::string& data) {
  std::ofstream stream(path, std::ios::binary);
  stream << data;
  if (!stream) throw std::runtime_error("fixture write failed");
}

const std::string duplicate = "I:" + std::string(220, 'd');
const std::string long_prefix = "I:" + std::string(70000, 'x') + u8"雪";
constexpr std::size_t occurrences = 400;

std::string numbered_key(int n) {
  return "I:prefix-shared-by-all-short-keys-" + std::to_string(100000 + n * 2);
}

void lookup_fixture(const fs::path& dir, const fs::path& fasta, int format) {
  fs::create_directory(dir);
  Fasta(fasta.string()).save((dir / "fasta.fai").string());
  write(dir / "core-version", "mss-preservation-" + std::to_string(format) + "\n");
  const std::string record = format == 3 ? "A\t1\t1\tA\tC\tr\t1\n"
                                        : "1\t1\tA\tC\tr\t1\tA\tC\n";
  std::string records;
  std::vector<std::string> rows;
  for (std::size_t i = 0; i < occurrences; ++i) {
    rows.push_back(duplicate + "\t" + std::to_string(records.size()));
    records += record;
  }
  for (int i = 0; i < 3000; ++i) rows.push_back(numbered_key(i) + "\t0");
  for (const auto& key : {std::string("I:a"), long_prefix + "0", long_prefix + "2"})
    rows.push_back(key + "\t0");
  std::sort(rows.begin(), rows.end());
  write(dir / "records", records);
  std::ofstream offsets(dir / "offsets", std::ios::binary);
  std::ofstream keys(dir / "keys", std::ios::binary);
  std::uint64_t pos = 0;
  for (const auto& row : rows) {
    write64(offsets, pos);
    keys << row << '\n';
    pos += row.size() + 1;
  }
  if (!offsets || !keys) throw std::runtime_error("fixture index write failed");
}

void fields(const Candidate& c) {
  for (const auto* s : {&c.chr,&c.ref,&c.alt,&c.id,&c.original_ref,&c.alleles})
    std::cout << ' ' << std::quoted(*s);
  std::cout << ' ' << c.pos << ' ' << c.original_pos << " ids";
  for (const auto& id : c.identifiers) std::cout << ' ' << std::quoted(id);
  std::cout << " positions";
  for (auto p : c.original_positions) std::cout << ' ' << p;
  std::cout << '\n';
}

void query(Reference& ref, const std::string& key, std::size_t cap, std::size_t work,
           const char* expected_error = nullptr) {
  std::cout << "cap " << cap << " initial " << work;
  try {
    const auto result = ref.lookup(key, cap, &work);
    assert(!expected_error);
    std::cout << " found " << result.size() << '\n';
    for (const auto& c : result) fields(c);
  } catch (const std::exception& e) {
    if (expected_error) {
      assert(std::strcmp(e.what(), expected_error) == 0);
      assert(work == 0);
    }
    std::cout << " error " << std::quoted(e.what()) << '\n';
  }
  std::cout << "left " << work << '\n';
}

void budgets(Reference& ref) {
  for (const auto& key : {std::string("I:a"), std::string("I:A"), long_prefix + "0",
                          long_prefix + "1", long_prefix + "2", std::string("I:zz")})
    query(ref, key, 1, 1);
  for (std::size_t cap : {occurrences, occurrences - 1, occurrences + 1})
    for (std::size_t work : {std::size_t(0), occurrences - 1, occurrences})
      query(ref, duplicate, cap, work);
  query(ref, duplicate, occurrences, occurrences);
  std::vector<int> indices(3000);
  for (int i = 0; i < 3000; ++i) indices[i] = i;
  std::mt19937 generator(257);
  std::shuffle(indices.begin(), indices.end(), generator);
  for (auto i : indices) {
    assert(ref.lookup(numbered_key(i)).size() == 1);
    assert(ref.lookup(numbered_key(i) + "z").empty());
  }
  query(ref, long_prefix + "0", 1, 1);
}

std::uint64_t bits(double value) {
  std::uint64_t raw;
  static_assert(sizeof(raw) == sizeof(value), "Double must be 64 bits");
  std::memcpy(&raw, &value, sizeof(raw));
  return raw;
}

void science(const fs::path& root, int format) {
  const auto fasta = root / "science.fa";
  const auto sites = root / "science.tsv";
  write(fasta, ">1\nCAAAAACGTACGT\n>2\nACGT\n>3\nCAAAAAAAAG\n");
  write(sites, "chr\tpos\tref\talt\tid\taliases\n"
    "1\t8\tG\tA,T\tr\told;shared\n1\t8\tG\tT\tbridge\t.\n"
    "1\t10\tA\tC\tbridge\tshared\n2\t1\tA\tG\tr\t.\n"
    "1\t4\tAA\tA\tdeletion\t.\n3\t5\tAA\tA\trepeat\told-repeat\n"
    "2\t2\tC\tA,G,T,CA,CAG,<DEL>\tdense\told-dense\n");
  const auto dir = root / ("science" + std::to_string(format));
  prepare(fasta.string(), sites.string(), dir.string(), 3, format);
  Reference ref(dir.string(), fasta.string());
  std::vector<Input> cases = {
    {"1",8,"old","G","T",0.125,0.2,"EAF"}, {"1",8,"old","T","G",-0.0,0.0,"EAF"},
    {"1",8,"old","A","T",0.25,0.5,"EAF"}, {"1",8,"old","C","A",0.75,1.0,"MAF"},
    {"",0,"shared","G","T",0.75,1.0,"EAF"}, {"1",10,"shared","A","C",-2,0.3,"EAF"},
    {"1",4,"deletion","AA","A",1.0,0.3,"EAF"}, {"3",8,"old-repeat","AA","A",0.5,0.9,"EAF"},
    {"",0,"old-repeat","AA","A",0.5,0.9,"EAF"}, {"2",2,"old-dense","C","A",0.5,0.9,"EAF"}
  };
  for (const auto& strand : {"+", "-", "?"}) {
    Input row = cases[1];
    row.source_contract = true; row.source_ref = "G"; row.source_alt = "T"; row.source_strand = strand;
    cases.push_back(row);
  }
  for (auto cap : {100, 1, 2, 4, 100}) for (const auto& row : cases) {
    try {
      const auto r = ref.harmonise(row, cap);
      std::cout << "result";
      for (const auto* s : {&r.status,&r.chr,&r.ref,&r.alt,&r.id}) std::cout << ' ' << std::quoted(*s);
      std::cout << ' ' << r.pos << ' ' << bits(r.effect) << ' ' << bits(r.eaf) << ' '
                << r.candidate_count << ' ' << r.swapped << ' ' << r.complemented << '\n';
    } catch (const std::exception& e) { std::cout << "science error " << std::quoted(e.what()) << '\n'; }
  }
}

void corruption(const fs::path& dir, const fs::path& fasta, int format) {
  fs::create_directory(dir);
  Fasta(fasta.string()).save((dir / "fasta.fai").string());
  write(dir / "core-version", "mss-preservation-" + std::to_string(format) + "\n");
  const std::string good = format == 3 ? "A\t1\t1\tA\tC\tr\t1\n" : "1\t1\tA\tC\tr\t1\tA\tC\n";
  write(dir / "records", good + "broken record\n");
  const std::string first = "I:r\t0\n";
  const auto second = "I:r\t" + std::to_string(good.size()) + "\n";
  write(dir / "keys", first + second + "I:s\t0\n");
  { std::ofstream o(dir / "offsets", std::ios::binary);
    write64(o,0); write64(o,first.size()); write64(o,first.size() + second.size()); }
  Reference ref(dir.string(), fasta.string());
  for (auto cap : {1U,2U,1U}) for (auto work : {0U,1U,2U}) query(ref,"I:r",cap,work);
  query(ref,"I:s",1,1);
  if (format == 3) {
    unsigned i = 0;
    for (const auto& bad : {"I:s\t01\n","I:s\t9999999\n","I:s\t1\n","I:s\t0","I:s\t0\tx\n"}) {
      const auto malformed = dir / std::to_string(i++);
      fs::create_directory(malformed);
      Fasta(fasta.string()).save((malformed / "fasta.fai").string());
      write(malformed / "core-version", "mss-preservation-3\n");
      write(malformed / "records", good); write(malformed / "keys",bad);
      { std::ofstream o(malformed / "offsets",std::ios::binary); write64(o,0); }
      Reference broken(malformed.string(),fasta.string());
      query(broken,"I:s",1,1); query(broken,"I:s",1,1);
    }

    const auto tail = dir / "record-tail";
    fs::create_directory(tail);
    Fasta(fasta.string()).save((tail / "fasta.fai").string());
    write(tail / "core-version", "mss-preservation-3\n");
    // Seven A fields, a valid eighth field, and an invalid UTF-8 ninth field.
    write(tail / "records", "A\t1\t1\tA\tC\tr\t1\textra\t\xff\n");
    write(tail / "keys", "I:r\t0\n");
    { std::ofstream o(tail / "offsets", std::ios::binary); write64(o, 0); }
    Reference malformed_tail(tail.string(), fasta.string());
    const char* error = "Typed reference field is not valid UTF-8/framing";
    std::cout << "record tail invalid UTF-8 in field 9\n";
    for (int attempt = 0; attempt < 2; ++attempt)
      query(malformed_tail, "I:r", 1, 1, error);
  }
}

int pressure(Reference& ref) {
  // A repeated absent key beyond the last rank has no candidate or boundary
  // record. After the first validated search no row decode/allocation is needed.
  std::string key = "I:zzzz-a-negative-key-just-beyond-the-last-record";
  assert(ref.lookup(key).empty());
  allocation_calls = 0;
  counting = true;
  for (int i=0;i<1024;++i) assert(ref.lookup(key).empty());
  counting = false;
  std::cout << "warm negative allocations " << allocation_calls << " for 1024 queries\n";
#ifdef MSS_LOOKUP_COUNTERS
  const auto warm=ref.lookup_counters();
  assert(warm.rank_hits==1024);
  std::vector<int> order(3000);
  for(int i=0;i<3000;++i)order[i]=i;
  std::mt19937 generator(257);
  std::shuffle(order.begin(),order.end(),generator);
  // The numbered keys occupy consecutive ranks after I:a and the duplicate
  // run. Replay baseline binary-search ordinal comparisons independently.
  std::uint64_t ordinary_probes=0;
  for(auto i:order) {
    std::size_t lo=0,hi=3000+occurrences+3;
    const auto target=1+occurrences+static_cast<std::size_t>(i);
    while(lo<hi) {
      ++ordinary_probes;
      const auto mid=lo+(hi-lo)/2;
      if(mid<target)lo=mid+1;else hi=mid;
    }
    assert(ref.lookup(numbered_key(i)).size()==1);
  }
  const auto distinct=ref.lookup_counters();
  const auto loaded=distinct.probe_entries-warm.probe_entries;
  const auto reused=distinct.pivot_hits-warm.pivot_hits;
  assert(loaded+reused==ordinary_probes);
  assert(loaded*2<ordinary_probes);
  assert(distinct.pivot_slots<=8191&&distinct.rank_slots<=4096&&distinct.rank_bytes<=1024*1024);
  std::cout<<"shuffled unique "<<order.size()<<" ordinary probes "<<ordinary_probes
           <<" loaded "<<loaded<<" fence hits "<<reused<<" key page requests "
           <<distinct.key_requests-warm.key_requests<<" reads "<<distinct.key_reads-warm.key_reads
           <<" offset page requests "<<distinct.offset_requests-warm.offset_requests
           <<" reads "<<distinct.offset_reads-warm.offset_reads<<" rank bytes "<<distinct.rank_bytes<<'\n';
#endif
  return allocation_calls > 64;
}
} // namespace

int main(int argc, char** argv) {
  if (argc != 3 || (std::string(argv[1]) != "parity" && std::string(argv[1]) != "pressure")) return 2;
  const fs::path root = argv[2];
  if (!fs::create_directory(root)) throw std::runtime_error("require a fresh fixture directory");
  const auto fasta = root / "tiny.fa";
  write(fasta, ">1\nACGT\n");
  int failed = 0;
  for (int format : {2,3}) {
    const auto dir = root / ("lookup" + std::to_string(format));
    lookup_fixture(dir, fasta, format);
    mss::Reference ref(dir.string(), fasta.string());
    std::cout << "format " << format << '\n';
    if (std::string(argv[1]) == "pressure") {
      failed |= pressure(ref);
    } else {
      budgets(ref); science(root,format);
      corruption(root / ("corrupt" + std::to_string(format)), fasta, format);
    }
  }
  if (std::string(argv[1]) == "parity") {
    const auto path = root / "blocks";
    write(path,"a\nbcdefghij\nkl\n");
    mss::BlockFile file(path.string(),3,1);
    for (int i=0;i<8;++i) {
      assert(file.strict_line(2)=="bcdefghij"); assert(file.strict_line(12)=="kl");
      assert(file.bytes(1,5)=="\nbcde"); assert(file.resident_bytes()<=3);
#ifdef MSS_LOOKUP_COUNTERS
      file.visit_line(2,true,[](std::string_view row){assert(row=="bcdefghij");});
      file.visit_line(12,true,[](std::string_view row){assert(row=="kl");});
      assert(file.little_endian_offset(2)==UINT64_C(0x6968676665646362));
#endif
    }
  }
  return failed;
}
