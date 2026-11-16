#include "preservation_core.h"
#include <iostream>

namespace {
using U64 = std::uint64_t;
void add(U64& total, U64 value) {
  if (value > std::numeric_limits<U64>::max() - total)
    throw std::runtime_error("Count-only byte/row arithmetic overflow");
  total += value;
}
class CountBuffer : public std::streambuf {
public:
  U64 bytes = 0;
  bool text_newlines;
  explicit CountBuffer(bool text = false) : text_newlines(text) {}
protected:
  std::streamsize xsputn(const char* data, std::streamsize n) override {
    if (n < 0) throw std::runtime_error("Negative stream length");
    add(bytes, static_cast<U64>(n));
#ifdef _WIN32
    if (text_newlines) for (std::streamsize i = 0; i < n; ++i)
      if (data[i] == '\n') add(bytes, 1);
#endif
    return n;
  }
  int_type overflow(int_type character) override {
    if (traits_type::eq_int_type(character, traits_type::eof()))
      return traits_type::not_eof(character);
    const char value = traits_type::to_char_type(character);
    xsputn(&value, 1);
    return character;
  }
};
class CountStream : public std::ostream {
public:
  CountBuffer buffer;
  explicit CountStream(bool text = false) : std::ostream(nullptr), buffer(text) {
    rdbuf(&buffer);
    exceptions(std::ios::badbit | std::ios::failbit);
  }
};
struct Counts {
  U64 source_rows = 0, source_bytes = 0, record_rows = 0;
  U64 site_records = 0, allele_records = 0, invalid_positions = 0;
  U64 raw_alt_tokens = 0, rejected_alts = 0;
  U64 raw_key_rows = 0, raw_sort_bytes = 0, final_key_rows = 0, keys_bytes = 0;
  U64 max_source_line_bytes = 0, max_key_row_bytes = 0, max_record_bytes = 0;
  U64 max_raw_keys_per_record = 0;
  void keys(U64 offset, const std::vector<std::string>& keys) {
    const auto pointer = std::to_string(offset);
    std::set<std::string> unique;
    max_raw_keys_per_record = std::max(max_raw_keys_per_record, static_cast<U64>(keys.size()));
    for (const auto& key : keys) {
      U64 length = key.size(); add(length, pointer.size()); add(length, 2);
      add(raw_sort_bytes, length); add(raw_key_rows, 1);
      max_key_row_bytes = std::max(max_key_row_bytes, length);
      if (unique.insert(key).second) { add(keys_bytes, length); add(final_key_rows, 1); }
    }
  }
};

void count(const std::string& fasta, const std::string& sites, int version) {
  if (version != 2 && version != 3) throw std::runtime_error("Unsupported reference format");
  mss::Fasta reference(fasta);
  std::ifstream input(sites);
  if (!input) throw std::runtime_error("Cannot open source sites");
  Counts result;
  result.source_bytes = std::filesystem::file_size(sites);
  CountStream records, rejected(true);
  rejected << "row\tstatus\n";
  std::string line;
  std::getline(input, line);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line != "chr\tpos\tref\talt\tid\taliases")
    throw std::runtime_error("Reference header must have the six declared TSV fields");
  U64 row = 1;
  while (std::getline(input, line)) {
    add(row, 1); add(result.source_rows, 1);
    result.max_source_line_bytes = std::max(result.max_source_line_bytes, static_cast<U64>(line.size()));
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto fields = mss::split(line, '\t');
    if (fields.size() != 6) throw std::runtime_error("Reference row must have six fields");
    if (version == 3) for (const auto& field : fields)
      if (!mss::typed_text(field)) throw std::runtime_error("Invalid typed source framing/UTF8");
    const auto alleles = mss::split(mss::upper(fields[3]), ',');
    add(result.raw_alt_tokens, alleles.size());
    std::int64_t position;
    try {
      std::size_t used;
      position = std::stoll(fields[1], &used);
      if (used != fields[1].size()) throw std::invalid_argument("position");
    } catch (...) {
      rejected << row << "\tinvalid_position\n";
      add(result.invalid_positions, 1);
      continue;
    }
    const auto original_ref = mss::upper(fields[2]);
    const auto original_alts = mss::upper(fields[3]);
    const auto source_position = std::to_string(position);
    auto offset = records.buffer.bytes;
    if (version == 2)
      records << fields[0] << '\t' << position << '\t' << original_ref << "\t\t" << fields[4]
              << '\t' << position << '\t' << original_ref << '\t' << original_alts << '\n';
    else records << "S\t" << fields[0] << '\t' << position << '\t' << original_ref << '\t' << original_alts << '\n';
    result.keys(offset, {"S:" + fields[0] + ":" + source_position});
    add(result.record_rows, 1); add(result.site_records, 1);
    result.max_record_bytes = std::max(result.max_record_bytes, records.buffer.bytes - offset);
    for (const auto& alt : alleles) {
      const auto normalized = mss::normalize(reference, fields[0], position, fields[2], alt);
      if (normalized.status != "ok") {
        rejected << row << '\t' << normalized.status << '\n';
        add(result.rejected_alts, 1);
        continue;
      }
      offset = records.buffer.bytes;
      if (version == 2)
        records << fields[0] << '\t' << normalized.pos << '\t' << normalized.ref << '\t'
                << normalized.alt << '\t' << fields[4] << '\t' << position << '\t'
                << original_ref << '\t' << original_alts << '\n';
      else records << "A\t" << fields[0] << '\t' << normalized.pos << '\t' << normalized.ref
                   << '\t' << normalized.alt << '\t' << fields[4] << '\t' << position << '\n';
      std::vector<std::string> keys = {"C:" + fields[0] + ":" + source_position,
                                      "C:" + fields[0] + ":" + std::to_string(normalized.pos)};
      if (!mss::missing(fields[4])) keys.push_back("I:" + fields[4]);
      for (const auto& alias : mss::split(fields[5], ';'))
        if (!mss::missing(alias)) keys.push_back("I:" + alias);
      result.keys(offset, keys);
      add(result.record_rows, 1); add(result.allele_records, 1);
      result.max_record_bytes = std::max(result.max_record_bytes, records.buffer.bytes - offset);
    }
  }
  if (input.bad()) throw std::runtime_error("Source read failed");
  U64 offsets_bytes = 0;
  if (result.final_key_rows > std::numeric_limits<U64>::max() / 8)
    throw std::runtime_error("Offset byte arithmetic overflow");
  offsets_bytes = 8 * result.final_key_rows;
  // Deliberately no admission boolean: small files, source copy, allocation,
  // concurrent work and actual free-space reserves are coordinator obligations.
  std::cout << "{\"format_version\":" << version
#define FIELD(name) << ",\"" #name "\":" << result.name
    FIELD(source_rows) FIELD(source_bytes) FIELD(record_rows) FIELD(site_records)
    FIELD(allele_records) FIELD(invalid_positions) FIELD(raw_alt_tokens) FIELD(rejected_alts)
    FIELD(raw_key_rows) FIELD(raw_sort_bytes) FIELD(final_key_rows) FIELD(keys_bytes)
    FIELD(max_source_line_bytes) FIELD(max_key_row_bytes) FIELD(max_record_bytes)
    FIELD(max_raw_keys_per_record)
#undef FIELD
    << ",\"records_bytes\":" << records.buffer.bytes
    << ",\"rejected_bytes\":" << rejected.buffer.bytes
    << ",\"offsets_bytes\":" << offsets_bytes << "}\n";
}
} // namespace

int main(int argc, char** argv) {
  if (argc != 4) { std::cerr << "Expected FASTA SITES FORMAT; emits counts only\n"; return 2; }
  try {
    const std::string version(argv[3]);
    if (version != "2" && version != "3") throw std::runtime_error("Expected exact format 2 or 3");
    count(argv[1], argv[2], version == "2" ? 2 : 3);
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
