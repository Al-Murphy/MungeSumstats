# Diagnostic copies only. Each complete input line must occur exactly once.
# Insert observations independently of whether the optimized helper is present.
BEGIN {
    after["inline Normalized normalize(Fasta& fa,const std::string& chr,std::int64_t pos,std::string ref,std::string alt) {"] = "  ::reuse_diag::normal(chr,pos,ref,alt);"
    before["      if(dna(other)&&dna(effect)&&other!=effect)for(int orientation=0;orientation<4;++orientation){"] = "      ::reuse_diag::phase=::reuse_diag::seed;"
    before["      for(auto& pos:positions) {"] = "      ::reuse_diag::phase=::reuse_diag::other;"
    before["    for(auto& c:candidates)for(int orientation=0;orientation<4;++orientation) {"] = "    ::reuse_diag::phase=::reuse_diag::hit;"
    before["    std::set<std::int64_t> resolved_anchors;"] = "    ::reuse_diag::phase=::reuse_diag::other;"
    before["      for(const auto& anchor:anchors)for(int orientation=0;orientation<4;++orientation) {"] = "      ::reuse_diag::phase=::reuse_diag::proof;"
    before["    std::size_t contexts=0;"] = "    ::reuse_diag::phase=::reuse_diag::other;"
    after["  std::vector<Candidate> lookup(const std::string& key,std::size_t cap=10000,std::size_t* context_work=nullptr) {"] = "    ::reuse_diag::lookup(key,cap,context_work);"
    before["        matched=true;"] = "        ::reuse_diag::occurrence(i,examined);"
    before["        auto c=decode_candidate(key,row);"] = "        ::reuse_diag::decode(offset,row);"
    after["inline void consume_context_work(std::size_t* remaining) {"] = "  ::reuse_diag::debit(remaining);"
    tail_marker = "      return std::tie(a.chr,a.original_pos,a.original_ref,a.alleles)<std::tie(b.chr,b.original_pos,b.original_ref,b.alleles);"
}
{
    if (in_lookup_tail == 2) {
        if ($0 != "    return out;") exit 2
        print "    ::reuse_diag::returned(out);"
        in_lookup_tail = 0
        returns++
    } else if (in_lookup_tail == 1) {
        if ($0 != "    });") exit 2
        in_lookup_tail = 2
    }
    if ($0 == tail_marker) { in_lookup_tail = 1; markers++ }
    if ($0 in before) { print before[$0]; seen[$0]++ }
    print
    if ($0 in after) { print after[$0]; seen[$0]++ }
}
END {
    if (markers != 1 || returns != 1 || in_lookup_tail) exit 2
    for (line in before) if (seen[line] != 1) exit 2
    for (line in after) if (seen[line] != 1) exit 2
}
