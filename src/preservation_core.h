#ifndef MSS_PRESERVATION_CORE_H
#define MSS_PRESERVATION_CORE_H
// Portable C++17 core. No R, htslib, shell commands, or whole-reference cache.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace mss {
inline std::vector<std::string> split(const std::string& s,char sep) {
  std::vector<std::string> out; std::string v; std::istringstream in(s);
  while(std::getline(in,v,sep)) out.push_back(v);
  if(!s.empty() && s.back()==sep) out.push_back(""); return out;
}
inline std::string upper(std::string s) { for(char& c:s) c=static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return s; }
inline bool dna(const std::string& s) { return !s.empty() && s.find_first_not_of("ACGT")==std::string::npos; }
inline bool missing(const std::string& s) {return s.empty() || s=="." || s=="NA";}
inline std::string complement(std::string s) {
  std::reverse(s.begin(),s.end()); for(char& c:s) { if(c=='A')c='T'; else if(c=='T')c='A'; else if(c=='C')c='G'; else if(c=='G')c='C'; } return s;
}
struct Contig {std::uint64_t size=0,offset=0,bases=0,bytes=0;};
class Fasta {
  std::ifstream file; std::map<std::string,Contig> contigs;
public:
  explicit Fasta(const std::string& path,const std::string& fai=""):file(path,std::ios::binary) {
    if(!file) throw std::runtime_error("Cannot open plain FASTA: "+path);
    char magic[2]={0,0}; file.read(magic,2); file.clear(); file.seekg(0);
    if(static_cast<unsigned char>(magic[0])==31 && static_cast<unsigned char>(magic[1])==139) throw std::runtime_error("Compressed FASTA is unsupported; provide plain FASTA");
    if(!fai.empty()) {
      std::ifstream in(fai); std::string line;
      if(!in) throw std::runtime_error("Missing FASTA index");
      while(std::getline(in,line)) {auto v=split(line,'\t'); if(v.size()!=5)throw std::runtime_error("Invalid FASTA index"); contigs[v[0]]={std::stoull(v[1]),std::stoull(v[2]),std::stoull(v[3]),std::stoull(v[4])};}
      return;
    }
    std::string line,name; std::uint64_t previous=0;
    while(file) {
      auto offset=file.tellg(); if(!std::getline(file,line))break;
      auto end=file.tellg(); auto bytes=end<0 ? static_cast<std::uint64_t>(line.size()) : static_cast<std::uint64_t>(end-offset);
      if(!line.empty() && line.back()=='\r')line.pop_back();
      if(!line.empty() && line[0]=='>') {
        name=line.substr(1,line.find_first_of(" \t",1)-1);
        if(name.empty() || contigs.count(name))throw std::runtime_error("Empty or duplicate FASTA contig");
        contigs[name]={}; previous=0;
      } else {
        if(name.empty() || line.empty())throw std::runtime_error("Invalid FASTA layout");
        auto& c=contigs[name];
        if(c.size && previous!=c.bases)throw std::runtime_error("FASTA must have fixed-width sequence lines except final line");
        if(!c.size) {c.offset=static_cast<std::uint64_t>(offset);c.bases=line.size();c.bytes=bytes;}
        if(line.size()>c.bases)throw std::runtime_error("FASTA line too long");
        c.size+=line.size(); previous=line.size();
      }
    }
    file.clear();
    if(contigs.empty())throw std::runtime_error("Empty FASTA");
  }
  void save(const std::string& path)const {std::ofstream out(path);for(const auto& p:contigs)out<<p.first<<'\t'<<p.second.size<<'\t'<<p.second.offset<<'\t'<<p.second.bases<<'\t'<<p.second.bytes<<'\n';if(!out)throw std::runtime_error("Cannot write FASTA index");}
  std::string sequence(const std::string& chr,std::int64_t pos,std::size_t n) {
    auto it=contigs.find(chr);if(it==contigs.end() || pos<1)return "";
    const auto& c=it->second; if(static_cast<std::uint64_t>(pos-1)+n>c.size)return "";
    std::string result; result.reserve(n); auto zero=static_cast<std::uint64_t>(pos-1);
    while(n) {auto take=std::min<std::uint64_t>(n,c.bases-zero%c.bases);std::string part(take,' ');file.clear();file.seekg(c.offset+(zero/c.bases)*c.bytes+zero%c.bases);file.read(&part[0],take);if(!file)return "";result+=part;zero+=take;n-=take;} return upper(result);
  }
};
struct Normalized {std::string status,chr,ref,alt;std::int64_t pos=0;};
inline Normalized normalize(Fasta& fa,const std::string& chr,std::int64_t pos,std::string ref,std::string alt) {
  ref=upper(ref);alt=upper(alt); Normalized out{"unsupported_allele",chr,ref,alt,pos};
  if(!dna(ref)||!dna(alt)||ref==alt)return out;
  if(fa.sequence(chr,pos,ref.size())!=ref) {out.status="reference_mismatch";return out;}
  // Right trim; an empty allele is re-anchored one base to the left. This also
  // rotates repeat indels to the leftmost equivalent reference representation.
  while(ref.back()==alt.back()) {
    if((ref.size()==1 || alt.size()==1) && pos==1)break;
    ref.pop_back();alt.pop_back();
    if(ref.empty()||alt.empty()) {auto base=fa.sequence(chr,pos-1,1);if(base.empty())break;ref=base+ref;alt=base+alt;--pos;}
  }
  while(ref.size()>1 && alt.size()>1 && ref.front()==alt.front()) {ref.erase(0,1);alt.erase(0,1);++pos;}
  return {"ok",chr,ref,alt,pos};
}
inline void write64(std::ostream& out,std::uint64_t x) {for(int i=0;i<8;++i)out.put(static_cast<char>((x>>(8*i))&255));}
inline std::uint64_t read64(std::istream& in) {std::uint64_t x=0;for(int i=0;i<8;++i){int c=in.get();if(c<0)throw std::runtime_error("Truncated disk index");x|=std::uint64_t(c)<<(8*i);}return x;}
// External sorting has a bounded row buffer and at most 32 open merge inputs.
class Sorter {
  std::string dir;std::size_t limit,serial=0;std::vector<std::string> buffer;std::vector<std::vector<std::string>> levels;
  void flush() {if(buffer.empty())return;std::sort(buffer.begin(),buffer.end());std::string p=dir+"/run"+std::to_string(serial++);std::ofstream out(p,std::ios::binary);for(auto& s:buffer)out<<s<<'\n';if(!out)throw std::runtime_error("Index run write failed");out.close();buffer.clear();store_run(p,0);}
  std::string merge(const std::vector<std::string>& paths) {
    std::string p=dir+"/run"+std::to_string(serial++);std::ofstream out(p,std::ios::binary);
    std::vector<std::ifstream> inputs(paths.size());using Item=std::pair<std::string,std::size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>> heap;
    for(std::size_t i=0;i<paths.size();++i){inputs[i].open(paths[i],std::ios::binary);std::string s;if(std::getline(inputs[i],s))heap.push({s,i});}
    while(!heap.empty()){auto item=heap.top();heap.pop();out<<item.first<<'\n';std::string s;if(std::getline(inputs[item.second],s))heap.push({s,item.second});}
    if(!out)throw std::runtime_error("Index merge failed");out.close();for(auto& input:inputs)input.close();for(auto& s:paths)std::filesystem::remove(s);return p;
  }
  void store_run(const std::string& path,std::size_t level) {
    if(levels.size()<=level)levels.resize(level+1);
    levels[level].push_back(path);
    if(levels[level].size()==32){auto paths=levels[level];levels[level].clear();auto merged=merge(paths);store_run(merged,level+1);}
  }
public:
  Sorter(std::string path,std::size_t n):dir(std::move(path)),limit(n){if(!n)throw std::runtime_error("batch_size must be positive");}
  void add(const std::string& key,std::uint64_t offset){buffer.push_back(key+'\t'+std::to_string(offset));if(buffer.size()>=limit)flush();}
  void finish() {
    flush();std::vector<std::string> runs;for(auto& level:levels)runs.insert(runs.end(),level.begin(),level.end());while(runs.size()>1){std::vector<std::string> next;for(std::size_t i=0;i<runs.size();i+=32){auto end=std::min(i+32,runs.size());std::vector<std::string> group(runs.begin()+i,runs.begin()+end);next.push_back(merge(group));}runs.swap(next);}
    std::ofstream index(dir+"/keys",std::ios::binary),offsets(dir+"/offsets",std::ios::binary);std::uint64_t pos=0;std::string previous,line;
    if(!runs.empty()){std::ifstream in(runs[0],std::ios::binary);while(std::getline(in,line)){if(line==previous)continue;write64(offsets,pos);index<<line<<'\n';pos+=line.size()+1;previous=line;}in.close();std::filesystem::remove(runs[0]);}
    if(!index||!offsets)throw std::runtime_error("Index final write failed");
  }
};
struct Candidate {
  std::string chr,ref,alt,id,original_ref,alleles;std::int64_t pos=0,original_pos=0;std::set<std::string> identifiers;
  std::string key()const{return chr+":"+std::to_string(pos)+":"+ref+":"+alt;}
};
inline void prepare(const std::string& fasta,const std::string& sites,const std::string& dir,std::size_t batch_size) {
  if(std::filesystem::exists(dir))throw std::runtime_error("Reference directory already exists; use a new immutable version path");
  std::filesystem::create_directories(dir);Fasta fa(fasta);fa.save(dir+"/fasta.fai");
  std::ifstream in(sites);std::ofstream data(dir+"/records",std::ios::binary),rejected(dir+"/rejected.tsv");if(!in||!data||!rejected)throw std::runtime_error("Cannot open reference input/output");
  std::string line;std::getline(in,line);if(!line.empty()&&line.back()=='\r')line.pop_back();if(line!="chr\tpos\tref\talt\tid\taliases")throw std::runtime_error("Reference header must be chr,pos,ref,alt,id,aliases (tab-separated)");
  Sorter sorter(dir,batch_size);std::uint64_t row=1;rejected<<"row\tstatus\n";
  while(std::getline(in,line)) {
    ++row;if(!line.empty()&&line.back()=='\r')line.pop_back();auto v=split(line,'\t');if(v.size()!=6)throw std::runtime_error("Reference row must have six fields");std::int64_t pos;
    try{std::size_t used;pos=std::stoll(v[1],&used);if(used!=v[1].size())throw std::invalid_argument("position");}catch(...){rejected<<row<<"\tinvalid_position\n";continue;}
    // Site context exists independently of sequence-resolved candidates. Even
    // symbolic-only source rows are evidence against a biallelic complement.
    auto context_offset=static_cast<std::uint64_t>(data.tellp());
    data<<v[0]<<'\t'<<pos<<'\t'<<upper(v[2])<<"\t\t"<<v[4]<<'\t'<<pos<<'\t'<<upper(v[2])<<'\t'<<upper(v[3])<<'\n';
    sorter.add("S:"+v[0]+":"+std::to_string(pos),context_offset);
    auto alleles=split(upper(v[3]),',');
    for(auto& alt:alleles) {
      auto n=normalize(fa,v[0],pos,v[2],alt);if(n.status!="ok"){rejected<<row<<'\t'<<n.status<<'\n';continue;}
      auto offset=static_cast<std::uint64_t>(data.tellp());
      data<<v[0]<<'\t'<<n.pos<<'\t'<<n.ref<<'\t'<<n.alt<<'\t'<<v[4]<<'\t'<<pos<<'\t'<<upper(v[2])<<'\t'<<upper(v[3])<<'\n';
      sorter.add("C:"+v[0]+":"+std::to_string(pos),offset);sorter.add("C:"+v[0]+":"+std::to_string(n.pos),offset);
      if(!missing(v[4]))sorter.add("I:"+v[4],offset);
      for(auto& alias:split(v[5],';'))if(!missing(alias))sorter.add("I:"+alias,offset);
    }
  }
  if(!data||!rejected)throw std::runtime_error("Reference data write failed");data.close();sorter.finish();
  std::ofstream meta(dir+"/core-version",std::ios::binary);meta<<"mss-preservation-2\n";
}
struct Input {std::string chr;std::int64_t pos;std::string id,other,effect_allele;double effect,eaf;std::string frequency_type;};
struct Result {std::string status="unmapped",chr,ref,alt,id;std::int64_t pos=0;double effect=std::numeric_limits<double>::quiet_NaN(),eaf=std::numeric_limits<double>::quiet_NaN();std::size_t candidate_count=0;bool swapped=false,complemented=false;};
class Reference {
  std::ifstream keys,offsets,records;std::uint64_t count;Fasta fa;
  std::string entry(std::uint64_t i){offsets.clear();offsets.seekg(i*8);auto p=read64(offsets);keys.clear();keys.seekg(p);std::string s;std::getline(keys,s);return s;}
public:
  Reference(const std::string& dir,const std::string& fasta):keys(dir+"/keys",std::ios::binary),offsets(dir+"/offsets",std::ios::binary),records(dir+"/records",std::ios::binary),fa(fasta,dir+"/fasta.fai") {
    std::ifstream meta(dir+"/core-version");std::string version;std::getline(meta,version);if(version!="mss-preservation-2"||!keys||!offsets||!records)throw std::runtime_error("Incomplete or unsupported reference index");offsets.seekg(0,std::ios::end);auto size=offsets.tellg();if(size<0||size%8)throw std::runtime_error("Invalid offsets file");count=static_cast<std::uint64_t>(size)/8;
  }
  std::vector<Candidate> lookup(const std::string& key,std::size_t cap=10000) {
    std::uint64_t low=0,high=count;while(low<high){auto mid=low+(high-low)/2;auto s=entry(mid);if(s.substr(0,s.find('\t'))<key)low=mid+1;else high=mid;}
    std::vector<Candidate> out;std::map<std::string,std::size_t> seen;std::size_t examined=0;
    for(auto i=low;i<count;++i){auto s=entry(i);auto tab=s.find('\t');if(s.substr(0,tab)!=key)break;if(++examined>cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");auto offset=std::stoull(s.substr(tab+1));records.clear();records.seekg(offset);std::getline(records,s);auto v=split(s,'\t');if(v.size()!=8)throw std::runtime_error("Corrupt reference record");Candidate c{v[0],v[2],v[3],v[4],v[6],v[7],std::stoll(v[1]),std::stoll(v[5]),{}};
      if(!missing(c.id))c.identifiers.insert(c.id);
      auto found=seen.find(c.key());
      if(found==seen.end()){seen[c.key()]=out.size();out.push_back(c);}
      else {
        auto& existing=out[found->second];existing.identifiers.insert(c.identifiers.begin(),c.identifiers.end());
        existing.id=existing.identifiers.size()==1?*existing.identifiers.begin():"";
        std::set<std::string> alleles;for(auto& allele:split(existing.alleles,','))alleles.insert(allele);for(auto& allele:split(c.alleles,','))alleles.insert(allele);
        existing.alleles="";for(auto& allele:alleles){if(!existing.alleles.empty())existing.alleles+=",";existing.alleles+=allele;}
      }
    }
    return out;
  }
  Result harmonise(const Input& input,std::size_t cap=10000) {
    Result result;std::vector<Candidate> coordinate,ids,candidates;
    if(!missing(input.chr)&&input.pos>0) {
      std::set<std::string> positions={std::to_string(input.pos)},seen;
      // Query representations can start at any equivalent repeat anchor, not
      // just the reference input anchor or its canonical anchor.
      auto other=upper(input.other),effect=upper(input.effect_allele);
      if(dna(other)&&dna(effect)&&other!=effect)for(int orientation=0;orientation<4;++orientation){
        auto a=orientation>=2?complement(other):other;auto b=orientation>=2?complement(effect):effect;
        if(orientation%2)std::swap(a,b);
        auto n=normalize(fa,input.chr,input.pos,a,b);if(n.status=="ok")positions.insert(std::to_string(n.pos));
      }
      for(auto& pos:positions)for(auto& candidate:lookup("C:"+input.chr+":"+pos,cap))
        if(seen.insert(candidate.key()).second){if(coordinate.size()>=cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");coordinate.push_back(candidate);}
    }
    if(!missing(input.id))ids=lookup("I:"+input.id,cap);
    if(!coordinate.empty()&&!ids.empty()) {std::set<std::string> wanted;for(auto& c:coordinate)wanted.insert(c.key());for(auto& c:ids)if(wanted.count(c.key()))candidates.push_back(c);if(candidates.empty()){result.status="coordinate_id_conflict";return result;}}
    else if(!ids.empty()) {if(!missing(input.chr)&&input.pos>0){result.status="coordinate_id_conflict";return result;}candidates=ids;}
    else candidates=coordinate;
    result.candidate_count=candidates.size();if(candidates.empty())return result;
    auto other=upper(input.other),effect=upper(input.effect_allele);
    if(!dna(other)||!dna(effect)||other==effect){result.status="unsupported_allele";return result;}
    struct Hit{Candidate candidate;bool swap,comp;};std::vector<Hit> hits;std::set<std::string> hitkeys;
    for(auto& c:candidates)for(int orientation=0;orientation<4;++orientation) {
      bool swap=orientation%2,comp=orientation>=2;auto a=comp?complement(other):other;auto b=comp?complement(effect):effect;if(swap)std::swap(a,b);
      auto chr=missing(input.chr)?c.chr:input.chr;
      std::set<std::int64_t> anchors=input.pos>0?std::set<std::int64_t>{input.pos}:std::set<std::int64_t>{c.original_pos,c.pos};
      for(auto pos:anchors){auto n=normalize(fa,chr,pos,a,b);
        if(n.status=="ok"&&n.chr==c.chr&&n.pos==c.pos&&n.ref==c.ref&&n.alt==c.alt){auto key=c.key()+":"+(swap?"1":"0")+":"+(comp?"1":"0");if(hitkeys.insert(key).second)hits.push_back({c,swap,comp});}}

    }
    if(hits.size()>1){result.status="ambiguous";return result;}
    if(hits.empty()) {
      for(auto& c:candidates){
        std::set<std::string> alleles;
        for(auto& allele:split(c.alleles,','))alleles.insert(allele);
        for(auto& site:lookup("S:"+c.chr+":"+std::to_string(c.original_pos),cap))
          if(site.original_pos==c.original_pos)for(auto& allele:split(site.alleles,','))alleles.insert(allele);
        for(int orientation=0;orientation<2;++orientation){auto a=orientation?complement(other):other;auto b=orientation?complement(effect):effect;
          if(a!=c.original_ref&&b!=c.original_ref&&alleles.count(a)&&alleles.count(b)){result.status="alt_alt_contrast";return result;}}
      }
      result.status="allele_mismatch";return result;
    }
    auto h=hits[0];auto c=h.candidate;result.status=h.comp?(h.swap?"swapped_complement":"complement"):(h.swap?"swap":"match");result.chr=c.chr;result.pos=c.pos;result.ref=c.ref;result.alt=c.alt;result.id=c.id;result.swapped=h.swap;result.complemented=h.comp;result.effect=h.swap?-input.effect:input.effect;
    if(input.frequency_type=="EAF"&&std::isfinite(input.eaf)&&input.eaf>=0&&input.eaf<=1){
      if(!h.swap)result.eaf=input.eaf;
      else {
        // Include other IDs at this site: split multiallelic reference rows must
        // not turn an effect-allele frequency into 1-p in the presence of a third allele.
        std::set<std::string> known;known.insert(c.original_ref);
        for(auto& allele:split(c.alleles,','))known.insert(allele);
        for(auto& site:lookup("S:"+c.chr+":"+std::to_string(c.original_pos),cap))
          if(site.original_pos==c.original_pos){known.insert(site.original_ref);for(auto& allele:split(site.alleles,','))known.insert(allele);}
        if(known.size()==2)result.eaf=1-input.eaf;
      }
    }
    return result;
  }
};
}
#endif
