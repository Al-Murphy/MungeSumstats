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
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <utility>
namespace mss {
// A fixed-size portable page cache bounds native reference memory independently
// of genome/reference size and avoids an OS seek for every binary-search step.
class BlockFile {
  std::ifstream file;std::uint64_t length=0,reads=0;std::size_t width,capacity;
  struct Page {std::string data;std::list<std::uint64_t>::iterator recent;};
  std::unordered_map<std::uint64_t,Page> pages;std::list<std::uint64_t> recent;
  const std::string& page(std::uint64_t number) {
    auto found=pages.find(number);
    if(found!=pages.end()){recent.splice(recent.begin(),recent,found->second.recent);return found->second.data;}
    if(pages.size()==capacity){pages.erase(recent.back());recent.pop_back();}
    std::string data(width,'\0');file.clear();file.seekg(number*width);file.read(&data[0],width);
    data.resize(static_cast<std::size_t>(file.gcount()));++reads;
    if(file.bad())throw std::runtime_error("Reference block read failed");
    recent.push_front(number);return pages.emplace(number,Page{std::move(data),recent.begin()}).first->second.data;
  }
public:
  explicit BlockFile(const std::string& path,std::size_t block_size=65536,std::size_t max_blocks=128):file(path,std::ios::binary),width(block_size),capacity(max_blocks) {
    if(!file||!width||!capacity)throw std::runtime_error("Cannot open bounded reference file: "+path);
    file.seekg(0,std::ios::end);auto end=file.tellg();if(end<0)throw std::runtime_error("Cannot measure reference file");length=static_cast<std::uint64_t>(end);
  }
  std::uint64_t size()const{return length;}
  std::uint64_t disk_reads()const{return reads;}
  std::size_t resident_bytes()const{std::size_t n=0;for(auto& p:pages)n+=p.second.data.size();return n;}
  std::string bytes(std::uint64_t pos,std::size_t n) {
    std::string out;
    while(n&&pos<length){const auto& data=page(pos/width);auto start=pos%width;
      if(start>=data.size())throw std::runtime_error("Truncated reference file");
      auto take=std::min(n,data.size()-static_cast<std::size_t>(start));out.append(data,start,take);pos+=take;n-=take;}
    return out;
  }
  std::string line(std::uint64_t pos) {
    std::string out;
    while(pos<length){const auto& data=page(pos/width);auto start=static_cast<std::size_t>(pos%width);
      if(start>=data.size())throw std::runtime_error("Truncated reference file");
      auto end=data.find('\n',start);auto take=(end==std::string::npos?data.size():end)-start;
      out.append(data,start,take);if(end!=std::string::npos)break;pos+=take;}
    return out;
  }
};
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
  std::ifstream file;BlockFile blocks; std::map<std::string,Contig> contigs;
public:
  explicit Fasta(const std::string& path,const std::string& fai=""):file(path,std::ios::binary),blocks(path,65536,64) {
    if(!file) throw std::runtime_error("Cannot open plain FASTA: "+path);
    char magic[2]={0,0}; file.read(magic,2); file.clear(); file.seekg(0);
    if(static_cast<unsigned char>(magic[0])==31 && static_cast<unsigned char>(magic[1])==139) throw std::runtime_error("Compressed FASTA is unsupported; provide plain FASTA");
    if(!fai.empty()) {
      std::ifstream in(fai); std::string line;
      if(!in) throw std::runtime_error("Missing FASTA index");
      while(std::getline(in,line)) {auto v=split(line,'\t'); if(v.size()!=5)throw std::runtime_error("Invalid FASTA index"); contigs[v[0]]={std::stoull(v[1]),std::stoull(v[2]),std::stoull(v[3]),std::stoull(v[4])};}
      return;
    }
    std::string line,name; std::uint64_t previous=0,offset=0;
    while(file) {
      if(!std::getline(file,line))break;
      // getline removes LF but retains CR. Count physical bytes before stripping
      // CR; eof means the final line had no delimiter. Avoid two seeks per line.
      const auto line_offset=offset;
      const auto bytes=static_cast<std::uint64_t>(line.size())+(file.eof()?0:1);
      offset+=bytes;
      if(!line.empty() && line.back()=='\r')line.pop_back();
      if(!line.empty() && line[0]=='>') {
        name=line.substr(1,line.find_first_of(" \t",1)-1);
        if(name.empty() || contigs.count(name))throw std::runtime_error("Empty or duplicate FASTA contig");
        contigs[name]={}; previous=0;
      } else {
        if(name.empty() || line.empty())throw std::runtime_error("Invalid FASTA layout");
        auto& c=contigs[name];
        if(c.size && previous!=c.bases)throw std::runtime_error("FASTA must have fixed-width sequence lines except final line");
        if(!c.size) {c.offset=line_offset;c.bases=line.size();c.bytes=bytes;}
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
    while(n) {auto take=std::min<std::uint64_t>(n,c.bases-zero%c.bases);auto part=blocks.bytes(c.offset+(zero/c.bases)*c.bytes+zero%c.bases,take);if(part.size()!=take)return "";result+=part;zero+=take;n-=take;} return upper(result);
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
  std::set<std::int64_t> original_positions;
  std::string key()const{return chr+":"+std::to_string(pos)+":"+ref+":"+alt;}
};
inline void merge_candidate_context(Candidate& existing,const Candidate& candidate,bool identifiers=true) {
  existing.original_positions.insert(candidate.original_positions.begin(),candidate.original_positions.end());
  if(identifiers){existing.identifiers.insert(candidate.identifiers.begin(),candidate.identifiers.end());existing.id=existing.identifiers.size()==1?*existing.identifiers.begin():"";}
  std::set<std::string> alleles;for(auto& allele:split(existing.alleles,','))alleles.insert(allele);for(auto& allele:split(candidate.alleles,','))alleles.insert(allele);
  existing.alleles="";for(auto& allele:alleles){if(!existing.alleles.empty())existing.alleles+=",";existing.alleles+=allele;}
}
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
struct Input {
  std::string chr;std::int64_t pos;std::string id,other,effect_allele;double effect,eaf;std::string frequency_type;
  std::string source_ref,source_alt,source_strand;bool source_contract=false;
  Input(std::string chromosome,std::int64_t position,std::string identifier,std::string other_allele,
        std::string measured_allele,double beta,double frequency,std::string frequency_kind):
    chr(std::move(chromosome)),pos(position),id(std::move(identifier)),other(std::move(other_allele)),
    effect_allele(std::move(measured_allele)),effect(beta),eaf(frequency),frequency_type(std::move(frequency_kind)){}
};
struct Result {std::string status="unmapped",chr,ref,alt,id;std::int64_t pos=0;double effect=std::numeric_limits<double>::quiet_NaN(),eaf=std::numeric_limits<double>::quiet_NaN();std::size_t candidate_count=0;bool swapped=false,complemented=false;};
class Reference {
  BlockFile keys,offsets,records;std::uint64_t count;Fasta fa;
  std::string entry(std::uint64_t i){auto bytes=offsets.bytes(i*8,8);if(bytes.size()!=8)throw std::runtime_error("Truncated disk index");std::uint64_t p=0;for(int j=0;j<8;++j)p|=std::uint64_t(static_cast<unsigned char>(bytes[j]))<<(8*j);return keys.line(p);}
public:
  Reference(const std::string& dir,const std::string& fasta):keys(dir+"/keys"),offsets(dir+"/offsets"),records(dir+"/records"),fa(fasta,dir+"/fasta.fai") {
    std::ifstream meta(dir+"/core-version");std::string version;std::getline(meta,version);if(version!="mss-preservation-2")throw std::runtime_error("Incomplete or unsupported reference index");auto size=offsets.size();if(size%8)throw std::runtime_error("Invalid offsets file");count=size/8;
  }
  std::vector<Candidate> lookup(const std::string& key,std::size_t cap=10000) {
    std::uint64_t low=0,high=count;while(low<high){auto mid=low+(high-low)/2;auto s=entry(mid);if(s.substr(0,s.find('\t'))<key)low=mid+1;else high=mid;}
    std::vector<Candidate> out;std::map<std::string,std::size_t> seen;std::size_t examined=0;
    for(auto i=low;i<count;++i){auto s=entry(i);auto tab=s.find('\t');if(s.substr(0,tab)!=key)break;if(++examined>cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");auto offset=std::stoull(s.substr(tab+1));s=records.line(offset);auto v=split(s,'\t');if(v.size()!=8)throw std::runtime_error("Corrupt reference record");Candidate c{v[0],v[2],v[3],v[4],v[6],v[7],std::stoll(v[1]),std::stoll(v[5]),{},{}};
      c.original_positions.insert(c.original_pos);
      if(!missing(c.id))c.identifiers.insert(c.id);
      auto found=seen.find(c.key());
      if(found==seen.end()){seen[c.key()]=out.size();out.push_back(c);}
      else {
        merge_candidate_context(out[found->second],c);
      }
    }
    return out;
  }
  Result harmonise(const Input& input,std::size_t cap=10000) {
    if(input.source_contract)return harmonise_declared(input,cap);
    Result result;std::vector<Candidate> coordinate,ids,candidates;
    std::set<std::string> indel_events;
    if(!missing(input.chr)&&input.pos>0) {
      std::set<std::string> positions={std::to_string(input.pos)};std::map<std::string,std::size_t> seen;
      // Query representations can start at any equivalent repeat anchor, not
      // just the reference input anchor or its canonical anchor.
      auto other=upper(input.other),effect=upper(input.effect_allele);
      if(dna(other)&&dna(effect)&&other!=effect)for(int orientation=0;orientation<4;++orientation){
        auto a=orientation>=2?complement(other):other;auto b=orientation>=2?complement(effect):effect;
        if(orientation%2)std::swap(a,b);
        auto n=normalize(fa,input.chr,input.pos,a,b);if(n.status=="ok"){positions.insert(std::to_string(n.pos));if(other.size()!=effect.size())indel_events.insert(n.chr+":"+std::to_string(n.pos)+":"+n.ref+":"+n.alt);}
      }
      for(auto& pos:positions)for(auto& candidate:lookup("C:"+input.chr+":"+pos,cap)) {
        auto found=seen.find(candidate.key());
        if(found==seen.end()){if(coordinate.size()>=cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");seen[candidate.key()]=coordinate.size();coordinate.push_back(candidate);}
        else merge_candidate_context(coordinate[found->second],candidate);
      }
    }
    if(!missing(input.id))ids=lookup("I:"+input.id,cap);
    if(!coordinate.empty()&&!ids.empty()) {std::map<std::string,std::size_t> wanted;for(std::size_t i=0;i<coordinate.size();++i)wanted[coordinate[i].key()]=i;for(auto& c:ids){auto found=wanted.find(c.key());if(found!=wanted.end()){merge_candidate_context(c,coordinate[found->second],false);candidates.push_back(c);}}if(candidates.empty()){result.status="coordinate_id_conflict";return result;}}
    else if(!ids.empty()) {if(!missing(input.chr)&&input.pos>0){result.status="coordinate_id_conflict";return result;}candidates=ids;}
    else candidates=coordinate;
    result.candidate_count=candidates.size();if(candidates.empty())return result;
    if(ids.empty()&&indel_events.size()>1){result.status="ambiguous";return result;}
    auto other=upper(input.other),effect=upper(input.effect_allele);
    if(!dna(other)||!dna(effect)||other==effect){result.status="unsupported_allele";return result;}
    struct Hit{Candidate candidate;bool swap,comp;std::int64_t anchor;};std::vector<Hit> hits;std::set<std::string> hitkeys;
    for(auto& c:candidates)for(int orientation=0;orientation<4;++orientation) {
      bool swap=orientation%2,comp=orientation>=2;auto a=comp?complement(other):other;auto b=comp?complement(effect):effect;if(swap)std::swap(a,b);
      auto chr=missing(input.chr)?c.chr:input.chr;
      std::set<std::int64_t> anchors=input.pos>0?std::set<std::int64_t>{input.pos}:c.original_positions;
      if(input.pos<=0)anchors.insert(c.pos);
      for(auto pos:anchors){auto n=normalize(fa,chr,pos,a,b);
        if(n.status=="ok"&&n.chr==c.chr&&n.pos==c.pos&&n.ref==c.ref&&n.alt==c.alt){auto key=c.key()+":"+(swap?"1":"0")+":"+(comp?"1":"0");if(hitkeys.insert(key).second)hits.push_back({c,swap,comp,pos});}}

    }
    if(hits.size()>1){result.status="ambiguous";return result;}
    // A unique REF-ALT hit does not prove the source strand: the opposite
    // strand may instead measure two alternate alleles at the same site.
    // Neither a site identifier nor EAF selects between these explanations.
    std::set<std::int64_t> resolved_anchors;
    if(hits.size()==1){resolved_anchors.insert(hits[0].anchor);resolved_anchors.insert(hits[0].candidate.pos);}
    const bool alt_alt=alt_alt_context(input,candidates,other,effect,cap,resolved_anchors);
    if(hits.empty()) {
      if(alt_alt){result.status="alt_alt_contrast";return result;}
      result.status="allele_mismatch";return result;
    }
    if(alt_alt){result.status="ambiguous";return result;}
    auto h=hits[0];auto c=h.candidate;result.status=h.comp?(h.swap?"swapped_complement":"complement"):(h.swap?"swap":"match");result.chr=c.chr;result.pos=c.pos;result.ref=c.ref;result.alt=c.alt;result.id=c.id;result.swapped=h.swap;result.complemented=h.comp;result.effect=h.swap?-input.effect:input.effect;
    if(input.frequency_type=="EAF"&&std::isfinite(input.eaf)&&input.eaf>=0&&input.eaf<=1){
      if(!h.swap)result.eaf=input.eaf;
      else if(biallelic_context(c.key(),c.chr,{c},input.pos,cap)) {
        result.eaf=1-input.eaf;
      }
    }
    return result;
  }
private:
  bool alt_alt_context(const Input& input,const std::vector<Candidate>& candidates,
                       const std::string& other,const std::string& effect,
                       std::size_t cap,const std::set<std::int64_t>& resolved_anchors) {
    std::set<std::pair<std::string,std::int64_t>> visited;
    for(const auto& c:candidates) {
      if(!missing(input.chr)&&input.chr!=c.chr)continue;
      for(auto position:c.original_positions) {
        if(!visited.emplace(c.chr,position).second)continue;
        for(const auto& site:lookup("S:"+c.chr+":"+std::to_string(position),cap)) {
          if(site.original_pos!=position||site.chr!=c.chr||!dna(site.original_ref)||
             fa.sequence(site.chr,site.original_pos,site.original_ref.size())!=site.original_ref)continue;
          // Own all strings until the row-scoped signature maps are destroyed.
          // Views avoid copying a long allele once for every possible trim.
          const auto alleles=split(site.alleles,',');
          std::vector<std::size_t> valid;
          for(std::size_t i=0;i<alleles.size();++i)
            if(dna(alleles[i])&&alleles[i]!=site.original_ref)valid.push_back(i);
          if(valid.size()<2)continue;
          for(int orientation=0;orientation<2;++orientation) {
            const std::string observed[2]={orientation?complement(other):other,
                                           orientation?complement(effect):effect};
            // One linear scan of ALT strings joins compatible common flanks;
            // no Cartesian product of ALT pairs is enumerated. A shared flank
            // may itself differ from REF, so retain it as part of the contrast.
            using Signature=std::tuple<std::int64_t,std::int64_t,std::string_view,std::string_view>;
            std::map<Signature,unsigned> signatures;
            for(auto i:valid)for(unsigned which=0;which<2;++which) {
              const auto& allele=alleles[i];const auto& query=observed[which];
              const auto length=static_cast<std::int64_t>(site.original_ref.size())+
                                static_cast<std::int64_t>(query.size())-static_cast<std::int64_t>(allele.size());
              if(length<1||query.size()>allele.size())continue;
              std::size_t offset=0;
              if(input.pos>0) {
                if(input.pos<site.original_pos)continue;
                offset=static_cast<std::size_t>(input.pos-site.original_pos);
              } else offset=allele.find(query);
              while(offset!=std::string::npos&&offset<=allele.size()-query.size()) {
                if(allele.compare(offset,query.size(),query)==0) {
                  const auto anchor=site.original_pos+static_cast<std::int64_t>(offset);
                  const auto view=std::string_view(allele);
                  auto& mask=signatures[{anchor,length,view.substr(0,offset),view.substr(offset+query.size())}];
                  mask|=1U<<which;if(mask==3)return true;
                }
                if(input.pos>0)break;
                offset=allele.find(query,offset+1);
              }
            }
            // Reference-shared padding and repeat-equivalent indel anchors can
            // also place the observed pair outside its original literal span.
            // Both ALT events must match under ONE inferred source REF span.
            // At most the original anchor and two actual single-hit anchors
            // are examined for unknown input position, independently of ALT count.
            std::set<std::int64_t> anchors=input.pos>0?std::set<std::int64_t>{input.pos}:resolved_anchors;
            if(input.pos<=0)anchors.insert(site.original_pos);
            for(auto anchor:anchors) {
              if(anchor<1)continue;
              std::map<std::int64_t,unsigned> matched;
              std::map<std::int64_t,std::pair<Normalized,Normalized>> query_events;
              for(auto i:valid) {
                const auto& allele=alleles[i];
                const auto event=normalize(fa,site.chr,site.original_pos,site.original_ref,allele);
                if(event.status!="ok")continue;
                for(unsigned which=0;which<2;++which) {
                  const auto length=static_cast<std::int64_t>(site.original_ref.size())+
                                    static_cast<std::int64_t>(observed[which].size())-static_cast<std::int64_t>(allele.size());
                  if(length<1)continue;
                  auto found=query_events.find(length);
                  if(found==query_events.end()) {
                    const auto ref=fa.sequence(site.chr,anchor,static_cast<std::size_t>(length));
                    found=query_events.emplace(length,std::make_pair(
                      normalize(fa,site.chr,anchor,ref,observed[0]),
                      normalize(fa,site.chr,anchor,ref,observed[1]))).first;
                  }
                  const auto& query=which?found->second.second:found->second.first;
                  if(query.status=="ok"&&query.chr==event.chr&&query.pos==event.pos&&
                     query.ref==event.ref&&query.alt==event.alt) {
                    auto& mask=matched[length];mask|=1U<<which;if(mask==3)return true;
                  }
                }
              }
            }
          }
        }
      }
    }
    return false;
  }

  // Count canonical alternate events, not raw allele strings. The same deletion
  // can be spelled GT>G, TG>G or TT>T at different anchors. Keep each site's
  // own REF/position when interpreting its ALT list; merged candidate strings
  // no longer retain that pairing and must not be used as frequency evidence.
  bool biallelic_context(const std::string& event,const std::string& chr,
                         const std::vector<Candidate>& matching,
                         std::int64_t input_pos,std::size_t cap) {
    if(matching.empty())return false; // FASTA alone does not prove two alleles.
    std::set<std::int64_t> positions;
    if(input_pos>0)positions.insert(input_pos);
    for(const auto& candidate:matching) {
      if(candidate.key()!=event)return false;
      positions.insert(candidate.pos);
      positions.insert(candidate.original_positions.begin(),candidate.original_positions.end());
    }
    bool witnessed=false;
    for(auto position:positions) {
      // S records are grouped by original position and original REF. Their ALT
      // lists therefore remain correctly paired even when lookup deduplicates.
      for(const auto& site:lookup("S:"+chr+":"+std::to_string(position),cap)) {
        const auto alleles=split(site.alleles,',');
        if(alleles.empty())return false;
        for(const auto& allele:alleles) {
          auto normalized=normalize(fa,chr,site.original_pos,site.original_ref,allele);
          if(normalized.status!="ok")return false; // Includes symbolic context.
          auto key=normalized.chr+":"+std::to_string(normalized.pos)+":"+normalized.ref+":"+normalized.alt;
          if(key!=event)return false; // A genuine additional alternate event.
          witnessed=true;
        }
      }
    }
    return witnessed;
  }

  Result harmonise_declared(const Input& input,std::size_t cap) {
    Result result;
    if(input.source_strand!="+"&&input.source_strand!="-"){result.status="source_strand_unknown";return result;}
    if(missing(input.chr)||input.pos<1){result.status="source_coordinates_missing";return result;}
    auto ref=upper(input.source_ref),alt=upper(input.source_alt);
    if(missing(ref)||missing(alt)){result.status="source_reference_missing";return result;}
    if(!dna(ref)||!dna(alt)||ref==alt){result.status="source_allele_unsupported";return result;}
    auto other=upper(input.other),effect=upper(input.effect_allele);
    bool swap=other==alt&&effect==ref;
    if(!(other==ref&&effect==alt)&&!swap){result.status="source_allele_conflict";return result;}
    bool comp=input.source_strand=="-";
    if(comp){ref=complement(ref);alt=complement(alt);}
    auto n=normalize(fa,input.chr,input.pos,ref,alt);
    if(n.status!="ok"){result.status="source_reference_mismatch";return result;}
    auto expected=n.chr+":"+std::to_string(n.pos)+":"+n.ref+":"+n.alt;
    std::vector<Candidate> ids;
    if(!missing(input.id))ids=lookup("I:"+input.id,cap);
    if(!ids.empty()&&std::none_of(ids.begin(),ids.end(),[&](const Candidate& c){return c.key()==expected;})) {
      result.status="source_id_conflict";return result;
    }
    std::vector<Candidate> matching;
    for(auto pos:std::set<std::int64_t>{input.pos,n.pos})for(auto& c:lookup("C:"+n.chr+":"+std::to_string(pos),cap))
      if(c.key()==expected)matching.push_back(c);
    for(auto& c:ids)if(c.key()==expected)matching.push_back(c);
    std::set<std::string> identifiers;
    for(auto& c:matching)identifiers.insert(c.identifiers.begin(),c.identifiers.end());
    result.status=matching.empty()?"sequence_only":(comp?(swap?"swapped_complement":"complement"):(swap?"swap":"match"));
    result.chr=n.chr;result.pos=n.pos;result.ref=n.ref;result.alt=n.alt;
    if(identifiers.size()==1)result.id=*identifiers.begin();
    result.candidate_count=matching.empty()?0:1;result.swapped=swap;result.complemented=comp;
    result.effect=swap?-input.effect:input.effect;
    if(input.frequency_type=="EAF"&&std::isfinite(input.eaf)&&input.eaf>=0&&input.eaf<=1){
      if(!swap)result.eaf=input.eaf;
      else if(biallelic_context(expected,n.chr,matching,input.pos,cap)) {
        result.eaf=1-input.eaf;
      }
    }
    return result;
  }
};
}
#endif
