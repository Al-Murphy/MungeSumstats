#ifndef MSS_PRESERVATION_CORE_H
#define MSS_PRESERVATION_CORE_H
// Portable C++17 core. No R, htslib, shell commands, or whole-reference cache.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <new>
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
  const Page* last_page=nullptr;
  std::uint64_t last_number=0;
#ifdef MSS_LOOKUP_COUNTERS
  std::uint64_t requests=0;
#endif
  const std::string& page(std::uint64_t number) {
#ifdef MSS_LOOKUP_COUNTERS
    ++requests;
#endif
    if(last_page&&last_number==number)return last_page->data;
    auto found=pages.find(number);
    if(found!=pages.end()) {
      recent.splice(recent.begin(),recent,found->second.recent);
    } else {
      std::string data;
      if(pages.size()==capacity) {
        auto victim=pages.find(recent.back());
        last_page=nullptr;
        data=std::move(victim->second.data);
        pages.erase(victim);
        recent.pop_back();
      }
      data.resize(width);
      file.clear();file.seekg(number*width);file.read(&data[0],width);
      data.resize(static_cast<std::size_t>(file.gcount()));++reads;
      if(file.bad())throw std::runtime_error("Reference block read failed");
      recent.push_front(number);
      try {
        found=pages.emplace(number,Page{std::move(data),recent.begin()}).first;
      } catch(...) {
        recent.pop_front();throw;
      }
    }
    last_number=number;last_page=&found->second;
    return last_page->data;
  }
public:
  explicit BlockFile(const std::string& path,std::size_t block_size=65536,std::size_t max_blocks=128):file(path,std::ios::binary),width(block_size),capacity(max_blocks) {
    if(!file||!width||!capacity)throw std::runtime_error("Cannot open bounded reference file: "+path);
    file.seekg(0,std::ios::end);auto end=file.tellg();if(end<0)throw std::runtime_error("Cannot measure reference file");length=static_cast<std::uint64_t>(end);
  }
  std::uint64_t size()const{return length;}
  std::uint64_t disk_reads()const{return reads;}
#ifdef MSS_LOOKUP_COUNTERS
  std::uint64_t page_requests()const{return requests;}
#endif
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
  std::string strict_line(std::uint64_t pos) {
    if(pos>=length||(pos&&bytes(pos-1,1)!="\n"))throw std::runtime_error("Typed reference offset is not a record boundary");
    auto out=line(pos);
    if(out.size()>=length-pos||bytes(pos+out.size(),1)!="\n")throw std::runtime_error("Typed reference record is truncated");
    return out;
  }
  std::uint64_t little_endian_offset(std::uint64_t pos) {
    if(pos>length||length-pos<8)throw std::runtime_error("Truncated disk index");
    const auto& data=page(pos/width);
    const auto start=static_cast<std::size_t>(pos%width);
    std::string crossing;
    const char* digits;
    if(start<=data.size()&&data.size()-start>=8)digits=data.data()+start;
    else {
      crossing=bytes(pos,8);
      if(crossing.size()!=8)throw std::runtime_error("Truncated disk index");
      digits=crossing.data();
    }
    std::uint64_t value=0;
    for(unsigned i=0;i<8;++i)value|=std::uint64_t(static_cast<unsigned char>(digits[i]))<<(8*i);
    return value;
  }
  // A visitor must consume the view before it next reads this BlockFile.
  template<class Visitor>
  void visit_line(std::uint64_t pos,bool strict,const Visitor& visit) {
    if(strict) {
      if(pos>=length)throw std::runtime_error("Typed reference offset is not a record boundary");
      if(pos) {
        const auto& before=page((pos-1)/width);
        const auto index=static_cast<std::size_t>((pos-1)%width);
        if(index>=before.size())throw std::runtime_error("Truncated reference file");
        if(before[index]!='\n')throw std::runtime_error("Typed reference offset is not a record boundary");
      }
    }
    if(pos>=length){visit(std::string_view());return;}
    const auto& data=page(pos/width);
    const auto start=static_cast<std::size_t>(pos%width);
    if(start>=data.size())throw std::runtime_error("Truncated reference file");
    const auto end=data.find('\n',start);
    if(end!=std::string::npos){visit(std::string_view(data).substr(start,end-start));return;}
    const auto owned=strict?strict_line(pos):line(pos);
    visit(std::string_view(owned));
  }
  // Advisory search scratch only: never materialize an unbounded crossing row.
  // false means the caller must discard speculation and use ordinary lookup.
  bool bounded_line(std::uint64_t pos,bool strict,char* output,std::size_t capacity,
                    std::size_t& used) {
    used=0;
    if(strict) {
      if(pos>=length)throw std::runtime_error("Typed reference offset is not a record boundary");
      if(pos) {
        const auto& before=page((pos-1)/width);
        const auto at=static_cast<std::size_t>((pos-1)%width);
        if(at>=before.size())throw std::runtime_error("Truncated reference file");
        if(before[at]!='\n')throw std::runtime_error("Typed reference offset is not a record boundary");
      }
    }
    while(pos<length) {
      const auto& data=page(pos/width);
      const auto start=static_cast<std::size_t>(pos%width);
      if(start>=data.size())throw std::runtime_error("Truncated reference file");
      const auto end=data.find('\n',start);
      const auto take=(end==std::string::npos?data.size():end)-start;
      if(take>capacity-used)return false;
      std::copy_n(data.data()+start,take,output+used);used+=take;
      if(end!=std::string::npos)return true;
      pos+=take;
    }
    if(strict)throw std::runtime_error("Typed reference record is truncated");
    return true;
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
// Format 3 has complete typed records, never an unchecked skipped tail.
inline bool typed_text(std::string_view text) {
  for(std::size_t i=0;i<text.size();) {
    const auto first=static_cast<unsigned char>(text[i]);
    if(first<128){if(!first||first=='\t'||first=='\n')return false;++i;continue;}
    std::size_t more;unsigned char low=128,high=191;
    if(first>=194&&first<=223)more=1;
    else if(first>=224&&first<=239){more=2;if(first==224)low=160;if(first==237)high=159;}
    else if(first>=240&&first<=244){more=3;if(first==240)low=144;if(first==244)high=143;}
    else return false;
    if(more>=text.size()-i)return false;
    const auto second=static_cast<unsigned char>(text[i+1]);if(second<low||second>high)return false;
    for(std::size_t j=2;j<=more;++j){const auto byte=static_cast<unsigned char>(text[i+j]);if(byte<128||byte>191)return false;}
    i+=more+1;
  }
  return true;
}
inline std::uint64_t typed_uint(std::string_view text) {
  if(text.empty()||(text.size()>1&&text[0]=='0'))throw std::runtime_error("Typed reference has an invalid unsigned integer");
  std::uint64_t value=0;
  for(char c:text){if(c<'0'||c>'9'||value>(std::numeric_limits<std::uint64_t>::max()-static_cast<unsigned>(c-'0'))/10)throw std::runtime_error("Typed reference has an invalid unsigned integer");value=value*10+static_cast<unsigned>(c-'0');}
  return value;
}
inline std::int64_t typed_int(std::string_view text,bool positive) {
  if(text.empty())throw std::runtime_error("Typed reference has an invalid signed integer");
  const bool negative=text[0]=='-';const auto start=negative?1U:0U;
  if(start==text.size())throw std::runtime_error("Typed reference has an invalid signed integer");
  const auto magnitude=typed_uint(text.substr(start));
  const auto maximum=static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
  if(magnitude>maximum+(negative?1U:0U))throw std::runtime_error("Typed reference signed integer overflows");
  const auto value=negative?(magnitude?-(static_cast<std::int64_t>(magnitude-1))-1:0):static_cast<std::int64_t>(magnitude);
  if((positive&&value<1)||(negative&&!magnitude))throw std::runtime_error("Typed reference has a noncanonical signed integer");
  return value;
}
inline void typed_key(std::string_view key) {
  if(key.size()<2||key[1]!=':'||!typed_text(key))throw std::runtime_error("Typed reference has an invalid lookup key");
  if(key[0]=='I'){if(key.size()==2)throw std::runtime_error("Typed reference has an empty identifier key");return;}
  if(key[0]!='C'&&key[0]!='S')throw std::runtime_error("Typed reference has an unknown key type");
  const auto colon=key.rfind(':');if(colon<2||colon==key.size()-1||(key[0]=='C'&&colon==2))throw std::runtime_error("Typed reference has an invalid coordinate key");
  typed_int(key.substr(colon+1),key[0]=='C');
}
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
      file.seekg(0,std::ios::end);const auto end=file.tellg();
      if(end<0)throw std::runtime_error("Cannot determine FASTA byte bounds");
      const auto fasta_bytes=static_cast<std::uint64_t>(end);file.seekg(0);
      auto number=[](const std::string& value) {
        if(value.empty())throw std::runtime_error("Invalid FASTA index number");
        std::uint64_t result=0;
        for(unsigned char c:value) {
          if(c<'0'||c>'9'||result>(std::numeric_limits<std::uint64_t>::max()-(c-'0'))/10)
            throw std::runtime_error("Invalid FASTA index number");
          result=result*10+(c-'0');
        }
        return result;
      };
      while(std::getline(in,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        auto v=split(line,'\t');
        if(v.size()!=5||v[0].empty()||contigs.count(v[0])||
           std::any_of(v[0].begin(),v[0].end(),[](unsigned char c){return c==0||std::isspace(c);}))
          throw std::runtime_error("Invalid FASTA index name or row");
        Contig c{number(v[1]),number(v[2]),number(v[3]),number(v[4])};
        if(!c.size) {
          // Preserve the native writer's empty-contig representation. No
          // positive-length query can enter line arithmetic for this contig.
          if(c.bases||c.bytes||c.offset>fasta_bytes)throw std::runtime_error("Invalid FASTA index empty contig");
        } else {
          if(!c.bases||c.bases>c.size||c.bytes<c.bases||c.offset>=fasta_bytes)
            throw std::runtime_error("Invalid FASTA index layout");
          const auto zero=c.size-1,rows=zero/c.bases,tail=zero%c.bases;
          const auto available=fasta_bytes-1-c.offset;
          // Division before multiplication prevents overflow even for malicious
          // uint64 field values; the last addressed base must exist in this file.
          if(rows>available/c.bytes||tail>available-rows*c.bytes)
            throw std::runtime_error("Invalid FASTA index byte bounds");
        }
        contigs.emplace(v[0],c);
      }
      if(in.bad()||contigs.empty())throw std::runtime_error("Invalid FASTA index read or empty index");
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
    const auto& c=it->second;auto zero=static_cast<std::uint64_t>(pos-1);
    if(zero>c.size||n>c.size-zero)return "";
    std::string result; result.reserve(n);
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
inline void prepare(const std::string& fasta,const std::string& sites,const std::string& dir,std::size_t batch_size,int format_version=2) {
  if(format_version!=2&&format_version!=3)throw std::runtime_error("Unsupported reference writer format");
  if(std::filesystem::exists(dir))throw std::runtime_error("Reference directory already exists; use a new immutable version path");
  std::filesystem::create_directories(dir);Fasta fa(fasta);fa.save(dir+"/fasta.fai");
  std::ifstream in(sites);std::ofstream data(dir+"/records",std::ios::binary),rejected(dir+"/rejected.tsv");if(!in||!data||!rejected)throw std::runtime_error("Cannot open reference input/output");
  std::string line;std::getline(in,line);if(!line.empty()&&line.back()=='\r')line.pop_back();if(line!="chr\tpos\tref\talt\tid\taliases")throw std::runtime_error("Reference header must be chr,pos,ref,alt,id,aliases (tab-separated)");
  Sorter sorter(dir,batch_size);std::uint64_t row=1;rejected<<"row\tstatus\n";
  while(std::getline(in,line)) {
    ++row;if(!line.empty()&&line.back()=='\r')line.pop_back();auto v=split(line,'\t');if(v.size()!=6)throw std::runtime_error("Reference row must have six fields");std::int64_t pos;
    if(format_version==3)for(const auto& field:v)if(!typed_text(field))throw std::runtime_error("Typed reference source text is not valid UTF-8/framing");
    try{std::size_t used;pos=std::stoll(v[1],&used);if(used!=v[1].size())throw std::invalid_argument("position");}catch(...){rejected<<row<<"\tinvalid_position\n";continue;}
    // Site context exists independently of sequence-resolved candidates. Even
    // symbolic-only source rows are evidence against a biallelic complement.
    auto context_offset=static_cast<std::uint64_t>(data.tellp());
    if(format_version==2)data<<v[0]<<'\t'<<pos<<'\t'<<upper(v[2])<<"\t\t"<<v[4]<<'\t'<<pos<<'\t'<<upper(v[2])<<'\t'<<upper(v[3])<<'\n';
    else data<<"S\t"<<v[0]<<'\t'<<pos<<'\t'<<upper(v[2])<<'\t'<<upper(v[3])<<'\n';
    sorter.add("S:"+v[0]+":"+std::to_string(pos),context_offset);
    auto alleles=split(upper(v[3]),',');
    for(auto& alt:alleles) {
      auto n=normalize(fa,v[0],pos,v[2],alt);if(n.status!="ok"){rejected<<row<<'\t'<<n.status<<'\n';continue;}
      auto offset=static_cast<std::uint64_t>(data.tellp());
      if(format_version==2)data<<v[0]<<'\t'<<n.pos<<'\t'<<n.ref<<'\t'<<n.alt<<'\t'<<v[4]<<'\t'<<pos<<'\t'<<upper(v[2])<<'\t'<<upper(v[3])<<'\n';
      else data<<"A\t"<<v[0]<<'\t'<<n.pos<<'\t'<<n.ref<<'\t'<<n.alt<<'\t'<<v[4]<<'\t'<<pos<<'\n';
      sorter.add("C:"+v[0]+":"+std::to_string(pos),offset);sorter.add("C:"+v[0]+":"+std::to_string(n.pos),offset);
      if(!missing(v[4]))sorter.add("I:"+v[4],offset);
      for(auto& alias:split(v[5],';'))if(!missing(alias))sorter.add("I:"+alias,offset);
    }
  }
  if(!data||!rejected)throw std::runtime_error("Reference data write failed");data.close();sorter.finish();
  std::ofstream meta(dir+"/core-version",std::ios::binary);meta<<"mss-preservation-"<<format_version<<"\n";
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
inline void consume_context_work(std::size_t* remaining) {
  if(!remaining)return;
  if(!*remaining)throw std::runtime_error("Reference context closure limit exceeded; no arbitrary context selected");
  --*remaining;
}
class Reference {
  BlockFile keys,offsets,records;std::uint64_t count;Fasta fa;int data_format=2;
  // Node indices correspond to the actual midpoint tree of [0,count). Keep
  // validated complete keys only; an oversized pivot always uses disk lookup.
  struct Pivot {
    char key[56];
    unsigned char stored_size=0; // one plus length, including empty format-2 keys
  };
  struct Rank {
    std::uint64_t rank_plus_one=0;
    unsigned char size=0;
    char key[224];
  };
  static_assert(sizeof(Pivot)<=64,"Fence directory slot exceeds its budget");
  std::unique_ptr<Pivot[]> pivots;
  std::unique_ptr<Rank[]> ranks;
  std::size_t pivot_count=0,rank_count=0;
  static_assert(sizeof(Rank)*4096+sizeof(ranks)+sizeof(rank_count)<1024*1024,"Rank cache exceeds one MiB");
  struct BatchRank {std::uint64_t rank;std::uint32_t text;std::uint16_t size;};
  struct BatchTask {std::uint64_t low,high;std::uint32_t first,last;std::size_t node;};
  struct BatchDelete {void operator()(unsigned char* p)const noexcept {::operator delete(p);}};
  using BatchStorage=std::unique_ptr<unsigned char,BatchDelete>;
  BatchStorage batch_storage;
  BatchRank* batch_ranks=nullptr;
  std::size_t batch_count=0,batch_bytes=0,batch_peak=0;
  static constexpr std::size_t batch_reserve=65536,batch_max=32U*1024*1024;
  static constexpr std::size_t batch_tasks=64,batch_scratch=256;
  static_assert(sizeof(BatchTask)*batch_tasks+batch_scratch<batch_reserve,"Batch fixed workspace exceeds reserve");
  std::string_view batch_key(const BatchRank& entry)const {
    return {reinterpret_cast<const char*>(batch_storage.get()+entry.text),entry.size};
  }
#ifdef MSS_LOOKUP_COUNTERS
  std::uint64_t rank_hits=0,pivot_hits=0,probes=0;
#endif
  template<class Visitor>
  void visit_entry(std::uint64_t i,const Visitor& visit) {
    const auto pos=offsets.little_endian_offset(i*8);
    keys.visit_line(pos,data_format==3,[&](std::string_view row) {
      const auto tab=row.find('\t');
      const auto key=row.substr(0,tab);
      // Format 2 historically passes the complete row to stoull when tab is
      // absent. Preserve conversion order and exceptions for that format.
      const auto offset=row.substr(tab+1);
      std::uint64_t parsed=0;
      if(data_format==3) {
        if(tab==std::string_view::npos||offset.find('\t')!=std::string_view::npos)
          throw std::runtime_error("Typed reference key row must have two fields");
        typed_key(key);
        parsed=typed_uint(offset);
        if(parsed>=records.size())throw std::runtime_error("Typed reference record offset is out of range");
      }
      visit(key,offset,parsed);
    });
  }
  std::uint64_t lower_bound_rank(const std::string& key) {
    if(batch_count) {
      std::size_t low=0,high=batch_count;
      const std::string_view wanted(key);
      while(low<high) {
        const auto mid=low+(high-low)/2;
        if(batch_key(batch_ranks[mid])<wanted)low=mid+1;else high=mid;
      }
      if(low<batch_count&&batch_key(batch_ranks[low])==wanted) {
#ifdef MSS_LOOKUP_COUNTERS
        ++rank_hits;
#endif
        return batch_ranks[low].rank;
      }
    }
    if(!count)return 0;
    Rank* rank=nullptr;
    if(key.size()<=sizeof(Rank::key)) {
      rank=&ranks[std::hash<std::string>{}(key)%rank_count];
      if(rank->rank_plus_one&&std::string_view(rank->key,rank->size)==std::string_view(key)) {
#ifdef MSS_LOOKUP_COUNTERS
        ++rank_hits;
#endif
        return rank->rank_plus_one-1;
      }
    }
    std::size_t node=0;
    std::uint64_t low=0,high=count;
    while(low<high) {
      const auto mid=low+(high-low)/2;
      auto* pivot=node<pivot_count?&pivots[node]:nullptr;
      bool less=false;
      if(pivot&&pivot->stored_size) {
        less=std::string_view(pivot->key,pivot->stored_size-1)<std::string_view(key);
#ifdef MSS_LOOKUP_COUNTERS
        ++pivot_hits;
#endif
      } else {
#ifdef MSS_LOOKUP_COUNTERS
        ++probes;
#endif
        visit_entry(mid,[&](std::string_view field,std::string_view,std::uint64_t) {
          less=field < std::string_view(key);
          if(pivot&&field.size()<=sizeof(pivot->key)) {
            std::copy(field.begin(),field.end(),pivot->key);
            pivot->stored_size=static_cast<unsigned char>(field.size()+1);
          }
        });
      }
      if(less)low=mid+1;else high=mid;
      if(pivot)node=node*2+(less?2:1);
    }
    if(rank) {
      std::copy(key.begin(),key.end(),rank->key);
      rank->size=static_cast<unsigned char>(key.size());
      rank->rank_plus_one=low+1;
    }
    return low;
  }
  // References have at most eight fields. Count and (for format 3) validate
  // the entire row even when an extra tail will cause rejection.
  struct RecordFields {
    std::string_view values[8];
    std::size_t size=0;
    RecordFields(std::string_view row,bool typed) {
      if(row.empty())return;
      std::size_t start=0;
      while(true) {
        const auto end=row.find('\t',start);
        const auto field=row.substr(start,end-start);
        if(typed&&!typed_text(field))throw std::runtime_error("Typed reference field is not valid UTF-8/framing");
        if(size<8)values[size]=field;
        ++size;
        if(end==std::string_view::npos)break;
        start=end+1;
      }
    }
  };
  Candidate decode_candidate(const std::string& key,std::string_view row) {
    const RecordFields fields(row,data_format==3);
    const auto& values=fields.values;
    Candidate out;
    if(data_format==2) {
      if(fields.size!=8)throw std::runtime_error("Corrupt reference record");
      out.chr=values[0];out.ref=values[2];out.alt=values[3];out.id=values[4];
      out.original_ref=values[6];out.alleles=values[7];
      out.pos=std::stoll(std::string(values[1]));out.original_pos=std::stoll(std::string(values[5]));
      if(!missing(out.id))out.identifiers.insert(out.id);
    } else if(fields.size==7&&values[0]=="A") {
      if(key[0]!='C'&&key[0]!='I')throw std::runtime_error("Typed reference A record has the wrong key type");
      out.chr=values[1];out.pos=typed_int(values[2],true);out.ref=values[3];out.alt=values[4];out.id=values[5];out.original_pos=typed_int(values[6],true);
      if(out.chr.empty()||!dna(out.ref)||!dna(out.alt)||out.ref==out.alt)throw std::runtime_error("Typed reference A record has invalid canonical fields");
      if(key[0]=='C'&&key!="C:"+out.chr+":"+std::to_string(out.pos)&&key!="C:"+out.chr+":"+std::to_string(out.original_pos))throw std::runtime_error("Typed reference A coordinate key disagrees with its record");
      if(!missing(out.id))out.identifiers.insert(out.id);
    } else if(fields.size==5&&values[0]=="S") {
      if(key[0]!='S')throw std::runtime_error("Typed reference S record has the wrong key type");
      out.chr=values[1];out.pos=out.original_pos=typed_int(values[2],false);out.ref=out.original_ref=values[3];out.alleles=values[4];
      if(key!="S:"+out.chr+":"+std::to_string(out.original_pos))throw std::runtime_error("Typed reference S key disagrees with its record");
    } else throw std::runtime_error("Typed reference record has an unknown tag or field count");
    out.original_positions.insert(out.original_pos);return out;
  }
public:
  void clear_prefetched_ranks()noexcept {
    batch_count=0;batch_ranks=nullptr;batch_bytes=0;batch_storage.reset();
  }
  std::size_t prefetched_rank_count()const noexcept{return batch_count;}
  // Conservative incremental bound, including fixed state and stack reserve.
  std::size_t prefetch_owned_bytes()const noexcept{return batch_storage?batch_bytes+batch_reserve:0;}
  std::size_t prefetch_peak_bytes()const noexcept{return batch_peak;}
  bool prefetch_rank_keys(const std::vector<std::string>& requested,std::size_t budget) {
    clear_prefetched_ranks();batch_peak=0;
    if(budget<=batch_reserve||budget>batch_max)return false;
    const auto available=budget-batch_reserve;
    // Task/scratch and all entries/owned key bytes live in one allocation.
    constexpr auto entry_start=(sizeof(BatchTask)*batch_tasks+batch_scratch+
                               alignof(BatchRank)-1)/alignof(BatchRank)*alignof(BatchRank);
    if(available<entry_start)return false;
    std::size_t selected=0,needed=entry_start;
    try {
      for(const auto& key:requested) {
        if(key.size()<2||key[1]!=':'||(key[0]!='C'&&key[0]!='I')||key.size()>224)continue;
        if(data_format==3)typed_key(key);
        if(sizeof(BatchRank)>available-needed||key.size()>available-needed-sizeof(BatchRank))return false;
        needed+=sizeof(BatchRank)+key.size();++selected;
      }
      if(!selected)return false;
      BatchStorage storage(static_cast<unsigned char*>(::operator new(needed)));
      batch_peak=needed+batch_reserve;
      auto* tasks=reinterpret_cast<BatchTask*>(storage.get());
      for(std::size_t i=0;i<batch_tasks;++i)::new(static_cast<void*>(tasks+i)) BatchTask{};
      auto* scratch=reinterpret_cast<char*>(storage.get()+sizeof(BatchTask)*batch_tasks);
      auto* entries=reinterpret_cast<BatchRank*>(storage.get()+entry_start);
      std::size_t text=entry_start+selected*sizeof(BatchRank),at=0;
      for(const auto& key:requested) {
        if(key.size()<2||key[1]!=':'||(key[0]!='C'&&key[0]!='I')||key.size()>224)continue;
        ::new(static_cast<void*>(entries+at++)) BatchRank{0,static_cast<std::uint32_t>(text),static_cast<std::uint16_t>(key.size())};
        std::copy(key.begin(),key.end(),reinterpret_cast<char*>(storage.get()+text));text+=key.size();
      }
      auto view=[&](const BatchRank& entry) {
        return std::string_view(reinterpret_cast<const char*>(storage.get()+entry.text),entry.size);
      };
      // Iterative heapsort: no implementation-dependent recursive sort stack.
      auto sift=[&](std::size_t root,std::size_t size) {
        while(root<size/2) {
          auto child=root*2+1;
          if(child+1<size&&view(entries[child])<view(entries[child+1]))++child;
          if(!(view(entries[root])<view(entries[child])))break;
          std::swap(entries[root],entries[child]);root=child;
        }
      };
      for(auto start=selected/2;start;--start)sift(start-1,selected);
      for(auto end=selected;end>1;--end){std::swap(entries[0],entries[end-1]);sift(0,end-1);}
      std::size_t unique=0;
      for(std::size_t i=0;i<selected;++i)
        if(!unique||view(entries[i])!=view(entries[unique-1]))entries[unique++]=entries[i];
      std::size_t pending=1;
      tasks[0]={0,count,0,static_cast<std::uint32_t>(unique),0};
      while(pending) {
        const auto task=tasks[--pending];
        if(task.first==task.last)continue;
        if(task.low==task.high) {
          for(auto i=task.first;i<task.last;++i)entries[i].rank=task.low;
          continue;
        }
        const auto mid=task.low+(task.high-task.low)/2;
        const auto* pivot=task.node<pivot_count?&pivots[task.node]:nullptr;
        std::string_view field;
        if(pivot&&pivot->stored_size) {
#ifdef MSS_LOOKUP_COUNTERS
          ++pivot_hits;
#endif
          field={pivot->key,static_cast<std::size_t>(pivot->stored_size-1)};
        } else {
#ifdef MSS_LOOKUP_COUNTERS
          ++probes;
#endif
          std::size_t size=0;
          if(!keys.bounded_line(offsets.little_endian_offset(mid*8),data_format==3,
                                scratch,batch_scratch,size))return false;
          const std::string_view row(scratch,size);
          const auto tab=row.find('\t');field=row.substr(0,tab);
          if(data_format==3) {
            const auto offset=row.substr(tab+1);
            if(tab==std::string_view::npos||offset.find('\t')!=std::string_view::npos)
              throw std::runtime_error("Typed reference key row must have two fields");
            typed_key(field);
            if(typed_uint(offset)>=records.size())throw std::runtime_error("Typed reference record offset is out of range");
          }
        }
        // Queries <= midpoint follow left; greater queries follow right.
        auto low=task.first,high=task.last;
        while(low<high) {
          const auto q=low+(high-low)/2;
          if(field<view(entries[q]))high=q;else low=q+1;
        }
        const auto split=low;
        if(pending+2>batch_tasks)return false;
        const auto left=pivot?task.node*2+1:pivot_count;
        const auto right=pivot?task.node*2+2:pivot_count;
        if(split<task.last)tasks[pending++]={mid+1,task.high,split,task.last,right};
        if(task.first<split)tasks[pending++]={task.low,mid,task.first,split,left};
      }
      // Publish only complete results. Original pivots/ranks were never written.
      batch_storage=std::move(storage);batch_ranks=entries;batch_bytes=needed;batch_count=unique;
      return true;
    } catch(const std::bad_alloc&) {return false;}
      catch(const std::runtime_error&) {return false;}
  }
  Reference(const std::string& dir,const std::string& fasta):keys(dir+"/keys",16384,512),offsets(dir+"/offsets"),records(dir+"/records"),fa(fasta,dir+"/fasta.fai") {
    std::ifstream meta(dir+"/core-version");std::string version;std::getline(meta,version);
    if(version!="mss-preservation-2"&&version!="mss-preservation-3")throw std::runtime_error("Incomplete or unsupported reference index");
    if(version=="mss-preservation-3") {
      std::ifstream exact(dir+"/core-version",std::ios::binary);const std::string marker((std::istreambuf_iterator<char>(exact)),std::istreambuf_iterator<char>());
      if(marker!="mss-preservation-3\n")throw std::runtime_error("Typed reference marker is malformed");data_format=3;
    }
    auto size=offsets.size();if(size%8)throw std::runtime_error("Invalid offsets file");count=size/8;
    pivot_count=static_cast<std::size_t>(std::min<std::uint64_t>(8191,count));
    rank_count=static_cast<std::size_t>(std::min<std::uint64_t>(4096,count));
    if(count){pivots=std::make_unique<Pivot[]>(pivot_count);ranks=std::make_unique<Rank[]>(rank_count);}
  }
#ifdef MSS_LOOKUP_COUNTERS
  struct LookupCounters {
    std::uint64_t rank_hits,pivot_hits,probe_entries,key_requests,offset_requests,key_reads,offset_reads;
    std::size_t pivot_slots,rank_slots,rank_bytes;
  };
  LookupCounters lookup_counters()const {
    return {rank_hits,pivot_hits,probes,keys.page_requests(),offsets.page_requests(),
            keys.disk_reads(),offsets.disk_reads(),pivot_count,rank_count,
            sizeof(ranks)+sizeof(rank_count)+rank_count*sizeof(Rank)};
  }
#endif
  std::vector<Candidate> lookup(const std::string& key,std::size_t cap=10000,std::size_t* context_work=nullptr) {
    if(data_format==3)typed_key(key);
    const auto low=lower_bound_rank(key);
    std::vector<Candidate> out;std::map<std::string,std::size_t> seen;std::size_t examined=0;
    for(auto i=low;i<count;++i) {
      bool matched=false;
      std::uint64_t offset=0;
      visit_entry(i,[&](std::string_view field,std::string_view offset_text,std::uint64_t parsed) {
        if(field!=std::string_view(key))return;
        matched=true;
        if(++examined>cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");
        consume_context_work(context_work);
        offset=data_format==3?parsed:std::stoull(std::string(offset_text));
      });
      if(!matched)break;
      records.visit_line(offset,data_format==3,[&](std::string_view row) {
        auto c=decode_candidate(key,row);
        if(out.empty()){out.push_back(std::move(c));return;}
        if(seen.empty())seen.emplace(out.front().key(),0);
        const auto inserted=seen.emplace(c.key(),out.size());
        if(inserted.second)out.push_back(std::move(c));
        else merge_candidate_context(out[inserted.first->second],c);
      });
    }
    // S lookup has already merged identical (chr, original position, REF) spans.
    // These leading fields uniquely order contexts; ALT token occurrence/order
    // remains unchanged for the scientific work budget, including duplicates.
    if(data_format==3&&key[0]=='S')std::sort(out.begin(),out.end(),[](const Candidate& a,const Candidate& b){
      return std::tie(a.chr,a.original_pos,a.original_ref,a.alleles)<std::tie(b.chr,b.original_pos,b.original_ref,b.alleles);
    });
    return out;
  }
  Result harmonise(const Input& input,std::size_t cap=10000) {
    if(input.source_contract)return harmonise_declared(input,cap);
    SourceNormalizations source(input);
    Result result;std::vector<Candidate> coordinate,ids,candidates;
    std::set<std::string> loaded_coordinates;
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
        source.retain(orientation,n);
      }
      for(auto& pos:positions) {
        const auto key="C:"+input.chr+":"+pos;loaded_coordinates.insert(key);
        for(auto& candidate:lookup(key,cap)) {
          auto found=seen.find(candidate.key());
          if(found==seen.end()){if(coordinate.size()>=cap)throw std::runtime_error("Reference candidate limit exceeded; no arbitrary candidate selected");seen[candidate.key()]=coordinate.size();coordinate.push_back(candidate);}
          else merge_candidate_context(coordinate[found->second],candidate);
        }
      }
    }
    if(!missing(input.id))ids=lookup("I:"+input.id,cap);
    if(!coordinate.empty()&&!ids.empty()) {std::map<std::string,std::size_t> wanted;for(std::size_t i=0;i<coordinate.size();++i)wanted[coordinate[i].key()]=i;for(auto& c:ids){auto found=wanted.find(c.key());if(found!=wanted.end()){merge_candidate_context(c,coordinate[found->second],false);candidates.push_back(c);}}if(candidates.empty()){result.status="coordinate_id_conflict";return result;}}
    else if(!ids.empty()) {if(!missing(input.chr)&&input.pos>0){result.status="coordinate_id_conflict";return result;}candidates=ids;}
    else candidates=coordinate;
    result.candidate_count=candidates.size();
    auto other=upper(input.other),effect=upper(input.effect_allele);
    // An identifier does not establish a source anchor for complex alleles.
    // Even one canonical event may admit opposite valid effect directions.
    if(input.pos<=0&&dna(other)&&dna(effect)&&other!=effect&&
       (other.size()!=1||effect.size()!=1)) {
      result.status="source_coordinates_missing";return result;
    }
    if(candidates.empty())return result;
    if(ids.empty()&&indel_events.size()>1){result.status="ambiguous";return result;}
    if(!dna(other)||!dna(effect)||other==effect){result.status="unsupported_allele";return result;}
    struct Hit{Candidate candidate;bool swap,comp;std::int64_t anchor;};std::vector<Hit> hits;std::set<std::string> hitkeys;
    for(auto& c:candidates)for(int orientation=0;orientation<4;++orientation) {
      bool swap=orientation%2,comp=orientation>=2;auto a=comp?complement(other):other;auto b=comp?complement(effect):effect;if(swap)std::swap(a,b);
      auto chr=missing(input.chr)?c.chr:input.chr;
      std::set<std::int64_t> anchors=input.pos>0?std::set<std::int64_t>{input.pos}:c.original_positions;
      if(input.pos<=0)anchors.insert(c.pos);
      for(auto pos:anchors) {
        const bool matches=source.compare(fa,chr,pos,orientation,a,b,[&](const NormalizationView& n) {
          return n.status=="ok"&&n.chr==c.chr&&n.pos==c.pos&&n.ref==c.ref&&n.alt==c.alt;
        });
        if(matches){auto key=c.key()+":"+(swap?"1":"0")+":"+(comp?"1":"0");if(hitkeys.insert(key).second)hits.push_back({c,swap,comp,pos});}
      }

    }
    if(hits.size()>1){result.status="ambiguous";return result;}
    // A unique REF-ALT hit does not prove the source strand: the opposite
    // strand may instead measure two alternate alleles at the same site.
    // Neither a site identifier nor EAF selects between these explanations.
    std::set<std::int64_t> resolved_anchors;
    if(hits.size()==1){resolved_anchors.insert(hits[0].anchor);resolved_anchors.insert(hits[0].candidate.pos);}
    const auto selected=hits.empty()?std::string():hits[0].candidate.key()+":"+(hits[0].swap?"1":"0")+":"+(hits[0].comp?"1":"0");
    bool competing_ref_alt=false,context_required=false;
    const bool alt_alt=alt_alt_context(input,coordinate,other,effect,source,cap,resolved_anchors,
                                      &ids,&loaded_coordinates,selected,&competing_ref_alt,&context_required);
    if(hits.empty()) {
      if(alt_alt){result.status="alt_alt_contrast";return result;}
      result.status="allele_mismatch";return result;
    }
    if(alt_alt||competing_ref_alt){result.status="ambiguous";return result;}
    if(context_required){result.status="source_allele_context_required";return result;}
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
  struct NormalizationView {
    std::string_view status;
    const std::string& chr;
    std::int64_t pos;
    std::string_view ref,alt;
  };
  // Only this row's coordinate loop writes slots. All three loops derive
  // orientation i from the same immutable other/effect alleles in the same order.
  class SourceNormalizations {
    struct Slot {
      std::int64_t pos;
      std::uint16_t status,ref,alt,end;
    };
    const Input& source;
    char bytes[256];
    Slot slots[4];
    std::size_t used=0;
    unsigned present=0;

    std::uint16_t append(std::string_view s) {
      const auto start=static_cast<std::uint16_t>(used);
      std::copy(s.begin(),s.end(),bytes+used);
      used+=s.size();
      return start;
    }
    std::string_view slice(std::uint16_t begin,std::uint16_t end) const {
      return {bytes+begin,static_cast<std::size_t>(end-begin)};
    }
  public:
    explicit SourceNormalizations(const Input& row):source(row) {}
    void retain(int orientation,const Normalized& n) {
      auto room=sizeof(bytes)-used;
      for(std::string_view s:{std::string_view(n.status),std::string_view(n.ref),std::string_view(n.alt)}) {
        if(s.size()>room)return;
        room-=s.size();
      }
      const auto status=append(n.status);
      const auto ref=append(n.ref);
      const auto alt=append(n.alt);
      slots[orientation]={n.pos,status,ref,alt,static_cast<std::uint16_t>(used)};
      present|=1U<<orientation;
    }
    template<class Compare>
    bool compare(Fasta& fa,const std::string& chr,std::int64_t pos,int orientation,
                 const std::string& ref,const std::string& alt,const Compare& compare_value) const {
      if((present&(1U<<orientation))&&pos==source.pos&&chr==source.chr) {
        const auto& s=slots[orientation];
        return compare_value(NormalizationView{slice(s.status,s.ref),source.chr,s.pos,
                             slice(s.ref,s.alt),slice(s.alt,s.end)});
      }
      const auto n=normalize(fa,chr,pos,ref,alt);
      return compare_value(NormalizationView{n.status,n.chr,n.pos,n.ref,n.alt});
    }
  };
  static_assert(sizeof(SourceNormalizations)<=576,"Row normalization storage exceeds its bound");
  struct ProofEvent {std::string chr,ref,alt;std::int64_t pos;};
  struct ContextClosure {
    std::set<std::pair<std::string,std::int64_t>> locations;
    std::map<std::string,ProofEvent> events;
    std::size_t remaining_work;
  };
  ContextClosure context_closure(const Input& input,const std::vector<Candidate>& coordinate,
                                 const std::vector<Candidate>* ids,std::size_t cap,
                                 const std::set<std::string>* loaded_coordinates) {
    if(cap>std::numeric_limits<std::size_t>::max()/16)
      throw std::runtime_error("Reference context closure limit is not representable");
    ContextClosure out;out.remaining_work=16*cap;
    std::set<std::string> current_ids,pending_ids,expanded_ids,pending_events,bridged_events;
    if(!missing(input.id))expanded_ids.insert(input.id); // Its I lookup is already complete.
    auto limit=[&](std::size_t size){if(size>cap)throw std::runtime_error("Reference context closure limit exceeded; no arbitrary context selected");};
    auto consider=[&](const Candidate& c) {
      // Never expand unrelated chromosomes through a shared identifier.
      if(!missing(input.chr)&&input.chr!=c.chr)return;
      for(auto position:c.original_positions) {
        out.locations.emplace(c.chr,position);limit(out.locations.size());
      }
      for(const auto& id:c.identifiers)if(current_ids.insert(id).second) {
        limit(current_ids.size());pending_ids.insert(id);
      }
      const auto key=c.key();
      if(out.events.emplace(key,ProofEvent{c.chr,c.ref,c.alt,c.pos}).second) {
        limit(out.events.size());pending_events.insert(key);
      }
    };
    for(const auto* pool:{&coordinate,ids})if(pool)for(const auto& c:*pool)consider(c);
    if(loaded_coordinates)for(const auto& entry:out.events) {
      const auto& event=entry.second;
      if(loaded_coordinates->count("C:"+event.chr+":"+std::to_string(event.pos)))bridged_events.insert(entry.first);
    }
    // Fixed-point current-ID closure. Only exact canonical events bridge from
    // C lookups; unrelated overlapping haplotypes are outside this universe.
    while(!pending_ids.empty()||!pending_events.empty()) {
      if(!pending_ids.empty()) {
        const auto id=*pending_ids.begin();pending_ids.erase(pending_ids.begin());
        if(!expanded_ids.insert(id).second)continue;
        const auto batch=lookup("I:"+id,cap,&out.remaining_work);
        for(const auto& c:batch)consider(c);
      } else {
        const auto key=*pending_events.begin();pending_events.erase(pending_events.begin());
        if(bridged_events.count(key))continue;
        const auto event=out.events.at(key);
        const auto batch=lookup("C:"+event.chr+":"+std::to_string(event.pos),cap,&out.remaining_work);
        for(const auto& c:batch) {
          const auto found=out.events.find(c.key());
          if(c.chr!=event.chr||c.pos!=event.pos||found==out.events.end())continue;
          bridged_events.insert(found->first);consider(c);
        }
        if(!bridged_events.count(key))throw std::runtime_error("Reference context closure lacks its canonical event");
      }
    }
    return out;
  }
  bool alt_alt_context(const Input& input,const std::vector<Candidate>& candidates,
                       const std::string& other,const std::string& effect,
                       const SourceNormalizations& source,std::size_t cap,const std::set<std::int64_t>& resolved_anchors,
                       const std::vector<Candidate>* ids=nullptr,
                       const std::set<std::string>* loaded_coordinates=nullptr,
                       const std::string& selected=std::string(),bool* competing_ref_alt=nullptr,
                       bool* context_required=nullptr) {
    auto closure=context_closure(input,candidates,ids,cap,loaded_coordinates);
    // Proof-only mappings may reject a selected mapping, but never replace it.
    // Complex missing-position rows have already been withheld. Single bases
    // cannot trim or shift, so canonical positions exhaust their possible anchors.
    if(!selected.empty()&&competing_ref_alt) {
      std::set<std::pair<std::string,std::int64_t>> anchors;
      for(const auto& entry:closure.events)anchors.emplace(entry.second.chr,input.pos>0?input.pos:entry.second.pos);
      for(const auto& anchor:anchors)for(int orientation=0;orientation<4;++orientation) {
        consume_context_work(&closure.remaining_work);
        const bool swap=orientation%2,comp=orientation>=2;
        auto a=comp?complement(other):other,b=comp?complement(effect):effect;if(swap)std::swap(a,b);
        const bool competes=source.compare(fa,anchor.first,anchor.second,orientation,a,b,[&](const NormalizationView& n) {
          if(n.status!="ok")return false;
          auto key=n.chr+":"+std::to_string(n.pos)+":";
          key.append(n.ref).append(":").append(n.alt);
          return closure.events.count(key)&&key+":"+(swap?"1":"0")+":"+(comp?"1":"0")!=selected;
        });
        if(competes) {
          *competing_ref_alt=true;return false;
        }
      }
    }
    std::size_t contexts=0;
    for(const auto& position:closure.locations) {
      for(const auto& site:lookup("S:"+position.first+":"+std::to_string(position.second),cap,&closure.remaining_work)) {
          if(++contexts>cap)throw std::runtime_error("Reference context closure limit exceeded; no arbitrary context selected");
          if(site.original_pos!=position.second||site.chr!=position.first||!dna(site.original_ref)||
             fa.sequence(site.chr,site.original_pos,site.original_ref.size())!=site.original_ref)continue;
          // Own all strings until the row-scoped signature maps are destroyed.
          // Views avoid copying a long allele once for every possible trim.
          const auto alleles=split(site.alleles,',');
          std::vector<std::size_t> valid;std::set<std::string_view> distinct;
          bool unmodelled=false;
          for(std::size_t i=0;i<alleles.size();++i) {
            consume_context_work(&closure.remaining_work);
            if(alleles[i]==site.original_ref||!distinct.insert(alleles[i]).second)continue;
            if(dna(alleles[i]))valid.push_back(i);else unmodelled=true;
          }
          if(context_required&&!selected.empty()&&unmodelled&&distinct.size()>=2)*context_required=true;
          if(valid.size()<2)continue;
          // A jointly padded pure SNP context has one proved genomic anchor.
          // Every shared flank must agree with original REF, not just the ALTs.
          if(other.size()==1&&effect.size()==1&&input.pos>0&&!unmodelled) {
            const auto width=site.original_ref.size();std::size_t variable=width;bool pure=true;
            for(auto i:valid) {
              consume_context_work(&closure.remaining_work);
              if(alleles[i].size()!=width){pure=false;break;}
              for(std::size_t j=0;j<width;++j)if(alleles[i][j]!=site.original_ref[j]) {
                if(variable==width)variable=j;else if(variable!=j){pure=false;break;}
              }
              if(!pure)break;
            }
            if(pure&&variable!=width&&input.pos!=site.original_pos+static_cast<std::int64_t>(variable))continue;
          }
          // Necessary conditions are not full haplotype proofs. Shared alternate
          // backgrounds can change both coordinate offset and REF-span length.
          // Defer compatible complex contrasts when complete context is absent.
          if(context_required&&!selected.empty()&&(other.size()!=1||effect.size()!=1)) {
            const auto delta=other.size()>effect.size()?other.size()-effect.size():effect.size()-other.size();
            std::set<std::size_t> lengths;
            for(auto i:valid) {
              consume_context_work(&closure.remaining_work);const auto length=alleles[i].size();
              if((length>=delta&&lengths.count(length-delta))||
                 (length<=std::numeric_limits<std::size_t>::max()-delta&&lengths.count(length+delta))) {
                *context_required=true;break;
              }
              lengths.insert(length);
            }
          }
          // For single bases, a necessary ALT pair must differ at exactly one
          // base with the observed pair (or its complement). Join shared flanks
          // without inventing a genomic anchor or a positive inferred REF span.
          if(context_required&&!selected.empty()&&other.size()==1&&effect.size()==1) {
            for(int orientation=0;orientation<2&&!*context_required;++orientation) {
              const std::string observed[2]={orientation?complement(other):other,
                                             orientation?complement(effect):effect};
              using Flanks=std::tuple<std::size_t,std::size_t,std::string_view,std::string_view>;
              std::map<Flanks,unsigned> pairs;
              for(auto i:valid) {
                const auto& allele=alleles[i];const auto view=std::string_view(allele);
                for(unsigned which=0;which<2;++which) {
                  auto offset=allele.find(observed[which]);
                  while(offset!=std::string::npos) {
                    consume_context_work(&closure.remaining_work);
                    auto& mask=pairs[{allele.size(),offset,view.substr(0,offset),view.substr(offset+1)}];
                    mask|=1U<<which;
                    if(mask==3){*context_required=true;break;}
                    offset=allele.find(observed[which],offset+1);
                  }
                  if(*context_required)break;
                }
                if(*context_required)break;
              }
            }
          }
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
                consume_context_work(&closure.remaining_work);
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
                consume_context_work(&closure.remaining_work);
                const auto event=normalize(fa,site.chr,site.original_pos,site.original_ref,allele);
                if(event.status!="ok")continue;
                for(unsigned which=0;which<2;++which) {
                  const auto length=static_cast<std::int64_t>(site.original_ref.size())+
                                    static_cast<std::int64_t>(observed[which].size())-static_cast<std::int64_t>(allele.size());
                  if(length<1)continue;
                  auto found=query_events.find(length);
                  if(found==query_events.end()) {
                    const auto ref=fa.sequence(site.chr,anchor,static_cast<std::size_t>(length));
                    consume_context_work(&closure.remaining_work);consume_context_work(&closure.remaining_work);
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
