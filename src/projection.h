#ifndef MSS_PROJECTION_H
#define MSS_PROJECTION_H
#include "preservation_core.h"
namespace mss {
struct ChainBlock {std::string chr,id;std::int64_t start,end,qstart,qsize;bool reverse;std::uint64_t group;};
struct ChainIndex {std::vector<ChainBlock> blocks;std::vector<std::int64_t> maxend;};
class Projection {
  std::map<std::string,ChainIndex> index, spans;
public:
  Projection(const std::string& path,const std::map<std::string,std::string>& source_names,
             const std::map<std::string,std::string>& target_names,
             const std::map<std::string,std::int64_t>& source_sizes,const std::map<std::string,std::int64_t>& target_sizes) {
    std::ifstream input(path);if(!input)throw std::runtime_error("Cannot open chain");
    std::string line,tname,qname,strand,id;std::int64_t t=0,q=0,tsize=0,qsize=0,tend=0,qend=0;bool active=false;
    std::set<std::pair<std::string,std::string>> ids;std::uint64_t group=0;
    while(std::getline(input,line)) {
      if(!line.empty()&&line.back()=='\r')line.pop_back();
      if(line.empty()||line[0]=='#')continue;
      std::istringstream fields(line);std::string tag;fields>>tag;
      if(tag=="chain") {
        if(active)throw std::runtime_error("Chain missing terminal block");
        std::string score,ts,extra;
        if(!(fields>>score>>tname>>tsize>>ts>>t>>tend>>qname>>qsize>>strand>>q>>qend>>id)||(fields>>extra))throw std::runtime_error("Malformed chain header");
        if(ts!="+"||(strand!="+"&&strand!="-")||t<0||q<0||tend<=t||qend<=q||tend>tsize||qend>qsize||!ids.insert({tname,id}).second)throw std::runtime_error("Invalid chain coordinates, strands or duplicate ID");
        if(source_names.count(tname))tname=source_names.at(tname);
        if(target_names.count(qname))qname=target_names.at(qname);
        if(source_sizes.count(tname)&&source_sizes.at(tname)!=tsize)throw std::runtime_error("Source chain/reference contig size mismatch");
        if(target_sizes.count(qname)&&target_sizes.at(qname)!=qsize)throw std::runtime_error("Target chain/reference contig size mismatch");
        ++group;spans[tname].blocks.push_back({qname,id,t,tend,q,qsize,strand=="-",group});
        active=true;
      } else {
        if(!active)throw std::runtime_error("Chain block without header");
        fields.clear();fields.str(line);std::vector<std::int64_t> v;std::int64_t x;
        while(fields>>x)v.push_back(x);
        if(!fields.eof()||(v.size()!=1&&v.size()!=3)||v[0]<=0||(v.size()==3&&(v[1]<0||v[2]<0))||v[0]>tend-t||v[0]>qend-q)throw std::runtime_error("Malformed chain block");
        index[tname].blocks.push_back({qname,id,t,t+v[0],q,qsize,strand=="-",group});
        t+=v[0];q+=v[0];
        if(v.size()==1){if(t!=tend||q!=qend)throw std::runtime_error("Chain endpoints disagree with blocks");active=false;}
        else {t+=v[1];q+=v[2];if(t>=tend||q>=qend)throw std::runtime_error("Chain gap outside endpoints");}
      }
    }
    if(active||index.empty())throw std::runtime_error("Truncated or empty chain");
    for(auto* collection:{&index,&spans})for(auto& entry:*collection){auto& a=entry.second;std::sort(a.blocks.begin(),a.blocks.end(),[](const ChainBlock& x,const ChainBlock& y){return x.start<y.start;});std::int64_t end=0;for(auto& b:a.blocks){end=std::max(end,b.end);a.maxend.push_back(end);}}
  }
  struct Mapping {std::string status="unmapped",chr,id;std::int64_t pos=0;bool reverse=false;};
  Mapping map(const std::string& chr,std::int64_t pos,std::size_t width) const {
    Mapping out;auto entry=index.find(chr);if(entry==index.end())return out;
    const auto& a=entry->second;auto start=pos-1,end=start+static_cast<std::int64_t>(width);
    auto hi=std::lower_bound(a.blocks.begin(),a.blocks.end(),end,[](const ChainBlock& b,std::int64_t e){return b.start<e;})-a.blocks.begin();
    std::vector<const ChainBlock*> complete;std::set<std::uint64_t> overlapping;
    while(hi>0&&a.maxend[hi-1]>start){const auto& b=a.blocks[--hi];if(b.end<=start)continue;overlapping.insert(b.group);if(b.start<=start&&b.end>=end)complete.push_back(&b);}
    if(overlapping.size()>1||complete.size()>1){out.status="multimap";return out;}
    if(complete.empty()){
      const auto& ranges=spans.at(chr);auto j=std::lower_bound(ranges.blocks.begin(),ranges.blocks.end(),end,[](const ChainBlock& b,std::int64_t e){return b.start<e;})-ranges.blocks.begin();
      while(j>0&&ranges.maxend[j-1]>start){const auto& b=ranges.blocks[--j];if(b.end>start)overlapping.insert(b.group);}
      if(!overlapping.empty())out.status=overlapping.size()>1?"multimap":"chain_gap";return out;
    }
    const auto& b=*complete[0];out.status="mapped";out.chr=b.chr;out.id=b.id;out.reverse=b.reverse;
    out.pos=b.reverse?b.qsize-(b.qstart+(start-b.start)+width)+1:b.qstart+(start-b.start)+1;
    return out;
  }
};
}
#endif
