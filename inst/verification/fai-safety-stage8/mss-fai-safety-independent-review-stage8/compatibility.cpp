#include "preservation_core.h"
#include <iostream>
#include <limits>

static void write_file(const std::string& path,const std::string& value) {
  std::ofstream out(path,std::ios::binary);out.write(value.data(),value.size());
  if(!out)throw std::runtime_error("review compatibility fixture write failed");
}
int main(int argc,char**argv) {
  if(argc!=2)return 2;
  const std::string root=argv[1];std::filesystem::create_directories(root);
  const auto fasta=root+"/sequence.fa",index=root+"/sequence.fai";
  std::size_t controls=0,queries=0,failures=0;
  for(std::size_t length:{std::size_t(0),std::size_t(1),std::size_t(13),std::size_t(65540)})
    for(std::size_t width:{std::size_t(1),std::size_t(2),std::size_t(7),std::size_t(32767),std::size_t(65536)})
      for(const std::string term:{std::string("\n"),std::string("\r\n"),std::string("\v"),std::string("\x1e\x1f")})
        for(bool terminated:{false,true}) {
          std::string seq;seq.reserve(length);
          for(std::size_t i=0;i<length;++i)seq.push_back("ACGT"[i%4]);
          const auto bases=std::min(width,length);
          const std::string header=">雪 description\n";
          std::string data=header;
          for(std::size_t i=0;i<length;i+=bases){auto n=std::min(bases,length-i);data+=seq.substr(i,n);if(i+n<length||terminated)data+=term;}
          write_file(fasta,data);
          const auto offset=length?header.size():0;
          const auto stride=length?bases+term.size():0;
          write_file(index,"雪\t"+std::to_string(length)+"\t"+std::to_string(offset)+"\t"+std::to_string(bases)+"\t"+std::to_string(stride)+(terminated?"\r\n":""));
          mss::Fasta reader(fasta,index);
          for(std::int64_t pos:{-1LL,0LL,1LL,2LL,13LL,65535LL,65536LL,65540LL,65541LL,std::numeric_limits<std::int64_t>::max()})
            for(std::size_t n:{std::size_t(0),std::size_t(1),std::size_t(2),std::size_t(7),std::size_t(65536),std::numeric_limits<std::size_t>::max()}) {
              ++queries;
              const bool within=pos>0&&static_cast<__uint128_t>(pos-1)+n<=length;
              const auto expected=within?seq.substr(static_cast<std::size_t>(pos-1),n):"";
              if(reader.sequence("雪",pos,n)!=expected){++failures;if(failures<4)std::cerr<<"compatibility mismatch length="<<length<<" width="<<width<<'\n';}
            }
          ++controls;
        }
  std::cout<<"valid_controls="<<controls<<" exact_range_queries="<<queries<<" mismatches="<<failures<<'\n';
  return failures?1:0;
}
