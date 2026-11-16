#include "preservation_core.h"
#include <iostream>
#include <limits>

static void write_file(const std::string& path,const std::string& value) {
  std::ofstream out(path,std::ios::binary);out.write(value.data(),value.size());
  if(!out)throw std::runtime_error("review fixture write failed");
}
static std::string unhex(const std::string& input) {
  auto val=[](char c){return c<='9'?c-'0':c-'a'+10;};
  std::string out;
  for(std::size_t i=0;i<input.size();i+=2)out.push_back(static_cast<char>((val(input[i])<<4)|val(input[i+1])));
  return out;
}
int main(int argc,char**argv) {
  if(argc!=3)return 2;
  const std::string root=argv[1];std::filesystem::create_directories(root);
  const auto fasta=root+"/tiny.fa",index=root+"/tiny.fai";
  const std::string content=">x\n"+std::string("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")+"\n";
  write_file(fasta,content);
  std::ifstream cases(argv[2]);std::string line;
  std::size_t total=0,accepted=0,refused=0,mismatch=0,queries=0,query_errors=0;
  while(std::getline(cases,line)) {
    auto fields=mss::split(line,'\t');if(fields.size()!=4)return 3;
    const bool expected=fields[2]=="1";write_file(index,unhex(fields[1]));
    bool actual=false;
    try {
      mss::Fasta ref(fasta,index);actual=true;++accepted;
      // Query only oracle-safe layouts, never malformed baseline acceptances.
      if(expected) {
        const auto expected_queries=mss::split(fields[3],';');std::size_t qi=0;
        if(expected_queries.size()!=35)return 4;
        for(std::int64_t pos: {-1LL,0LL,1LL,2LL,63LL,64LL,std::numeric_limits<std::int64_t>::max()})
          for(std::size_t n:{std::size_t(0),std::size_t(1),std::size_t(2),std::size_t(64),std::numeric_limits<std::size_t>::max()}) {
            ++queries;
            try {auto seq=ref.sequence("x",pos,n);
              if(seq!=unhex(expected_queries[qi]))++query_errors;
            } catch(const std::exception& e) {++query_errors;if(query_errors<6)std::cerr<<"query exception "<<fields[0]<<" "<<e.what()<<'\n';}
            ++qi;
          }
      }
    }catch(const std::exception&){++refused;}
    if(expected!=actual){++mismatch;if(mismatch<12)std::cerr<<"constructor mismatch "<<fields[0]<<" expected="<<expected<<" actual="<<actual<<'\n';}
    ++total;
  }
  std::cout<<"total="<<total<<" accepted="<<accepted<<" refused="<<refused<<" constructor_mismatches="<<mismatch
           <<" range_queries="<<queries<<" query_errors="<<query_errors<<'\n';
  return mismatch||query_errors?1:0;
}
