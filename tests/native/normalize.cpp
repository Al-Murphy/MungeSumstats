// Test-only normalization adapter for optional independent oracle comparisons.
#include "../../src/preservation_core.h"
#include <iostream>
int main(int argc,char** argv){
  if(argc!=2)return 2;
  mss::Fasta fasta(argv[1]);std::string line;
  while(std::getline(std::cin,line)){
    auto fields=mss::split(line,'\t');if(fields.size()!=5)return 3;
    auto result=mss::normalize(fasta,fields[0],std::stoll(fields[1]),fields[2],fields[3]);
    std::cout<<fields[4]<<'\t'<<result.status<<'\t'<<result.pos<<'\t'<<result.ref<<'\t'<<result.alt<<'\n';
  }
}
