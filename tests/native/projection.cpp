#include "../../src/projection.h"
#include <cassert>
#include <iostream>
int main(int argc,char** argv) {
  assert(argc==2);std::string file=argv[1];
  {std::ofstream f(file);f<<"chain 100 1 20 + 0 20 1 20 - 0 20 1\n20\n";}
  mss::Projection p(file,{}, {},{{"1",20}},{{"1",20}});
  for(int repeat=0;repeat<100;++repeat)for(int pos=1;pos<=20;++pos)for(int width=1;width<=10;++width){
    auto m=p.map("1",pos,width);
    if(pos+width-1<=20){assert(m.status=="mapped");assert(m.pos==22-pos-width);assert(m.reverse);}
    else assert(m.status=="chain_gap");
  }
  assert(p.map("2",1,1).status=="unmapped");
  {std::ofstream f(file);f<<"chain 100 1 20 + 0 20 1 20 + 0 20 1\n20\n\nchain 100 chr1 20 + 0 20 1 20 + 0 20 1\n1 1 1\n18\n";}
  mss::Projection aliases(file,{{"chr1","1"}}, {},{{"1",20}},{{"1",20}});
  assert(aliases.map("1",1,2).status=="multimap");
  std::filesystem::remove(file);
  std::cout<<"20000 span checks and alias-overlap regression passed\n";
}
