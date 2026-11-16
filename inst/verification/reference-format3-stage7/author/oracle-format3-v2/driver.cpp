#include "preservation_core.h"
#include <iostream>
int main(int argc,char** argv) {
  if(argc!=5)return 2;
  try {
    if(std::string(argv[1])=="prepare") {mss::prepare(argv[2],argv[3],argv[4],7,3);return 0;}
    mss::Reference r(argv[4],argv[2]);std::ifstream in(argv[3]);std::string line;
    while(std::getline(in,line)) {
      auto v=mss::split(line,'	');if(v.size()!=9)return 3;
      mss::Input q(v[0],std::stoll(v[1]),v[2],v[3],v[4],std::stod(v[5]),std::stod(v[6]),v[7]);
      try {
        auto x=r.harmonise(q,std::stoull(v[8]));
        std::cout<<"ok\t"<<x.status<<'	'<<x.chr<<'	'<<x.pos<<'	'<<x.ref<<'	'<<x.alt<<'	'<<x.id
          <<'	'<<std::hexfloat<<x.effect<<'	'<<x.eaf<<'	'<<x.candidate_count<<'	'<<x.swapped<<'	'<<x.complemented<<'\n';
      }catch(const std::exception& e){std::cout<<"error\t"<<e.what()<<'\n';}
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
