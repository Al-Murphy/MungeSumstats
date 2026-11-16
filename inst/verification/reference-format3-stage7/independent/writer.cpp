#include "preservation_core.h"
#include <iostream>
int main(int argc,char** argv) {
 if(argc!=6)return 2;
 try {
#ifdef LEGACY_SCI7
 mss::prepare(argv[1],argv[2],argv[3],std::stoull(argv[5]));
#else
 mss::prepare(argv[1],argv[2],argv[3],std::stoull(argv[5]),std::stoi(argv[4]));
#endif
 }
 catch(const std::exception& e){ std::cerr<<e.what()<<'\n';return 1; }
}
