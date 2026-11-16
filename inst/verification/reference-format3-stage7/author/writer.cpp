#include "preservation_core.h"
#include <iostream>
int main(int argc,char** argv){if(argc!=6)return 2;try{mss::prepare(argv[1],argv[2],argv[3],std::stoull(argv[5]),std::stoi(argv[4]));}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
