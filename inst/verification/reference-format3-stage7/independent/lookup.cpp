#include "preservation_core.h"
#include <iostream>
static std::string hex(const std::string& s){const char* d="0123456789abcdef";std::string out;for(unsigned char c:s){out+=d[c>>4];out+=d[c&15];}return out;}
static std::string unhex(const std::string& s){std::string out;for(std::size_t i=0;i<s.size();i+=2)out+=static_cast<char>(std::stoul(s.substr(i,2),nullptr,16));return out;}
int main(int argc,char** argv){
 if(argc!=4)return 2;
 try{
  mss::Reference ref(argv[2],argv[1]);std::ifstream in(argv[3]);std::string line;
  while(std::getline(in,line)){
   auto fields=mss::split(line,'\t');if(fields.size()!=2)return 3;
   try{
    auto values=ref.lookup(unhex(fields[0]),std::stoull(fields[1]));std::cout<<"ok\t"<<values.size();
    for(const auto& c:values){
     std::cout<<'\t'<<hex(c.chr)<<'|'<<c.pos<<'|'<<hex(c.ref)<<'|'<<hex(c.alt)<<'|'<<hex(c.id)<<'|'<<c.original_pos<<'|'<<hex(c.original_ref)<<'|'<<hex(c.alleles)<<'|';
     for(auto pos:c.original_positions)std::cout<<pos<<',';std::cout<<'|';for(const auto& id:c.identifiers)std::cout<<hex(id)<<',';
    }
    std::cout<<'\n';
   }catch(const std::exception& e){std::cout<<"error\t"<<e.what()<<'\n';}
  }
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
