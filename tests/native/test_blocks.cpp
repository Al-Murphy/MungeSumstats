#include "../../src/preservation_core.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
int main(int,char** argv) {
  std::string path=std::string(argv[1])+"/pages.txt";
  // Deliberately put a newline and a numeric offset across page boundaries.
  {std::ofstream out(path,std::ios::binary);out<<"abc\ndefghijklmnop\nlast";}
  mss::BlockFile file(path,4,2);
  assert(file.line(0)=="abc");
  assert(file.line(4)=="defghijklmnop");
  assert(file.line(18)=="last");
  assert(file.bytes(2,8)=="c\ndefghi");
  assert(file.bytes(21,5)=="t");
  assert(file.bytes(22,5).empty());
  for(int i=0;i<100;++i) {
    assert(file.bytes(0,3)=="abc");
    assert(file.bytes(18,4)=="last");
    assert(file.resident_bytes()<=8);
  }
  // Repeated reads in resident pages must avoid additional OS reads.
  auto before=file.disk_reads();
  assert(file.bytes(18,4)=="last");
  assert(file.disk_reads()==before);
  assert(file.size()==22);
  std::cout<<"bounded block I/O fixtures passed\n";
}
