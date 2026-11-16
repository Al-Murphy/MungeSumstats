#include "../../src/preservation_core.h"
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

static void write(const std::string& path,const std::string& data) {
  std::ofstream out(path,std::ios::binary);out<<data;
  if(!out)throw std::runtime_error("Fixture write failed");
}
int main(int,char** argv) {
  const std::string root=std::string(argv[1])+"/fai-safety";
  std::filesystem::create_directories(root);
  const auto fasta=root+"/reference.fa",index=root+"/reference.fai";
  write(fasta,">one\nACGT\nAC\n>two\nT\n");
  const std::string valid="one\t6\t5\t4\t5\ntwo\t1\t18\t1\t2\n";
  std::vector<std::pair<std::string,std::string>> bad={
    {"zero bases","one\t6\t5\t0\t5\n"},
    {"zero bytes","one\t6\t5\t4\t0\n"},
    {"narrow bytes","one\t6\t5\t4\t3\n"},
    {"bases exceed size","one\t6\t5\t7\t8\n"},
    {"offset beyond file","one\t6\t999\t4\t5\n"},
    {"last base beyond file","one\t100\t5\t4\t5\n"},
    {"offset overflow","one\t6\t18446744073709551615\t4\t5\n"},
    {"multiplication overflow","one\t18446744073709551615\t5\t1\t18446744073709551615\n"},
    {"duplicate name",valid+"one\t6\t5\t4\t5\n"},
    {"empty name","\t6\t5\t4\t5\n"},
    {"whitespace name","one bad\t6\t5\t4\t5\n"},
    {"empty index",""},
    {"empty numeric","one\t\t5\t4\t5\n"},
    {"negative","one\t-1\t5\t4\t5\n"},
    {"plus sign","one\t+6\t5\t4\t5\n"},
    {"suffix junk","one\t6x\t5\t4\t5\n"},
    {"leading whitespace","one\t 6\t5\t4\t5\n"},
    {"uint overflow","one\t18446744073709551616\t5\t4\t5\n"},
    {"extra field","one\t6\t5\t4\t5\t0\n"},
    {"inconsistent empty","empty\t0\t0\t1\t1\n"},
  };
  std::size_t failures=0,refused=0;
  for(const auto& item:bad) {
    write(index,item.second);bool rejected=false;
    try{mss::Fasta f(fasta,index);}catch(const std::exception&){rejected=true;}
    if(!rejected){++failures;std::cerr<<"UNSAFE_ACCEPTANCE "<<item.first<<'\n';}
    else ++refused;
  }
  auto control=[&](const std::string& data,const std::string& fai,const std::string& chr,const std::string& expected) {
    write(fasta,data);write(index,fai);mss::Fasta ref(fasta,index);
    if(ref.sequence(chr,1,expected.size())!=expected)throw std::runtime_error("Valid indexed sequence changed");
    if(ref.sequence(chr,1,std::numeric_limits<std::size_t>::max())!="")throw std::runtime_error("Oversized query accepted");
    if(ref.sequence(chr,std::numeric_limits<std::int64_t>::max(),2)!="")throw std::runtime_error("Out-of-bounds query accepted");
  };
  control(">one\nACGT\nAC\n>two\nT\n",valid,"one","ACGTAC");
  control(">one\r\nACGT\r\nAC\r\n>two\r\nT\r\n","one\t6\t6\t4\t6\r\ntwo\t1\t22\t1\t3\r\n","one","ACGTAC");
  control(">one\nACGTAC","one\t6\t5\t6\t6\n","one","ACGTAC");
  control(">one\nACGTAC\n","one\t0006\t0005\t0006\t0007\n","one","ACGTAC");
  control(">雪\nAC\n","雪\t2\t5\t2\t3\n","雪","AC");
  write(fasta,">empty\n>one\nAC\n");mss::Fasta generated(fasta);generated.save(index);mss::Fasta reopened(fasta,index);
  if(reopened.sequence("empty",1,1)!=""||reopened.sequence("one",1,2)!="AC")throw std::runtime_error("Empty-contig writer compatibility changed");
  std::cout<<refused<<" invalid indexes refused; 6 valid controls; range boundaries passed\n";
  return failures?1:0;
}
