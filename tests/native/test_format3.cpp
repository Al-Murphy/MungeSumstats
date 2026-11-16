#include "../../src/preservation_core.h"
#include <cassert>
#include <filesystem>
#include <iostream>
#include <functional>
using namespace mss;
static void write(const std::string& path,const std::string& text){std::ofstream out(path,std::ios::binary);out<<text;}
static void fails(const std::function<void()>& action){bool failed=false;try{action();}catch(const std::exception&){failed=true;}assert(failed);}
int main(int,char** argv){
  const std::string root=std::string(argv[1])+"/typed";std::filesystem::create_directories(root);
  const auto fasta=root+"/tiny.fa",sites=root+"/sites.tsv";
  write(fasta,">1\nCAAAAACGTACGT\n");write(sites,"chr\tpos\tref\talt\tid\taliases\n1\t8\tG\tA,T\tr\told;雪\n1\t8\tG\tA,T\tr\told\n1\t0\tA\tG\tz\t.\n");
  prepare(fasta,sites,root+"/two",2);prepare(fasta,sites,root+"/three",2,3);
  Reference two(root+"/two",fasta),three(root+"/three",fasta);
  assert(three.lookup("S:1:0").size()==1);assert(three.lookup("I:雪").size()==2);
  fails([&]{three.lookup("I:r",3);});assert(three.lookup("I:r",4).size()==2);
  for(const auto& other:{"A","C","G","T"})for(const auto& effect:{"A","C","G","T"}){
    Input q("1",8,"old",other,effect,.3,.2,"EAF");const auto a=two.harmonise(q),b=three.harmonise(q);
    assert(a.status==b.status&&a.chr==b.chr&&a.pos==b.pos&&a.ref==b.ref&&a.alt==b.alt&&a.id==b.id&&a.candidate_count==b.candidate_count&&a.swapped==b.swapped&&a.complemented==b.complemented);
    assert(a.effect==b.effect||(std::isnan(a.effect)&&std::isnan(b.effect)));assert(a.eaf==b.eaf||(std::isnan(a.eaf)&&std::isnan(b.eaf)));
  }
  assert(typed_text("opaque-é-雪"));assert(!typed_text(std::string("x\0y",3)));assert(!typed_text("\xc0\x80"));assert(!typed_text("\xed\xa0\x80"));assert(!typed_text("\xf4\x90\x80\x80"));
  for(const auto& bad:{"","01","-1","1x","18446744073709551616"})fails([&]{typed_uint(bad);});
  for(const auto& bad:{"-0","+1","9223372036854775808","-9223372036854775809"})fails([&]{typed_int(bad,false);});
  assert(typed_int("-9223372036854775808",false)==std::numeric_limits<std::int64_t>::min());
  auto mini=[&](const std::string& body,const std::string& key="I:r",const std::string& offset="0"){
    const auto path=root+"/mini";std::filesystem::remove_all(path);std::filesystem::create_directories(path);Fasta fa(fasta);fa.save(path+"/fasta.fai");
    write(path+"/core-version","mss-preservation-3\n");write(path+"/records",body);write(path+"/keys",key+"\t"+offset+"\n");std::ofstream offsets(path+"/offsets",std::ios::binary);write64(offsets,0);offsets.close();
    Reference ref(path,fasta);return ref.lookup(key);
  };
  const std::string good="A\t1\t8\tG\tA\tr\t8\n";assert(mini(good).size()==1);
  for(const auto& bad:{"X\t1\t8\tG\tA\tr\t8\n","A\t1\t0\tG\tA\tr\t8\n","A\t1\t8x\tG\tA\tr\t8\n","A\t1\t8\tG\tG\tr\t8\n","A\t1\t8\tG\tA\tr\t8","A\t1\t8\tG\tA\tr\t8\textra\n","S\t1\t8\tG\tA\n"})fails([&]{mini(bad);});
  for(const auto& bad:{"01","1","-1","18446744073709551616"})fails([&]{mini(good,"I:r",bad);});
  fails([&]{mini(good,"C:1:7");});fails([&]{mini("S\t1\t8\tG\tA\n","S:1:9");});
  const std::string longid(70000,'x');assert(mini("A\t1\t8\tG\tA\t"+longid+"\t8\n","I:"+longid).size()==1);
  std::cout<<"typed format3 structural, occurrence, cap and scientific parity fixtures passed\n";
}
