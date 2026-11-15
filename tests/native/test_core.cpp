#include "../../src/preservation_core.h"
#include <cassert>
#include <filesystem>
#include <iostream>
using namespace mss;
int main(int argc, char** argv) {
  std::string dir = argc > 1 ? argv[1] : "/tmp/mss-native-fixture";
  std::filesystem::remove_all(dir); std::filesystem::create_directories(dir);
  { std::ofstream f(dir+"/tiny.fa",std::ios::binary); f << ">1\nCAAAAA\nCGTACG\nT\n>2\nACGT\n"; }
  Fasta fa(dir+"/tiny.fa");
  assert(fa.sequence("1",6,4)=="ACGT");
  auto snp=normalize(fa,"1",7,"CG","CT");
  assert(snp.status=="ok" && snp.pos==8 && snp.ref=="G" && snp.alt=="T");
  auto del=normalize(fa,"1",4,"AA","A");
  assert(del.status=="ok" && del.pos==1 && del.ref=="CA" && del.alt=="C");
  auto ins=normalize(fa,"1",4,"A","AA");
  assert(ins.status=="ok" && ins.pos==1 && ins.ref=="C" && ins.alt=="CA");
  auto again=normalize(fa,"1",ins.pos,ins.ref,ins.alt);
  assert(again.pos==ins.pos && again.ref==ins.ref && again.alt==ins.alt);
  assert(normalize(fa,"1",8,"A","T").status=="reference_mismatch");
  assert(normalize(fa,"1",8,"G","<DEL>").status=="unsupported_allele");
  { std::ofstream f(dir+"/sites.tsv");
    f << "chr\tpos\tref\talt\tid\taliases\n"
      << "1\t8\tG\tA,T\trs8\told8;merged8\n"
      << "1\t4\tAA\tA\trsDel\t.\n"
      << "1\t10\tA\tC\trsDup\tshared\n"
      << "1\t12\tG\tA\trsDup\tshared\n";
  }
  prepare(dir+"/tiny.fa",dir+"/sites.tsv",dir+"/index",2);
  Reference index(dir+"/index",dir+"/tiny.fa");
  assert(index.lookup("I:merged8").size()==2);
  assert(index.lookup("I:rsDup").size()==2);
  assert(index.lookup("C:1:1").size()==1);
  assert(index.lookup("I:absent").empty());
  auto match=index.harmonise({"1",8,"old8","G","T",2.0,0.0,"EAF"});
  assert(match.status=="match" && match.effect==2.0 && match.eaf==0.0 && match.alt=="T");
  auto swap=index.harmonise({"1",8,"rs8","T","G",2.0,0.2,"EAF"});
  assert(swap.status=="swap" && swap.effect==-2.0 && std::isnan(swap.eaf));
  auto comp=index.harmonise({"1",8,"rs8","C","A",2.0,0.2,"EAF"});
  assert(comp.status=="complement" && comp.effect==2.0 && comp.eaf==0.2);
  auto sc=index.harmonise({"1",8,"rs8","A","C",2.0,0.2,"EAF"});
  assert(sc.status=="swapped_complement" && sc.effect==-2.0 && std::isnan(sc.eaf));
  auto aa=index.harmonise({"1",8,"rs8","A","T",2.0,0.2,"EAF"});
  assert(aa.status=="alt_alt_contrast" && std::isnan(aa.effect));
  auto noid=index.harmonise({"1",8,"","G","T",2.0,0.0,"EAF"});
  assert(noid.status=="match");
  auto indel=index.harmonise({"1",4,"","AA","A",2.0,0.2,"EAF"});
  assert(indel.status=="match" && indel.pos==1 && indel.ref=="CA");
  auto shifted=index.harmonise({"1",5,"","AA","A",2.0,0.2,"EAF"});
  assert(shifted.status=="match" && shifted.pos==1 && shifted.ref=="CA");
  auto id_canonical=index.harmonise({"",0,"rsDel","CA","C",2.0,0.2,"EAF"});
  assert(id_canonical.status=="match" && id_canonical.pos==1);
  auto zero_swap=index.harmonise({"1",4,"rsDel","A","AA",2.0,0.0,"EAF"});
  assert(zero_swap.status=="swap" && zero_swap.eaf==1.0);
  auto bad_eaf=index.harmonise({"1",8,"rs8","G","T",2.0,1.1,"EAF"});
  assert(bad_eaf.status=="match" && std::isnan(bad_eaf.eaf));
  auto maf=index.harmonise({"1",8,"rs8","G","T",2.0,0.1,"MAF"});
  assert(maf.status=="match" && std::isnan(maf.eaf));
  auto ambiguous=index.harmonise({"",0,"rsDup","G","A",2.0,0.2,"EAF"});
  assert(ambiguous.status=="match" && ambiguous.pos==12);
  auto unknown=index.harmonise({"1",9,"new","T","C",2.0,0.2,"EAF"});
  assert(unknown.status=="unmapped" && std::isnan(unknown.effect));
  auto conflict=index.harmonise({"1",10,"rs8","A","C",2.0,0.2,"EAF"});
  assert(conflict.status=="coordinate_id_conflict");
  { std::ofstream f(dir+"/many.tsv");f<<"chr\tpos\tref\talt\tid\taliases\n";
    for(int i=0;i<250;++i)f<<"2\t1\tA\tC\trs"<<i<<"\talias"<<i<<"\n";
  }
  prepare(dir+"/tiny.fa",dir+"/many.tsv",dir+"/many-index",1);
  Reference many(dir+"/many-index",dir+"/tiny.fa");
  for(int i=0;i<250;++i){auto candidates=many.lookup("I:alias"+std::to_string(i));assert(candidates.size()==1 && candidates[0].id=="rs"+std::to_string(i));}
  assert(many.lookup("I:rs250").empty());
  std::cout << "native fixtures passed\n";
}
