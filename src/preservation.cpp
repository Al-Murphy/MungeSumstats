#include <Rcpp.h>
#include <memory>
#include "preservation_core.h"
// [[Rcpp::plugins(cpp17)]]
using namespace Rcpp;
// [[Rcpp::export]]
void mss_prepare_cpp(std::string fasta, std::string sites, std::string path, int batch_size) {
  if(batch_size<1)stop("batch_size must be positive");
  mss::prepare(fasta,sites,path,static_cast<std::size_t>(batch_size));
}
// [[Rcpp::export]]
SEXP mss_open_reference_cpp(std::string path, std::string fasta) {
  XPtr<mss::Reference> pointer(new mss::Reference(path,fasta),true);
  return pointer;
}
// [[Rcpp::export]]
DataFrame mss_harmonise_cpp(std::string path, std::string fasta,
                          CharacterVector chr, NumericVector pos, CharacterVector id,
                          CharacterVector other, CharacterVector effect_allele,
                          NumericVector effect, NumericVector eaf,
                          std::string frequency_type, int candidate_limit, SEXP reference_handle = R_NilValue,
                          Nullable<CharacterVector> source_ref = R_NilValue,
                          Nullable<CharacterVector> source_alt = R_NilValue,
                          Nullable<CharacterVector> source_strand = R_NilValue) {
  R_xlen_t n=chr.size();
  if(pos.size()!=n || id.size()!=n || other.size()!=n || effect_allele.size()!=n || effect.size()!=n || eaf.size()!=n)stop("Native input columns differ in length");
  if(candidate_limit<1)stop("candidate_limit must be positive");
  bool contract=source_ref.isNotNull()||source_alt.isNotNull()||source_strand.isNotNull();
  CharacterVector source_refs,source_alts,source_strands;
  if(contract){
    if(source_ref.isNull()||source_alt.isNull()||source_strand.isNull())stop("Source REF, ALT and strand are required together");
    source_refs=source_ref.get();source_alts=source_alt.get();source_strands=source_strand.get();
    if(source_refs.size()!=n||source_alts.size()!=n||source_strands.size()!=n)stop("Source contract columns differ in length");
  }
  std::unique_ptr<mss::Reference> local;
  mss::Reference* reference;
  if(reference_handle==R_NilValue){local.reset(new mss::Reference(path,fasta));reference=local.get();}
  else {XPtr<mss::Reference> pointer(reference_handle);reference=pointer.get();if(!reference)stop("Reference handle is no longer live; reopen its reference path");}
  CharacterVector status(n),out_chr(n,NA_STRING),ref(n,NA_STRING),alt(n,NA_STRING),out_id(n,NA_STRING);
  NumericVector out_pos(n,NA_REAL),out_effect(n,NA_REAL),out_eaf(n,NA_REAL),candidates(n,0.0);
  LogicalVector swapped(n,NA_LOGICAL),complemented(n,NA_LOGICAL);
  auto str=[](CharacterVector& x,R_xlen_t i){return CharacterVector::is_na(x[i])?std::string(""):as<std::string>(x[i]);};
  for(R_xlen_t i=0;i<n;++i){
    if(i%1000==0)checkUserInterrupt();
    bool valid_pos=NumericVector::is_na(pos[i]) || (std::isfinite(pos[i]) && pos[i]>=1 && pos[i]<=9007199254740991.0 && std::floor(pos[i])==pos[i]);
    if(!valid_pos){status[i]="invalid_position";continue;}
    mss::Input input{str(chr,i),NumericVector::is_na(pos[i])?0:static_cast<std::int64_t>(pos[i]),str(id,i),str(other,i),str(effect_allele,i),effect[i],eaf[i],frequency_type};
    if(contract){input.source_contract=true;input.source_ref=str(source_refs,i);input.source_alt=str(source_alts,i);input.source_strand=str(source_strands,i);}
    auto result=reference->harmonise(input,static_cast<std::size_t>(candidate_limit));
    status[i]=result.status;candidates[i]=result.candidate_count;
    if(result.pos>0){out_chr[i]=result.chr;out_pos[i]=result.pos;ref[i]=result.ref;alt[i]=result.alt;if(!mss::missing(result.id))out_id[i]=result.id;swapped[i]=result.swapped;complemented[i]=result.complemented;}
    if(std::isfinite(result.effect))out_effect[i]=result.effect;
    if(std::isfinite(result.eaf))out_eaf[i]=result.eaf;
  }
  return DataFrame::create(_["status"]=status,_["chr"]=out_chr,_["pos"]=out_pos,_["ref"]=ref,_["alt"]=alt,_["id"]=out_id,_["effect"]=out_effect,_["eaf"]=out_eaf,_["swapped"]=swapped,_["complemented"]=complemented,_["candidate_count"]=candidates);
}
