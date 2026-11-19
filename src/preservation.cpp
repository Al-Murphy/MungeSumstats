#include <Rcpp.h>
#include <memory>
#include "preservation_core.h"
// [[Rcpp::plugins(cpp17)]]
using namespace Rcpp;
// [[Rcpp::export]]
void mss_prepare_cpp(std::string fasta, std::string sites, std::string path, int batch_size, int format_version=2) {
  if(batch_size<1)stop("batch_size must be positive");
  mss::prepare(fasta,sites,path,static_cast<std::size_t>(batch_size),format_version);
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
#include "preservation_stream.h"

// [[Rcpp::export]]
SEXP mss_parser_open_cpp(std::string delimiter,int max_record_bytes) {
  if(delimiter.size()!=1||max_record_bytes<1)stop("Invalid parser limits");
  return XPtr<MssDelimitedParser>(new MssDelimitedParser(delimiter[0],max_record_bytes),true);
}
// [[Rcpp::export]]
List mss_parser_feed_cpp(SEXP parser,RawVector input,bool eof=false) {
  XPtr<MssDelimitedParser> p(parser);return p->feed(input,eof);
}

// [[Rcpp::export]]
List mss_parser_batch_cpp(SEXP parser, RawVector input, int offset, bool eof,
                          int max_rows, int batch_bytes, int field_count, IntegerVector keep) {
    XPtr<MssDelimitedParser> p(parser);
    return p->batch(input, offset, eof, max_rows, batch_bytes, field_count, keep);
}

// [[Rcpp::export]]
void mss_write_normalized_cpp(DataFrame table, SEXP connection, bool header, int block_bytes = 262144) {
    MssNormalizedWriter writer(connection, block_bytes);
    writer.write(table, header);
}

#include "projection.h"
// [[Rcpp::export]]
SEXP mss_projection_open_cpp(std::string chain,CharacterVector source_names,CharacterVector target_names,NumericVector source_sizes,NumericVector target_sizes) {
  auto aliases=[](CharacterVector x){std::map<std::string,std::string> out;if(x.size()){CharacterVector names=x.names();for(R_xlen_t i=0;i<x.size();++i)out[as<std::string>(names[i])]=as<std::string>(x[i]);}return out;};
  auto sizes=[](NumericVector x){std::map<std::string,std::int64_t> out;CharacterVector names=x.names();for(R_xlen_t i=0;i<x.size();++i)out[as<std::string>(names[i])]=static_cast<std::int64_t>(x[i]);return out;};
  return XPtr<mss::Projection>(new mss::Projection(chain,aliases(source_names),aliases(target_names),sizes(source_sizes),sizes(target_sizes)),true);
}
// [[Rcpp::export]]
DataFrame mss_project_cpp(SEXP projection,SEXP target_handle,CharacterVector chr,NumericVector pos,
                         CharacterVector ref,CharacterVector alt,NumericVector effect,NumericVector eaf,
                         int candidate_limit,LogicalVector source_swapped) {
  XPtr<mss::Projection> p(projection);XPtr<mss::Reference> target(target_handle);
  if(!p.get()||!target.get())stop("Projection handles are no longer live; reopen immutable inputs");
  R_xlen_t n=chr.size();CharacterVector status(n,"native_unresolved"),oc(n,NA_STRING),orr(n,NA_STRING),oa(n,NA_STRING),id(n,NA_STRING),chain_id(n,NA_STRING),strand(n,NA_STRING);
  NumericVector op(n,NA_REAL),oe(n,NA_REAL),of(n,NA_REAL);
  for(R_xlen_t i=0;i<n;++i){
    if(i%1000==0)checkUserInterrupt();
    if(NumericVector::is_na(pos[i])||CharacterVector::is_na(ref[i])||CharacterVector::is_na(alt[i]))continue;
    auto r=as<std::string>(ref[i]),a=as<std::string>(alt[i]);
    auto m=p->map(as<std::string>(chr[i]),static_cast<std::int64_t>(pos[i]),r.size());status[i]=m.status;
    if(m.status!="mapped")continue;
    chain_id[i]=m.id;strand[i]=m.reverse?"-":"+";
    if(m.reverse){r=mss::complement(r);a=mss::complement(a);}
    mss::Input in{m.chr,m.pos,"",r,a,effect[i],eaf[i],"EAF"};
    in.source_contract=true;in.source_ref=r;in.source_alt=a;in.source_strand="+";
    auto result=target->harmonise(in,static_cast<std::size_t>(candidate_limit));
    if(result.pos<1){status[i]=result.status=="source_reference_mismatch"?"target_ref_conflict":result.status;continue;}
    status[i]="projected";oc[i]=result.chr;op[i]=result.pos;orr[i]=result.ref;oa[i]=result.alt;
    if(!mss::missing(result.id))id[i]=result.id;
    if(std::isfinite(result.effect))oe[i]=result.effect;
    // A source REF-effect complement additionally requires target biallelic
    // evidence. Unchanged ALT frequencies need no two-allele assumption.
    if(source_swapped[i]==TRUE && std::isfinite(result.eaf)) {
      auto guard=in;guard.other=a;guard.effect_allele=r;guard.eaf=0.5;
      if(!std::isfinite(target->harmonise(guard,static_cast<std::size_t>(candidate_limit)).eaf))result.eaf=std::numeric_limits<double>::quiet_NaN();
    }
    if(std::isfinite(result.eaf))of[i]=result.eaf;
  }
  return DataFrame::create(_["projection_status"]=status,_["projected_chr"]=oc,_["projected_pos"]=op,_["projected_ref"]=orr,_["projected_alt"]=oa,_["projected_id"]=id,_["projected_effect"]=oe,_["projected_eaf"]=of,_["projection_chain_id"]=chain_id,_["projection_strand"]=strand);
}
