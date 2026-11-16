#include "preservation_core.h"
#include <iostream>
int main(int argc,char** argv) {
    if(argc<5||argc>6)return 2;
    try {
        if(std::string(argv[1])=="prepare") {
            mss::prepare(argv[2],argv[3],argv[4],100000,argc==6?std::stoi(argv[5]):2);
            return 0;
        }
        mss::Reference reference(argv[4],argv[2]);
        if(std::string(argv[1])=="lookup") {
            for(const auto& c:reference.lookup(argv[3],argc==6?std::stoull(argv[5]):10000))std::cout<<c.key()<<"\t"<<c.alleles<<"\t"<<c.id<<"\n";return 0;
        }
        std::ifstream queries(argv[3]);std::string line;
        while(std::getline(queries,line)) {
            auto v=mss::split(line,'\t');if(v.size()!=13)return 3;
            mss::Input input(v[0],std::stoll(v[1]),v[2],v[3],v[4],std::stod(v[5]),std::stod(v[6]),v[7]);
            input.source_contract=v[8]=="declared";
            input.source_ref=v[9];input.source_alt=v[10];input.source_strand=v[11];
            try {
                auto r=reference.harmonise(input,std::stoull(v[12]));
                std::cout<<"ok\t"<<r.status<<'\t'<<r.chr<<'\t'<<r.pos<<'\t'<<r.ref<<'\t'<<r.alt<<'\t'<<r.id
                         <<'\t'<<std::hexfloat<<r.effect<<'\t'<<r.eaf<<'\t'<<r.candidate_count
                         <<'\t'<<r.swapped<<'\t'<<r.complemented<<'\n';
            }catch(const std::exception& e){std::cout<<"error\t"<<e.what()<<'\n';}
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
