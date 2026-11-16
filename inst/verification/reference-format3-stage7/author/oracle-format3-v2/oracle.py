"""Independent finite SNP haplotype oracle, invariant to common VCF padding."""
from collections import Counter
from pathlib import Path
import hashlib, itertools, json, subprocess

ROOT=Path(__file__).parent
DNA='ACGT'
RC=dict(zip(DNA,'TGCA'))
REPS=('plain','left','right','both')
MODES=('coordinate','id','primary','alias')
seq={}; sites={rep:[] for rep in REPS}; metadata=[]
for ref in DNA:
    others=[a for a in DNA if a!=ref]
    for mask in range(1,8):
        alts=[a for i,a in enumerate(others) if mask>>i&1]
        chrom=f'c_{ref}_{mask}'
        seq[chrom]='G'+ref+'CAAAAA'
        identifier='id_'+chrom
        for rep in REPS:
            left='G' if rep in ('left','both') else ''
            right='C' if rep in ('right','both') else ''
            sites[rep].append((chrom,1 if left else 2,left+ref+right,
                              ','.join(left+a+right for a in alts),identifier,'old_'+identifier))
        for a,b in itertools.permutations(DNA,2):
            hits=[]; altalt=False
            for complement in (False,True):
                x,y=(RC[a],RC[b]) if complement else (a,b)
                if x==ref and y in alts:hits.append((False,complement,y))
                if y==ref and x in alts:hits.append((True,complement,x))
                if x in alts and y in alts:altalt=True
            if len(hits)>1 or (hits and altalt):expected='ambiguous'
            elif hits:
                swap,comp,_=hits[0]
                expected=('swapped_complement' if swap else 'complement') if comp else ('swap' if swap else 'match')
            else:expected='alt_alt_contrast' if altalt else 'allele_mismatch'
            for frequency,mode in itertools.product(('EAF','MAF','unknown'),MODES):
                c,p,idvalue=chrom,2,identifier
                if mode=='coordinate':idvalue='.'
                if mode=='id':c,p='.',0
                if mode=='alias':idvalue='old_'+identifier
                query=(c,p,idvalue,a,b,-.625,.137,frequency,1000)
                metadata.append(dict(query=query,expected=expected,ref=ref,alts=alts,chr=chrom,id=identifier,hits=hits))

fasta=ROOT/'oracle.fa';fasta.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in seq.items()))
queryfile=ROOT/'oracle-queries.tsv';queryfile.write_text(''.join('\t'.join(map(str,m['query']))+'\n' for m in metadata))
results={}; mismatch_examples=[]; totals=Counter()
for rep in REPS:
    source=ROOT/(rep+'-oracle-sites.tsv')
    source.write_text('chr\tpos\tref\talt\tid\taliases\n'+''.join('\t'.join(map(str,r))+'\n' for r in sites[rep]))
    index=ROOT/(rep+'-oracle-index')
    subprocess.run([str(ROOT/'core'),'prepare',str(fasta),str(source),str(index)],check=True)
    output=subprocess.check_output([str(ROOT/'core'),'query',str(fasta),str(queryfile),str(index)],text=True)
    (ROOT/(rep+'-oracle-output.tsv')).write_text(output)
    counts=Counter();rows=output.splitlines()
    assert len(rows)==len(metadata)
    for m,line in zip(metadata,rows):
        fields=line.split('\t');q=m['query'];actual=fields[1]
        if fields[0]!='ok' or actual!=m['expected']:
            counts['wrong_status']+=1
            counts[m['expected']+'->'+actual]+=1
            if len(mismatch_examples)<24:mismatch_examples.append(dict(rep=rep,case=m,actual=fields))
            continue
        counts['correct_status']+=1
        if actual in ('match','swap','complement','swapped_complement'):
            swap,comp,alt=m['hits'][0]
            assert fields[2:7]==[m['chr'],'2',m['ref'],alt,m['id']],(rep,m,fields)
            assert float.fromhex(fields[7])==(.625 if swap else -.625),(rep,m,fields)
            expected_eaf=(1-.137 if swap else .137) if q[7]=='EAF' and (not swap or len(m['alts'])==1) else None
            if expected_eaf is None:assert fields[8]=='nan',(rep,m,fields)
            else:assert float.fromhex(fields[8])==expected_eaf,(rep,m,fields)
            assert fields[10:]==[str(int(swap)),str(int(comp))]
            counts['correct_resolved_payload']+=1
        else:
            assert fields[2:9]==['','0','','','','nan','nan'],(rep,m,fields)
            counts['correct_empty_unresolved_payload']+=1
    results[rep]=dict(counts)
    totals.update(counts)
summary=dict(source='5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e',queries_per_representation=len(metadata),
             representations=REPS,total_queries=len(metadata)*len(REPS),
             mathematical_oracle='Enumerate forward and reverse-complement source pairs over one finite SNP allele set; unchanged under common REF/ALT prefix/suffix padding.',
             results=results,totals=dict(totals),examples=mismatch_examples)
(ROOT/'oracle-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('examples','mathematical_oracle')}))
