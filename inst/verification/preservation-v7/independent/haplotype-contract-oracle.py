"""Independent bounded whole-haplotype and pair-flank witnesses; no native oracle."""
from pathlib import Path
from collections import Counter
import itertools,json,random,subprocess

ROOT=Path(__file__).parent
rng=random.Random(6161987)
DNA='ACGT'
def rc(s):return s.translate(str.maketrans(DNA,'TGCA'))[::-1]
words=[''.join(w) for length in range(1,5) for w in itertools.product(DNA,repeat=length)]
contexts=[]
manual=[('AAAAAAAAAA',1,'AA',['GC','GT','AG']),
        ('ACAAAAAAAA',1,'AC',['GC','GT','AT']),
        ('CAAAAAAAAA',1,'CAA',['CA','CAAA','CAAAA']),
        ('AAAAAAAAAA',1,'AA',['ATA','TA']),
        ('GACAAAAAAA',1,'GAC',['GAT','GGT'])]
for i in range(120):
    if i<len(manual):seq,pos,ref,alts=manual[i]
    else:
        seq=rng.choice(['AAAAAAAAAAAA','CAAAAAAAAAAA','ACACACACACAC','ATATATATATAT',
                        ''.join(rng.choices(DNA,k=12))])
        pos=rng.randint(1,4);length=rng.randint(1,4);ref=seq[pos-1:pos-1+length]
        pool=[w for w in words if w!=ref]
        alts=rng.sample(pool,3)
    contexts.append(dict(chr=f'h{i}',id=f'event{i}',seq=seq,pos=pos,ref=ref,alts=alts))

def trimmed_pairs(a,b,ref_length,position):
    # Cartesian enumeration belongs only to the independent tiny oracle.
    for prefix in range(min(len(a),len(b))):
        if a[:prefix]!=b[:prefix]:continue
        for suffix in range(min(len(a)-prefix,len(b)-prefix)):
            if prefix+suffix>=ref_length:continue
            if suffix and a[-suffix:]!=b[-suffix:]:continue
            yield position+prefix,a[prefix:len(a)-suffix if suffix else len(a)],b[prefix:len(b)-suffix if suffix else len(b)]

def altalt(c,pos,a,b):
    seq=c['seq'];start=c['pos']-1
    haplotypes={seq[:start]+alt+seq[start+len(c['ref']):] for alt in c['alts']}
    for x,y in [(a,b),(rc(a),rc(b))]:
        for first,second in itertools.permutations(c['alts'],2):
            for anchor,q0,q1 in trimmed_pairs(first,second,len(c['ref']),c['pos']):
                if (pos==0 or pos==anchor) and (x,y)==(q0,q1):return 'paired_shared_flank'
        for anchor in ([pos] if pos else range(1,len(seq)+1)):
            if anchor<1 or anchor>len(seq):continue
            for length in range(1,len(seq)-anchor+2):
                left=seq[:anchor-1];right=seq[anchor-1+length:]
                if left+x+right in haplotypes and left+y+right in haplotypes:
                    return 'whole_chromosome_haplotype'
    return None

def refalt_hypotheses(c,pos,a,b):
    seq=c['seq'];start=c['pos']-1
    haplotypes={seq[:start]+alt+seq[start+len(c['ref']):] for alt in c['alts']}
    hypotheses=set()
    for comp,swap in itertools.product((False,True),repeat=2):
        x,y=(rc(a),rc(b)) if comp else (a,b)
        if swap:x,y=y,x
        for anchor in ([pos] if pos else range(1,len(seq)+1)):
            if seq[anchor-1:anchor-1+len(x)]!=x:continue
            hap=seq[:anchor-1]+y+seq[anchor-1+len(x):]
            if hap in haplotypes:hypotheses.add((hap,swap,comp))
    return hypotheses

def potential(c,pos,a,b):
    alts=[x for x in set(c['alts']) if x!=c['ref']]
    if len(a)!=1 or len(b)!=1:
        delta=abs(len(a)-len(b))
        return any(abs(len(x)-len(y))==delta for x,y in itertools.combinations(alts,2))
    # Independently derive the pure-SNP anchor from all genomic differences.
    differences=set()
    pure=all(len(x)==len(c['ref']) for x in alts)
    if pure:
        for x in alts:differences.update(i for i,(r,v) in enumerate(zip(c['ref'],x)) if r!=v)
        pure=len(differences)==1
    if pure and pos and pos!=c['pos']+next(iter(differences)):return False
    pairs={(a,b),(b,a),(rc(a),rc(b)),(rc(b),rc(a))}
    for x,y in itertools.combinations(alts,2):
        if len(x)!=len(y):continue
        differences=[(u,v) for u,v in zip(x,y) if u!=v]
        if len(differences)==1 and differences[0] in pairs:return True
    return False

queries=[];meta=[]
for c in contexts:
    local=set()
    for a,b in itertools.permutations(c['alts'],2):
        local.update(trimmed_pairs(a,b,len(c['ref']),c['pos']))
    for alt in c['alts']:
        local.update(trimmed_pairs(c['ref'],alt,len(c['ref']),c['pos']))
    for _ in range(28):
        a,b=rng.sample(words,2);local.add((rng.randint(1,len(c['seq'])),a,b))
    for anchor,a,b in sorted(local):
        if a==b:continue
        for comp,swap in itertools.product((False,True),repeat=2):
            q0,q1=(rc(a),rc(b)) if comp else (a,b)
            if swap:q0,q1=q1,q0
            for mode in ('known','id_only','chr_only','pos_only'):
                pos=0 if mode in ('id_only','chr_only') else anchor
                chrvalue='.' if mode in ('id_only','pos_only') else c['chr']
                freq=('EAF','MAF','unknown')[len(queries)%3]
                q=(chrvalue,pos,c['id'],q0,q1,-.625,.137,freq,1000)
                queries.append(q);meta.append((c,pos,q0,q1,mode))
fasta=ROOT/'haplotypes.fa';fasta.write_text(''.join('>'+c['chr']+'\n'+c['seq']+'\n' for c in contexts))
sites=ROOT/'haplotype-sites.tsv';sites.write_text('chr\tpos\tref\talt\tid\taliases\n'+''.join(f"{c['chr']}\t{c['pos']}\t{c['ref']}\t{','.join(c['alts'])}\t{c['id']}\told_{c['id']}\n" for c in contexts))
index=ROOT/'haplotype-index'
if not index.exists():subprocess.run([str(ROOT/'core'),'prepare',str(fasta),str(sites),str(index)],check=True)
qfile=ROOT/'haplotype-queries.tsv';qfile.write_text(''.join('\t'.join(map(str,q))+'\n' for q in queries))
output=subprocess.check_output([str(ROOT/'core'),'query',str(fasta),str(qfile),str(index)],text=True)
(ROOT/'haplotype-output.tsv').write_text(output)
rows=output.splitlines();assert len(rows)==len(queries)
counts=Counter();violations=[];cache={};refcache={}
for i,(q,m,line) in enumerate(zip(queries,meta,rows)):
    fields=line.split('\t');c,pos,a,b,mode=m
    assert fields[0]=='ok',(i,q,fields)
    resolved=bool(fields[2])
    key=(c['chr'],pos,a,b)
    if key not in cache:cache[key]=altalt(c,pos,a,b)
    if key not in refcache:refcache[key]=refalt_hypotheses(c,pos,a,b)
    witness=cache[key]
    hypotheses=refcache[key]
    possible=potential(c,pos,a,b)
    complex_missing=pos==0 and (len(a)!=1 or len(b)!=1)
    if complex_missing:assert fields[1]=='source_coordinates_missing',(q,fields)
    if fields[1]=='source_coordinates_missing':
        assert complex_missing,(q,fields)
        counts['coordinate_deferrals_verified']+=1
    if fields[1]=='source_allele_context_required':
        assert possible and not complex_missing,(q,fields,'unjustified necessary-pair deferral')
        counts['necessary_pair_deferrals_verified']+=1
    counts['queries']+=1;counts['status_'+fields[1]]+=1
    if witness:
        counts['independent_'+witness+'_witness']+=1
        if resolved:
            violations.append(dict(kind='unresolved_ALTALT_explanation',query_index=i,context=c,query=q,witness=witness,actual=fields))
            counts['false_resolved_with_ALTALT_witness']+=1
        else:counts['withheld_with_independent_ALTALT_witness']+=1
    if fields[1]=='alt_alt_contrast':
        assert witness,(q,fields,'unwitnessed ALTALT')
        counts['all_ALTALT_diagnostics_independently_witnessed']+=1
    if fields[1]=='ambiguous':
        assert witness or len(hypotheses)>1,(q,fields,'false ambiguity')
        counts['all_ambiguity_outputs_independently_witnessed']+=1
    if resolved:
        if possible:
            violations.append(dict(kind='resolved_despite_necessary_ALTALT_pair',query_index=i,context=c,query=q,actual=fields))
            counts['false_resolved_despite_necessary_ALTALT_pair']+=1
        if len(hypotheses)!=1:
            violations.append(dict(kind='nonunique_global_REFALT_identity',query_index=i,context=c,query=q,
                hypotheses=sorted(hypotheses),actual=fields))
            counts['false_resolved_with_nonunique_REFALT_identity']+=1
        chrvalue,anchor,ref,alt=fields[2:6];anchor=int(anchor);seq=c['seq']
        assert chrvalue==c['chr'] and seq[anchor-1:anchor-1+len(ref)]==ref
        hap=seq[:anchor-1]+alt+seq[anchor-1+len(ref):]
        swap,comp=map(int,fields[10:12]);x,y=(rc(a),rc(b)) if comp else (a,b)
        if swap:x,y=y,x
        anchors=[pos] if pos else range(1,len(seq)+1)
        assert any(seq[p-1:p-1+len(x)]==x and seq[:p-1]+y+seq[p-1+len(x):]==hap for p in anchors),(q,fields)
        assert float.fromhex(fields[7])==(.625 if swap else -.625)
        if q[7]=='EAF' and not swap:assert float.fromhex(fields[8])==.137
        else:assert fields[8]=='nan'
        counts['resolved_payload_independently_validated']+=1
    else:assert fields[2:9]==['','0','','','','nan','nan'],(q,fields)
result=dict(source='5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e',contexts=len(contexts),counts=dict(counts),
            oracle='Full alternate chromosome sequences compared after replacing one inferred source REF span; independently enumerated tiny ALT-pair common-flank trims including shared backgrounds absent from REF.',
            violations=violations[:30],all_violation_count=len(violations))
(ROOT/'haplotype-summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('violations','oracle')}))
