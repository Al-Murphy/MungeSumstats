"""Independent explicit contract fixtures; expected facts predate candidate outputs."""
from pathlib import Path
from collections import Counter
import itertools,json,subprocess

ROOT=Path(__file__).parent
OLD=Path('/private/tmp/gwas-harvest-stage2/review-mss-multiallelic-reanchored-stage6')
counts=Counter();cases=[]
def query(chr,pos,id,a,b,cap=1000):return (chr,pos,id,a,b,.625,.137,'EAF',cap)
def check(fields,status,payload=None,count=None,frequency='EAF'):
    if status=='error':
        assert fields[0]=='error' and ('candidate limit exceeded' in fields[1] or 'context closure limit exceeded' in fields[1]),fields
        counts['fail_closed_resource_refusals']+=1;return
    assert fields[0]=='ok' and fields[1]==status,(status,fields)
    if count is not None:assert int(fields[9])==count,(count,fields)
    counts['status_'+status]+=1
    if payload is None:
        assert fields[2:9]==['','0','','','','nan','nan'],fields
        counts['empty_uncertainty_payloads']+=1
    else:
        chr,pos,ref,alt,id,swap,comp,biallelic=payload
        assert fields[2:7]==[chr,str(pos),ref,alt,id],fields
        assert float.fromhex(fields[7])==(-.625 if swap else .625),fields
        if frequency=='EAF' and (not swap or biallelic):assert float.fromhex(fields[8])==(1-.137 if swap else .137),fields
        else:assert fields[8]=='nan',fields
        assert fields[10:12]==[str(int(swap)),str(int(comp))],fields
        counts['verified_resolved_payloads']+=1

def run(name,fasta,index,items,frequencies=('EAF','MAF','unknown')):
    folder=ROOT/name;folder.mkdir(exist_ok=True)
    qs=[];expect=[]
    for q,status,payload,count in items:
        for freq in frequencies:
            q=list(q);q[7]=freq;qs.append(q);expect.append((status,payload,count,freq))
    qfile=folder/'queries.tsv';qfile.write_text(''.join('\t'.join(map(str,q))+'\n' for q in qs))
    output=subprocess.check_output([str(ROOT/'core'),'query',str(fasta),str(qfile),str(index)],text=True)
    (folder/'output.tsv').write_text(output);rows=output.splitlines();assert len(rows)==len(expect)
    for fields,e in zip(rows,expect):check(fields.split('\t'),*e)
    counts['queries']+=len(rows);cases.append(dict(name=name,queries=len(rows)))
    return output

def prepare(name,sequences,rows):
    folder=ROOT/name;folder.mkdir(exist_ok=True)
    fasta=folder/'tiny.fa';fasta.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in sequences.items()))
    sites=folder/'sites.tsv';sites.write_text('chr\tpos\tref\talt\tid\taliases\n'+'\n'.join(rows)+'\n')
    index=folder/'index';subprocess.run([str(ROOT/'core'),'prepare',str(fasta),str(sites),str(index)],check=True)
    return fasta,index

# Distinct identities and opposite signs at known positions are positive controls;
# neither a unique ID event nor a biallelic context supplies a missing anchor.
d=OLD/'unknown-anchor'
items=[]
for chr,pos in [('h95',3),('h95',4),('.',0),('h95',0),('.',3),('.',4)]:
    missing=pos==0
    payload=None if missing else ('h95',1 if pos==3 else 3,'C' if pos==3 else 'A','CCAC' if pos==3 else 'AACC','event95',False,False,False)
    items.append((query(chr,pos,'event95','A','ACCA'),'source_coordinates_missing' if missing else 'match',payload,3))
items.append((query('.',0,'absent','A','AA'),'source_coordinates_missing',None,0))
run('counter-unknown',d/'tiny.fa',d/'index',items)
d=OLD/'boundary'
run('counter-sign-and-pool',d/'tiny.fa',d/'index',[
    (query('pool',2,'old_r','C','T'),'ambiguous',None,1),
    (query('pool',2,'.','C','T'),'ambiguous',None,1),
    (query('.',0,'old_r','C','T'),'ambiguous',None,3),
    (query('unique',1,'old_u','CCAG','CCTG'),'match',('unique',3,'A','T','u',False,False,True),1),
    (query('unique',2,'old_u','CCAG','CCTG'),'swapped_complement',('unique',3,'A','T','u',True,True,True),1),
    (query('.',0,'old_u','CCAG','CCTG'),'source_coordinates_missing',None,1),
    (query('unique',0,'old_u','CCAG','CCTG'),'source_coordinates_missing',None,1)])

for name,pairs,pos in [('alternate-background-repeat',[('A','AAA'),('AAA','A'),('T','TTT'),('TTT','T')],4),
                       ('singlebase-alternate-indel-background',[('A','G'),('G','A'),('T','C'),('C','T')],2)]:
    d=OLD/name
    run('counter-'+name,d/'tiny.fa',d/'index',[(query('c',pos,'old_r',a,b),'source_allele_context_required',None,3 if pos==4 else 1) for a,b in pairs])

# Closed ID/event hypotheses reject the main selection but cannot replace it.
base_rows=['c\t1\tA\tG\tr\told_r','c\t1\tA\tG\tq\told_q','c\t3\tA\tG\tq\told_q']
order_outputs=[]
for i,order in enumerate(itertools.permutations(base_rows)):
    f,index=prepare('closure-order-'+str(i),{'c':'AAAAA'},list(order))
    out=run('closure-order-'+str(i),f,index,[
        (query('.',0,'r','A','G'),'ambiguous',None,1),
        (query('c',0,'old_r','A','G'),'ambiguous',None,1),
        (query('.',0,'old_q','A','G'),'ambiguous',None,2),
        (query('c',1,'r','A','G'),'match',('c',1,'A','G','r',False,False,True),1),
        (query('c',3,'r','A','G'),'coordinate_id_conflict',None,0)])
    if order_outputs:assert out==order_outputs[0]
    order_outputs.append(out)

# Unrelated overlapping C events must NOT bridge a new ID. Shared r/q events may.
f,index=prepare('unrelated-coordinate',{'c':'AAAAA'},['c\t1\tA\tG\tr\told_r','c\t1\tA\tC\tq\t.','c\t3\tA\tG\tq\t.'])
run('unrelated-coordinate',f,index,[(query('.',0,'old_r','A','G'),'match',('c',1,'A','G','r',False,False,False),1)])

# Position exclusion is proved jointly with REF/all ALTs, not by dropping padding.
for i,other in enumerate(['c\t3\tA\tT,C\tr\t.','c\t2\tAA\tAT,AC\tr\t.','c\t3\tAA\tTA,CA\tr\t.','c\t2\tAAA\tATA,ACA\tr\t.']):
    for order in (False,True):
        rows=['c\t1\tA\tG\tr\told_r',other]
        if order:rows.reverse()
        name=f'pure-position-{i}-{order}';f,index=prepare(name,{'c':'AAAAAAA'},rows)
        run(name,f,index,[(query('c',1,'old_r',a,b),status,('c',1,'A','G','r',swap,comp,True),1)
            for a,b,status,swap,comp in [('A','G','match',False,False),('G','A','swap',True,False),('T','C','complement',False,True),('C','T','swapped_complement',True,True)]])

# Valid length/base exclusions and duplicate tokens retain verified coverage.
for name,seq,row,a,b,ref,alt in [
    ('biallelic-indel','AAAAAAA','c\t1\tA\tAT\tr\told_r','A','AT','A','AT'),
    ('noncompeting-lengths','AAAAAAA','c\t1\tA\tAT,ATTT\tr\told_r','A','AT','A','AT'),
    ('duplicate-ALT','AAAAAAA','c\t1\tA\tC,C\tr\told_r','A','C','A','C')]:
    f,index=prepare(name,{'c':seq},[row])
    run(name,f,index,[(query('c',1,'old_r',a,b),'match',('c',1,ref,alt,'r',False,False,True),1 if name!='noncompeting-lengths' else 2)])

# Unknowns require separate conservative deferral; one unknown-only different
# original REF span cannot be pooled with another span's single concrete ALT.
for i,marker in enumerate(['<DEL>','.','N']):
    name=f'unknown-marker-{i}';f,index=prepare(name,{'c':'AAAAAAA'},[f'c\t1\tA\tC,{marker}\tr\told_r'])
    run(name,f,index,[(query('c',1,'old_r','A','C'),'source_allele_context_required',None,1)])
f,index=prepare('single-unknown-separate-span',{'c':'AAAAAAA'},['c\t1\tA\tC\tr\told_r','c\t1\tAA\t<DEL>\tother\t.'])
run('single-unknown-separate-span',f,index,[(query('c',1,'old_r','A','C'),'match',('c',1,'A','C','r',False,False,False),1)])

# Current-ID cycles, raw caps, context counts and raw work all fail closed.
f,index=prepare('id-cycle',{'c':'AAAAA'},['c\t1\tA\tC\tr\tq','c\t1\tA\tC\tq\tr'])
run('id-cycle',f,index,[(query('c',1,'r','A','C',2),'match',('c',1,'A','C','',False,False,True),1)])
f,index=prepare('raw-token-work',{'c':'AAAAA'},['c\t1\tA\tG,'+','.join(['N']*70)+'\tr\t.'])
run('raw-token-work',f,index,[(query('c',1,'r','A','G',2),'error',None,None),
                            (query('c',1,'r','A','G',10),'source_allele_context_required',None,1)])
f,index=prepare('signature-work',{'c':'AAAAA'},['c\t1\tA\tG,C'+'A'*100+'\tr\t.'])
run('signature-work',f,index,[(query('c',1,'r','A','G',2),'error',None,None),
                            (query('c',1,'r','A','G',100),'match',('c',1,'A','G','r',False,False,False),2)])
f,index=prepare('examined-context-cap',{'c':'AAAAA'},['c\t1\tA\tG\tr\t.','c\t1\tN\tC\tx\t.',
                                                'c\t3\tA\tC\tr\t.','c\t3\tN\tT\tx\t.'])
run('examined-context-cap',f,index,[(query('c',1,'r','A','G',2),'error',None,None),
                                 (query('c',1,'r','A','G',10),'match',('c',1,'A','G','r',False,False,False),1)])

summary=dict(source='5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e',counts=dict(counts),cases=cases,
             six_order_closed_hypotheses_full_fields_exact=True)
(ROOT/'contract-gold-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='cases'}))
