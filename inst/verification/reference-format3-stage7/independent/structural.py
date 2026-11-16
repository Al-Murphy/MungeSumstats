"""Independent tiny structural, occurrence and framing controls; no throughput runs."""
from pathlib import Path
from collections import Counter
import hashlib,itertools,json,shutil,struct,subprocess

ROOT=Path(__file__).parent
DATA=ROOT/'structural';DATA.mkdir()
fasta=DATA/'tiny.fa';fasta.write_text('>x\n'+'A'*128+'\n')
sites=DATA/'base.tsv';sites.write_text('chr\tpos\tref\talt\tid\taliases\nx\t1\tA\tC\tr\told\n')
base=DATA/'base';subprocess.run([str(ROOT/'core'),'prepare3',str(fasta),str(sites),str(base)],check=True)
fai=(base/'fasta.fai').read_bytes()
GOOD=b'A\tx\t1\tA\tC\tr\t1\n'
results=[]

def mini(name,record=GOOD,key=b'I:r',value=b'0',marker=b'mss-preservation-3\n',
         keys=None,offsets=None,query=None,success=False):
    d=DATA/name;d.mkdir()
    (d/'fasta.fai').write_bytes(fai);(d/'core-version').write_bytes(marker)
    (d/'records').write_bytes(record)
    (d/'keys').write_bytes(key+b'\t'+value+b'\n' if keys is None else keys)
    (d/'offsets').write_bytes(struct.pack('<Q',0) if offsets is None else offsets)
    q=DATA/(name+'-query.tsv');q.write_text((query or key).hex()+'\t1000\n')
    out=subprocess.run([str(ROOT/'lookup'),str(fasta),str(d),str(q)],capture_output=True)
    (DATA/(name+'-stdout.txt')).write_bytes(out.stdout);(DATA/(name+'-stderr.txt')).write_bytes(out.stderr)
    ok=out.returncode==0 and out.stdout.startswith(b'ok\t')
    assert ok==success,(name,out.returncode,out.stdout[:300],out.stderr[:300])
    results.append(dict(name=name,expected='accepted' if success else 'refused',exit=out.returncode,
        outcome=out.stdout.decode('utf8','replace').strip()[:180] or out.stderr.decode('utf8','replace').strip()[:180]))
    return out.stdout

mini('good-A',success=True)
mini('good-opaque-alias',key='I:ancien_é_雪'.encode(),success=True)
mini('good-UTF8-id',record='A\tx\t1\tA\tC\trs雪\t1\n'.encode(),key='I:同義'.encode(),success=True)
mini('good-S-unknown',record=b'S\tx\t1\tA\tC,<DEL>,.,N\n',key=b'S:x:1',success=True)
mini('good-S-empty-ALT',record=b'S\tx\t1\tA\t\n',key=b'S:x:1',success=True)
mini('good-S-zero',record=b'S\tx\t0\tA\tC\n',key=b'S:x:0',success=True)
mini('good-S-negative',record=b'S\tx\t-1\tA\tC\n',key=b'S:x:-1',success=True)
mini('good-S-minimum',record=b'S\tx\t-9223372036854775808\tA\tC\n',key=b'S:x:-9223372036854775808',success=True)
mini('good-canonical-C-key',record=b'A\tx\t2\tA\tC\tr\t1\n',key=b'C:x:2',success=True)
mini('good-original-C-key',record=b'A\tx\t2\tA\tC\tr\t1\n',key=b'C:x:1',success=True)

bad_records={
 'unknown-tag':GOOD.replace(b'A\t',b'X\t',1),'short-A':b'A\tx\t1\tA\tC\tr\n',
 'extra-A-field':GOOD[:-1]+b'\textra\n','wrong-key-type':b'S\tx\t1\tA\tC\n',
 'no-record-LF':GOOD[:-1],'NUL-id':GOOD.replace(b'\tr\t',b'\tr\x00bad\t'),
 'empty-chr':GOOD.replace(b'\tx\t',b'\t\t',1),'nonDNA-ref':GOOD.replace(b'\tA\tC\t',b'\tN\tC\t'),
 'nonDNA-alt':GOOD.replace(b'\tA\tC\t',b'\tA\tN\t'),'same-alleles':GOOD.replace(b'\tA\tC\t',b'\tA\tA\t'),
 'lowercase-canonical':GOOD.replace(b'\tA\tC\t',b'\ta\tC\t'),
 'original-zero':GOOD[:-2]+b'0\n','original-negative':GOOD[:-2]+b'-1\n',
 'canonical-zero':GOOD.replace(b'\tx\t1\t',b'\tx\t0\t'),
 'canonical-negative':GOOD.replace(b'\tx\t1\t',b'\tx\t-1\t'),
 'canonical-plus':GOOD.replace(b'\tx\t1\t',b'\tx\t+1\t'),
 'canonical-leading-zero':GOOD.replace(b'\tx\t1\t',b'\tx\t01\t'),
 'canonical-trailing-junk':GOOD.replace(b'\tx\t1\t',b'\tx\t1x\t'),
 'canonical-overflow':GOOD.replace(b'\tx\t1\t',b'\tx\t9223372036854775808\t'),
 'original-overflow':GOOD[:-2]+b'9223372036854775808\n',
}
for label,invalid in [('overlong',b'\xc0\x80'),('surrogate',b'\xed\xa0\x80'),('above-Unicode',b'\xf4\x90\x80\x80'),('truncated-UTF8',b'\xe2'),('continuation',b'\x80')]:
    bad_records[label]=GOOD.replace(b'\tr\t',b'\tr'+invalid+b'\t')
for name,record in bad_records.items():mini(name,record=record)
for name,key in [('unknown-key-type',b'X:r'),('empty-ID-key',b'I:'),('empty-C-chr',b'C::1'),('C-zero',b'C:x:0'),('C-plus',b'C:x:+1'),('C-leading-zero',b'C:x:01'),('NUL-key',b'I:r\x00'),('invalid-UTF8-key',b'I:\xff')]:
    mini(name,key=key)
mini('C-record-chr-disagreement',key=b'C:y:1')
mini('C-record-coordinate-disagreement',key=b'C:x:3')
mini('S-record-coordinate-disagreement',record=b'S\tx\t1\tA\tC\n',key=b'S:x:2')
mini('S-as-C',record=b'S\tx\t1\tA\tC\n',key=b'C:x:1')
mini('A-as-S',key=b'S:x:1')
for value in [b'',b'00',b'-1',b'+0',b'1',b'99999',b'18446744073709551616']:
    mini('bad-record-offset-'+str(len(results)),value=value)
mini('key-no-LF',keys=b'I:r\t0')
mini('key-no-separator',keys=b'I:r0\n')
mini('key-extra-field',keys=b'I:r\t0\textra\n')
mini('key-offset-middle',offsets=struct.pack('<Q',1))
mini('key-offset-out-of-range',offsets=struct.pack('<Q',999999))
mini('short-binary-offset-table',offsets=b'\0'*7)
for i,marker in enumerate([b'mss-preservation-4\n',b'mss-preservation-3',b'mss-preservation-3\nextra',b'mss-preservation-3\r\n']):mini('bad-marker-'+str(i),marker=marker)

# Page-spanning selected records/keys and UTF-8 sequences at the byte boundary.
for position in [65535,65536,65537]:
    prefix=b'Z'*(position-1)+b'\n'
    mini('page-start-'+str(position),record=prefix+GOOD,value=str(position).encode(),success=True)
    key_prefix=b'Z'*(position-1)+b'\n'
    mini('key-page-start-'+str(position),keys=key_prefix+b'I:r\t0\n',offsets=struct.pack('<Q',position),success=True)
head=b'A\tx\t1\tA\tC\t';identifier=b'x'*(65535-len(head))+'雪'.encode()+b'y'*100
mini('UTF8-spans-page',record=head+identifier+b'\t1\n',key=b'I:opaque',success=True)
longid=b'x'*70000+'é'.encode();mini('long-record-and-key',record=head+longid+b'\t1\n',key=b'I:'+longid,success=True)

# Physical occurrence/raw-key cardinality and S reconstruction stay independent
# of omitted S IDs. Distinct Unicode aliases remain byte-distinct.
fixture=DATA/'occurrences';fixture.mkdir()
fa=fixture/'tiny.fa';fa.write_text('>1\nCAAAAACGTACGT\n')
rows=['1\t8\tG\tA,T\tr\told;同義;é;é;old',
      '1\t8\tG\tA,T\tr\told',
      '1\t8\tG\t<DEL>\tunsupported-only\tunsupported-alias',
      '1\t4\tAA\tA\tdel\told_del',
      '1\t0\tA\tG\tzero\t.',
      '1\t-1\tA\tG\tnegative\t.',
      '1\tbad\tA\tG\tbad_position\t.']
s=fixture/'sites.tsv';s.write_text('chr\tpos\tref\talt\tid\taliases\n'+'\n'.join(rows)+'\n')
for fmt in [2,3]:subprocess.run([str(ROOT/'core'),'prepare'+str(fmt),str(fa),str(s),str(fixture/f'index{fmt}')],check=True)
counts=[]
for fmt in [2,3]:
    keys=(fixture/f'index{fmt}'/'keys').read_bytes().splitlines();counts.append(Counter(x.rsplit(b'\t',1)[0] for x in keys))
assert counts[0]==counts[1]
assert counts[1][b'I:r']==4 and counts[1][b'I:old']==4 and counts[1]['I:同義'.encode()]==2
assert b'I:unsupported-alias' not in counts[1] and counts[1][b'S:1:8']==3
records=(fixture/'index3/records').read_bytes().splitlines();assert sum(x.startswith(b'S\t') for x in records)==6 and sum(x.startswith(b'A\t') for x in records)==5
assert all(len(x.split(b'\t'))==(7 if x.startswith(b'A\t') else 5) for x in records)
queries=[(b'I:r',3),(b'I:r',4),(b'S:1:8',2),(b'S:1:8',3),('I:é'.encode(),2),('I:é'.encode(),2),('I:同義'.encode(),2)]
q=fixture/'queries.tsv';q.write_text(''.join(k.hex()+'\t'+str(cap)+'\n' for k,cap in queries))
outputs={}
for fmt in [2,3]:
    out=subprocess.check_output([str(ROOT/'lookup'),str(fa),str(fixture/f'index{fmt}'),str(q)])
    (fixture/f'lookup{fmt}.tsv').write_bytes(out);outputs[fmt]=out.splitlines()
assert outputs[2][0].startswith(b'error\tReference candidate limit exceeded') and outputs[3][0].startswith(b'error\tReference candidate limit exceeded')
assert outputs[2][2].startswith(b'error\tReference candidate limit exceeded') and outputs[3][2].startswith(b'error\tReference candidate limit exceeded')
assert all(outputs[3][i].startswith(b'ok\t') for i in [1,3,4,5,6])
result=dict(selected_record_cases=len(results),accepted=sum(r['expected']=='accepted' for r in results),
            refused=sum(r['expected']=='refused' for r in results),cases=results,
            occurrence_payload_tags={'S':6,'A':5},raw_key_multiplicity_exact=True,unicode_aliases_distinct=True,
            raw_candidate_cap_boundaries_exact=True)
(ROOT/'structural-summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
