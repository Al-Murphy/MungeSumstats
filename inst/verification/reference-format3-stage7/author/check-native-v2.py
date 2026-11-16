from pathlib import Path
import subprocess,shutil,struct,json,hashlib
BASE=Path(__file__).resolve().parent
ROOT=BASE/'native-small-v2';ROOT.mkdir(exist_ok=False)
DRIVER=BASE/'driver-v2'
fa=ROOT/'reference.fa';fa.write_text('>1\nCAAAAACGTACGT\n')
sites=ROOT/'sites.tsv';sites.write_text('chr\tpos\tref\talt\tid\taliases\n1\t8\tG\tA,T\tr雪\té;同義\n1\t4\tAA\tA\td\t.\n1\t8\tG\tA,T\tr雪\té\n1\t10\tN\tG\tu\t.\n1\t0\tA\tG\tzero\t.\n1\t-1\tA\tG\tnegative\t.\n')
def run(args,ok=True):
 p=subprocess.run([str(DRIVER),*map(str,args)],capture_output=True)
 if ok:assert p.returncode==0,(args,p.stderr)
 else:assert p.returncode!=0,(args,p.stdout)
 return p
refs={}
for fmt in (2,3):
 ref=ROOT/f'ref-{fmt}';run(['prepare',fa,sites,ref,fmt]);refs[fmt]=ref
rec=refs[3].joinpath('records').read_bytes();lines=rec.splitlines();assert sum(x.startswith(b'S\t') for x in lines)==6;assert sum(x.startswith(b'A\t') for x in lines)==5
assert run(['lookup',fa,'S:1:0',refs[3]]).stdout.startswith(b'1:0:A:')
assert run(['lookup',fa,'S:1:-1',refs[3]]).stdout.startswith(b'1:-1:A:')
assert run(['lookup',fa,'I:同義',refs[3]]).stdout
# Raw occurrence accounting: two equal raw rows are not globally deduplicated.
for fmt in (2,3):run(['lookup',fa,'I:r雪',refs[fmt],3],False)
assert run(['lookup',fa,'I:r雪',refs[3],4]).stdout.count(b'\n')==2
cases=[]
def corrupt(name,mutate,key='I:r雪'):
 dst=ROOT/name;shutil.copytree(refs[3],dst);mutate(dst);p=run(['lookup',fa,key,dst],False);cases.append({'case':name,'error':p.stderr.decode(errors='replace').strip()})
def records_replace(old,new):
 def f(d):
  p=d/'records';p.write_bytes(p.read_bytes().replace(old,new,1))
 return f
# First A starts after first S; equal-length changes keep offsets intact.
for name,old,new in [('wrong-A-tag',b'A\t1\t8\tG\tA',b'X\t1\t8\tG\tA'),('wrong-A-pos',b'A\t1\t8\tG\tA',b'A\t1\t0\tG\tA'),('A-ref-N',b'A\t1\t8\tG\tA',b'A\t1\t8\tN\tA'),('A-same-allele',b'A\t1\t8\tG\tA',b'A\t1\t8\tG\tG'),('A-invalidUTF8',b'r\xe9\x9b\xaa',b'r\xff\x9b\xaa'),('A-NUL',b'A\t1\t8\tG\tA',b'A\t1\t8\tG\t\0')]:corrupt(name,records_replace(old,new))
corrupt('S-wrong-tag',records_replace(b'S\t1\t8',b'A\t1\t8'),'S:1:8')
corrupt('S-wrong-position',records_replace(b'S\t1\t8',b'S\t1\t9'),'S:1:8')
corrupt('S-invalidUTF8',records_replace(b'S\t1\t8\tG\tA,T',b'S\t1\t8\tG\tA,\xff'),'S:1:8')
corrupt('marker-extra',lambda d:(d/'core-version').write_bytes(b'mss-preservation-3\nextra\n'))
corrupt('offsets-truncated',lambda d:(d/'offsets').write_bytes((d/'offsets').read_bytes()[:-1]))
corrupt('keyoffset-interior',lambda d:(d/'offsets').write_bytes(struct.pack('<Q',1)*((d/'offsets').stat().st_size//8)))
corrupt('keyoffset-beyond',lambda d:(d/'offsets').write_bytes(struct.pack('<Q',2**64-1)*((d/'offsets').stat().st_size//8)))
# Direct minimal indexes let each malformed number/field/line be exercised without stale offsets.
def mini(name,key,body,offset='0',keyrow=None):
 d=ROOT/name;d.mkdir();shutil.copy(refs[3]/'fasta.fai',d/'fasta.fai');(d/'core-version').write_bytes(b'mss-preservation-3\n');(d/'records').write_bytes(body);(d/'keys').write_bytes(keyrow if keyrow is not None else key.encode()+b'\t'+offset.encode()+b'\n');(d/'offsets').write_bytes(struct.pack('<Q',0));return d
A=b'A\t1\t8\tG\tA\tr\t8\n'
for i,bad in enumerate([b'A\t1\t8\tG\tA\tr\n',b'A\t1\t8\tG\tA\tr\t8\textra\n',A[:-1],A.replace(b'\t8\t',b'\t8x\t',1),A.replace(b'\t8\t',b'\t08\t',1),A.replace(b'\t8\t',b'\t+8\t',1),A.replace(b'\t8\t',b'\t9223372036854775808\t',1),A.replace(b'\t8\t',b'\t-9223372036854775809\t',1),A.replace(b'\t8\t',b'\t-0\t',1)]):
 d=mini('badrecord-'+str(i),'I:r',bad);p=run(['lookup',fa,'I:r',d],False);cases.append({'case':d.name,'error':p.stderr.decode().strip()})
for i,bad in enumerate(['','-1','01','1x','18446744073709551616',str(len(A)),'1']):
 d=mini('badoffset-'+str(i),'I:r',A,offset=bad);p=run(['lookup',fa,'I:r',d],False);cases.append({'case':d.name,'error':p.stderr.decode().strip()})
for i,keyrow in enumerate([b'I:r\t0\textra\n',b'I:r\t0',b'I:\t0\n',b'X:r\t0\n',b'I:\xff\t0\n',b'C:1:08\t0\n']):
 d=mini('badkey-'+str(i),'I:r',A,keyrow=keyrow);p=run(['lookup',fa,'I:r',d],False);cases.append({'case':d.name,'error':p.stderr.decode().strip()})
# Whole records and key lines may cross cache pages; valid UTF-8 remains opaque.
longid='雪'*23000
body=('A\t1\t8\tG\tA\t'+longid+'\t8\n').encode();d=mini('cross-page-valid','I:'+longid,body);assert run(['lookup',fa,'I:'+longid,d]).stdout
(d/'records').write_bytes(body.replace('雪'.encode(),b'\xf0\x80\x80',1));p=run(['lookup',fa,'I:'+longid,d],False);cases.append({'case':'cross-page-invalidUTF8','error':p.stderr.decode().strip()})
# The old frozen native reader must reject new data rather than misinterpret it.
query=ROOT/'q.tsv';query.write_text('1\t8\tr雪\tG\tA\t.3\t.2\tEAF\tgeneric\t\t\t\t10000\n')
old=BASE.parent/'mss-preservation-v7/compare-candidate-final';p=subprocess.run([str(old),'query',str(fa),str(query),str(refs[3])],capture_output=True);assert p.returncode!=0;cases.append({'case':'old-reader-rejects3','error':p.stderr.decode().strip()})
(BASE/'native-small-summary-v2.json').write_text(json.dumps({'passed':True,'corruption_cases':cases,'raw_S':6,'raw_A':5,'source_sha256':hashlib.sha256(sites.read_bytes()).hexdigest()},indent=2)+'\n')
print('PASS',len(cases),'corruption/compatibility checks; raw S6 A5; UTF8/cross-page/cap controls')
