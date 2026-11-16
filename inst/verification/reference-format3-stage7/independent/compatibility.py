import collections,hashlib,json,pathlib,shutil,subprocess
ROOT=pathlib.Path(__file__).resolve().parent
OLD=ROOT.parent/'review-mss-preservation-v7-stage7'
D=ROOT/'compatibility-v2';D.mkdir()
PAYLOAD=['core-version','fasta.fai','keys','offsets','records','rejected.tsv']
for tag,inc in [('7',OLD/'source/src'),('3',ROOT/'source/src')]:
 p=subprocess.run(['clang++','-std=c++17','-O2']+(['-DLEGACY_SCI7'] if tag=='7' else [])+['-I'+str(inc),str(ROOT/'writer.cpp'),'-o',str(D/('writer'+tag))],capture_output=True,text=True)
 assert p.returncode==0,p.stderr
fasta=D/'tiny.fa';fasta.write_bytes((OLD/'writer-integrity/tiny.fa').read_bytes())
counter=[]
for order in range(6):
 source=D/f'sites-{order}.tsv'
 lines=(OLD/f'writer-integrity/sites-{order}.tsv').read_text().splitlines()
 fields=lines[1].split('\t');fields[5]+=';'+';'.join('alias%03d'%i for i in range(300));lines[1]='\t'.join(fields)
 source.write_text('\n'.join(lines)+'\n')
 for batch in [1,2,3,7,31,100]:
  outputs={}
  for code,fmt in [('7',2),('3',2),('3',3)]:
   index=D/f'index-{order}-{batch}-{code}-{fmt}'
   p=subprocess.run([str(D/('writer'+code)),str(fasta),str(source),str(index),str(fmt),str(batch)],capture_output=True)
   assert p.returncode==0,p.stderr
   outputs[(code,fmt)]=index
  for name in PAYLOAD:
   assert (outputs['7',2]/name).read_bytes()==(outputs['3',2]/name).read_bytes(),(order,batch,name)
  def keys(index):return collections.Counter(l.rsplit(b'\t',1)[0] for l in (index/'keys').read_bytes().splitlines())
  assert keys(outputs['7',2])==keys(outputs['3',3]),(order,batch,'raw key multiplicity')
  for name in ['fasta.fai','rejected.tsv']:
   assert (outputs['7',2]/name).read_bytes()==(outputs['3',3]/name).read_bytes()
  counter.append({'order':order,'batch':batch,'default2_six_payloads_exact':True,'typed3_raw_key_multiplicity_exact':True})
work=[]
W=OLD/'work-order-qualification'
for fmt in [2,3]:
 for order in [0,1]:
  index=D/f'work-{fmt}-{order}'
  p=subprocess.run([str(D/'writer3'),str(W/'tiny.fa'),str(W/f'sites-{order}.tsv'),str(index),str(fmt),'7'],capture_output=True);assert p.returncode==0,p.stderr
  for cap in [5,200]:
   p=subprocess.run([str(ROOT/'core'),'query',str(W/'tiny.fa'),str(W/f'query-{order}-{cap}.tsv'),str(index)],capture_output=True,text=True);assert p.returncode==0,p.stderr
   output=p.stdout; (D/f'work-{fmt}-{order}-{cap}.tsv').write_text(output)
   expected='error\tReference context closure limit exceeded' if (fmt,order,cap)==(2,1,5) else 'ok\tambiguous'
   assert output.startswith(expected),(fmt,order,cap,output)
   work.append({'format':fmt,'order':order,'cap':cap,'output':output.rstrip()})
parity=[]
for name in ['plain-oracle-output.tsv','left-oracle-output.tsv','right-oracle-output.tsv','both-oracle-output.tsv','haplotype-output.tsv']:
 a=(ROOT/'format2'/name).read_bytes();b=(ROOT/'format3'/name).read_bytes();assert a==b,name
 parity.append({'name':name,'rows':a.count(b'\n'),'sha256':hashlib.sha256(a).hexdigest()})
log=(ROOT/'science.log').read_text(); sections={}
for label in ['SCIENCE FORMAT 2','SCIENCE FORMAT 3','NEW API TESTS']:
 chunk=log.split('INDEPENDENT '+label+'\n',1)[1].split('INDEPENDENT ',1)[0]
 import re
 counts=[int(x) for x in re.findall(r'Test passed with (\d+) successes',chunk)]
 singles=len(re.findall(r'Test passed with 1 success',chunk))
 assert not re.search(r'Failure|Error|Test failed',chunk)
 sections[label]={'blocks':len(counts)+singles,'assertions':sum(counts)+singles}
assert sections=={'SCIENCE FORMAT 2':{'blocks':65,'assertions':1170},'SCIENCE FORMAT 3':{'blocks':65,'assertions':1170},'NEW API TESTS':{'blocks':7,'assertions':65}},sections
summary={'writer_matrix':counter,'work_order':work,'independent_oracle_full_payload_parity':parity,'R_science_sections':sections}
(ROOT/'compatibility-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'writer_matrix':len(counter),'work_order_cases':len(work),'oracle_exact_rows':sum(x['rows'] for x in parity),'R_science_sections':sections}))
