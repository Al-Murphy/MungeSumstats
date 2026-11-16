from pathlib import Path
import itertools,json,hashlib
B=Path(__file__).resolve().parent
pairs=[('actual-prefix-1M',B/'prefix-1M-v1/format2/index',B/'prefix-1M-v1/format3/index')]
pairs += [('ALT'+str(n),B/f'highALT-whole-query-v1/ALT{n}/format2/index',B/f'highALT-whole-query-v1/ALT{n}/format3/index') for n in (256,2048,10000)]
out=[]
for label,two,three in pairs:
 counts={'S':0,'A':0};context=None
 with open(two/'records','rb') as a,open(three/'records','rb') as z:
  for row,(old,new) in enumerate(itertools.zip_longest(a,z),1):
   assert old is not None and new is not None,(label,row,'rowcount')
   v=old.rstrip(b'\n').split(b'\t');w=new.rstrip(b'\n').split(b'\t');assert len(v)==8,(label,row)
   if v[3]==b'':
    expected=[b'S',v[0],v[5],v[6],v[7]];assert v[1]==v[5] and v[2]==v[6]
    context=v;counts['S']+=1
   else:
    expected=[b'A',*v[:6]];counts['A']+=1
    assert context is not None and v[0]==context[0] and v[5:8]==context[5:8],(label,row,'source pairing')
   assert w==expected,(label,row,'typedfield mismatch')
 for name in ('source.tsv','fasta.fai','rejected.tsv'):
  assert hashlib.file_digest(open(two/name,'rb'),'sha256').hexdigest()==hashlib.file_digest(open(three/name,'rb'),'sha256').hexdigest(),(label,name)
 if label.startswith('ALT'):assert counts=={'S':1,'A':int(label[3:])}
 else:assert counts['S']==1000000
 out.append(dict(label=label,counts=counts,all_typed_fields_and_occurrence_order_exact=True,original_REF_and_full_ALT_pairing_exact=True,source_FAI_rejected_exact=True))
(B/'record-parity-final.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
