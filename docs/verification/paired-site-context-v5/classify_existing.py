"""Classify the already-completed bounded old/new outputs, including cap errors."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

ROOT=Path(__file__).parent/'resolved-comparison'
transitions=Counter();error_transitions=Counter();query_classes=Counter();resolved=unchanged=errors=0;rows_seen=0;examples=[]
with (ROOT/'baseline.tsv').open() as old,(ROOT/'candidate.tsv').open() as new,(ROOT/'queries.tsv').open() as queries:
 for ordinal,(left,right,query) in enumerate(itertools.zip_longest(old,new,queries),1):
  assert left is not None and right is not None and query is not None
  a=left.rstrip('\n').split('\t');b=right.rstrip('\n').split('\t');rows_seen+=1
  q=query.rstrip('\n').split('\t')
  query_classes['declared' if q[8]=='declared' else ('bounded_cap' if q[12]!='10000' else 'generic')]+=1
  is_resolved=lambda r:r[0]=='ok' and (r[2]!='' or r[3]!='0' or r[4]!='' or r[5]!='')
  if is_resolved(a) or is_resolved(b):
   resolved+=1
   assert left==right,('resolved scientific output changed',ordinal,left,right)
  if left==right:
   unchanged+=1;errors+=int(a[0]=='error');continue
  if a[0]!='ok' or b[0]!='ok':
   error_transitions[a[0]+' -> '+b[0]]+=1
   assert a[0]=='error' and b[0]=='ok',('new or changed error',ordinal,left,right)
   assert not is_resolved(a) and not is_resolved(b)
   if len(examples)<30:examples.append(dict(row=ordinal,query=query.rstrip('\n'),old=left.rstrip('\n'),new=right.rstrip('\n')))
   continue
  assert a[2:]==b[2:],('non-status output changed',ordinal,left,right)
  assert {a[1],b[1]}=={'alt_alt_contrast','allele_mismatch'},('unexpected status change',ordinal,left,right)
  transitions[a[1]+' -> '+b[1]]+=1
  if len(examples)<20:examples.append(dict(row=ordinal,query=query.rstrip('\n'),old=a[1],new=b[1]))
result=dict(scope='Bounded synthetic/legacy-fixture functional comparison, not performance or full-genome validation',
 source_frozen_commit='5cf1e23c907ad9c923141faf0728964951be7c66',candidate_commit='ba91c55',
 queries=rows_seen,query_classes=dict(query_classes),resolved_scientific_rows_exact=resolved,complete_rows_exact=unchanged,
 identical_errors=errors,status_transitions=dict(transitions),error_transitions=dict(error_transitions),
 every_resolved_result_exact=True,all_successful_changed_rows_have_nonstatus_fields_exact=True,
 examples=examples,
 sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['reference.fa','sites.tsv','queries.tsv','baseline.tsv','candidate.tsv']})
(ROOT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['examples','sha256']},indent=2))
