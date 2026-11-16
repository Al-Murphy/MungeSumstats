"""Streaming independent SNP-panel oracle; no source/output row is sampled.

This performance workload requires unique biallelic non-palindromic SNP loci
and unique ID events. It refuses broader contexts rather than interpreting
multiallelic or complex contexts through a narrow performance oracle. The
separately reviewed scientific suites cover those cases.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math
from pathlib import Path

FIELDS=['status','chr','pos','ref','alt','id','effect','eaf','swapped','complemented',
        'candidate_count','entity_type','entity_id','normalized_build','source_row',
        'grch38_chr','grch38_pos','grch38_ref','grch38_alt']
RESOLVED={'match','swap','complement','swapped_complement','sequence_only'}
NA='NA'
COMP=str.maketrans('ACGT','TGCA')

def sha256(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()

def numeric(value):
 return None if value==NA else float(value)

class Fasta:
 def __init__(self,path,fai):
  self.file=Path(path).open('rb');self.contigs={};self.pages=collections.OrderedDict()
  for line in Path(fai).read_text().splitlines():
   c,n,o,b,w=line.split('\t');self.contigs[c]=tuple(map(int,(n,o,b,w)))
 def close(self):
  if not self.file.closed:self.file.close()
 def __del__(self):
  file=getattr(self,'file',None)
  if file is not None and not file.closed:file.close()
 def base(self,chrom,pos):
  n,offset,bases,width=self.contigs[chrom]
  if not 1<=pos<=n:raise AssertionError('Reference coordinate out of FASTA')
  absolute=offset+(pos-1)//bases*width+(pos-1)%bases;page=absolute//65536
  if page not in self.pages:
   self.file.seek(page*65536);self.pages[page]=self.file.read(65536)
   if len(self.pages)>64:self.pages.popitem(last=False)
  self.pages.move_to_end(page)
  return chr(self.pages[page][absolute%65536]).upper()

class Gold:
 def __init__(self,source,fasta,fai,build):
  self.by_id={};self.loci={};self.rows=0;self.panel=[];self.build=build;self.cache={};self.fa=Fasta(fasta,fai)
  with Path(source).open() as f:
   reader=csv.DictReader(f,delimiter='\t')
   if reader.fieldnames!=['chr','pos','ref','alt','id','aliases']:raise AssertionError('Reference source schema changed')
   for r in reader:
    self.rows+=1;chrom,pos,ref=r['chr'],int(r['pos']),r['ref'].upper()
    alts=r['alt'].upper().split(',')
    if len(ref)!=1 or ref not in 'ACGT' or any(len(a)!=1 or a not in 'ACGT' or a==ref for a in alts):
     raise AssertionError('Whole-file oracle requires concrete SNP reference; complex/symbolic context is not assumed biallelic')
    if self.fa.base(chrom,pos)!=ref:raise AssertionError('Reference REF/FASTA conflict')
    locus=(chrom,pos);site=self.loci.setdefault(locus,{'ref':ref,'alts':set()})
    if site['ref']!=ref:raise AssertionError('Incoherent original SNP REF context')
    site['alts'].update(alts)
    if len(alts)==1 and ref+alts[0] not in ('AT','TA','CG','GC') and r['id'].startswith('rs'):
     self.panel.append((r['id'],chrom,pos,ref,alts[0]))
    names=[r['id']]+r['aliases'].split(';')
    for alt in alts:
     event=(chrom,pos,ref,alt)
     for name in names:
      if name not in ('','.','NA'):
       ids=self.by_id.setdefault(name,{}).setdefault(event,set())
       if r['id'] not in ('','.','NA'):ids.add(r['id'])
  # This performance workload deliberately requires unique biallelic loci/IDs.
  # Refuse a broader context rather than mislabeling a narrow oracle general.
  if len(self.panel)!=self.rows or len(self.loci)!=self.rows or len({p[0] for p in self.panel})!=self.rows or any(len(v['alts'])!=1 for v in self.loci.values()) or any(len(v)!=1 for v in self.by_id.values()):
   raise AssertionError('Performance oracle requires unique biallelic non-palindromic SNP loci and identifiers')
  self.inventory=dict(source_rows=self.rows,loci=len(self.loci),ids=len(self.by_id),eligible_panel_records=len(self.panel),
      multiallelic_loci=sum(len(v['alts'])>1 for v in self.loci.values()),
      quadallelic_loci=sum(len(v['alts'])==3 for v in self.loci.values()),
      ids_with_multiple_events=sum(len(v)>1 for v in self.by_id.values()),
      reference_scope='All concrete SNP REF/ALT records included; no assumption of one ALT or one locus per ID')
 def shape(self,r,version):
  key=(r['CHR'],int(r['BP']),r['SNP'],r['A1'].upper(),r['A2'].upper(),version)
  if key in self.cache:return self.cache[key]
  chrom,pos,identifier,other,effect=key[:5]
  if len(other)!=1 or len(effect)!=1 or other not in 'ACGT' or effect not in 'ACGT' or other==effect:
   raise AssertionError('Pinned benchmark input is not its disclosed SNP workload')
  site=self.loci[(chrom,pos)];id_candidates=self.by_id[identifier]
  candidates={event:ids for event,ids in id_candidates.items() if event[:2]==(chrom,pos)}
  if not candidates:raise AssertionError('Pinned input coordinate/ID contract mismatch')
  hits=[];altalt=False
  for comp in (False,True):
   a,b=(other.translate(COMP),effect.translate(COMP)) if comp else (other,effect)
   if a in site['alts'] and b in site['alts']:altalt=True
   for event,ids in candidates.items():
    if (a,b)==event[2:]:hits.append((event,False,comp,ids))
    if (b,a)==event[2:]:hits.append((event,True,comp,ids))
  out={field:NA for field in FIELDS}
  out.update(candidate_count=str(len(candidates)),entity_type='small_variant',entity_id=identifier,normalized_build=self.build)
  if len(hits)>1 or (version=='mss7' and hits and altalt):out['status']='ambiguous'
  elif not hits:out['status']='alt_alt_contrast' if altalt else 'allele_mismatch'
  else:
   event,swap,comp,ids=hits[0]
   out.update(status=('swapped_complement' if swap else 'complement') if comp else ('swap' if swap else 'match'),
      chr=event[0],pos=str(event[1]),ref=event[2],alt=event[3],id=next(iter(ids)) if len(ids)==1 else NA,
      swapped='TRUE' if swap else 'FALSE',complemented='TRUE' if comp else 'FALSE',
      effect_factor=-1 if swap else 1,eaf_rule='complement' if swap and len(site['alts'])==1 else ('missing' if swap else 'same'))
  if len(self.cache)>=4*self.rows:raise AssertionError('Gold cache exceeded its bounded panel domain')
  self.cache[key]=out;return out
 def expected(self,r,version,ordinal):
  out=self.shape(r,version).copy();out['source_row']=str(ordinal)
  beta,eaf=float(r['BETA']),float(r['FRQ'])
  if out['status'] in RESOLVED:
   out['effect']=beta*out.pop('effect_factor')
   rule=out.pop('eaf_rule');out['eaf']=None if rule=='missing' else (1-eaf if rule=='complement' else eaf)
  else:out['effect']=None;out['eaf']=None
  return out

def check_row(row,want,ordinal):
 if list(row)!=FIELDS or any(type(value) is not str for value in row.values()):raise AssertionError((ordinal,'Unexpected normalized row width'))
 for field in FIELDS:
  if field in ('effect','eaf'):
   if numeric(row[field])!=want[field]:raise AssertionError((ordinal,field,row[field],want[field]))
  elif field in ('pos','source_row','candidate_count') and want[field]!=NA:
   if numeric(row[field])!=int(want[field]):raise AssertionError((ordinal,field,row[field],want[field]))
  elif row[field]!=want[field]:raise AssertionError((ordinal,field,row[field],want[field]))

def verify(profile,stage,editions,receipt):
 handles={};g=None
 try:
  p=json.loads(Path(profile).read_text());s=p['stages'][stage];ref=p['reference'];g=Gold(ref['source'],ref['fasta'],ref['fai'],ref['build'])
  versions=list(editions);counters={v:collections.Counter() for v in versions};statuses={v:collections.Counter() for v in versions};deltas=collections.Counter()
  handles={};readers={}
  for version,folder in editions.items():
   folder=Path(folder)
   for name in ['source.tsv','source-artifact']:
    path=folder/name
    if not path.is_file() or path.is_symlink() or sha256(path)!=s['input']['sha256']:raise AssertionError((version,name,'raw source bytes differ'))
   handles[version]=(folder/'normalized.tsv').open();readers[version]=csv.DictReader(handles[version],delimiter='\t')
   if readers[version].fieldnames!=FIELDS:raise AssertionError((version,'normalized output schema changed'))
  with Path(s['input']['path']).open() as src:
   source=csv.DictReader(src,delimiter='\t')
   if source.fieldnames!=['SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P']:raise AssertionError('Input schema changed')
   n=0
   for n,r in enumerate(source,1):
    panel=g.panel[(n-1)%len(g.panel)];swap=bool(((n-1)//len(g.panel))%2)
    expected_source=[panel[0],panel[1],str(panel[2]),panel[4] if swap else panel[3],panel[3] if swap else panel[4],
                     '-.125' if swap else '.125','.75' if swap else '.25','.05','.01']
    if list(r.values())!=expected_source:raise AssertionError((n,'Source deviates from pinned replicated-panel generation contract'))
    outputs={};expected={}
    for version in versions:
     row=next(readers[version],None)
     if row is None:raise AssertionError((version,'output shorter than source',n))
     want=g.expected(r,version,n);check_row(row,want,n);outputs[version]=row;expected[version]=want
     count=counters[version];count['source_rows']+=1;count['normalized_rows']+=1;count['all_19_fields_correct_rows']+=1
     count['expected_resolved_rows']+=want['status'] in RESOLVED;count['observed_resolved_rows']+=row['status'] in RESOLVED
     count['effect_nonmissing_denominator']+=want['effect'] is not None;count['eaf_nonmissing_denominator']+=want['eaf'] is not None
     count['effect_missing_correct_rows']+=want['effect'] is None;count['eaf_missing_correct_rows']+=want['eaf'] is None
     statuses[version][row['status']]+=1
    if set(versions)=={'mss4','mss7'}:
     if outputs['mss4']==outputs['mss7']:deltas['all_19_fields_byte_exact_rows']+=1
     else:
      if expected['mss4']['status'] not in RESOLVED or expected['mss7']['status']!='ambiguous':raise AssertionError((n,'Unaccounted 4/7 delta'))
      deltas['gold_confirmed_resolved_to_ambiguous_rows']+=1
   if n!=s['rows']:raise AssertionError(('Source rows',n,s['rows']))
  for version in versions:
   if next(readers[version],None) is not None:raise AssertionError((version,'output longer than source'))
   handles[version].close()
  result=dict(schema=1,stage=stage,rows_expected=s['rows'],passed=True,
    all_rows_examined=True,source_bytes_sha256=s['input']['sha256'],source_copies_exact_each=2,
    reference_inventory=g.inventory,versions={v:dict(counts=dict(counters[v]),statuses=dict(statuses[v])) for v in versions},
    deltas=dict(deltas),unaccounted_deltas=0,
    scope='Repeated real chr22 unique biallelic non-palindromic SNP panel; every source and normalized row/all19fields verified; indels/multiallelic contexts are validated in separate scientific suites, not timed here')
  with Path(receipt).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
  print(json.dumps(result));return result
 finally:
  for handle in handles.values():handle.close()
  if g is not None:g.fa.close()

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--profile',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--edition',action='append',required=True,help='mss4=/new/edition or mss7=/new/edition');ap.add_argument('--receipt',required=True);a=ap.parse_args()
 verify(a.profile,a.stage,dict(item.split('=',1) for item in a.edition),a.receipt)
