"""Bounded functional matrix; not a speed, throughput or genome benchmark."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess

ROOT=Path(__file__).parent/'resolved-comparison'
ROOT.mkdir()
TOOLS=ROOT.parent
RNG=random.Random(651902)
fasta=ROOT/'reference.fa'
fasta.write_text('>1\nCAAAAACGTACGT\n>2\nACGT\n>C\nCAAAC\n>G\nGGGGG\n')
sites=ROOT/'sites.tsv'
rows=[
 '1\t8\tG\tA,T\trs8\told8;merged8',
 '1\t4\tAA\tA\trsDel\t.',
 '1\t10\tA\tC\trsDup\tshared',
 '1\t12\tG\tA\trsDup\tshared',
 '2\t2\tC\tG\tpalindrome\t.',
 '2\t1\tA\tC\tduplicate\t.',
 '1\t10\tA\tC\tduplicate\t.',
 '2\t3\tG\tA\tsplitThird\t.',
 '2\t3\tG\tC\tsplitThird2\t.',
 'C\t1\tC\tT\tpad\t.',
 'C\t2\tAA\tAAA\tctx\told_ctx',
 'C\t4\tA\tAA\tctx\told_ctx',
 'G\t3\tG\tA\tspan\told_span',
 'G\t3\tGG\tCC\tspan\told_span',
 'G\t2\tG\tA,T,<DEL>\tsymbolic\told_symbolic',
 'G\t4\tG\tC\tvalid\t.',
 'G\t4\tN\tA,T\t.\t.',
]
sites.write_text('chr\tpos\tref\talt\tid\taliases\n'+'\n'.join(rows)+'\n')
index=ROOT/'frozen-format2-index'
subprocess.run([str(TOOLS/'compare-baseline'),'prepare',str(fasta),str(sites),str(index)],check=True)
alleles=['A','C','G','T','AA','AAA','AC','GG','<DEL>','N','.','*']
locations=[('1',p) for p in [0,1,2,4,6,8,10,12,13,14]]
locations += [('2',p) for p in [0,1,2,3,4,5]]
locations += [('C',p) for p in [0,1,2,4,5]]
locations += [('G',p) for p in [0,1,2,3,4,5]]
locations += [('.',0),('unknown',3)]
ids=['.','rs8','old8','merged8','rsDel','rsDup','shared','duplicate','palindrome','ctx','old_ctx','old_span','symbolic','valid']
queries=ROOT/'queries.tsv'
counts=Counter()
with queries.open('w') as out:
 for (chrom,pos),identifier,other,effect,frequency in itertools.product(locations,ids,alleles,alleles,['EAF','MAF','unknown']):
  values=[chrom,str(pos),identifier,other,effect,'-0.25','0.3',frequency,'generic','.','.','unknown','10000']
  out.write('\t'.join(values)+'\n');counts['generic']+=1
 for _ in range(12000):
  row=RNG.choice(rows).split('\t')
  source_ref=row[2];source_alt=RNG.choice(row[3].split(','))
  other,effect=(source_ref,source_alt) if RNG.randrange(2) else (source_alt,source_ref)
  values=[row[0],row[1],RNG.choice([row[4],'.',row[5].split(';')[0]]),other,effect,
          RNG.choice(['0','-0.25','2','nan']),RNG.choice(['0','0.3','1','nan']),
          RNG.choice(['EAF','MAF','unknown']),'declared',source_ref,source_alt,
          RNG.choice(['+','-','unknown']),'10000']
  out.write('\t'.join(values)+'\n');counts['declared']+=1
 for _ in range(3000):
  chrom,pos=RNG.choice(locations)
  values=[chrom,str(pos),RNG.choice(ids),RNG.choice(alleles),RNG.choice(alleles),'0.5','0.3',
          'EAF','generic','.','.','unknown',str(RNG.choice([1,2,3]))]
  out.write('\t'.join(values)+'\n');counts['bounded_cap']+=1
for mode in ('baseline','candidate'):
 with (ROOT/(mode+'.tsv')).open('wb') as out:
  subprocess.run([str(TOOLS/('compare-'+mode)),'query',str(fasta),str(queries),str(index)],stdout=out,check=True)
# The classifier asserts exact resolved payloads and rejects new errors while
# allowing documented diagnostic corrections and avoided foreign-context cap errors.
import runpy
runpy.run_path(str(TOOLS/'classify_existing.py'), run_name='__main__')
