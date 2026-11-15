#!/usr/bin/env python3
"""Derive auditable actual-association and separately synthetic gold panels."""
import collections,csv,json,pathlib
base=pathlib.Path('/private/tmp/gwas-harvest-stage2/assets')
out=base/'mss'/'functionality';out.mkdir(exist_ok=True)
def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def write(path,rows,fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
sites=read(base/'reference/ensembl104.chr22.reference.tsv')
aliases=collections.defaultdict(list)
for r in read(base/'reference/dbSNP151.chr22.merged_aliases.tsv'):aliases[r['current_id']].append(r['old_id'])
for r in sites:r['aliases']=';'.join(aliases[r['id']]) or '.'
write(out/'real-reference.tsv',sites,['chr','pos','ref','alt','id','aliases'])
lookup=collections.defaultdict(list)
for r in sites:lookup[r['chr'],r['pos'],r['ref'],r['alt']].append(r['id'])
association=read(base/'panel/panukb_biomarker30600_chr22_16050001_18626005.tsv')
fields=['BENCH_ROW','SNP','CHR','BP','A1','A2','BETA','FRQ','SE','LOG10P','reported_neglog10_pval_EUR','id_origin']
real=[]
for i,r in enumerate(association,1):
    ids=lookup[r['chr'],r['pos'],r['ref'],r['alt']]
    ident=ids[0] if len(ids)==1 else '.'
    real.append(dict(zip(fields,[str(i),ident,r['chr'],r['pos'],r['ref'],r['alt'],r['beta_EUR'],r['af_EUR'],r['se_EUR'],r['neglog10_pval_EUR'],r['neglog10_pval_EUR'],'derived_unique_reference_join' if len(ids)==1 else 'absent_or_nonunique_reference_join'])))
write(out/'real-input.tsv',real,fields)
write(out/'real-comparable-input.tsv',[r for r in real if r['SNP']!='.'],fields)
# Actual historical merged aliases: retain the same association measurement,
# replacing only its derived current rsID with a documented old rsID.
merged=[]
for r in real:
    if aliases[r['SNP']]:
        x=r.copy();x['BENCH_ROW']='merge-'+x['BENCH_ROW'];x['SNP']=aliases[r['SNP']][0];x['id_origin']='derived_dbSNP151_historical_merge';merged.append(x)
write(out/'real-merged-input.tsv',merged,fields)
(out/'synthetic.fa').write_text('>1\nCAAAAACGTACGT\n>2\nACGT\n')
(out/'synthetic-reference.tsv').write_text('chr\tpos\tref\talt\tid\taliases\n1\t8\tG\tA,T\trs900001\trs900101\n1\t4\tAA\tA\trs900002\t.\n1\t10\tA\tC\trs900003\trs900103\n1\t12\tG\tA\trs900004\t.\n2\t2\tC\tG\trs900005\t.\n2\t3\tG\tA\trs900006\t.\n2\t3\tG\tC\trs900007\t.\n')
challenge=[];gold=[]
def add(name,chrom,pos,ident,other,effect,beta,frq,status,canon=None,expected_effect=None,expected_eaf=None,p='.01',se='.05'):
    challenge.append(dict(BENCH_ROW=name,SNP=ident,CHR=chrom,BP=pos,A1=other,A2=effect,BETA=beta,FRQ=frq,SE=se,P=p))
    gold.append(dict(case=name,status=status,canonical=canon,effect=expected_effect,eaf=expected_eaf,source_p=p))
c=['1',10,'A','C'];indel=['1',1,'CA','C'];multi=['1',8,'G','T']
add('snp_match','1',10,'rs900003','A','C','.125','.25','match',c,.125,.25)
add('snp_swap','1',10,'rs900003','C','A','-.125','.75','swap',c,.125,.25)
add('snp_complement','1',10,'rs900003','T','G','.125','.25','complement',c,.125,.25)
add('snp_swapped_complement','1',10,'rs900003','G','T','-.125','.75','swapped_complement',c,.125,.25)
add('multiallelic_ref_alt','1',8,'rs900001','G','T','.125','.25','match',multi,.125,.25)
add('multiallelic_swap','1',8,'rs900001','T','G','-.125','.75','swap',multi,.125,None)
add('multiallelic_alt_alt','1',8,'rs900001','A','T','.125','.25','alt_alt_contrast')
add('indel_left_align','1',4,'rs900002','AA','A','.125','.25','match',indel,.125,.25)
add('indel_swap','1',4,'rs900002','A','AA','-.125','.75','swap',indel,.125,.25)
add('indel_shifted','1',5,'rs900002','AA','A','.125','.25','match',indel,.125,.25)
add('synthetic_alias_swap','1',10,'rs900103','C','A','-.125','.75','swap',c,.125,.25)
add('unknown_id_known_coordinate','1',10,'rs999999991','A','C','.125','.25','match',c,.125,.25)
add('unknown_id_unknown_coordinate','1',9,'rs999999992','T','C','.125','.25','unmapped')
add('missing_beta','1',10,'rs900003','A','C','NA','.25','match',c,None,.25)
add('invalid_beta','1',10,'rs900003','A','C','bad','.25','match',c,None,.25)
add('missing_eaf','1',10,'rs900003','A','C','.125','NA','match',c,.125,None)
add('invalid_eaf','1',10,'rs900003','A','C','.125','1.2','match',c,.125,None)
add('palindrome','2',2,'rs900005','C','G','.125','.25','ambiguous')
add('p_underflow','1',10,'rs900003','A','C','.125','.25','match',c,.125,.25,p='1e-400')
add('missing_optional_se','1',10,'rs900003','A','C','.125','.25','match',c,.125,.25,se='NA')
add('id_coordinate_conflict','1',10,'rs900004','A','C','.125','.25','coordinate_id_conflict')
add('unsupported_symbolic','1',10,'rs900003','A','<DEL>','.125','.25','unsupported_allele')
write(out/'synthetic-input.tsv',challenge,list(challenge[0]))
(out/'synthetic-gold.json').write_text(json.dumps(gold,indent=2)+'\n')
receipt=dict(real_source_rows=len(real),real_unique_reference_ids=sum(r['SNP']!='.' for r in real),
             real_merged_alias_rows=len(merged),synthetic_challenge_rows=len(challenge),
             negative_log10_semantics='Actual neglog10_pval_EUR header and official per-phenotype-files docs; raw token retained; no source-header rewriting',
             effect_semantics='PanUKB beta_EUR/eaf_EUR declared relative to source ALT; synthetic challenge truth hand-checked',
             historical_alias_release='NCBI dbSNP151 RsMergeArch (2018), transitive aliases; not current alias database')
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
