#!/usr/bin/env python3
"""Full gold/oracle verification, with explicit comparable denominators."""
import csv,hashlib,json,math,pathlib,subprocess
base=pathlib.Path('/private/tmp/gwas-harvest-stage2')
assets=base/'assets/mss/functionality';results=base/'results/mss/functionality'
def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def md5(path):return hashlib.md5(path.read_bytes()).hexdigest()
def number(s):
    try:
        value=float(s)
        return value if math.isfinite(value) else None
    except (ValueError,TypeError):return None
def canonical(row,new):
    names=['chr','pos','ref','alt'] if new else ['CHR','BP','A1','A2']
    values=[row[n] for n in names]
    if values[1] in ('NA','.',''):return None
    return [values[0],int(float(values[1])),values[2],values[3]]
summary={}
gold={r['case']:r for r in json.loads((assets/'synthetic-gold.json').read_text())}
inputs=read(assets/'synthetic-input.tsv')
for mode in ('new','legacy'):
    directory=results/f'{mode}-synthetic'
    if mode=='legacy':directory=results/'legacy-synthetic-v2'
    path=directory/'edition/normalized.tsv' if mode=='new' else directory/'legacy.tsv'
    if not path.exists():
        summary[f'{mode}_synthetic']={'pipeline_failed':True};continue
    output=read(path)
    keyed={(inputs[int(r['source_row'])-1]['BENCH_ROW'] if mode=='new' else r['BENCH_ROW']):r for r in output}
    count=dict(source_rows=len(inputs),retained_rows=len(output),effect_gold_eligible=sum(g['effect'] is not None for g in gold.values()),
               signed_effect_correct=0,signed_effect_and_canonical_correct=0,
               canonical_gold_eligible=sum(g['canonical'] is not None for g in gold.values()),
               canonical_correct=0,eaf_gold_eligible=sum(g['eaf'] is not None for g in gold.values()),eaf_correct=0,eaf_and_canonical_correct=0)
    details=[]
    for name,g in gold.items():
        r=keyed.get(name)
        effect=number(r['effect' if mode=='new' else 'BETA']) if r else None
        eaf=number(r['eaf'] if mode=='new' else r.get('FRQ',r.get('MAJOR_ALLELE_FRQ'))) if r else None
        canon=canonical(r,mode=='new') if r else None
        if g['effect'] is not None:
            count['signed_effect_correct']+=bool(r and effect==g['effect'])
            count['signed_effect_and_canonical_correct']+=bool(r and effect==g['effect'] and canon==g['canonical'])
        if g['eaf'] is not None:
            count['eaf_correct']+=bool(r and eaf==g['eaf'])
            count['eaf_and_canonical_correct']+=bool(r and eaf==g['eaf'] and canon==g['canonical'])
        if g['canonical'] is not None:count['canonical_correct']+=bool(r and canon==g['canonical'])
        details.append(dict(case=name,retained=r is not None,status=r.get('status') if r else None,canonical=canon,effect=effect,eaf=eaf))
        if mode=='new':
            assert r and r['status']==g['status'],(name,r,g)
            assert canon==g['canonical'] and effect==g['effect'] and eaf==g['eaf'],(name,r,g)
    if mode=='new':
        assert md5(assets/'synthetic-input.tsv')==md5(directory/'edition/source-artifact')
        assert (assets/'synthetic-input.tsv').read_bytes()==(directory/'edition/source.tsv').read_bytes()
        count['exact_source_preservation']=True
    count['underflow_source_token_retained']=mode=='new' or bool(keyed.get('p_underflow') and keyed['p_underflow']['P']=='1e-400')
    summary[f'{mode}_synthetic']=count
    (directory/'gold-details.json').write_text(json.dumps(details,indent=2)+'\n')
# Independent bcftools oracle for every actual source REF/ALT, not merely those
# retained/normalized by either path. No P or effect conversion in the oracle.
real=read(assets/'real-input.tsv')
vcf=results/'real-source.vcf'
with vcf.open('w') as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=22,length=51304566>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for r in real:f.write(f"{r['CHR']}\t{r['BP']}\t{r['BENCH_ROW']}\t{r['A1']}\t{r['A2']}\t.\t.\t.\n")
oracle=subprocess.run(['bcftools','norm','-f',str(base/'assets/reference/GRCh37.ensembl104.chr22.fa'),'-c','s','-Ov',str(vcf)],capture_output=True,text=True)
(results/'real-oracle.log').write_text(oracle.stderr)
assert oracle.returncode==0,oracle.stderr
oracle_rows={v[2]:[v[0],int(v[1]),v[3],v[4]] for line in oracle.stdout.splitlines() if not line.startswith('#') for v in [line.split('\t')]}
(results/'real-normalized.vcf').write_text(oracle.stdout)
for panel in ('real','real-comparable','real-merged'):
    source=read(assets/f'{panel}-input.tsv');source_by_id={r['BENCH_ROW']:r for r in source}
    for mode in ('new','legacy'):
        if mode=='legacy' and panel=='real':continue
        directory=results/f'{mode}-{panel}'
        path=directory/'edition/normalized.tsv' if mode=='new' else directory/'legacy.tsv'
        if not path.exists():
            summary[f'{mode}_{panel}']={'pipeline_failed':True};continue
        rows=read(path)
        counts=dict(source_rows=len(source),retained_rows=len(rows),canonical_calls=0,canonical_oracle_agree=0,
                    signed_effect_calls=0,signed_effect_source_agree=0,eaf_calls=0,eaf_source_agree=0,statuses={},discordant_cases=[])
        for r in rows:
            original=source[int(r['source_row'])-1] if mode=='new' else source_by_id[r['BENCH_ROW']]
            source_case=original['BENCH_ROW'].removeprefix('merge-')
            canon=canonical(r,mode=='new')
            if canon:
                counts['canonical_calls']+=1
                counts['canonical_oracle_agree']+=canon==oracle_rows[source_case]
                if canon!=oracle_rows[source_case]:counts['discordant_cases'].append(original['BENCH_ROW'])
            effect=number(r['effect' if mode=='new' else 'BETA']);eaf=number(r['eaf'] if mode=='new' else r.get('FRQ',r.get('MAJOR_ALLELE_FRQ')))
            # Correct the six single-base source REF/ALT inversions flagged by
            # bcftools, keeping the original association measurement unchanged.
            swapped=(len(original['A1'])==len(original['A2'])==1 and oracle_rows[source_case][2]==original['A2'])
            expected_effect=number(original['BETA'])
            if swapped and expected_effect is not None:expected_effect=-expected_effect
            expected_eaf=number(original['FRQ'])
            if swapped and expected_eaf is not None:expected_eaf=1-expected_eaf
            if effect is not None:
                counts['signed_effect_calls']+=1;counts['signed_effect_source_agree']+=effect==expected_effect and canon==oracle_rows[source_case]
            if eaf is not None:
                counts['eaf_calls']+=1;counts['eaf_source_agree']+=eaf==expected_eaf and canon==oracle_rows[source_case]
            if mode=='new':counts['statuses'][r['status']]=counts['statuses'].get(r['status'],0)+1
        if mode=='new':
            counts['exact_source_preservation']=md5(assets/f'{panel}-input.tsv')==md5(directory/'edition/source-artifact')
            assert counts['exact_source_preservation'] and counts['retained_rows']==counts['source_rows']
            counts['all_scientific_calls_oracle_agree']=(counts['canonical_calls']==counts['canonical_oracle_agree'] and counts['signed_effect_calls']==counts['signed_effect_source_agree'] and counts['eaf_calls']==counts['eaf_source_agree'])
        summary[f'{mode}_{panel}']=counts
summary['bcftools_version']=subprocess.check_output(['bcftools','--version'],text=True).splitlines()[0]
summary['real_oracle_source_rows']=len(oracle_rows)
summary['real_denominator_note']='Full actual source14428; comparable subset13419unique reference-joined IDs, same derived input both paths; merged3178old IDs from historical dbSNP151 map, same association measurements.'
(results/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
