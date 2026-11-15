#!/usr/bin/env python3
"""Before/after real-source scientific audit with independent FASTA and oracle."""
import collections,csv,hashlib,json,math,pathlib
base=pathlib.Path('/private/tmp/gwas-harvest-stage2')
assets=base/'assets/mss/functionality';results=base/'results/mss'
def read(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def number(s):
    try:
        v=float(s)
        return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
fasta=''.join(line.strip() for line in (base/'assets/reference/GRCh37.ensembl104.chr22.fa').open() if not line.startswith('>'))
oracle={v[2]:[v[0],int(v[1]),v[3],v[4]] for line in (results/'functionality/real-normalized.vcf').read_text().splitlines() if not line.startswith('#') for v in [line.split('\t')]}
summary={}
for panel in ('real','real-comparable','real-merged'):
    source=read(assets/f'{panel}-input.tsv')
    valid=[fasta[int(r['BP'])-1:int(r['BP'])-1+len(r['A1'])]==r['A1'] for r in source]
    for version,mode in (('functionality','new'),('functionality-v2','new'),('functionality-v2','trusted'),
                         ('functionality-v3','new'),('functionality-v3','trusted')):
        directory=results/version/f'{mode}-{panel}'/'edition'
        if not directory.exists():continue
        rows=read(directory/'normalized.tsv')
        counts=dict(source_rows=len(source),retained_rows=len(rows),forward_ref_valid=sum(valid),forward_ref_mismatch=len(source)-sum(valid),
                    effect_eligible=sum(ok and number(r['BETA']) is not None for r,ok in zip(source,valid)),
                    eaf_eligible=sum(ok and number(r['FRQ']) is not None for r,ok in zip(source,valid)),
                    eligible_unresolved=0,canonical_correct=0,canonical_wrong=0,signed_effect_correct=0,wrong_sign=0,
                    eaf_correct=0,eaf_wrong=0,source_ref_mismatch_unresolved=0,source_ref_mismatch_called=0,
                    sequence_only_calls=0,catalog_calls=0,statuses={})
        details=[]
        for r in rows:
            index=int(r['source_row'])-1;s=source[index]
            counts['statuses'][r['status']]=counts['statuses'].get(r['status'],0)+1
            called=r['pos']!='NA'
            if not valid[index]:
                counts['source_ref_mismatch_called' if called else 'source_ref_mismatch_unresolved']+=1
                if mode=='trusted':assert not called and r['status']=='source_reference_mismatch'
                continue
            if not called:counts['eligible_unresolved']+=1;continue
            if r['status']=='sequence_only':
                counts['sequence_only_calls']+=1;assert r['id']=='NA'
            else:counts['catalog_calls']+=1
            canon=[r['chr'],int(float(r['pos'])),r['ref'],r['alt']]
            case=s['BENCH_ROW'].removeprefix('merge-')
            correct=canon==oracle[case]
            counts['canonical_correct' if correct else 'canonical_wrong']+=1
            effect=number(r['effect']);expected=number(s['BETA'])
            if effect is not None and expected is not None:
                counts['signed_effect_correct']+=correct and effect==expected
                counts['wrong_sign']+=expected!=0 and effect==-expected
            eaf=number(r['eaf']);expected_eaf=number(s['FRQ'])
            if eaf is not None and expected_eaf is not None:counts['eaf_correct' if correct and eaf==expected_eaf else 'eaf_wrong']+=1
            if not correct:details.append(dict(case=case,status=r['status'],observed=canon,oracle=oracle[case],effect=effect,source_effect=expected))
        counts['exact_source_preservation']=(assets/f'{panel}-input.tsv').read_bytes()==(directory/'source-artifact').read_bytes()
        assert counts['retained_rows']==counts['source_rows'] and counts['exact_source_preservation']
        if version in ('functionality-v2','functionality-v3'):assert counts['canonical_wrong']==counts['wrong_sign']==counts['eaf_wrong']==0
        if mode=='trusted':
            assert counts['eligible_unresolved']==0
            assert counts['canonical_correct']==counts['forward_ref_valid']
            assert counts['signed_effect_correct']==counts['effect_eligible']
            assert counts['eaf_correct']==counts['eaf_eligible']
        summary[f'{version}:{mode}:{panel}']=counts
        (directory.parent/'source-contract-discordances.json').write_text(json.dumps(details,indent=2)+'\n')
summary['contract']='Trusted source_ref=A1,source_alt=A2,file-level source_strand=+; same original input bytes. Generic input does not infer REF from other/effect names.'
summary['oracle']='bcftools1.23.1 norm; six single-base reported REF mismatches independently identified against FASTA, not counted as trusted-forward eligible.'
destination=results/('functionality-v3' if (results/'functionality-v3').exists() else 'functionality-v2')/'source-contract-summary.json'
destination.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
