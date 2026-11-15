#!/usr/bin/env python3
"""Make disclosed replicated real-variant scale inputs and verify outputs fully."""
import argparse, collections, csv, json, pathlib

p = argparse.ArgumentParser()
p.add_argument('action',choices=['generate','verify'])
p.add_argument('--panel',required=True,type=pathlib.Path)
p.add_argument('--input',required=True,type=pathlib.Path)
p.add_argument('--rows',required=True,type=int)
p.add_argument('--output',type=pathlib.Path)
p.add_argument('--mode',choices=['new','legacy'])
p.add_argument('--receipt',required=True,type=pathlib.Path)
a = p.parse_args()
with a.panel.open() as f:
    records = list(csv.DictReader(f,delimiter='\t'))
panel = [r for r in records if len(r['ref'])==len(r['alt'])==1 and
         r['ref']+r['alt'] not in ('AT','TA','CG','GC') and r['id'].startswith('rs')]
if not panel: raise SystemExit('No nonpalindromic biallelic rsID SNPs')
complement = str.maketrans('ACGT','TGCA')
if a.action=='generate':
    with a.input.open('x') as f:
        f.write('SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\tSE\tP\n')
        for i in range(a.rows):
            r = panel[i%len(panel)]
            ref,alt = r['ref'],r['alt']
            swap,comp = bool((i//len(panel))%2),False
            other,effect = (alt,ref) if swap else (ref,alt)
            if comp: other,effect = other.translate(complement),effect.translate(complement)
            f.write(f"{r['id']}\t{r['chr']}\t{r['pos']}\t{other}\t{effect}\t{'-.125' if swap else '.125'}\t{'.75' if swap else '.25'}\t.05\t.01\n")
    result = dict(action=a.action,rows=a.rows,unique_panel_records=len(panel),
                  description='Real variant coordinates/rsIDs/alleles; synthetic exact effects/EAF; repeated panel, match/swap orientations in full cycles',
                  expected_beta=.125,expected_eaf=.25)
else:
    lookup = {r['id']:(r['chr'],int(r['pos']),r['ref'],r['alt']) for r in panel}
    counts = dict(rows=0,alleles_correct=0,effects_correct=0,eaf_correct=0,identity_correct=0,fully_correct=0)
    statuses = {}
    expected_counts=collections.Counter()
    with a.input.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['SNP'] not in lookup:raise AssertionError('Input ID absent from matched panel')
            expected_key=lookup[r['SNP']]
            if (r['CHR'],int(r['BP']))!=expected_key[:2]:raise AssertionError('Input coordinate/ID mismatch')
            pair=(r['A1'],r['A2'])
            if pair==expected_key[2:]:swap=False
            elif pair==expected_key[2:][::-1]:swap=True
            else:raise AssertionError('Input alleles differ from disclosed matched workload')
            beta=float(r['BETA'])*(-1 if swap else 1)
            eaf=1-float(r['FRQ']) if swap else float(r['FRQ'])
            if beta!=.125 or eaf!=.25:raise AssertionError('Input differs from exact disclosed synthetic truth')
            expected_counts[r['SNP']]+=1
    if sum(expected_counts.values())!=a.rows:raise AssertionError('Input count differs from acquisition receipt')
    source_handle=a.input.open()
    source_rows=csv.DictReader(source_handle,delimiter='\t')
    output = a.output/'normalized.tsv' if a.mode=='new' else a.output
    with output.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        typed=a.mode=='new' and 'entity_type' in reader.fieldnames
        if typed:
            if 'entity_id' not in reader.fieldnames:raise AssertionError('Incomplete typed entity output')
            counts['typed_entity_correct']=0
        for r in reader:
            counts['rows'] += 1
            if a.mode=='new':
                key = (r['chr'],int(float(r['pos'])) if r['pos']!='NA' else -1,r['ref'],r['alt'])
                effect,eaf = r['effect'],r['eaf']
                original=next(source_rows,None)
                if original is None:raise AssertionError('Output exceeds input')
                expected_id=original['SNP']
                expected_key=lookup[expected_id]
                expected_status='match' if (original['A1'],original['A2'])==expected_key[2:] else 'swap'
                identity=(r['id']==expected_id and key==expected_key and r['status']==expected_status)
                statuses[r['status']] = statuses.get(r['status'],0)+1
                if int(r['source_row'])!=counts['rows']: raise AssertionError('Broken source row association')
            else:
                key = (r['CHR'],int(r['BP']),r['A1'],r['A2'])
                effect,eaf = r['BETA'],r['FRQ']
                expected_id=r['SNP']
                expected_key=lookup.get(expected_id)
                identity=(expected_key==key and expected_counts[expected_id]>0)
            expected_counts[expected_id]-=1
            ok = (key==expected_key, effect!='NA' and float(effect)==.125, eaf!='NA' and float(eaf)==.25)
            for name,value in zip(('alleles_correct','effects_correct','eaf_correct'),ok): counts[name]+=value
            counts['identity_correct']+=identity
            typed_ok=not typed or (r['entity_type']=='small_variant' and r['entity_id']==expected_id)
            if typed:counts['typed_entity_correct']+=typed_ok
            counts['fully_correct'] += all(ok) and identity and typed_ok
    source_handle.close()
    result = dict(action=a.action,mode=a.mode,expected_rows=a.rows,counts=counts,statuses=statuses,
                  identity_scope='Exact source_row key+canonical rsID+status for new; full canonical identity multiset for reordered legacy',
                  typed_entity_columns_verified=typed,
                  passed=all(value==a.rows for value in counts.values()) and all(n==0 for n in expected_counts.values()))
a.receipt.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
if a.action=='verify' and not result['passed']: raise SystemExit(1)
