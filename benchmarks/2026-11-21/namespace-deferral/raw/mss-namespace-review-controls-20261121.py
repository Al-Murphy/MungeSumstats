import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

T = Path('/private/tmp')
R = T/'mss-fullfile-20261121/runtime'
B = T/'mss-namespace-builds-20261121'
V = T/'mss-namespace-benchmark-v2-20261121'
D = V/'diagnostic'
O = T/'mss-namespace-oracle-20261121'
OLD = T/'mss-namespace-benchmark-20261121'
H = T/'mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade'
P = T/'mss-namespace-preparation-v2-20261121/metadata'
commits = {'baseline':'14ea7af51a19272c7792fbe582a3d42b291b4a1a','candidate':'f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b'}
expected = set('AnnotationDbi DBI GenomicFeatures KEGGREST RSQLite VariantAnnotation bit bit64 blob cachem fastmap memoise otel png'.split())
rows_expected = {'hard-cases-f2':32,'hard-cases-f3':32,'declared-f2':9,'declared-f3':9,'matched-1M':1000000}
fields = set('source_row status chr pos ref alt id effect eaf swapped complemented candidate_count entity_type entity_id normalized_build grch38_chr grch38_pos grch38_ref grch38_alt'.split())
def read(p):
    p=Path(p)
    assert p.stat().st_size <= 600000, p
    return p.read_bytes()
def j(p):
    return json.loads(read(p))
def sha(p):
    return hashlib.sha256(read(p)).hexdigest()
def pin(rec):
    assert sha(rec['path'])==rec['sha256'], rec
def raw(p):
    return [json.loads(row) for row in read(p).splitlines()]
def terminal(t, exit_code=0):
    assert t['raw_exit']==exit_code and t['closure']=='confirmed' and t['leader_reaped'] and not t['cleanup_errors']
def obs(p, digest, exit_code=0):
    assert sha(p)==digest
    data=j(p)
    assert data['actual_exit_code']==exit_code and not data['live_members'] and not data['inspection_errors']
    if 'guard' in data:
        pin(data['guard'])
        t=raw(data['guard']['path'])[-1]
        terminal(t,exit_code)
        assert data['terminal']==t
    return data
def table(p, header):
    records=list(csv.reader(read(p).decode().splitlines(),delimiter='\t',strict=True))
    assert records[0]==header and len({a[0] for a in records[1:]})==len(records)-1
    return {a[0]: a[1:] for a in records[1:]}
def kv(p):
    return {k:v[0] for k,v in table(p,['key','value']).items()}
bo=obs(T/'mss-namespace-builds-observation-20261121.json','c5177c03c405857cddcf7f6a56521ec98309b16ce4651bb606222e4bd2ec66a9')
pin(bo['producer'])
pin(bo['operator_result'])
op=j(B/'operator-result.json')
assert op['commits']==commits and op['status']=='completed' and len(bo['phases'])==len(op['phases'])==9
assert {p['group_session'] for p in bo['phases']}==set(bo['absent_group_sessions'])
for o,p in zip(bo['phases'],op['phases']):
    pin(o['guard'])
    records=raw(o['guard']['path'])
    assert o['label']==p['label']==records[0]['label'] and p['argv']==records[0]['argv']
    assert o['closure']=='confirmed' and o['leader_reaped'] and o['raw_exit']==0
    terminal(records[-1])
    print('BUILD', o['label'],p['argv'])
roles={r:j(B/r/'build.json') for r in commits}
for r,m in roles.items():
    pin(op['roles'][r])
    assert m['source_commit']==commits[r] and m['library']==str(B/r/'library')
    assert m['archive']['path']==str(B/r/'MungeSumstats.tar') and re.fullmatch('[a-f0-9]{64}',m['archive']['sha256'])
    assert m['dll_sha256']==m['dll']['sha256'] and m['compile_and_link_argv']
    for p in m['install_transcripts']: pin(p)
    assert '* DONE (MungeSumstats)' in read(B/r/'install.stderr').decode()
    terminal(m['terminal'])
assert roles['baseline']['toolchain']==roles['candidate']['toolchain']
bart=obs(T/'mss-namespace-benchmark-v2-observation-20261121.json','3a3359c916164758a3b9586c7f03cc40cc6e2634120fa5421cbdd7e82da2359f')
for k in ('summary','namespace_policy'): pin(bart[k])
summary=j(D/'summary.json')
policy=j(D/'namespace-policy.json')
admission=j(V/'admission.json')
for k in ('build_observation','configuration','preparation_binding','preparation_observation'): pin(admission[k])
assert bart['guard']['path']==str(V/'diagnostic.jsonl') and raw(bart['guard']['path'])[0]['argv']==admission['argv']
assert not read(V/'diagnostic.stderr') and not read(V/'diagnostic.stdout')
prep=obs(T/'mss-namespace-preparation-v2-observation-20261121.json','f9d5e0e01c7fa918bb7b09f9f53cbdc9fe10221ffe995a7c2370e3ac9bc602ed')
for k in ('binding','preparation','configuration'): pin(policy[k])
override=policy['namespace_policy_override']
for k in ('wrapper','audit'): pin(override[k])
assert override['added']==['otel'] or override['added']=='otel'
assert override['wrapper']['sha256']=='c5071a4503b583c51874127425aafb52aa2d0014f7e7fd3be33adadc2e9885b3'
assert override['audit']['sha256']=='1baeba297218e4d7fca22cda8706019e9093b20f9530ba9f62997e393d603e79'
assert policy['source_commits']==commits and set(policy['required_baseline_only'])==expected
assert set(override['revised_eager_only'])==expected and set(override['prior_eager_only'])==expected-{'otel'}
inv = read(D/'namespace-inputs-before.json')
assert inv==read(D/'namespace-inputs-after.json')==read(P/'preparation-inventory.json')
ns_inv=json.loads(inv)
libs=ns_inv['installed_dependency_tree_sha256']
assert sorted(map(len,libs.values())) == [30,96]
assert all(ns_inv['role_trees'][r]['tree_sha256']==m['package_tree_sha256'] for r,m in roles.items())
pre=j(D/'pins.json')
post=j(D/'post-pins.json')
assert pre['assets']==post['assets'] and pre['installed_packages']==post['installed_packages']==ns_inv['role_trees']
assert summary['fresh_processes'] and summary['suite']=='diagnostic' and summary['profile']=='none'
assert summary['target_status']=='not_applicable' and summary['processing_speed_and_rss_target_met'] is None
assert set(rows_expected)==set(summary['cases'])==set(policy['cases'])
raw_pins={(w['case'], w['worker']):w for w in policy['workers']}
observed={(w['case'],w['worker']):w for w in bart['workers']}
assert len(raw_pins)==len(policy['workers'])==len(observed)==len(bart['workers'])==30
assert len({w['pid'] for w in observed.values()})==30
all_role_state={}
for case, n in rows_expected.items():
    expected_order=['pair-01-baseline','pair-01-candidate','pair-02-candidate','pair-02-baseline','pair-03-baseline','pair-03-candidate']
    assert [w['worker'] for w in bart['workers'] if w['case']==case]==expected_order
    s=summary['cases'][case]
    assert s['rows']==n and [Path(w['directory']).name for w in s['observations']]==expected_order
    cert=policy['cases'][case]
    assert set(cert['baseline_only'])==expected and not cert['candidate_only']
    verification=j(D/case/'verification.json')
    assert verification==s['verification'] and verification['normalized_rows_checked']==n and verification['edition_count']==6
    assert verification['stable_byte_equality'] and verification['stable_rds_manifest_equality']
    assert set(verification['normalized_columns'])==fields and sum(verification['status_counts'].values())==n
    assert verification['oracle_rows_checked']==(0 if case=='matched-1M' else n)
    assert set(verification['oracle_fields'] or [])==(fields if case!='matched-1M' else set())
    if case!='matched-1M':
        cfg=[x for x in pre['configuration']['cases'] if x.get('oracle') and x['rows']==n]
        oracle_paths={x['oracle']['path']:x['oracle']['sha256'] for x in cfg}
        assert len(oracle_paths)==1
        gold_path, gold_sha=next(iter(oracle_paths.items()))
        assert sha(gold_path)==gold_sha
        with open(gold_path,encoding='utf-8',newline='') as h:
            gold=list(csv.DictReader(h,delimiter='\t',strict=True))
        for worker in expected_order:
            with open(D/case/worker/'edition/normalized.tsv',encoding='utf-8',newline='') as h:
                actual=list(csv.DictReader(h,delimiter='\t',strict=True))
            assert len(gold)==len(actual)==n
            assert all(set(a)==fields and set(g)==fields and a==g and a['source_row']==str(i)
                       for i,(a,g) in enumerate(zip(actual,gold),1))
    mverify=j(D/case/'manifest-verification.json')
    assert mverify['exit_code']==0 and read(D/case/'stable-manifests.txt').strip()==b'pass'
    assert mverify['command'][-6:]==[str(D/case/w/'edition/manifest.rds') for w in expected_order]
    artifacts=j(D/case/'artifact-digests.json')
    assert set(artifacts)=={str(D/case/w) for w in expected_order}
    first=next(iter(artifacts.values()))
    assert first['source-artifact']==verification['source_artifact'] and first['source.tsv']==verification['decoded_source']
    byrole={r:[] for r in commits}
    for worker, w in zip(expected_order, s['observations']):
        work=D/case/worker
        role=worker.split('-')[-1]
        m=j(work/'process.json')
        metrics=kv(work/'metrics.tsv')
        settings=kv(work/'settings.tsv')
        assert m==w['process'] and metrics==w['metrics'] and w['pair']==int(worker[5:7])
        assert m['case']==case and m['role']==role and m['pid']==observed[case,worker]['pid']
        assert m['exit_code']==observed[case,worker]['raw_exit']==0 and m['rss_complete'] and m['rss_source']=='wait4 ru_maxrss'
        assert m['command']==[pre['configuration']['r_executable']['path'],'--vanilla',str(H/'worker.R'),'run',str(work/'settings.tsv')]
        assert metrics['threads']=='1' and float(metrics['source_rows'])==float(metrics['normalized_rows'])==n
        assert metrics['implementation']=='mss-preservation-7'
        dll=roles[role]['dll']['path']
        buildfiles=ns_inv['role_trees'][role]['files']
        assert metrics['dll']==settings['dll']==dll and metrics['dll_md5']==settings['dll_md5']==buildfiles['libs/MungeSumstats.so']['md5']
        assert buildfiles['libs/MungeSumstats.so']['sha256']==roles[role]['dll']['sha256']
        assert settings['library.0']==roles[role]['library'] and settings['library.1'] in libs
        assert settings['mode']=='diagnostic' and settings['profile']=='none'
        pin(raw_pins[case,worker]['raw_inventory'])
        assert raw_pins[case,worker]['role']==role and raw_pins[case,worker]['raw_inventory']['path']==str(work/'namespaces.tsv')
        state=table(work/'namespaces.tsv',['package','version','path'])
        assert state==cert['actual_role_states'][role]
        assert state==all_role_state.setdefault(role,state)
        assert state['MungeSumstats']==['1.21.4',str(B/role/'library/MungeSumstats')]
        for name,(version, path) in state.items():
            pp=Path(path)
            assert pp.resolve(strict=True)==pp and pp.name==name
            if name!='MungeSumstats': assert name in libs[str(pp.parent)]
            desc=read(pp/'DESCRIPTION').decode('latin-1')
            assert re.search(r'^Version:[ \t]*([^\s]+)',desc,re.M).group(1).replace('-','.')==version, (name,version)
        entries=artifacts[str(work)]
        assert entries.keys()==first.keys()
        for name, entry in entries.items():
            if name!='manifest.rds': assert entry==first[name]
            p=work/'edition'/name
            assert p.stat().st_size==entry['bytes']
            if entry['bytes']<=16384:
                data=read(p)
                assert hashlib.sha256(data).hexdigest()==entry['sha256'] and hashlib.md5(data).hexdigest()==entry['md5']
        byrole[role].append((metrics,m))
    a,b=[all_role_state[r] for r in commits]
    assert (len(a),len(b))==(85,71) and set(a)-set(b)==expected and not set(b)-set(a)
    assert all(a[x]==b[x] for x in b if x!='MungeSumstats')
    assert cert['common_comparison_state']=={**b,'MungeSumstats':['1.21.4','<selected installation>']}
    assert len({entries['manifest.rds']['sha256'] for entries in artifacts.values()})==6
    med={}
    for role, obsrows in byrole.items():
        med[role]={}
        for phase in ('package_setup','reference_open','processing'):
            med[role][phase]={}
            for part in ('elapsed','user.self','sys.self','user.child','sys.child'):
                med[role][phase][part]=statistics.median(float(x[phase+'.'+part]) for x,_ in obsrows)
        for key,s_key in [('wall_seconds','process_wall_median_seconds'),('user_seconds','process_user_median_seconds'),('system_seconds','process_system_median_seconds'),('peak_rss_bytes','median_peak_rss_bytes')]:
            med[role][s_key]=statistics.median(p[key] for _,p in obsrows)
        med[role]['max_peak_rss_bytes']=max(p['peak_rss_bytes'] for _,p in obsrows)
        med[role]['rss_complete']=True
        assert med[role]==s['medians'][role]
    pair=[float(a[0]['processing.elapsed'])/float(b[0]['processing.elapsed']) for a,b in zip(*byrole.values())]
    assert pair==s['baseline_over_candidate_by_pair']
    a,b=med.values()
    assert a['processing']['elapsed']/b['processing']['elapsed']==s['ratio_of_median_processing_seconds']
    print('CASE',case,n,verification['status_counts'])
    for label,av,bv in [('setup',a['package_setup']['elapsed'],b['package_setup']['elapsed']),('wall',a['process_wall_median_seconds'],b['process_wall_median_seconds']),('rss',a['median_peak_rss_bytes'],b['median_peak_rss_bytes']),('processing',a['processing']['elapsed'],b['processing']['elapsed'])]:
        print(label, av,bv,'ratio',av/bv,'less pct',100*(1-bv/av), 'absolute',av-bv)
    print('pairs processing',pair)
print('INVENTORY',hashlib.sha256(inv).hexdigest())
print('PIN COUNTS',len(pre['assets']),{p:len(x) for p,x in libs.items()})
print('BASELINE ONLY VERSIONS', {p:all_role_state['baseline'][p][0] for p in sorted(expected)})
failure=T/'mss-namespace-benchmark-failure-observation-20261121.json'
fo=obs(failure,sha(failure),1)
assert len(fo['completed_workers'])==2 and 'otel' in j(OLD/'diagnostic/failure.json')['error']
for w in fo['completed_workers']:
    pin(w['process'])
    old=j(w['process']['path'])
    assert old['exit_code']==0 and old['pid'] not in {x['pid'] for x in observed.values()}
for role in commits:
    assert table(OLD/'diagnostic/hard-cases-f2'/('pair-01-'+role)/'namespaces.tsv',['package','version','path'])==all_role_state[role]
alignment=R/'NAMESPACE-DEFERRAL-SOURCE-ALIGNMENT.json'
assert sha(alignment)=='fe349f66af494a00a67d9284970b0c420813cb47f44da181dee83bceb7df3494'
assert j(alignment)['roles']['eager_ann_fixed']['commit']==commits['baseline'] and j(alignment)['roles']['deferred']['commit']==commits['candidate']
ore=obs(T/'mss-namespace-oracle-observation-20261121.json','565acb8b46d8be55be6c43bdc4c6ae6820f63b5ea27b7fc1a7e91405f31f760b')
for name in ('oracle','input_window'):
    pin(ore[name])
oa=j(O/'admission.json')
ow=j(O/'input-window.json')
o=j(O/'snp-oracle-1M.json')
assert oa['source_commit']==commits['candidate'] and oa['argv']==raw(ore['guard']['path'])[0]['argv']
assert ow['equal'] and ow['inputs_before']==ow['inputs_after']==oa['inputs_before'] and len(ow['inputs_before'])==5
assert ow['oracle']==ore['oracle'] and ow['oracle']['sha256']=='cefb1f6d1308ccdbc71f18d0dea872eba2d8b4a6298be32653668e5792415b50'
pin(oa['profile'])
assert sha(oa['argv'][4])==oa['oracle_sha256']
assert oa['argv'][1:4]==['-I','-B','-S']
assert oa['argv'][oa['argv'].index('--edition')+1]==f'mss7={D}/matched-1M/pair-01-candidate/edition'
edition=j(D/'matched-1M/artifact-digests.json')[str(D/'matched-1M/pair-01-candidate')]
assert ow['inputs_before'][str(D/'matched-1M/pair-01-candidate/edition/normalized.tsv')]==edition['normalized.tsv']
assert o['source_bytes_sha256']==edition['source-artifact']['sha256']==edition['source.tsv']['sha256']
assert o['rows_expected']==ore['all_19_fields_correct_rows']==o['versions']['mss7']['counts']['all_19_fields_correct_rows']==1000000
assert o['passed'] and o['all_rows_examined'] and not o['deltas'] and o['unaccounted_deltas']==0
assert o['source_copies_exact_each']==ore['source_copies_exact_each']==2
assert o['versions']['mss7']['statuses']==summary['cases']['matched-1M']['verification']['status_counts']
assert not read(O/'oracle.stderr')
print('COMMANDS',admission['argv'], oa['argv'])
for p in [B/'operator-result.json',*(B/r/'build.json' for r in commits),D/'summary.json',D/'namespace-policy.json',D/'pins.json',D/'post-pins.json',D/'namespace-inputs-before.json',D/'namespace-inputs-after.json',*(D/c/n for c in rows_expected for n in ('artifact-digests.json','verification.json','manifest-verification.json')),R/'NAMESPACE-BENCH-V2-REVIEW.md',R/'NAMESPACE-DEFERRAL-RUNTIME-REVIEW.md',failure]:
    print('PIN',sha(p),p)
print('PASS retained controls')
