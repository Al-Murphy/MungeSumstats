#!/usr/bin/env python3
"""Sequential benchmark supervision, immutable scripts, OS RSS and disk receipts.

macOS implementation: external time is authoritative for per-process peak RSS.
The process-tree sampler is a host-protection guard, not the RSS estimator.
"""
import argparse, csv, hashlib, json, os, pathlib, re, shutil, signal, subprocess, time

p=argparse.ArgumentParser()
p.add_argument('--base',type=pathlib.Path,default=pathlib.Path('/private/tmp/gwas-harvest-stage2'))
p.add_argument('--tag',default='matrix-v1')
p.add_argument('--memory-cap-gib',type=float,default=32)
p.add_argument('--cases',default='unoptimized:unique:cold,new:unique:cold,new:1M:cold,legacy:1M:cold,new:1M:warm,legacy:1M:warm,new:10M:cold,legacy:10M:cold,new:10M:warm,legacy:10M:warm,new:40M:cold,legacy:40M:cold,new:40M:warm,legacy:40M:warm')
a=p.parse_args()
repo=pathlib.Path(__file__).resolve().parents[1]
root=a.base/'results'/'mss'/a.tag
root.mkdir(exist_ok=False)
runner=root/'runner.R'
shutil.copyfile(repo/'scripts/benchmark-matched.R',runner)
shutil.copyfile(repo/'scripts/benchmark-scale.py',root/'verify.py')
def digest(path):
    h=hashlib.md5()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()
def disk_bytes(path):
    return sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
def tree_rss(pid):
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True)
    rows=[tuple(map(int,line.split())) for line in text.splitlines() if line.strip()]
    descendants={pid}
    for _ in range(8):
        descendants.update(child for child,parent,rss in rows if parent in descendants)
    return sum(rss*1024 for child,parent,rss in rows if child in descendants),descendants
for impl,size,cache in (case.split(':') for case in a.cases.split(',')):
    label=f'{impl}-{size}-{cache}'
    dest=root/label
    dest.mkdir()
    source=a.base/'assets'/'mss'/f'matched-{size}.tsv'
    expected={'unique':437458,'1M':1000000,'10M':10000000,'40M':40000000}[size]
    output=dest/('legacy.tsv' if impl=='legacy' else 'edition')
    metrics=dest/'metrics.tsv'
    env=os.environ.copy()
    env['MSS_EXPECTED_ROWS']=str(expected)
    if impl=='unoptimized':env['MSS_NEW_RLIB']=str(a.base/'Rlib'/'unoptimized')
    tmp=dest/'tmp';tmp.mkdir();env['TMPDIR']=str(tmp)
    cmd=['/usr/bin/time','-l','Rscript',str(runner),'legacy' if impl=='legacy' else 'new',
         str(source),str(a.base/'assets/reference/matched-index'),str(output),str(metrics),cache,'1']
    receipt=dict(label=label,command=cmd,expected_rows=expected,runner_md5=digest(runner),
                 started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                 memory_guard_bytes=int(a.memory_cap_gib*1024**3),os_cache='uncontrolled, never flushed',
                 application_state=cache,source_scope='437458 real SNP records repeated with synthetic exact BETA/EAF',
                 bounded_page_cache_bytes=0 if impl in ('unoptimized','legacy') else 28*1024**2)
    with (dest/'process.log').open('w') as log:
        begin=time.monotonic();proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
        receipt['time_wrapper_pid']=proc.pid
        (root/'live.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print('START',label,'PID',proc.pid,flush=True)
        peak_disk=0;peak_sample=0;guard=False
        while proc.poll() is None:
            used,children=tree_rss(proc.pid)
            peak_sample=max(peak_sample,used);peak_disk=max(peak_disk,disk_bytes(dest))
            if used>receipt['memory_guard_bytes']:
                guard=True
                # Child first lets time reap and print its resource receipt.
                for child in children-{proc.pid}:
                    try:os.kill(child,signal.SIGKILL)
                    except ProcessLookupError:pass
                print('MEMORY_GUARD',label,used,flush=True)
            time.sleep(2)
        receipt.update(returncode=proc.returncode,supervisor_elapsed=time.monotonic()-begin,
                       peak_output_and_temp_bytes=max(peak_disk,disk_bytes(dest)),
                       sampled_peak_process_tree_rss_bytes=peak_sample,memory_guard_triggered=guard,
                       completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    logtext=(dest/'process.log').read_text()
    for pattern,name,conv in [(r'(\d+)\s+maximum resident set size','external_peak_rss_bytes',int),
             (r'([0-9.]+) real','external_wall_seconds',float),
             (r'([0-9.]+) user','external_user_seconds',float),
             (r'([0-9.]+) sys','external_system_seconds',float)]:
        found=re.search(pattern,logtext)
        if found:receipt[name]=conv(found.group(1))
    if metrics.exists():
        with metrics.open() as f:receipt['metrics']={r['key']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
    if proc.returncode==0:
        verify=['python3',str(root/'verify.py'),'verify','--panel',str(a.base/'assets/reference/matched-snps.tsv'),
                '--input',str(source),'--rows',str(expected),'--mode','legacy' if impl=='legacy' else 'new',
                '--output',str(output),'--receipt',str(dest/'correctness.json')]
        check=subprocess.run(verify,capture_output=True,text=True)
        (dest/'verify.log').write_text(check.stdout+check.stderr)
        receipt['correctness_returncode']=check.returncode
        if check.returncode:print('CORRECTNESS_FAILURE',label,flush=True)
    (dest/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('DONE',label,'rc',proc.returncode,'wall',receipt.get('external_wall_seconds'),
          'RSS',receipt.get('external_peak_rss_bytes'),'disk',receipt['peak_output_and_temp_bytes'],flush=True)
    # Free complete priming editions after retaining their timing/hash receipts.
    for prime in dest.glob('*.prime*'):
        if prime.is_dir():
            manifests=prime/'manifest.rds'
            if manifests.exists():shutil.copyfile(manifests,dest/'prime-manifest.rds')
            shutil.rmtree(prime)
        else:prime.unlink()
(root/'live.json').write_text(json.dumps({'status':'complete','completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})+'\n')
print('ALL_CASES_COMPLETE',flush=True)
