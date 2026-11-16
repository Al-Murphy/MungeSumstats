"""Small native fixture-only independent review; never runs prefix/full counts."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

os.umask(0o077)
ROOT=Path('/private/tmp/gwas-harvest-stage2')
SOURCE=ROOT/'reference-count-independent-stage9'
EVIDENCE=ROOT/'results/reference-count-independent-stage9'
AUTHOR=ROOT/'results/reference-count-only-stage9'
HEADER=ROOT/'mss-context-format3-release/src/preservation_core.h'
FROZEN='1bf1a08685488798496e9d7978b1792dbd290a9d'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*a):return subprocess.check_output(['git',*a],cwd=SOURCE,text=True).strip()
assert git('rev-parse','HEAD')==FROZEN and not git('status','--porcelain')
index=AUTHOR/'evidence-index.json'
assert sha(index)=='94b45796456938b95509fcb47cbf7a0fb6ec01ec626e5b71c0b81b7d6fdf84b4'
author=json.loads(index.read_text())
assert len(author['source_files'])==10
for name,pin in author['source_files'].items():assert sha(SOURCE/name)==pin['sha256']
for group in ['artifacts','generated_binaries']:
    for name,pin in author[group].items():
        p=ROOT/name
        assert sha(p)==pin['sha256'] and p.stat().st_size==pin['bytes']
assert sha(HEADER)=='4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc'
builds=[]
for name in ['count_reference','writer_fixture','arithmetic_fixture']:
    command=['c++','-std=c++17','-O2','-I',str(HEADER.parent),str(SOURCE/(name+'.cpp')),'-o',str(SOURCE/name)]
    p=subprocess.run(command,capture_output=True,text=True,timeout=60)
    (EVIDENCE/(name+'-build.log')).write_text(p.stdout+p.stderr)
    assert p.returncode==0
    builds.append({'command':command,'returncode':p.returncode,'binary_sha256':sha(SOURCE/name)})
command=[sys.executable,'-m','pytest',str(SOURCE/'test_count.py'),'-q']
p=subprocess.run(command,capture_output=True,text=True,timeout=60)
(EVIDENCE/'author-29-fresh.log').write_text(p.stdout+p.stderr)
assert p.returncode==0
p=subprocess.run([str(SOURCE/'arithmetic_fixture')],capture_output=True,timeout=30)
assert p.returncode==0 and not p.stderr
cases={
 'repeated-alt-events':b'1\t4\tA\tAA,AA,AAA\trs1\trs1;old;old\n',
 'offset-width-and-cross-record-keys':(b'1\t4\tA\tAA,C,G\tshared\tshared;old;old\n')*111,
 'signed-extremes':b'1\t-9223372036854775808\tA\tC\ti\t.\n1\t9223372036854775807\tA\tC\ti\t.\n',
 'interior-cr-alias':b'1\t1\tA\tC\ti\tx\ry\n',
 'nul-alias':b'1\t1\tA\tC\ti\tx\x00y\n',
 'invalid-utf8':b'1\t1\tA\tC\ti\tx\xc3\n',
 'overlong-utf8':b'1\t1\tA\tC\ti\t\xc0\xaf\n',
 'surrogate-utf8':b'1\t1\tA\tC\ti\t\xed\xa0\x80\n',
 'truncated-six-field-row':b'1\t1\tA\tC\ti\n',
 'extra-field':b'1\t1\tA\tC\ti\talias\textra\n',
 'empty-header-only':b'',
}
outcomes=[]
for version in [2,3]:
    for label,rows in cases.items():
        with tempfile.TemporaryDirectory(prefix='count-small-',dir=EVIDENCE) as td:
            d=Path(td);fa=d/'reference.fa';tsv=d/'sites.tsv';native=d/'native'
            fa.write_bytes(b'>1\r\nAAAAAAAAAAAAAAAAAAAA\r\n')
            tsv.write_bytes((b'' if label=='empty-header-only' else b'chr\tpos\tref\talt\tid\taliases\n')+rows)
            before={str(p.relative_to(d)):sha(p) for p in d.iterdir() if p.is_file()}
            count=subprocess.run([str(SOURCE/'count_reference'),str(fa),str(tsv),str(version)],capture_output=True,timeout=30)
            after={str(p.relative_to(d)):sha(p) for p in d.iterdir() if p.is_file()}
            assert before==after and not native.exists(), 'Counter mutated fixture input/created index'
            writer=subprocess.run([str(SOURCE/'writer_fixture'),str(fa),str(tsv),str(native),'1',str(version)],capture_output=True,timeout=30)
            assert bool(count.returncode)==bool(writer.returncode),(label,version,count.stderr,writer.stderr)
            if count.returncode:
                assert not count.stdout
                outcomes.append({'case':label,'format':version,'result':'both-refuse'})
            else:
                assert not count.stderr and not writer.stderr
                v=json.loads(count.stdout)
                for filename,key in [('records','records_bytes'),('keys','keys_bytes'),('offsets','offsets_bytes'),('rejected.tsv','rejected_bytes')]:
                    assert (native/filename).stat().st_size==v[key],(label,version,key)
                assert len((native/'offsets').read_bytes())//8==v['final_key_rows']
                outcomes.append({'case':label,'format':version,'result':'exact-file-byte-parity','counts':v})
(EVIDENCE/'independent-small-cases.json').write_text(json.dumps(outcomes,indent=2)+'\n')
assert not git('status','--porcelain')
verification={'source_commit':FROZEN,'source_tree':git('rev-parse','HEAD^{tree}'),'tracked_files_verified':10,
 'source_clean_before_and_after':True,'author_index_and_all_entries_verified':True,'header_sha256':sha(HEADER),
 'builds':builds,'author_cases':'29 passed','overflow_control':'passed with empty stderr',
 'additional_small_cases':len(outcomes),'decision':'no blocking finding for count-only diagnostic',
 'limits':'No repeat prefix1M/full-input count, index admission, actual MSS benchmark, cloud/provider/credential operation or native Windows qualification.'}
(EVIDENCE/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
print(json.dumps(verification,indent=2))
