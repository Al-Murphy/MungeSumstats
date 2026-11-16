"""Actual tiny producer, native writer, and unchanged process-tree monitor."""
import bz2
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from test_coordinator import api, fixture, spec, write, HERE, fixture_roots

BASE = Path('/private/tmp/gwas-harvest-stage2')
GUARD = BASE/'results/verification-stage7/guard_common.py'
HARNESS = GUARD.with_name('process_memory.py')
PRODUCER = BASE/'biodb/src/biodb/gwas_reference_stream.py'
HEADER = BASE/'mss-context-format3-release/src/preservation_core.h'


def actual_guard():
    sys.path.insert(0, str(GUARD.parent))
    module = importlib.util.spec_from_file_location('actual_guard_for_tiny_test', GUARD)
    value = importlib.util.module_from_spec(module); module.loader.exec_module(value)
    return value


def actual_profile(tmp_path):
    p = fixture(tmp_path)
    for key, path in [('guard', GUARD), ('harness', HARNESS), ('producer', PRODUCER), ('header', HEADER)]:
        p['artifacts'][key] = spec(path)
    return p


def produce_both(tmp_path):
    """Guard real stdlib producer on newly owned complete tiny sources."""
    p = actual_profile(tmp_path)
    runner = tmp_path/'tiny-producer.py'
    runner.write_text('import json,runpy,sys\nn=runpy.run_path(sys.argv[1])\n'
                      'n["prepare_reference"](json.load(open(sys.argv[2])),sys.argv[3])\n')
    p['builds'] = []
    guard = actual_guard()
    for build in ['GRCh37','GRCh38']:
        root = tmp_path/(build+'-sources'); root.mkdir()
        fasta = root/'original.fa'; fasta.write_bytes(b'>chr1\r\nAAAAAAAA\r\n')
        assembly = root/'assembly.txt'
        assembly.write_text('# Assembly name: '+build+'\n# RefSeq assembly accession: GCF_000001405.25\n'
                            '1\tassembled-molecule\t1\tChromosome\tCM1.1\t=\tNC_1.1\tPrimary Assembly\t8\tchr1\n')
        vcf = root/'source.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n##reference='+build+'\n#CHROM\tPOS\tID\tREF\tALT\n'
                       'chr1\t1\trs1;rs2\tA\tC,G,AA,<DEL>,*\nchr1\t2\trs3\tA\tAA,AAA\n')
        merged, withdrawn = root/'merged.json.bz2', root/'withdrawn.json.bz2'
        merged.write_bytes(bz2.compress((json.dumps({'refsnp_id':'9','merged_snapshot_data':{'merged_into':['1']}})+'\n').encode()))
        withdrawn.write_bytes(bz2.compress(b''))
        rights = root/'terms.txt'; rights.write_text('Owned synthetic test fixture, unrestricted test reuse.\n')
        def receipt(path):
            return dict(spec(path), md5=hashlib.md5(path.read_bytes()).hexdigest(),
                        build=build, release='157', scope='tiny owned complete fixture',
                        source_url='https://example.invalid/owned-fixture',
                        licence_url='https://example.invalid/terms', rights_basis='newly authored fixture',
                        rights_receipt=dict(spec(rights), md5=hashlib.md5(rights.read_bytes()).hexdigest(),
                                            source_url='https://example.invalid/terms'), completeness='complete_source')
        plan = dict(build=build, release='tiny-actual-v1', scope='tiny_fixture', assembly='GCF_000001405.25',
                    vcf_headers=['##reference='+build], sources={k:receipt(f) for k,f in
                    [('fasta',fasta),('assembly_report',assembly),('vcf',vcf),('merged',merged),('withdrawn',withdrawn)]})
        plan_pin = write(root/'source-plan.json',plan)
        final = tmp_path/(build+'-published')
        owned = tmp_path/(build+'-producer-evidence'); owned.mkdir()
        command = [sys.executable,str(runner),str(PRODUCER),plan_pin['path'],str(final)]
        guard.guarded(command,owned,owned/'guard.json',owned/'producer.log',
                      {'max_tree_memory_bytes':8*1024**3,'min_filesystem_free_bytes':160*1024**3,
                       'poll_seconds':.02,'max_case_disk_bytes':16*1024**2},time.monotonic()+30,os.environ.copy())
        r = json.loads((final/'receipts.json').read_text())
        p['builds'].append(dict(build=build,release=plan['release'],final_dir=str(final),source_plan=plan_pin,
                                receipt=spec(final/'receipts.json'),producer_guard=spec(owned/'guard.json'),
                                producer_command=command,preparation_version=r['preparation_version'],
                                expected_emitted_rows=r['counts']['emitted_rows']))
    return p


def run_profile(p, tmp_path, name, *, expected_publications):
    pin = write(tmp_path/(name+'-profile.json'),p)
    return api().launch(pin['path'],pin['sha256'],tmp_path/name,expected_publications=expected_publications)


def test_actual_producer_writer_and_guarded_count(tmp_path):
    p = produce_both(tmp_path)
    source_before = {b['build']: {x.name: spec(x) for x in Path(b['final_dir']).iterdir()} for b in p['builds']}
    result = run_profile(p,tmp_path,'guarded-count',expected_publications=fixture_roots(tmp_path,published=True))
    assert result['guard']['status']=='success'
    assert result['guard']['peak_process_tree']['processes'] >= 1
    assert result['full_reference_admitted'] is False and result['index_builds_executed']==0
    for seal, counted in zip(p['builds'],result['builds'],strict=True):
        final = Path(seal['final_dir'])
        assert counted['counts']['source_rows']==3  # 2 original contexts + one historical-alias context.
        baseline = Path(os.environ.get('MSS_BASELINE_COUNTER', HERE/'fix-evidence'/'original-counter'))
        original = subprocess.run([str(baseline),str(final/'reference.fa'),str(final/'variants.tsv'),'3'],capture_output=True)
        assert original.returncode==0 and not original.stderr
        (tmp_path/(seal['build']+'.original.stdout')).write_bytes(original.stdout)
        (tmp_path/(seal['build']+'.original.stderr')).write_bytes(original.stderr)
        original_counts=json.loads(original.stdout)
        assert len(original_counts)==20
        assert {k:v for k,v in counted['counts'].items() if k!='fai_upper_bytes'}==original_counts
        writer = tmp_path/(seal['build']+'-tiny-writer-oracle')
        call = subprocess.run([str(HERE/'writer_fixture'),str(final/'reference.fa'),str(final/'variants.tsv'),str(writer),'2','3'],capture_output=True)
        assert call.returncode==0 and not call.stderr
        for field,name in [('records_bytes','records'),('keys_bytes','keys'),('offsets_bytes','offsets'),('rejected_bytes','rejected.tsv')]:
            assert counted['counts'][field]==(writer/name).stat().st_size
        assert (writer/'fasta.fai').stat().st_size <= counted['counts']['fai_upper_bytes']
        assert not (final/'fasta.fai').exists()
    assert source_before == {b['build']: {x.name: spec(x) for x in Path(b['final_dir']).iterdir()} for b in p['builds']}
    assert sorted(x.name for x in (tmp_path/'guarded-count').iterdir())==['guard.json','summary.json','worker','worker.log']


@pytest.mark.parametrize('failure',['deadline','output'])
def test_actual_guard_censor_preserves_failure_evidence(tmp_path,failure):
    p=actual_profile(tmp_path)
    native=tmp_path/'slow-native'
    native.write_text('#!'+sys.executable+'\nimport os,time\n'+
                      ('time.sleep(10)\n' if failure=='deadline' else
                       'time.sleep(.2)\nfor i in range(256):\n os.write(1,b"x"*8192)\n time.sleep(.002)\n'))
    native.chmod(0o700);p['artifacts']['native']=spec(native)
    if failure=='deadline':p['limits']['deadline_seconds']=.5
    else:p['limits']['owned_case_bytes']=128*1024
    with pytest.raises((api().Refusal,RuntimeError)):
        run_profile(p,tmp_path,'censored',expected_publications=fixture_roots(tmp_path))
    case=tmp_path/'censored'
    assert not (case/'summary.json').exists()
    receipt=json.loads((case/'guard.json').read_text())
    assert receipt['status']=='censored' and receipt['natural_oom'] is False
    assert receipt['peak_process_tree']['processes']>=2
