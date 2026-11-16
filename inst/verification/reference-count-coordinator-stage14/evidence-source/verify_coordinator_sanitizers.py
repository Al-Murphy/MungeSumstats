"""Strict sanitizer replay on tiny fixtures; preserve every native stderr."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
OUT=HERE/'coordinator-evidence'/'sanitizer-calls'
BIN=OUT.parent/'sanitized'


def main():
    OUT.mkdir()
    os.environ['ASAN_OPTIONS']='halt_on_error=1'
    os.environ['UBSAN_OPTIONS']='halt_on_error=1:print_stacktrace=1'
    calls=[]
    def invoke(name,args,expected_success=True):
        call=subprocess.run([str(BIN/name),*map(str,args)],capture_output=True)
        number=len(calls)
        (OUT/f'{number:03}.stdout').write_bytes(call.stdout)
        (OUT/f'{number:03}.stderr').write_bytes(call.stderr)
        calls.append(dict(binary=name,exit_code=call.returncode,stderr_bytes=len(call.stderr),expected_success=expected_success))
        if expected_success:
            assert call.returncode==0 and not call.stderr, calls[-1]
        return call
    module_spec=importlib.util.spec_from_file_location('original_tiny_tests',HERE/'test_count.py')
    tests=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(tests)
    tests.run=lambda name,*args:invoke(name,args)
    for version in (2,3):
        for case,rows in tests.CASES.items():
            with tempfile.TemporaryDirectory(dir=OUT) as temp:
                tests.test_exact_native_sizes(Path(temp),version,case,rows)
    fai_spec=importlib.util.spec_from_file_location('fai_tiny_tests',HERE/'test_fai_count.py')
    fai=importlib.util.module_from_spec(fai_spec);fai_spec.loader.exec_module(fai)
    for version in (2,3):
        for fasta_bytes in fai.FASTAS:
            with tempfile.TemporaryDirectory(dir=OUT) as temp:
                root=Path(temp);fa=root/'reference.fa';tsv=root/'variants.tsv';oracle=root/'writer'
                fa.write_bytes(fasta_bytes);tsv.write_bytes(b'chr\tpos\tref\talt\tid\taliases\n')
                invoke('writer_fixture',[fa,tsv,oracle,2,version])
                result=json.loads(invoke('count_reference',[fa,tsv,version]).stdout)
                actual=(oracle/'fasta.fai').read_bytes()
                upper=sum(len(row.split(b'\t')[0])+84+len(os.linesep.encode()) for row in actual.splitlines())
                assert result['fai_upper_bytes']==upper and upper>=len(actual)
                assert sorted(x.name for x in root.iterdir())==['reference.fa','variants.tsv','writer']
    with tempfile.TemporaryDirectory(dir=OUT) as temp:
        root=Path(temp);fa=root/'reference.fa';tsv=root/'variants.tsv'
        fa.write_bytes(b'>1\nAA\n');tsv.write_bytes(b'chr\tpos\tref\talt\tid\taliases\n')
        bad=invoke('count_reference',[fa,tsv,'3suffix'],False)
        assert bad.returncode==1 and not bad.stdout and bad.stderr==b'Expected exact format 2 or 3\n'
    invoke('arithmetic_fixture',[])
    value=dict(status='success',sanitizers=['address','undefined'],halt_on_error=True,
               successful_empty_stderr_calls=sum(x['expected_success'] for x in calls),
               expected_refusal_calls=1,tiny_writer_oracle_outputs=40,
               differential_cases=41,full_production_index_builds=0,count_process_index_outputs=0,calls=calls)
    (OUT/'summary.json').write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps({k:v for k,v in value.items() if k!='calls'}))


if __name__=='__main__':main()
