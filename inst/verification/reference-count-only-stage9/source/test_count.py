"""Differential oracle uses actual native files, not planner formulae."""
import json
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
CASES = {
    'empty': [],
    'snps': ['1\t1\tA\tC\trs1\t.'],
    'duplicate_aliases': ['1\t1\tA\tC,G,T\trs1\trs1;rsOld;rsOld;NA;.;'],
    'repeated_contexts': ['1\t2\tA\tC,G\trs2\trsOld'] * 3,
    'indels': ['1\t5\tA\tAA,AAA,AT\trs3\toldA;oldB', '2\t4\tTA\tT,TAA\trs4\t.'],
    'unsupported': ['1\t1\tA\tN,<DEL>,*,A,\trs5\told'],
    'empty_alt': ['1\t1\tA\t\t.\t.'],
    'missing_ids': ['1\t1\tA\tC\tNA\tNA;.;'],
    'positions': [f'1\t{p}\tA\tC\trsPos\t.' for p in
                  ['0', '-1', '001', '+1', ' 1', 'bad', '1x', '9223372036854775808']],
    'mismatch_contig': ['missing\t1\tA\tC\trs6\t.', '1\t1\tC\tA\trs7\t.'],
    'unicode': ['染色体\t1\tA\tC,G\t識別子\t旧名;旧名'],
    'crlf': ['1\t1\tA\tC\trs8\trsOld'],
    'no_final_newline': ['1\t1\tA\tC\trs8\trsOld'],
    'many_aliases': ['1\t1\tA\tC,G,T\trs9\t' + ';'.join(['rs9', *('rsOld' + str(i) for i in range(100))])],
}


def run(binary, *args):
    return subprocess.run([str(HERE / binary), *map(str, args)], capture_output=True, check=True)


@pytest.mark.parametrize('version', [2, 3])
@pytest.mark.parametrize('case,rows', CASES.items())
def test_exact_native_sizes(tmp_path, version, case, rows):
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(('>1\n' + 'A' * 80 + '\n>2\nACGTAACGTACGTA\n>染色体\nAAAA\n').encode())
    sites = tmp_path / 'sites.tsv'
    separator = '\r\n' if case == 'crlf' else '\n'
    tail = '' if case == 'no_final_newline' else separator
    sites.write_bytes((separator.join(['chr\tpos\tref\talt\tid\taliases', *rows]) + tail).encode())
    index = tmp_path / 'native'
    run('writer_fixture', fasta, sites, index, 2, version)
    count = json.loads(run('count_reference', fasta, sites, version).stdout)
    expected = {
        'records_bytes': (index / 'records').stat().st_size,
        'rejected_bytes': (index / 'rejected.tsv').stat().st_size,
        'keys_bytes': (index / 'keys').stat().st_size,
        'offsets_bytes': (index / 'offsets').stat().st_size,
    }
    assert {key: count.get(key) for key in expected} == expected
    assert count['final_key_rows'] * 8 == expected['offsets_bytes']
    assert count['source_rows'] == len(rows)
    assert count['source_bytes'] == sites.stat().st_size
    # Raw sorter rows retain duplicate C/current/alias keys until final emission.
    records = (index / 'records').read_bytes().splitlines(keepends=True)
    assert count['record_rows'] == len(records)
    assert count['raw_key_rows'] >= count['final_key_rows']
    assert count['raw_sort_bytes'] >= count['keys_bytes']
    if case == 'duplicate_aliases':
        lines = (index / 'keys').read_bytes().splitlines(keepends=True)
        assert count['raw_key_rows'] == 19
        assert count['final_key_rows'] == 10
        assert count['raw_sort_bytes'] == sum(len(x) * (1 if x.startswith(b'S:') else 2) for x in lines)
    assert not (tmp_path / 'count-index').exists()


def test_format_suffix_refused(tmp_path):
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(b'>1\nAA\n')
    sites = tmp_path / 'sites.tsv'
    sites.write_bytes(b'chr\tpos\tref\talt\tid\taliases\n')
    result = subprocess.run([str(HERE / 'count_reference'), str(fasta), str(sites), '3suffix'], capture_output=True)
    assert result.returncode != 0
    assert not result.stdout
