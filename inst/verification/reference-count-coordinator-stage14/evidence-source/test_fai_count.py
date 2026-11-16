"""Compare diagnostic FAI bounds to files emitted by the actual writer."""
import json
import os
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
COUNTER = Path(os.environ.get('MSS_DIAGNOSTIC_COUNTER', HERE / 'count_reference'))
FASTAS = [
    b'>1\nAAAA\n',
    b'>1\r\nAAAA\r\n>2\r\nAAA\r\n',
    b'>1\nAAAA\nAA',
    '>染色体\nAAAA\n>空\n'.encode(),
    b'>long_contig_name description is excluded\nAAAA\nAAAA\nAAA\n',
    b'>empty\n>1\nA\n',
]


@pytest.mark.parametrize('version', [2, 3])
@pytest.mark.parametrize('fasta_bytes', FASTAS)
def test_fai_upper_bytes_cover_real_writer(tmp_path, version, fasta_bytes):
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(fasta_bytes)
    sites = tmp_path / 'sites.tsv'
    sites.write_bytes(b'chr\tpos\tref\talt\tid\taliases\n')
    index = tmp_path / 'actual-writer'
    subprocess.run([str(HERE / 'writer_fixture'), str(fasta), str(sites), str(index),
                    '2', str(version)], check=True, capture_output=True)
    before = set(tmp_path.iterdir())
    counted = subprocess.run([str(COUNTER), str(fasta), str(sites),
                              str(version)], check=True, capture_output=True)
    result = json.loads(counted.stdout)
    actual = (index / 'fasta.fai').read_bytes()
    # Each physical FAI row has four uint64 numbers of at most20 decimal bytes,
    # four TABs and the platform text newline. Name bytes remain exact.
    upper = sum(len(row.split(b'\t')[0]) + 84 + len(os.linesep.encode())
                for row in actual.splitlines())
    assert result['fai_upper_bytes'] == upper
    assert upper >= len(actual)
    assert set(tmp_path.iterdir()) == before
