"""One guarded count-only qualification against the existing frozen 1M prefix.

No full input, active preparation or new reference index is opened/created.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import time

BASE = Path('/private/tmp/gwas-harvest-stage2')
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / 'results/verification-stage7'))
from guard_common import guarded


def pin(path):
    state = path.stat()
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    after = path.stat()
    assert (state.st_size, state.st_mtime_ns, state.st_ino) == (
        after.st_size, after.st_mtime_ns, after.st_ino
    ), 'Input changed while hashing'
    return {'path': str(path), 'bytes': state.st_size, 'sha256': digest,
            'inode': state.st_ino, 'mtime_ns': state.st_mtime_ns}


def main():
    os.umask(0o077)
    out = HERE / 'prefix-1M-v1'
    out.mkdir()
    inputs = json.loads((BASE / 'results/mss-context-format3-stage7/prefix-1M-v1/inputs.json').read_text())
    fasta, sites = Path(inputs['fasta']['path']), Path(inputs['variants']['path'])
    header = BASE / 'mss-context-format3-release/src/preservation_core.h'
    paths = [fasta, sites, header, HERE / 'count_reference', HERE / 'count_reference.cpp',
             HERE / 'test_count.py', Path(__file__)]
    before = [pin(path) for path in paths]
    for actual, expected in zip(before[:2], [inputs['fasta'], inputs['variants']], strict=True):
        assert actual['sha256'] == expected['sha256'] and actual['bytes'] == expected['size_bytes']
    assert before[2]['sha256'] == '4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc'
    profile = {
        'scope': 'Existing incomplete 1M prefix only; count diagnostic, not production reference/admission',
        'mode': 'Count only; no writer/sorter/index output', 'format_version': 3,
        'pins_before': before,
        'existing_external_input_bytes': before[0]['bytes'] + before[1]['bytes'],
        'external_input_accounting': 'Inputs remain existing frozen files; their disk is included in actual free space, not counted as a new owned copy',
        'small_files_upper_bytes': 25906 + 19 + 1024**2,
        'small_files_explanation': 'Pinned FAI25,906 + format marker19 + explicit1MiB manifest allowance, no contig-alias file in this prefix; full input requires its own exact allowance',
        'limits': {'tree_memory_bytes': 8 * 1024**3, 'minimum_free_bytes': 160 * 1024**3,
                   'owned_case_bytes': 16 * 1024**2, 'deadline_seconds': 600},
        'expected': {'source_rows': 1000000, 'records_bytes': 63883154, 'raw_sort_bytes': 100705522,
                     'keys_bytes': 75230231, 'final_key_rows': 3484369, 'offsets_bytes': 27874952,
                     'raw_key_rows': 4726213, 'site_records': 1000000, 'allele_records': 1242071,
                     'raw_alt_tokens': 1242077, 'rejected_bytes': 167},
    }
    (out / 'profile.json').write_text(json.dumps(profile, indent=2) + '\n')
    guard = guarded([str(HERE / 'count_reference'), str(fasta), str(sites), '3'],
                    out, out / 'guard.json', out / 'counts.log',
                    {'max_tree_memory_bytes': 8 * 1024**3,
                     'min_filesystem_free_bytes': 160 * 1024**3,
                     'poll_seconds': .1, 'max_case_disk_bytes': 16 * 1024**2},
                    time.monotonic() + 600, os.environ.copy())
    assert guard['status'] == 'success', 'Count phase did not qualify'
    actual = json.loads((out / 'counts.log').read_text())
    assert all(actual[name] == value for name, value in profile['expected'].items()), actual
    after = [pin(path) for path in paths]
    assert after == before, 'Pinned input/source/binary changed'
    (out / 'pins-after.json').write_text(json.dumps(after, indent=2) + '\n')
    r, u, k, n, j, c, t = (actual['records_bytes'], actual['raw_sort_bytes'], actual['keys_bytes'],
                           actual['final_key_rows'], actual['rejected_bytes'],
                           profile['small_files_upper_bytes'], actual['source_bytes'])
    retained = t + r + k + 8*n + j + c
    envelope = max(retained, r + 2*u + j + c, r + u + k + 8*n + j + c)
    result = {'status': 'success', 'counts': actual, 'retained_logical_upper_bytes': retained,
              'phase_logical_upper_bytes': envelope,
              'not_an_admission': 'Allocation/filesystem allowance, concurrent reservations and170GiB planning reserve remain separate; both full source preparations remain unsealed',
              'guard': guard, 'all_input_source_binary_pins_unchanged': True}
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'success', 'source_rows': actual['source_rows'],
                      'retained_upper_bytes': retained, 'phase_upper_bytes': envelope,
                      'guard_seconds': guard['wall_seconds']}))


if __name__ == '__main__':
    main()
