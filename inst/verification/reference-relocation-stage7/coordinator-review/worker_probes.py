"""Independent actual-worker configuration/closure probes; tiny local fixtures."""
import csv
import json
import shutil
from pathlib import Path

import pytest
import test_mss_worker
import test_mss_worker_relocation as supplied

import gwas_harvest.mss_worker as module
from gwas_harvest.pilot import ReferenceBundle

live_reference = test_mss_worker.live_reference
relocation_profile = supplied.relocation_profile


def assert_no_scientific_job(workspace):
    assert not (workspace / 'worker-receipt.json').exists()
    assert not list(workspace.glob('job-*'))
    assert not list(workspace.glob('measure-*'))


def test_same_path_retains_old_api_and_operational_identity_scope(tmp_path, live_reference):
    source = supplied._input(tmp_path / 'source.tsv')
    workspace = tmp_path / 'same-worker'
    with supplied._worker(live_reference, workspace) as worker:
        result = supplied._run(worker, source, tmp_path / 'edition')
        assert result['reference_open_count'] == 1
        assert worker.jobs == 1
    with (workspace / 'worker.tsv').open() as handle:
        config = {r['key']: r['value'] for r in csv.DictReader(handle, delimiter='\t')}
    assert config['runtime_fasta'] == next(p.path for p in live_reference.files if p.label == 'fasta')
    with (tmp_path / 'edition/metadata.tsv').open() as handle:
        metadata = {r['key']: r['value'] for r in csv.DictReader(handle, delimiter='\t')}
    assert not any(key.startswith('reference_fasta_') for key in metadata)
    receipt = json.loads((workspace / 'worker-receipt.json').read_text())
    assert receipt['reference_content_identity'] == live_reference.content_identity
    assert receipt['operational_reference'] == live_reference.model_dump(mode='json')


def test_available_recorded_fasta_cannot_replace_the_distinct_verified_pin(
    tmp_path, live_reference, relocation_profile
):
    value = live_reference.model_dump()
    pin = next(p for p in value['files'] if p['label'] == 'fasta')
    original = Path(pin['path'])
    physical = tmp_path / 'verified é physical copy.fa'
    shutil.copyfile(original, physical)
    pin['path'] = str(physical)
    bundle = ReferenceBundle.model_validate(value)
    assert bundle.content_identity == live_reference.content_identity
    bundle.verify()
    workspace = tmp_path / 'distinct-worker'
    if relocation_profile == 'unsupported':
        with pytest.raises(RuntimeError, match='explicit runtime FASTA'):
            supplied._worker(bundle, workspace)
        assert_no_scientific_job(workspace)
    else:
        with supplied._worker(bundle, workspace) as worker:
            output = tmp_path / 'edition'
            supplied._run(worker, supplied._input(tmp_path / 'source.tsv'), output)
        with (workspace / 'reference-files.tsv').open() as handle:
            actual = {r['label']: r['path'] for r in csv.DictReader(handle, delimiter='\t')}
        assert actual['fasta'] == str(physical)
        with (output / 'metadata.tsv').open() as handle:
            metadata = {r['key']: r['value'] for r in csv.DictReader(handle, delimiter='\t')}
        assert metadata['reference_fasta_recorded_path'] == str(original)
        assert metadata['reference_fasta_runtime_path'] == str(physical)
        assert original.exists()


@pytest.mark.parametrize('configured', ['missing-old-path', 'other-identical-copy'])
def test_stale_or_mismatched_runtime_config_never_starts_a_scientific_job(
    tmp_path, live_reference, relocation_profile, monkeypatch, configured
):
    transferred = supplied._relocate(live_reference, tmp_path / 'transferred')
    pinned = Path(next(p.path for p in transferred.files if p.label == 'fasta'))
    if configured == 'missing-old-path':
        override = tmp_path / 'tiny.fa'
        assert not override.exists()
        expected = 'path|file|exist'
    else:
        override = tmp_path / 'unselected but identical.fa'
        shutil.copyfile(pinned, override)
        expected = 'membership' if relocation_profile == 'supported' else 'explicit runtime FASTA'
    original_write = module._write_tsv
    def wrong_config(path, values):
        if path.name == 'worker.tsv':
            values = dict(values, runtime_fasta=str(override))
        return original_write(path, values)
    monkeypatch.setattr(module, '_write_tsv', wrong_config)
    workspace = tmp_path / 'stale-worker'
    with pytest.raises((ValueError, RuntimeError), match=expected):
        supplied._worker(transferred, workspace)
    assert_no_scientific_job(workspace)


def test_byte_identical_nonfasta_pin_outside_live_index_fails_closure(tmp_path, live_reference):
    value = live_reference.model_dump()
    pin = next(p for p in value['files'] if p['label'] == 'records')
    outside = tmp_path / 'outside-records'
    shutil.copyfile(pin['path'], outside)
    pin['path'] = str(outside)
    bundle = ReferenceBundle.model_validate(value)
    assert bundle.content_identity == live_reference.content_identity
    bundle.verify()
    workspace = tmp_path / 'wrong-closure-worker'
    with pytest.raises(ValueError, match='membership'):
        supplied._worker(bundle, workspace)
    assert_no_scientific_job(workspace)
