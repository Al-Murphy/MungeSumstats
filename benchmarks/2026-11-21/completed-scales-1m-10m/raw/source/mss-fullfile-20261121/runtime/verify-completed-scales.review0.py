"""Verify the retained 1M and 10M pairs from the interrupted cold matrix."""

import argparse
from datetime import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import stat
import sys


SOURCE = Path('/private/tmp/mss-fullfile-20261121/runtime/verify-completed-scales.py')
WARM_SOURCE = SOURCE.with_name('warm-current.py')
WARM_SHA256 = 'c58b14cda2ad650cefa078e4a73cd2a2a4d5ab8ffc4f7477e5375e0aaa794220'
ORIGINAL = Path('/private/tmp/mss-fullmatrix-20261121')
COLD = ORIGINAL / 'cold-01'
MATRIX = COLD / 'matrix'
CONFIG = ORIGINAL / 'prep-01/matrix.json'
PROFILE = CONFIG.with_name('snp-oracle.json')
FAILURE = ORIGINAL / 'cold-failure-observation.json'
CONFIG_SHA256 = '156143e169744ead89eb7a01721ac433da0d0983be68924b417bcc72f00426c7'
PROFILE_SHA256 = 'b387920dbab4d53af316a8b2b320d577be463abdad7ce80dbf8b9de411e303d0'
FAILURE_SHA256 = '39a96694ac60f09bfad5e662c24ec2759fc88823149b2e94d2e3bd1628ecf995'
OPENING_SHA256 = 'a757eae6a2dc5f2738fb0b03f865a6bb31cb7b2a4540999bebb737e09fea6350'
OUTPUT = Path('/private/tmp/mss-completed-scales-20261121/verify-01')
SCOPED = {
    'scope': 'completed-1M-and-10M-only',
    'original_full_matrix_complete': False,
    'warm_eligible': False,
    'independent_oracles_pending': ['1M', '10M'],
}
SCALES = {'1M': 1_000_000, '10M': 10_000_000}
ORACLE_ARTIFACTS = ('source.tsv', 'source-artifact', 'normalized.tsv')
OUTPUT_FILES = ('current-verification-pins.json', 'scoped-summary.json', 'oracle-window-plan.json', 'RESULT.json')


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def load_frozen():
    require(Path(__file__) == SOURCE and sys.flags.isolated and sys.dont_write_bytecode and sys.flags.no_site,
            'Expected the physical source with isolated Python -I -B -S')
    require(WARM_SOURCE.resolve(strict=True) == WARM_SOURCE, 'Indirect helper source')
    fd = os.open(WARM_SOURCE, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    with os.fdopen(fd, 'rb') as stream:
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and 0 < info.st_size <= 1024**2, 'Invalid helper source')
        data = stream.read(1024**2 + 1)
    require(len(data) == info.st_size and hashlib.sha256(data).hexdigest() == WARM_SHA256,
            'Changed helper source')
    spec = importlib.util.spec_from_file_location('completed_scales_frozen', WARM_SOURCE)
    warm = importlib.util.module_from_spec(spec)
    exec(compile(data, str(WARM_SOURCE), 'exec'), warm.__dict__)
    warm.modules()
    return warm


def read_pinned(controls, path, expected):
    data = controls.data(path)
    require(controls.files[str(path)] == expected, f'Unexpected control digest: {path}')
    return json.loads(data)


def failure_evidence(warm, controls):
    failure = read_pinned(controls, FAILURE, FAILURE_SHA256)
    admission_path, receipt_path = (COLD / name for name in ('admission.json', 'matrix.jsonl'))
    admission = read_pinned(controls, admission_path, warm.ACCEPTED[admission_path])
    records = [json.loads(line) for line in controls.data(receipt_path).splitlines()]
    require(failure['schema'] == 'current-cold-failure-observation-v1'
            and failure['root'] == warm.root_identity(COLD)
            and failure['admission'] == controls.pin(admission_path)
            and failure['receipt'] == controls.pin(receipt_path)
            and failure['receipt_records'] == len(records) == 2
            and failure['terminal_receipt_present'] is False
            and failure['observed_supervisor_exit'] == 1, 'Original guard failure binding differs')
    planned, running = records
    outcome, absence = failure['returned_outcome'], failure['independent_absence']
    argv = [warm.COLD_PY, '-E', '-s', '-B', str(warm.BENCH / 'run.py'), 'run', str(CONFIG),
            str(MATRIX), '--suite', 'matrix', '--profile', 'none']
    require(planned['status'] == 'planned' and running['status'] == 'running'
            and planned['label'] == running['label'] == 'full-cold-matrix'
            and planned['argv'] == admission['argv'] == argv and planned['cwd'] == str(COLD)
            and planned['receipt'] == str(receipt_path) and planned['allowed_exits'] == [0, 2, 3]
            and admission['config_sha256'] == CONFIG_SHA256
            and admission['guard_sha256'] == warm.SOURCES[warm.GUARD]
            and planned['limits']['floor_bytes'] == admission['guard_floor_bytes']
            and running['uid'] == os.geteuid()
            and running['pid'] == running['expected_pgid'] == running['expected_sid']
            == absence['pgid'] == absence['sid'] > 1, 'Original launch binding differs')
    require(outcome['status'] == 'error' and outcome['closure'] == 'confirmed' and outcome['raw_exit'] == -15
            and outcome['primary']['reason'] == 'filesystem_floor'
            and outcome['receipt_error']['reason'] == 'receipt_reserve' and outcome['cleanup_errors'] == []
            and absence['group_absent'] is True and absence['session_absent'] is True
            and absence['matches'] == absence['permission_errors'] == [], 'Expected recorded disk-floor failure')
    for stream in ('stdout', 'stderr'):
        path = COLD / f'matrix.{stream}'
        require(planned[stream] == str(path), f'Unexpected original {stream} path')
        controls.pin(path)
    require(all(not os.path.lexists(MATRIX / name) for name in ('post-pins.json', 'summary.json')),
            'Unexpected original closing files')
    return failure, running


def select_cases(config, roster):
    cases = {case['id']: case for case in config['cases']}
    require(set(cases) == set(roster) and len(cases) == len(config['cases']) == 19, 'Original matrix roster differs')
    selected = {}
    for size, rows in SCALES.items():
        name = f'matched-{size}'
        case = cases[name]
        ref = config['references'][case['reference']]
        require((case['rows'], case['pairs'], ref['format']) == roster[name] == (rows, 3, 2)
                and case['target'] is True and case['settings']['delimiter'] == 'tab' and 'oracle' not in case,
                f'Unexpected selected case: {name}')
        selected[size] = case
    return selected


def worker_dirs(case, roles):
    work = MATRIX / case['id']
    runs = []
    for pair in range(1, case['pairs'] + 1):
        for role in roles if pair % 2 else roles[::-1]:
            dest = work / f'pair-{pair:02d}-{role}'
            for path in (dest, dest / 'edition'):
                require(path.is_dir() and path.resolve(strict=True) == path, f'Missing or indirect retained directory: {path}')
            runs.append((pair, role, dest))
    return runs


def live_inventory(warm, config, resolved):
    hashes, packages, builds = resolved
    assets = warm.r.asset_inventory(hashes)
    installed, toolchain = warm.r.read_provenance(config, packages, builds, assets)
    deps = warm.dependencies(config)
    require(installed == deps['role_trees'], 'Role trees differ from accepted dependency inventory')
    return {'assets': assets, 'installed_packages': installed, 'r_toolchain': toolchain, 'dependencies': deps}


def collect_workers(warm, controls, config, cases, runs, current, opening, failure):
    observations = {size: [] for size in cases}
    previous = datetime.fromisoformat(opening['started_utc'])
    end = datetime.fromisoformat(failure['observed_utc'])
    shared_namespaces = None
    for size, case in cases.items():
        for pair, role, dest in runs[size]:
            _, process, metrics, namespaces = warm.checked_run(
                controls, config, case, dest, role, current['assets'], current['dependencies'],
                current['r_toolchain'], warm.WORKER)
            started, finished = (datetime.fromisoformat(process[f'{label}_utc']) for label in ('started', 'finished'))
            require(previous <= started < finished <= end, f'Nonserial or out-of-order worker: {dest}')
            previous = finished
            require(shared_namespaces is None or namespaces == shared_namespaces, f'Namespace drift: {dest}')
            shared_namespaces = namespaces
            for phase in ('package_setup', 'reference_open', 'processing'):
                for stat_name in ('elapsed', 'user.self', 'sys.self', 'user.child', 'sys.child'):
                    value = float(metrics[f'{phase}.{stat_name}'])
                    require(math.isfinite(value) and value >= 0, f'Invalid retained timing: {dest}')
            for name in ('process.stdout.log', 'process.stderr.log'):
                controls.pin(dest / name)
            observations[size].append({'directory': str(dest), 'pair': pair, 'process': process, 'metrics': metrics})
    return observations, shared_namespaces


def status(checks):
    return 'miss' if False in checks else 'indeterminate' if None in checks else 'pass'


def case_summary(warm, case, observations):
    summary = warm.r.summarize_case(case, observations, suite='matrix')
    require(summary['target'] is True, f'Selected target missing: {case["id"]}')
    ratio = summary['ratio_of_median_processing_seconds']
    checks = (None if ratio is None else ratio >= 2, summary['median_peak_rss_nonregression'])
    baseline, candidate = (summary['medians'][role] for role in warm.r.ROLES)
    reductions = {}
    for key in ('median_peak_rss_bytes', 'max_peak_rss_bytes'):
        reductions[key] = {'bytes': baseline[key] - candidate[key],
                           'fraction': (baseline[key] - candidate[key]) / baseline[key]}
    absolute = warm.d.summarize_rss([row['process'] for row in observations if row['process']['role'] == 'candidate'])
    return {**summary, 'target_status': status(checks), 'rss_reduction': reductions,
            'candidate_8GiB_absolute': {**absolute, 'limit_bytes': warm.d.CANDIDATE_RSS}}


def oracle_window(profile, size, case, cfg_ref, edition, inventory, assets, digest_pin):
    ref = profile['reference']
    refdir = Path(cfg_ref['directory'])
    require(profile['stages'][size] == {'rows': case['rows'], 'input': case['input']}
            and cfg_ref['format'] == 2 and ref == {'build': 'GRCh37', 'source': str(refdir / 'source.tsv'),
                                                  'fai': str(refdir / 'fasta.fai'), 'fasta': cfg_ref['fasta']['path']}
            and all(assets[str(refdir / name)]['sha256'] == cfg_ref['files'][name] for name in ('source.tsv', 'fasta.fai'))
            and assets[ref['fasta']]['sha256'] == cfg_ref['fasta']['sha256'], f'Original oracle profile differs: {size}')
    original_input = case['input']['path']
    info = edition.stat()
    require(edition.resolve(strict=True) == edition and stat.S_ISDIR(info.st_mode), 'Oracle edition path differs')
    return {'stage': size, 'rows': case['rows'], 'selector': f'mss7={edition}', 'edition': str(edition),
            'edition_identity': {'dev': info.st_dev, 'ino': info.st_ino}, 'artifact_digests': digest_pin,
            'files': {original_input: assets[original_input], **{
                str(edition / name): inventory[str(edition.parent)][name] for name in ORACLE_ARTIFACTS}}}


def verify(output):
    require(output == OUTPUT, f'Expected verification output root: {OUTPUT}')
    warm = load_frozen()
    d, r, v = warm.d, warm.r, warm.v
    output_identity = warm.root_identity(output)
    matrix_identity = warm.root_identity(MATRIX)
    require(all(not os.path.lexists(output / name) for name in (*OUTPUT_FILES, *(f'matched-{s}' for s in SCALES))),
            'Verifier output already exists')
    original = warm.Controls()
    evidence = warm.Controls()
    sources = {str(path): sha for path, sha in warm.SOURCES.items()}
    sources[str(WARM_SOURCE)] = WARM_SHA256
    sources[str(SOURCE)] = evidence.pin(SOURCE)['sha256']
    config = read_pinned(original, CONFIG, CONFIG_SHA256)
    profile = read_pinned(original, PROFILE, PROFILE_SHA256)
    failure, running = failure_evidence(warm, original)
    opening_path = MATRIX / 'pins.json'
    opening = read_pinned(original, opening_path, OPENING_SHA256)
    require(opening['suite'] == 'matrix' and opening['profile'] == 'none'
            and opening['configuration_file'] == {'path': str(CONFIG), **r.digest(CONFIG)}
            and running['utc'] <= opening['started_utc'] <= failure['observed_utc'], 'Original scope or configuration differs')
    resolved = r.resolve_configuration(config, 'matrix')
    require(config == opening['configuration'], 'Original resolved configuration differs')
    cases = select_cases(config, d.expected_roster())
    runs = {size: worker_dirs(case, r.ROLES) for size, case in cases.items()}
    for pin in (config['dependency_inventory'], *(config[role]['build_receipt'] for role in r.ROLES)):
        read_pinned(original, Path(pin['path']), pin['sha256'])
    harness = {name: r.digest(warm.BENCH / name) for name in ('run.py', 'verify.py', 'worker.R', 'process_metrics.py')}
    require(opening['harness'] == harness and all(row['sha256'] == sources[str(warm.BENCH / name)]
                                                 for name, row in harness.items()), 'Frozen harness differs')
    current = live_inventory(warm, config, resolved)
    require(all(current[key] == opening[key] for key in ('assets', 'installed_packages', 'r_toolchain')),
            'Current assets, installations or toolchain differ from original opening')
    observations, namespaces = collect_workers(warm, original, config, cases, runs, current, opening, failure)
    cardinalities = {}
    for size, case in cases.items():
        pin = case['cardinality_receipt']
        cardinal = read_pinned(original, Path(pin['path']), pin['sha256'])
        require(cardinal == opening['cardinalities'][case['id']] == r.cardinality(case), 'Cardinality differs')
        cardinalities[size] = cardinal
    current_path, summary_path, oracle_path, result_path = (output / name for name in OUTPUT_FILES)
    r.write_json(current_path, {'schema': 'current-completed-scales-verification-pins-v1', **SCOPED,
                                'observed_utc': r.utc_now(), 'original_opening': original.pin(opening_path),
                                'configuration': original.pin(CONFIG), 'harness': harness,
                                'observed_worker_namespaces': namespaces, **current})
    current_pin = evidence.pin(current_path)
    summaries, windows, selected_evidence = {}, {}, {}
    for size, case in cases.items():
        work = output / case['id']
        work.mkdir(mode=0o700)
        directories = [dest for _, _, dest in runs[size]]
        ref = config['references'][case['reference']]
        assets = current['assets']
        verified = v.verify_editions(case, ref, directories, work, assets, config['edition_encoding'])
        r.verify_manifests(config, work, directories)
        verified['stable_rds_manifest_equality'] = True
        r.write_json(work / 'verification.json', verified)
        verified, inventory = warm.checked_artifacts(evidence, config, case, directories, work, assets, rehash=False)
        edition = MATRIX / case['id'] / 'pair-01-candidate/edition'
        windows[size] = oracle_window(profile, size, case, ref, edition, inventory, assets,
                                      evidence.pin(work / 'artifact-digests.json'))
        selected_evidence[case['id']] = {path.name: evidence.pin(path) for path in (
            work / 'verification.json', work / 'artifact-digests.json', work / 'stable-manifests.txt',
            work / 'manifest-verification.json', work / 'manifest-verification.stdout.log',
            work / 'manifest-verification.stderr.log')}
        summaries[case['id']] = {**case_summary(warm, case, observations[size]), 'rows': case['rows'],
                                 'cardinality': cardinalities[size], 'observations': observations[size],
                                 'verification': verified}

    require(live_inventory(warm, config, resolved) == current, 'Scientific or installation inventories changed')
    require(warm.root_identity(output) == output_identity and warm.root_identity(MATRIX) == matrix_identity,
            'Output or original matrix root changed')
    require(failure_evidence(warm, original) == (failure, running), 'Original guard failure changed')
    for path, sha in {**original.files, **evidence.files, **sources}.items():
        d.check_control(warm.pin(path, sha))

    failure_pin = original.pin(FAILURE)
    checks = {'pass': True, 'miss': False, 'indeterminate': None}
    target_status = status([checks[case['target_status']] for case in summaries.values()])
    absolute = d.summarize_rss([row['process'] for rows in observations.values() for row in rows
                                if row['process']['role'] == 'candidate'])
    r.write_json(summary_path, {'schema': 'current-completed-scales-summary-v1', **SCOPED,
                                'original_failure': failure_pin, 'suite': 'matrix', 'profile': 'none',
                                'timing_source': 'original serial fresh-worker receipts',
                                'target_rule': 'Processing ratio of medians >= 2 and median peak RSS nonregression.',
                                'target_status': target_status, 'processing_speed_and_rss_target_met': checks[target_status],
                                'candidate_8GiB_absolute': {**absolute, 'limit_bytes': d.CANDIDATE_RSS},
                                'cases': summaries})
    r.write_json(oracle_path, {'schema': 'current-completed-scales-oracle-window-plan-v1', **SCOPED,
                               'configuration': original.pin(CONFIG), 'profile': original.pin(PROFILE),
                               'original_failure': failure_pin, 'oracle_source': warm.pin(warm.ORACLE, sources[str(warm.ORACLE)]),
                               'windows': windows})
    r.write_json(result_path, {'schema': 'current-completed-scales-verification-result-v1', **SCOPED,
                                'status': 'retained-editions-verified-pending-oracles-and-observer',
                                'verified_scales': list(SCALES), 'verified_workers': sum(map(len, observations.values())),
                                'root_observer_pending': True,
                                'output_root': output_identity, 'observed_utc': r.utc_now(),
                                'source': warm.pin(SOURCE, sources[str(SOURCE)]), 'sources': sources,
                                'configuration': original.pin(CONFIG), 'original_failure': failure_pin,
                                'original_guard_outcome': failure['returned_outcome'],
                                'original_supervisor_exit': failure['observed_supervisor_exit'],
                                'original_terminal_receipt_present': failure['terminal_receipt_present'],
                                'original_opening': original.pin(opening_path), 'original_small_controls': original.files,
                                'current_verification_pins': current_pin, 'summary': evidence.pin(summary_path),
                                'oracle_window_plan': evidence.pin(oracle_path), 'selected_case_evidence': selected_evidence,
                                'current_inventories_equal_at_end': True, 'target_status': target_status})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='Private output root created by the guard operator')
    return verify(parser.parse_args().output)


if __name__ == '__main__':
    raise SystemExit(main())
