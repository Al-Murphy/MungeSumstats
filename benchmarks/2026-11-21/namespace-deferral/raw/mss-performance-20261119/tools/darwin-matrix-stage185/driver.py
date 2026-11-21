"""Guarded local orchestration for the frozen MSS matrix and warm companion."""
from datetime import datetime, timezone
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import stat
import sys
import time


HERE = Path(__file__).resolve().parent
TREE = HERE.parent.parent / 'repo'
BENCH = TREE / 'benchmarks/2026-11-19/performance-upgrade'
ORACLE = TREE / 'inst/verification/fai-safety-stage8/whole-file/harness/verify.py'
ASSETS = Path('/private/tmp/gwas-harvest-stage2/assets')
PENDING_SHA256 = '57b3c1db7909e13708b47faf7868c7e960f0dac0eb452daa41570b3f894feab9'
NATIVE_PROOF_SOURCE = {'path': str(HERE / 'native-proof.py'), 'bytes': 17672, 'sha256': '6af8d7650a64edf85ec699dba532f726fb548e26f13691cac0e89a030aad5b70'}
NATIVE_BINDING = {
    "externals": {
        "/opt/homebrew/bin/python3": "6cc907e4a022fefc843bcd93954eaf5a2cbfc77a095640bfadeeed96154190c5",
        "/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/process_metrics.py": "f402957540febe6e0f9e8c56422b7acc5f68bcc32afe4377705347de199a05d7"
    },
    "inventory": {
        "test_guard.AdmissionTests": 2,
        "test_guard.CreationBudgetTests": 2,
        "test_guard.DescendantObserverTests": 5,
        "test_guard.DeterministicTests": 12,
        "test_guard.OwnedProcessTests": 8,
        "test_guard.ReviewRegressionTests": 9,
        "test_guard.TerminalObserverTests": 4,
        "test_review_fix2.ReviewFix2": 4,
        "test_wait_status.WaitStatusTests": 12
    },
    "methods": [
        "test_guard.AdmissionTests.test_future_and_oversized_admissions_refuse_run",
        "test_guard.AdmissionTests.test_stale_admission_refuses_run",
        "test_guard.CreationBudgetTests.test_logical_cap_blocks_each_create_before_the_entry_exists",
        "test_guard.CreationBudgetTests.test_pending_entry_and_record_require_joint_reserve",
        "test_guard.DescendantObserverTests.test_archived_descendant_state_requires_fresh_strict_snapshot",
        "test_guard.DescendantObserverTests.test_cleanup_rechecks_a_previously_escaped_descendant",
        "test_guard.DescendantObserverTests.test_only_observed_descendant_can_trigger_strict_requery",
        "test_guard.DescendantObserverTests.test_parent_reaps_descendant_before_fresh_native_query",
        "test_guard.DescendantObserverTests.test_partial_rss_and_escaped_members_remain_failures",
        "test_guard.DeterministicTests.test_cadence_failure_requires_two_consecutive_misses",
        "test_guard.DeterministicTests.test_child_wait_and_reap_preserve_exact_status_once",
        "test_guard.DeterministicTests.test_environment_rejection_does_not_launch",
        "test_guard.DeterministicTests.test_group_parser_reports_unavailable_and_malformed_evidence",
        "test_guard.DeterministicTests.test_metadata_counts_symlink_without_following_and_tolerates_vanishing_entry",
        "test_guard.DeterministicTests.test_path_controls_refuse_aliases_collisions_and_wrong_volume",
        "test_guard.DeterministicTests.test_prelaunch_filesystem_boundaries_record_refusals",
        "test_guard.DeterministicTests.test_prelaunch_receipt_failure_preserves_both_errors_and_blocks",
        "test_guard.DeterministicTests.test_shared_expired_deadline_never_spawns",
        "test_guard.DeterministicTests.test_signal_esrch_requires_recheck_and_eperm_stays_explicit",
        "test_guard.DeterministicTests.test_snapshot_bounds_races_and_rejects_changed_sessions",
        "test_guard.DeterministicTests.test_unconfirmed_esrch_stays_unclosed",
        "test_guard.OwnedProcessTests.test_a_fast_exits_and_signal_preserve_raw_wait4_status",
        "test_guard.OwnedProcessTests.test_foreground_descendant_is_sampled_and_can_exit_cleanly",
        "test_guard.OwnedProcessTests.test_leftover_descendant_escalates_kill_and_decoy_is_left_alive",
        "test_guard.OwnedProcessTests.test_rss_guard_triggers_and_keeps_sample_separate",
        "test_guard.OwnedProcessTests.test_shutdown_during_phase_uses_cleanup_and_marks_next_pending",
        "test_guard.OwnedProcessTests.test_term_aware_child_closes_before_escalation",
        "test_guard.OwnedProcessTests.test_uncertain_closure_uses_only_direct_child_and_blocks_pending",
        "test_guard.OwnedProcessTests.test_wait4_metrics_separately_report_bytes_and_exact_signal",
        "test_guard.ReviewRegressionTests.test_final_receipt_cannot_cross_logical_allocated_or_free_space_bounds",
        "test_guard.ReviewRegressionTests.test_observer_timeout_preserves_secondary_helper_cleanup_failure",
        "test_guard.ReviewRegressionTests.test_oversized_plan_never_spawns_or_records_unbounded_input",
        "test_guard.ReviewRegressionTests.test_partial_terminal_write_preserves_primary_and_blocks",
        "test_guard.ReviewRegressionTests.test_receipt_close_failure_preserves_primary_and_blocks_next_phase",
        "test_guard.ReviewRegressionTests.test_receipt_open_failure_blocks_next_phase_and_preserves_existing_file",
        "test_guard.ReviewRegressionTests.test_reserve_prevents_spawn_but_allows_small_failure_receipt",
        "test_guard.ReviewRegressionTests.test_terminal_second_late_sample_fails_with_confirmed_raw_exit",
        "test_guard.ReviewRegressionTests.test_third_party_owner_rejected_even_for_0755_ancestor",
        "test_guard.TerminalObserverTests.test_recheck_preserves_real_group_and_observer_errors",
        "test_guard.TerminalObserverTests.test_rejected_identity_diagnostic_never_includes_arbitrary_text",
        "test_guard.TerminalObserverTests.test_transition_requires_complete_matching_single_row_and_waitable_anchor",
        "test_guard.TerminalObserverTests.test_transition_requires_fresh_native_zombie_and_exact_reap",
        "test_review_fix2.ReviewFix2.test_fresh_refusal_preserves_measured_rss_in_terminal_receipt",
        "test_review_fix2.ReviewFix2.test_malformed_companion_cannot_clear_escaped_sibling_at_cleanup",
        "test_review_fix2.ReviewFix2.test_partial_failed_suite_queries_only_recorded_identities",
        "test_review_fix2.ReviewFix2.test_suite_failure_and_timeout_preserve_status_and_execute_both_postchecks",
        "test_wait_status.WaitStatusTests.test_changed_pid_or_status_refuses",
        "test_wait_status.WaitStatusTests.test_cleanup_does_not_poll_or_use_popen_kill",
        "test_wait_status.WaitStatusTests.test_fixture_extent_mutations_refuse_before_write",
        "test_wait_status.WaitStatusTests.test_foreign_pid_never_reaped_or_signalled",
        "test_wait_status.WaitStatusTests.test_helper_positive_kernel_status",
        "test_wait_status.WaitStatusTests.test_pending_sleep_interrupt_preserves_cleanup",
        "test_wait_status.WaitStatusTests.test_plausible_output_with_echild_refuses",
        "test_wait_status.WaitStatusTests.test_scoped_absence_needs_observed_one_and_empty_streams",
        "test_wait_status.WaitStatusTests.test_timeout_preserves_cleanup_unknown",
        "test_wait_status.WaitStatusTests.test_workflow_observed_signal_and_nonzero_preserved",
        "test_wait_status.WaitStatusTests.test_workflow_synthetic_zero_refuses",
        "test_wait_status.WaitStatusTests.test_workflow_timeout_runs_postchecks"
    ],
    "phase_cases": {
        "SIGHUP": "test_shutdown_during_phase_uses_cleanup_and_marks_next_pending/SIGHUP.jsonl",
        "SIGINT": "test_shutdown_during_phase_uses_cleanup_and_marks_next_pending/SIGINT.jsonl",
        "SIGTERM": "test_shutdown_during_phase_uses_cleanup_and_marks_next_pending/SIGTERM.jsonl",
        "allowed-two": "test_a_fast_exits_and_signal_preserve_raw_wait4_status/allowed-two.jsonl",
        "descendant-transition": "test_parent_reaps_descendant_before_fresh_native_query/descendant-transition.jsonl",
        "fast": "test_a_fast_exits_and_signal_preserve_raw_wait4_status/fast.jsonl",
        "foreground": "test_foreground_descendant_is_sampled_and_can_exit_cleanly/foreground.jsonl",
        "leftover": "test_leftover_descendant_escalates_kill_and_decoy_is_left_alive/leftover.jsonl",
        "lost": "test_uncertain_closure_uses_only_direct_child_and_blocks_pending/lost.jsonl",
        "nonzero": "test_a_fast_exits_and_signal_preserve_raw_wait4_status/nonzero.jsonl",
        "rss": "test_rss_guard_triggers_and_keeps_sample_separate/rss.jsonl",
        "signal": "test_a_fast_exits_and_signal_preserve_raw_wait4_status/signal.jsonl",
        "term-aware": "test_term_aware_child_closes_before_escalation/term-aware.jsonl",
        "transition": "test_transition_requires_fresh_native_zombie_and_exact_reap/transition.jsonl"
    },
    "root": "/private/tmp/mss-performance-20261119/tools/darwin-full-stage185",
    "sources": {
        "full-native-suite.py": [
            1176,
            "acb6c3b7c418a05995200eeca2dabbfcbbf9b1b27fa26b8554ed2ad0c01de906"
        ],
        "guard.py": [
            59374,
            "58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2"
        ],
        "native-body.sh": [
            427,
            "3c878ad5368ec65d69aa1eef5f5424a1344c1d993a80a9ec10b2c2c66c7827d5"
        ],
        "native-pins.sha256": [
            1079,
            "bc3b70fcdd5deb06b8e82313e075bf4ed1578a86c79707d3977b31d7c2eb2323"
        ],
        "native-verify.py": [
            14528,
            "e0d9e37f4537c7a91b72004a7e14c299124af615af0f8ff3bbcdff7d87d5aac1"
        ],
        "pure_suite.py": [
            16374,
            "baae077349e2f8261c87fd9fa702f6c0768766a624feb51222df207af6db3266"
        ],
        "run-native.sh": [
            149,
            "e27ef8538ed028ba6f01fcd34a7ef6e42d1c6c3ed4052bf265a38517a17d8810"
        ],
        "run-tests.sh": [
            992,
            "2753cb8cc28c3c49c9e86ef3282bb13581ba7513b864d126508a281b79dcf07a"
        ],
        "test_guard.py": [
            64606,
            "0badfc4b6b8beadf3bdd9bbfecef739cd59d2ca0a9158135eaf09ccbc8415838"
        ],
        "test_review_fix2.py": [
            8250,
            "13988ab007f6d9273371887523258c5853e2c63d7469a9636ab4d6bea2ddeee8"
        ],
        "test_wait_status.py": [
            25753,
            "16ae8e8621fa302e9413f61d4195a934ad3d08794f9190e80f153a5e527d9893"
        ]
    }
}
PRIMARY = Path('/private/tmp/mss-matrix-primary-stage185')
WARM = Path('/private/tmp/mss-matrix-warm-stage185')
BASE = 'a86e8c94d753630eff8c6a695ed2f240d1123d0d'
CANDIDATE = '737ad24784ee9a2445a3eaf4f8262c87561a72f7'
CANDIDATE_TREE = '785ce7bb01ae6b43417432909681bb355cc8e599'
HARNESS = 'bd82dfb9f7830126ba385283b1f1ffb9f98d757a'
ROLES = ('baseline', 'candidate')
WARM_ORDER = ('1M', '10M', '40M')
SCALES = dict(zip(WARM_ORDER, ((1_000_000, 3), (10_000_000, 3), (40_000_000, 1)), strict=True))
SUPPLEMENTS = {'distinct-shuffled': 437458, 'coordinate-only': 437458, 'id-only': 437458,
               'negative': 437458, 'dense': 10000, 'mixed': 100000, 'hard-cases': 32, 'declared': 9}
FLOOR = 186 * 1024**3
DISK = 60_000_000_000
GROUP_RSS = 24_000_000_000
CANDIDATE_RSS = 8 * 1024**3
BUDGETS = {'primary': {'identity': 10800, 'prep-diagnostics': 7200,
                       'matrix': 25200, 'oracles': 14400},
           'warm': {'identity': 7200, 'warm': 14400, 'warm-verify': 14400}}
RUNTIME = ('r_executable', 'r_engine', 'runtime_files')
COMMON_TOOLCHAIN = ('compiler', 'compiler_flags', 'R_version', 'R_platform')
ORDERS = ('MungeSumstats', 'testthat')
CANDIDATE_REGRESSION_FILE = 'test-preservation-stream-performance.R'
# Match the exact named assertions in the pinned qualification producer.
CANDIDATE_REGRESSIONS = frozenset({
    'file errors respect raw block and deferred-row boundaries',
    'file publication no longer materializes legacy row lists',
    'column batches match row decoding and row/byte membership at every split',
    'collectors reject mixed use and invalid source records even in discarded columns',
    'bounded writer matches exact baseline bytes across block boundaries',
    'unknown converting connections and marked strings use the compatible formatter',
    'deferred character elements stay live across native block flushes',
    'file and resident batches retain source contracts, classes and edition bytes',
})
FUNCTIONAL_PHASES = frozenset({
    'inputs-before', 'inputs-after', 'compiler-version', 'compiler-sdk', 'stage-test-sources',
    'compare-lookup-parity', 'smoke-config', 'scientific-smoke', 'smoke-verification',
    *(f'{role}-{case}' for role in ROLES for case in ('native', 'preservation', 'registered',
                                                   'testthat', 'lookup-build', 'lookup')),
})
HARD_PIN_NAMES = (
    'python', 'cmp', 'guard', 'functional_runtime_guard', 'historical_reference_guard', 'installer',
    'reference_plan', 'baseline', 'candidate',
    'focused', 'candidate_compiler', 'candidate_sdk', 'candidate_stdout', 'candidate_stderr',
    'run', 'verify', 'process_metrics', 'worker', 'template', 'assets', 'build', 'science',
    'references', 'profile', 'oracle',
)


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def control_bytes(path):
    """A bounded control read. Payload hashing belongs to dispatched phases."""
    def identity(info):
        return (info.st_dev, info.st_ino, info.st_mode, info.st_nlink, info.st_uid,
                info.st_size, info.st_mtime_ns, info.st_ctime_ns)

    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    with os.fdopen(fd, 'rb') as stream:
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and info.st_size <= 1024**2, f'Invalid control file: {path}')
        data = stream.read(1024**2 + 1)
        require(identity(info) == identity(os.fstat(stream.fileno())) == identity(Path(path).lstat())
                and len(data) == info.st_size,
                f'Control changed during read: {path}')
    return data


def pin_control(path):
    return {'path': str(path), 'sha256': sha(control_bytes(path))}


def check_control(pin):
    path = Path(pin['path'])
    require(path.is_absolute() and path.resolve(strict=True) == path, f'Control path changed: {path}')
    data = control_bytes(path)
    require(sha(data) == pin['sha256'] and pin.get('bytes', len(data)) == len(data),
            f'Control pin differs: {path}')
    return data


def private_grant_bytes(path):
    path = Path(path)
    require(path.is_absolute() and path.resolve(strict=True) == path, 'Grant must have a physical path')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    with os.fdopen(fd, 'rb') as stream:
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and stat.S_IMODE(info.st_mode) == 0o600
                and info.st_uid == os.geteuid() and info.st_nlink == 1 and 0 < info.st_size <= 1024**2,
                'Expected private single-link owner grant')
        data = stream.read(info.st_size + 1)
        current = path.lstat()
        require(len(data) == info.st_size and current == os.fstat(stream.fileno())
                and (current.st_dev, current.st_ino, current.st_ctime_ns, current.st_mtime_ns)
                == (info.st_dev, info.st_ino, info.st_ctime_ns, info.st_mtime_ns), 'Grant changed during read')
    return data


def encode(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def module(pin, name):
    data = check_control(pin)
    spec = importlib.util.spec_from_file_location(name, pin['path'])
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    exec(compile(data, pin['path'], 'exec'), loaded.__dict__)
    return loaded


def expected_roster():
    return {**{f'matched-{size}': (rows, pairs, 2) for size, (rows, pairs) in SCALES.items()},
            **{f'{name}-f{fmt}': (rows, 1, fmt) for name, rows in SUPPLEMENTS.items() for fmt in (2, 3)}}


def require_functional(plan):
    proof = plan.get('functional')
    require(isinstance(proof, dict) and set(proof) == {'qualification', 'root_observation', 'lookup'},
            'Matrix disabled until actual functional qualification and independent observation are pinned')
    return proof


def phase_status(row, allowed):
    raw = row['raw_exit']
    require(row['status'] == 'completed' and row['closure'] == 'confirmed'
            and row['leader_reaped'] is True and row['primary'] is None
            and not row['cleanup_errors'] and row.get('receipt_error') is None and type(raw) is int and raw in allowed,
            f'Phase closure, ordinary status or receipt is incomplete: {raw}')
    return raw


def validate_dispatch(dispatch, rows, argv, process, now):
    require(len(rows) == 2, 'Expected a live dispatch receipt without a terminal record')
    planned, running = rows
    parent, pid, pgid, sid, uid = process
    require(planned['status'] == 'planned' and running['status'] == 'running'
            and planned['label'] == running['label'] == dispatch['label']
            and planned['argv'] == argv and planned['cwd'] == dispatch['root']
            and planned['allowed_exits'] == dispatch['allowed'] and now < planned['monotonic_deadline']
            and parent == dispatch['parent_pid'] > 1 and uid == running['uid']
            and pid == pgid == sid == running['pid'] == running['expected_pgid'] == running['expected_sid'] > 1,
            'Worker dispatch parent, group, argv or deadline differs')


def read_activation(path):
    template = json.loads(check_control({'path': str(HERE / 'plan.pending.json'), 'sha256': PENDING_SHA256}))
    path = Path(path)
    require(path.is_absolute() and path.resolve(strict=True) == path, 'Activation must have a physical path')
    data = control_bytes(path)
    active = json.loads(data)
    require(all(template[key] is None for key in ('functional', 'native', 'primary'))
            and {**active, 'functional': None, 'native': None, 'primary': None} == template
            and set(active) == set(template), 'Activation differs from the frozen pending plan')
    return active, sha(data)


def load_request(plan_path, grant_path):
    plan, plan_sha = read_activation(plan_path)
    require(plan['schema'] == 'mss-darwin-matrix-plan-v1' and plan['harness_commit'] == HARNESS
            and set(plan['pins']) == set(HARD_PIN_NAMES), 'Unrecognized matrix plan')
    require_functional(plan)
    module(NATIVE_PROOF_SOURCE, 'matrix_native_proof').validate_native_proof(
        plan['native'], NATIVE_BINDING, check_control, os.path.lexists)
    grant_data = private_grant_bytes(grant_path)
    grant = json.loads(grant_data)
    require(grant['schema'] == 'mss-darwin-matrix-admission-v1' and grant['plan_sha256'] == plan_sha,
            'Unreviewed matrix plan or admission')
    mode = grant['mode']
    root = {'primary': PRIMARY, 'warm': WARM}[mode]
    require(grant['root'] == str(root) and not Path(grant_path).is_relative_to(root), 'Admission root differs')
    require(grant['limits'] == {'disk_bytes': DISK, 'floor_bytes': FLOOR, 'rss_bytes': GROUP_RSS,
                                'wall_seconds': sum(BUDGETS[mode].values())}, 'Unexpected resource envelope')
    names = ('driver.py', 'warm.py', 'prime.R')
    require(grant['controls'] == {name: pin_control(HERE / name) for name in names}, 'Owned controls changed')
    require(Path(sys.executable).resolve(strict=True) == Path(plan['pins']['python']['path']),
            'Wrong Python executable')
    require(plan['primary'] is not None if mode == 'warm' else plan['primary'] is None,
            'Warm requires separately observed primary results')
    return plan, grant, root, sha(grant_data)


def closed_phase(pin, root, label, allowed):
    path = root / f'{label}.jsonl'
    require(Path(pin['path']) == path, f'Unexpected phase path: {label}')
    rows = [json.loads(line) for line in check_control(pin).splitlines()]
    require(len(rows) == 3, f'Incomplete phase history: {label}')
    planned, running, terminal = rows
    require(all(row['label'] == label for row in rows)
            and planned['status'] == 'planned' and planned['cwd'] == str(root)
            and planned['receipt'] == str(path) and running['status'] == 'running'
            and running['pid'] == running['expected_pgid'] == running['expected_sid'] > 1,
            f'Unbound planned/running phase identity: {label}')
    phase_status(terminal, allowed)
    return planned, running, terminal


def validate_functional_header(q, observer, q_pin, controls, grant_hash):
    require(q['status'] == 'darwin_functional_qualified' and q['schema'] == 'mss-darwin-functional-stage177-v1'
            and observer['schema'] == 'mss-stage177-functional-root-observation-v1'
            and observer['qualification'] == q_pin and type(observer['observed_supervisor_exit']) is int
            and observer['observed_supervisor_exit'] == 0 and observer['all_groups_independently_absent'] is True
            and q['controls'] == controls and q['grant_sha256'] == grant_hash and q['failures'] == [],
            'Functional qualification lacks an independently observed successful exit')
    phases = {row['label']: row for row in observer['phases']}
    require(len(phases) == len(observer['phases']) == 21 and set(phases) == set(q['phases']) == FUNCTIONAL_PHASES,
            'Incomplete functional phase coverage')
    return phases


def check_functional_science(q, producer_plan, root):
    historical = producer_plan['historical']
    cold = {f'{role}-{order}-first': historical[f'{key}_cold'][order]
            for role, key in (('baseline', 'baseline'), ('new', 'candidate')) for order in ORDERS}
    require(q['cold_identity'] == cold, 'Functional cold evidence differs from the approved plan')
    for pin in cold.values():
        check_control(pin)
    smoke = q['smoke']
    cases = {f'{case}-f{fmt}': n for case, n in (('hard-cases', 32), ('declared', 9)) for fmt in (2, 3)}
    require(smoke['editions'] == 8 and smoke['normalized_field_count'] == 19
            and set(smoke['verifications']) == set(cases) and set(smoke['normalized_cross_format']) == {
                f'{role}-{case}' for role in ROLES for case in ('hard-cases', 'declared')},
            'Functional smoke coverage differs')
    require(Path(smoke['receipt']['path']) == root / 'smoke-verification.stdout', 'Unexpected smoke report')
    require(json.loads(check_control(smoke['receipt'])) == {key: val for key, val in smoke.items() if key != 'receipt'},
            'Smoke verification receipt differs')
    for name, n in cases.items():
        pin = smoke['verifications'][name]
        require(Path(pin['path']) == root / 'smoke' / name / 'verification.json', 'Unexpected smoke case path')
        data = json.loads(check_control(pin))
        require(data['stable_byte_equality'] is True and data['stable_rds_manifest_equality'] is True
                and data['edition_count'] == 2 and data['oracle_rows_checked'] == n
                and len(set(data['oracle_fields'])) == 19 and set(data['normalized_columns']) == set(data['oracle_fields']),
                f'Incomplete functional scientific smoke: {name}')
    require(set(q['testthat']) == set(ROLES), 'Missing exact-pair installed suites')
    for role in ROLES:
        report = q['testthat'][role]
        files = report['receipts']
        require(set(files) == {'qualification', 'conditions', 'tests'}, 'Missing unfiltered suite receipts')
        tables = {}
        for name, pin in files.items():
            require(Path(pin['path']) == root / f'{role}-testthat' / f'{name}.tsv', 'Unexpected suite receipt path')
            tables[name] = list(csv.DictReader(check_control(pin).decode().splitlines(keepends=True),
                                               delimiter='\t', strict=True))
        counts, = tables['qualification']
        tests = tables['tests']
        require(counts == report['counts'] and counts['qualified'] == 'TRUE' and int(counts['passed']) > 0
                and all(int(counts[key]) == 0 for key in ('skipped', 'failed', 'errors'))
                and len(tests) == int(counts['tests'])
                and len({row['file'] for row in tests}) == report['discovered_files'] - len(report['comment_only']),
                f'Functional installed suite was incomplete: {role}')
        if role == 'candidate':
            regressions = [row for row in tests if row['file'] == CANDIDATE_REGRESSION_FILE]
            require(len(regressions) == len(CANDIDATE_REGRESSIONS)
                    and {row['test'] for row in regressions} == CANDIDATE_REGRESSIONS
                    and all(int(row['passed']) > 0 for row in regressions)
                    and report.get('regressions') == regressions,
                    'Candidate regressions differ from the pinned suite or require positive passed assertions')


def closed_qualification(plan):
    proof = require_functional(plan)
    producer = plan['producer']
    root = Path(producer['root'])
    controls = producer['controls']
    directory = Path(producer['directory'])
    q_pin = proof['qualification']
    require(Path(q_pin['path']) == root / 'qualification.json', 'Unexpected qualification root')
    q = json.loads(check_control(q_pin))
    observer = json.loads(check_control(proof['root_observation']))
    payloads = {name: check_control({'path': str(directory / name), 'sha256': digest})
                for name, digest in controls.items()}
    producer_plan = json.loads(payloads['plan.json'])
    grant_data = private_grant_bytes(directory / producer['grant'])
    grant = json.loads(grant_data)
    grant_hash = sha(grant_data)
    started = json.loads(control_bytes(root / 'started.json'))
    require(producer_plan['root'] == str(root) and grant['authorize'] == producer['authorize']
            and grant['bindings'] == controls and grant['limits'] == producer_plan['limits']
            and started['grant'] == grant and started['grant_sha256'] == grant_hash,
            'Unbound functional producer grant or controls')
    matrix_guard = plan['pins']['guard']
    functional_guard = plan['pins']['functional_runtime_guard']
    historical_guard = plan['pins']['historical_reference_guard']
    reference_plan = json.loads(check_control(plan['pins']['reference_plan']))
    require(matrix_guard not in (functional_guard, historical_guard) and functional_guard != historical_guard
            and producer_plan['runtime_guard'] == grant['runtime_guard'] == functional_guard
            and all(record['runtime_guard'] == functional_guard and record['historical_reference_guard'] == historical_guard
                    for record in (started, q, observer))
            and reference_plan['pins']['guard'] == historical_guard, 'Matrix, functional or historical guard binding differs')
    for pin in (matrix_guard, functional_guard, historical_guard):
        check_control(pin)
    phases = validate_functional_header(q, observer, q_pin, controls, grant_hash)
    for label, row in phases.items():
        require(type(row['raw_exit']) is int and row['raw_exit'] == 0 and row['group_observed_absent'] is True,
                f'Unbound functional closure: {label}')
        planned, running, terminal = closed_phase(row, root, label, (0,))
        require(terminal == q['phases'][label] and row.get('group', running['pid']) == running['pid'],
                f'Invalid functional terminal record: {label}')
        if label == 'candidate-lookup-build':
            require('-DMSS_LOOKUP_COUNTERS' in planned['argv'] and planned['argv'][-1] == proof['lookup']['path'],
                    'Functional lookup counter build does not bind the candidate binary')
    check_functional_science(q, producer_plan, root)
    return q


def compare_candidate_build(plan, original, build):
    pins = plan['pins']
    baseline = json.loads(check_control(original['baseline']['build_receipt']))
    toolchain = baseline['toolchain']
    observed_compiler = check_control(pins['candidate_compiler']).decode().splitlines()
    observed_sdk = check_control(pins['candidate_sdk']).decode().strip()
    physical_sdk = Path(build['toolchain']['sdk'])
    require(observed_compiler == toolchain['compiler'] and build['compiler_flags'] == toolchain['compiler_flags']
            and Path(observed_sdk).is_absolute() and physical_sdk.is_absolute()
            and Path(observed_sdk).resolve(strict=True) == physical_sdk.resolve(strict=True) == physical_sdk
            and build['runtime'] == {key: original[key] for key in RUNTIME},
            'Observed candidate and baseline compiler, flags, SDK or R selections differ')
    observed_r = json.loads(check_control(toolchain['original_runtime_receipt']))
    desc = check_control(build['description']).decode()
    built, = [line.removeprefix('Built: ').split('; ') for line in desc.splitlines() if line.startswith('Built: ')]
    observed_version = observed_r['r_version'][0].split(' -- ', 1)[0]
    require(observed_r['r_runtime_sha256'] == original['r_engine']['sha256']
            and observed_r['compiler_version'] == toolchain['compiler']
            and observed_r['r_config'] == toolchain['compiler_flags']
            and observed_version.startswith('R version ') and built[0] == 'R ' + observed_version.split()[2],
            'Original R runtime or candidate installed Built field differs')
    common = {'baseline': {key: toolchain[key] for key in COMMON_TOOLCHAIN},
              'candidate': {'compiler': observed_compiler, 'compiler_flags': build['compiler_flags'],
                            'R_version': observed_version, 'R_platform': built[1]}}
    require(common['baseline'] == common['candidate'], 'Projected normal toolchains differ')
    transcript = '\n'.join(check_control(pins[name]).decode()
                           for name in ('candidate_stdout', 'candidate_stderr'))
    def canonical(line):
        return [f'-I{Path(arg[2:]).resolve(strict=True)}' if arg.startswith('-I') else arg
                for arg in shlex.split(line)]

    candidate_commands = [line for line in transcript.splitlines() if line.startswith('clang++ ')]
    lines = [canonical(line) for line in candidate_commands]
    expect = [canonical(line) for line in toolchain['actual_compiler_and_link_commands']]
    require(lines == expect and all(line in transcript for line in toolchain['actual_compiler_and_sdk']),
            'Candidate normal install transcript and baseline commands differ')
    return baseline, common, candidate_commands, observed_sdk


def matrix_build_receipts(plan, original, build, candidate):
    pins = plan['pins']
    baseline, common, candidate_commands, observed_sdk = compare_candidate_build(plan, original, build)
    base_role = original['baseline']
    require(baseline['source_commit'] == base_role['commit']
            and baseline['dll_sha256'] == base_role['dll']['sha256']
            and baseline['package_tree_sha256'] == base_role['package_tree_sha256'], 'Baseline build pin differs')
    upstream = {
        'baseline': {key: baseline[key] for key in ('source_archive', 'source_pin', 'install_argv', 'install_receipt',
                                                    'install_stdout', 'install_stderr', 'toolchain')},
        'candidate': {'source': build['source'], 'install_argv': build['argv'], 'toolchain': build['toolchain'],
                      'original_runtime': build['runtime'], 'installed_description': build['description'],
                      'root_observation': pins['focused'], 'compiler_output': pins['candidate_compiler'],
                      'sdk_output': pins['candidate_sdk'], 'sdk_literal': observed_sdk,
                      'sdk_physical': build['toolchain']['sdk'], 'install_stdout': pins['candidate_stdout'],
                      'install_stderr': pins['candidate_stderr'], 'actual_compiler_and_link_commands': candidate_commands},
    }
    return {name: {'scope': 'matrix_derived_normal_build_provenance', 'source_commit': role['commit'],
                   'dll_sha256': role['dll']['sha256'], 'package_tree_sha256': role['package_tree_sha256'],
                   'toolchain': common[name], 'upstream_build': role['build_receipt'], 'upstream': upstream[name]}
            for name, role in (('baseline', base_role), ('candidate', candidate))}


def inspect_installations(install, plan, root, *, closing):
    """Executed only in a dispatched group. Revalidates scientific package trees."""
    pins = plan['pins']
    q = closed_qualification(plan)
    functional = Path(plan['producer']['root'])
    original = json.loads(check_control(pins['baseline']))
    build = json.loads(check_control(pins['candidate']))
    focused = json.loads(check_control(pins['focused']))
    require(focused['observed_supervisor_exit'] == 0 and focused['all_groups_independently_absent'] is True
            and focused['source'] == build['source'] and focused['dll'] == build['dll']
            and build['source']['commit'] == CANDIDATE and build['source']['tree'] == CANDIDATE_TREE
            and q['stage176_exit'] == pins['focused'] and q['frozen_failure'] == focused['result']
            and set(focused['comparisons']) == {'width', 'precedence', 'frozen-failure'}
            and all(focused['comparisons'][key]['byte_equal'] is True
                    and focused['comparisons'][key]['case_count'] == count
                    for key, count in (('width', 4), ('precedence', 16), ('frozen-failure', 14))),
            'Focused source/build proof differs')
    baseline = original['baseline']
    candidate = q['roles']['candidate']
    require(q['roles']['baseline'] == baseline and baseline['commit'] == BASE
            and q['runtime'] == {key: original[key] for key in RUNTIME}
            and candidate['commit'] == CANDIDATE and candidate['dll'] == build['dll']
            and candidate['package_tree_sha256'] == build['tree']['tree_sha256'] == focused['package_tree_sha256']
            and candidate['build_receipt'] == pins['candidate']
            and Path(candidate['library']) / 'MungeSumstats/libs/MungeSumstats.so' == Path(build['dll']['path']),
            'Functional roles, build or runtime pins differ')
    records = matrix_build_receipts(plan, original, build, candidate)
    for field in ('source_archive', 'source_pin', 'install_receipt', 'install_stdout', 'install_stderr'):
        install.check_pin(records['baseline']['upstream'][field])
    install.check_pin(build['source'])
    reference_plan = json.loads(check_control(pins['reference_plan']))
    install.check_runtime(reference_plan, original)
    require(json.loads(check_control(q['reference_installation']['installations'])) == original,
            'Qualified dependency overlay differs')
    overlay = json.loads(check_control(q['reference_installation']['overlay']))
    new = {package['name'] for package in reference_plan['packages']}
    install.verify_overlay(reference_plan, {name: overlay['targets'][name] for name in new})
    tree_receipts = q['reference_installation']['trees']
    require(set(tree_receipts) == {name for name in overlay['targets']
                                 if name.startswith(('BSgenome.Hsapiens.', 'SNPlocs.Hsapiens.'))}
            and len(tree_receipts) == 6, 'Functional reference tree coverage differs')
    for name, expected in tree_receipts.items():
        target = overlay['targets'][name]['target']
        require(target == expected['target'] and install.installed_tree(Path(target))['tree_sha256']
                == expected['tree_sha256'], f'Qualified reference tree changed: {name}')
    for role in (baseline, candidate):
        require(install.installed_tree(Path(role['library']) / 'MungeSumstats')['tree_sha256']
                == role['package_tree_sha256'], f"Qualified installed tree changed: {role['commit']}")
        install.check_pin(role['dll'])
    require([Path(item['path']) for item in q['input_revalidation']]
            == [functional / f'inputs-{suffix}.stdout' for suffix in ('before', 'after')],
            'Missing functional opening or closing input record')
    for item in q['input_revalidation']:
        observed = json.loads(check_control(item))
        require(observed['roles'] == q['roles'] and observed['runtime'] == q['runtime']
                and observed['runtime_guard'] == pins['functional_runtime_guard']
                and observed['reference_installation'] == q['reference_installation'],
                'Functional input receipts disagree')
        expected_rcpp = observed['rcpp_headers']['tree_sha256']
        require(expected_rcpp == '74916f9a0362deee1c30a37c72dc2547b2a9413ec7a626734d5dcf92c6fa07a2',
                'Functional Rcpp header pin differs')
    rcpp = Path(reference_plan['targets']['Rcpp']['target']) / 'include'
    require(install.installed_tree(rcpp)['tree_sha256'] == expected_rcpp, 'Qualified Rcpp headers changed')
    lookup_path = install.check_pin(require_functional(plan)['lookup'])
    require(lookup_path == functional / 'candidate-lookup', 'Unexpected lookup counter binary')

    result = dict(original)
    for name, role in (('baseline', baseline), ('candidate', candidate)):
        path = root / f'{name}-build-normalized.json'
        expected = encode(records[name])
        if closing:
            require(control_bytes(path) == expected, f'Derived build provenance changed: {name}')
        else:
            with path.open('xb') as stream:
                stream.write(expected)
        result[name] = {**role, 'build_receipt': pin_control(path)}
    return result


def import_runner():
    sys.path.insert(0, str(BENCH))
    import run
    return run


def check_roster(root, runner):
    path = root / 'assets/matrix.json'
    config = json.loads(control_bytes(path))
    runner.resolve_configuration(config, 'matrix')
    roster = {case['id']: (case['rows'], case['pairs'], config['references'][case['reference']]['format'])
              for case in config['cases']}
    require(roster == expected_roster() and len(config['cases']) == len(roster), 'Generated matrix roster differs')
    require(all(case.get('target', False) == case['id'].startswith('matched-') for case in config['cases']),
            'Generated primary target membership differs')


def summarize_rss(observations):
    measured = []
    for process in observations:
        value = process['peak_rss_bytes']
        positive = value if type(value) is int and value > 0 else None
        measured.append({'peak_rss_bytes': positive, 'complete': positive is not None and process['rss_complete'] is True})
    largest = max((obs['peak_rss_bytes'] for obs in measured if obs['peak_rss_bytes'] is not None), default=None)
    complete = bool(measured) and all(obs['complete'] for obs in measured)
    status = ('miss' if largest is not None and largest > CANDIDATE_RSS
              else 'pass' if complete else 'indeterminate')
    return {'peak_rss_bytes': largest, 'rss_complete': complete, 'observations': measured, 'limit_status': status}


def summarize_primary(root):
    summary_pin = pin_control(root / 'matrix/summary.json')
    summary = json.loads(check_control(summary_pin))
    require(set(summary['cases']) == set(expected_roster()) and summary['fresh_processes'] is True
            and summary['profile'] == 'none' and summary['suite'] == 'matrix', 'Unexpected matrix summary scope')
    maxima, baselines, rss_by_case = {}, {}, {}
    for name, (rows, repeats, _) in expected_roster().items():
        case = summary['cases'][name]
        all_runs = case['observations']
        order = [(i, role) for i in range(1, repeats + 1) for role in (ROLES if i % 2 else ROLES[::-1])]
        require(case['rows'] == rows and [(run['pair'], run['process']['role']) for run in all_runs] == order
                and case['verification']['stable_rds_manifest_equality'] is True
                and case['verification']['stable_byte_equality'] is True, f'Incomplete matrix evidence: {name}')
        measurements = {role: summarize_rss([run['process'] for run in all_runs if run['process']['role'] == role])
                        for role in ROLES}
        rss_by_case[name] = measurements
        maxima[name] = measurements['candidate']['peak_rss_bytes']
        baselines[name] = measurements['baseline']['peak_rss_bytes']
        # Baseline retains observations and completeness without a candidate target annotation.
        measurements['baseline'].pop('limit_status')
    statuses = {case['candidate']['limit_status'] for case in rss_by_case.values()}
    status = 'miss' if 'miss' in statuses else 'indeterminate' if 'indeterminate' in statuses else 'pass'
    receipts = {size: pin_control(root / f'snp-oracle-{size}.json') for size in WARM_ORDER}
    controls = {}
    for size, pinned in receipts.items():
        data = json.loads(check_control(pinned))
        expected, _ = SCALES[size]
        require(data['all_rows_examined'] is True and data['passed'] is True and data['rows_expected'] == expected
                and data['versions']['mss7']['counts']['all_19_fields_correct_rows'] == expected,
                f'Incomplete independent SNP oracle: {size}')
        case = root / 'matrix' / f'matched-{size}'
        controls[size] = {'artifact_digests': pin_control(case / 'artifact-digests.json')}
        for role in ROLES:
            work = case / f'pair-01-{role}'
            controls[size][role] = {name: pin_control(work / name)
                                    for name in ('settings.tsv', 'metrics.tsv', 'namespaces.tsv', 'process.json')}
    config_pin = pin_control(root / 'assets/matrix.json')
    matrix_pin = pin_control(root / 'matrix/pins.json')
    pre = json.loads(check_control(matrix_pin))
    original = pre['configuration_file']
    require(all(original[key] == config_pin[key] for key in ('path', 'sha256')),
            'Retained matrix configuration differs from the executed configuration')
    diagnostic_paths = ('legacy/summary.json', 'alloc/summary.json', 'cpu/summary.json', 'lookup-pressure.stdout')
    history = TREE / 'benchmarks/2026-11-15/matrix-receipts.json'
    return {'summary': summary_pin, 'runner_target_status': summary['target_status'],
            'candidate_8GiB_status': status, 'candidate_rss_max_bytes': maxima,
            'baseline_rss_max_bytes': baselines, 'rss_by_case': rss_by_case,
            'snp_oracles': receipts, 'scale_controls': controls,
            'diagnostics': {name: pin_control(root / name) for name in diagnostic_paths},
            'historical_legacy': {'receipt': pin_control(history), '40M-cold': 'memory-censored',
                                  '40M-warm': 'memory-censored', 'scope': 'legacy b4a2d21 historical private footprint'},
            'matrix_config': config_pin, 'matrix_pins': matrix_pin,
            'matrix_post_pins': pin_control(root / 'matrix/post-pins.json')}


def run_primary(phase, root, py, bench, oracle, historical, lookup):
    assets, matrix = root / 'assets', root / 'matrix'
    run_cmd = [*py, bench / 'run.py', 'run']
    phase('prepare-assets', 'prep-diagnostics', [*py, bench / 'operations/assets.py', 'full',
          root / 'installations.json', assets, '--inputs', historical, '--admission-note', root / 'admission.txt'])
    phase('check-roster', 'prep-diagnostics')
    raw = phase('matrix', 'matrix', [*run_cmd, assets / 'matrix.json', matrix, '--suite', 'matrix'], (0, 1, 2, 3))
    failed = []

    def check(label, category, argv, allowed=(0, 1)):
        if phase(label, category, argv, allowed):
            failed.append(label)

    for name in SUPPLEMENTS:
        for role in ROLES:
            tables = [matrix / f'{name}-f{fmt}' / f'pair-01-{role}/edition/normalized.tsv' for fmt in (2, 3)]
            check(f'cross-{name}-{role}', 'matrix', ['/usr/bin/cmp', *tables], (0, 1, 2))
    for size in WARM_ORDER:
        edition = matrix / f'matched-{size}/pair-01-candidate/edition'
        check(f'snp-{size}', 'oracles', [*py, oracle, '--profile', assets / 'snp-oracle.json',
              '--stage', size, '--edition', f'mss7={edition}', '--receipt', root / f'snp-oracle-{size}.json'])
    check('legacy', 'prep-diagnostics', [*run_cmd, assets / 'legacy.json', root / 'legacy', '--suite', 'legacy'])
    check('allocation', 'prep-diagnostics', [*run_cmd, assets / 'diagnostic.json', root / 'alloc',
          '--suite', 'diagnostic', '--profile', 'alloc', '--allocation-threshold', '65536'])
    check('cpu', 'prep-diagnostics', [*run_cmd, assets / 'diagnostic.json', root / 'cpu',
          '--suite', 'diagnostic', '--profile', 'cpu'])
    check('lookup-pressure', 'prep-diagnostics', [lookup, 'pressure', root / 'lookup-pressure-fixture'])
    phase('inputs-close', 'identity')
    if raw == 1:
        failed.append('matrix')
    return {'measurement_exit': raw, 'failed': failed}


def primary_phases():
    return {'inputs-open', 'prepare-assets', 'check-roster', 'matrix', 'legacy', 'allocation', 'cpu',
            'lookup-pressure', 'inputs-close', 'matrix-evidence',
            *(f'cross-{case}-{role}' for case in SUPPLEMENTS for role in ROLES),
            *(f'snp-{size}' for size in WARM_ORDER)}


def prior_primary(plan, controls):
    spec = plan['primary']
    primary = json.loads(check_control(spec['result']))
    observed = json.loads(check_control(spec['root_observation']))
    require(Path(spec['result']['path']) == PRIMARY / 'result.json'
            and observed['schema'] == 'mss-darwin-matrix-primary-root-observation-v1'
            and observed['result'] == spec['result'] and primary['outcome'] == 'scientifically_complete'
            and primary['mode'] == 'primary' and type(observed['observed_supervisor_exit']) is int
            and type(primary['outer_exit']) is int and observed['observed_supervisor_exit'] == primary['outer_exit']
            and primary['outer_exit'] in (0, 2, 3) and observed['all_groups_independently_absent'] is True
            and primary['functional'] == plan['functional'] and primary['failed'] == [],
            'No independently closed matching primary science')
    require(observed['phases'] == primary['phases'] and set(primary['phases']) == primary_phases(),
            'Primary phase inventory changed')
    binding_pin = primary['admission_binding']
    require(Path(binding_pin['path']) == PRIMARY / 'admission-binding.json', 'Invalid primary launch binding')
    binding = json.loads(check_control(binding_pin))
    admission = binding['grant']
    active, active_sha = read_activation(binding['plan']['path'])
    grant_bytes = private_grant_bytes(admission['path'])
    grant = json.loads(grant_bytes)
    info = PRIMARY.lstat()
    require(sha(grant_bytes) == admission['sha256'] and binding['plan']['sha256'] == active_sha == grant['plan_sha256']
            and active == {**plan, 'primary': None}
            and grant['controls'] == binding['controls'] == controls and grant['mode'] == binding['mode'] == 'primary'
            and grant['root'] == str(PRIMARY) and grant['root_id'] == [info.st_dev, info.st_ino],
            'Primary launch controls, grant or frozen template differ')
    for label, pin in observed['phases'].items():
        _, _, terminal = closed_phase(pin, PRIMARY, label, (0, 2, 3) if label == 'matrix' else (0,))
        if label == 'matrix':
            require(terminal['raw_exit'] == primary['measurement_exit'], 'Primary runner raw exit differs')
    return primary['science']


def worker(label, digest, plan_path, grant_path):
    plan, grant, root, grant_sha = load_request(plan_path, grant_path)
    dispatch = control_bytes(root / f'{label}.dispatch.json')
    data = json.loads(dispatch)
    identity = root.lstat()
    require(sha(dispatch) == digest and data['label'] == label and data['root'] == str(root)
            and data['grant_sha256'] == grant_sha and data['root_id'] == grant['root_id']
            == [identity.st_dev, identity.st_ino] and Path.cwd() == root and data['parent_pid'] == os.getppid(),
            'Wrong guarded dispatch context')
    process = (os.getppid(), os.getpid(), os.getpgrp(), os.getsid(0), os.geteuid())
    require(process[1] == process[2] == process[3] > 1, 'Wrong guarded dispatch group')
    end = time.monotonic() + 5
    while True:
        try:
            blob = control_bytes(root / f'{label}.jsonl')
        except ValueError as error:
            if not str(error).startswith('Control changed during read:'):
                raise
            blob = b''  # The guard may have appended its running record during the read.
        if blob.endswith(b'\n') and blob.count(b'\n') >= 2:
            break
        require(time.monotonic() < end, 'Missing running dispatch record')
        time.sleep(.01)
    process = (os.getppid(), *process[1:])
    validate_dispatch(data, [json.loads(line) for line in blob.splitlines()],
                      [*python_argv(plan), str(HERE / 'driver.py'), '--worker', plan_path, grant_path, label, digest],
                      process, time.monotonic())

    pins = plan['pins']
    if label in ('inputs-open', 'inputs-close'):
        install = module(pins['installer'], 'matrix_reference_install')
        closing = label == 'inputs-close'
        result = inspect_installations(install, plan, root, closing=closing)
        for pin in pins.values():
            install.check_pin(pin)
        if closing:
            require(json.loads(control_bytes(root / 'installations.json')) == result, 'Installations changed')
        else:
            with (root / 'installations.json').open('xb') as stream:
                stream.write(encode(result))
        if grant['mode'] == 'warm':
            prior_primary(plan, grant['controls'])
        return 0
    runner = import_runner()
    if label == 'check-roster' and grant['mode'] == 'primary':
        check_roster(root, runner)
        return 0
    if label == 'matrix-evidence' and grant['mode'] == 'primary':
        print(encode(summarize_primary(root)).decode(), end='')
        return 0
    if grant['mode'] == 'warm':
        for prefix, action in (('warm-run-', 'measure'), ('warm-verify-', 'verify')):
            if label.startswith(prefix) and (size := label.removeprefix(prefix)) in SCALES:
                warm = module(grant['controls']['warm.py'], 'matrix_warm')
                return getattr(warm, action)(runner, prior_primary(plan, grant['controls']), root, PRIMARY, size,
                                             Path(pins['worker']['path']), HERE / 'prime.R')
    raise ValueError(f'Unexpected dispatch action: {label}')


def python_argv(plan):
    return [plan['pins']['python']['path'], '-I', '-B']


def environment(root):
    env = {'PATH': '/usr/bin:/bin:/usr/sbin:/sbin', 'HOME': str(root / 'home'),
           'TZ': 'UTC', 'LANG': 'en_US.UTF-8', 'LC_ALL': 'en_US.UTF-8',
           'R_ENVIRON_USER': '/dev/null', 'R_PROFILE_USER': '/dev/null',
           'R_LIBS_USER': '/dev/null', 'R_LIBS_SITE': '/dev/null'}
    env.update({key: str(root / 'tmp') for key in ('TMPDIR', 'TMP', 'TEMP')})
    env.update({name: '1' for name in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
                'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                'R_DATATABLE_NUM_THREADS')})
    return {**env, 'OMP_DYNAMIC': 'FALSE', 'MKL_DYNAMIC': 'FALSE'}


def run_warm(phase):
    failed = []
    for size in WARM_ORDER:
        for prefix, category in (('warm-run-', 'warm'), ('warm-verify-', 'warm-verify')):
            label = prefix + size
            if phase(label, category, allowed=(0, 1)):
                failed.append(label)
    phase('inputs-close', 'identity')
    return {'failed': failed}


def publish(install, home, guard, result):
    draft = '.result.pending'
    final = 'result.json'
    payload = encode(result)
    install.check_control_space(home, guard, (payload, payload))
    home.write(draft, payload)
    install.check_control_space(home, guard, (payload,))
    expected = install.directory_identity(home.stat(draft))
    linked = False
    try:
        os.link(draft, final, src_dir_fd=home.fd, dst_dir_fd=home.fd, follow_symlinks=False)
        linked = True
        require(install.directory_identity(home.stat(final)) == expected, 'Completion inode changed')
        install.check_control_space(home, guard, ())
        os.unlink(draft, dir_fd=home.fd)
        install.check_control_space(home, guard, ())
    except BaseException:
        if linked:
            home.withdraw(final, expected)
        raise


def execute(plan_path, grant_path):
    require(all(Path(path).is_absolute() for path in (plan_path, grant_path)), 'Use absolute control paths')
    plan, grant, root, grant_sha = load_request(plan_path, grant_path)
    mode = grant['mode']
    require(grant['all_other_heavy_work_excluded'] is True and grant['GRCh38_window_excluded'] is True
            and grant['retention_and_continuing_reservations_charged'] is True
            and isinstance(grant['contention_note'], str) and 0 < len(grant['contention_note'].strip()) <= 4096,
            'Missing fresh resource and contention observations')
    admitted = datetime.fromisoformat(grant['issued_utc'])
    now = datetime.now(timezone.utc)
    require(0 <= (now - admitted).total_seconds() <= 120, 'Admission is stale or future-dated')
    install = module(plan['pins']['installer'], 'matrix_reference_install')
    accepted = module(plan['pins']['guard'], 'matrix_darwin_guard')
    with install.AdmittedRoot(root, tuple(grant['root_id'])) as home:
        require(set(os.listdir(home.fd)) == {'home', 'tmp'}, 'Expected fresh root with only private home and tmp')
        for name in ('home', 'tmp'):
            info = home.stat(name)
            require(stat.S_ISDIR(info.st_mode) and stat.S_IMODE(info.st_mode) == 0o700
                    and info.st_uid == os.geteuid() and info.st_dev == grant['root_id'][0]
                    and not any((root / name).iterdir()), f'Unexpected private directory: {name}')
        fs = os.fstatvfs(home.fd)
        require(fs.f_bavail * fs.f_frsize >= FLOOR + DISK, 'Fresh capacity does not fund the entire admitted pool')
        require(ASSETS.stat().st_dev == grant['root_id'][0]
                and (mode == 'primary' or PRIMARY.stat().st_dev == grant['root_id'][0]),
                'Historical inputs, retained primary evidence and fresh root must share the admitted filesystem')
        guard = home.guard(accepted, accepted.GuardConfig(root, admitted, **grant['limits']))
        remaining = BUDGETS[mode].copy()
        phases = {}
        py = python_argv(plan)
        env = environment(root)

        def save(name, payload):
            install.check_control_space(home, guard, (payload,))
            home.write(name, payload)
            guard._stop_check()

        def phase(label, budget, argv=None, allowed=(0,)):
            home.check()
            guard._stop_check()
            started = time.monotonic()
            whole_deadline = guard.deadline
            guard.deadline = min(whole_deadline, started + remaining[budget])
            try:
                if argv is None:
                    dispatch = encode({'label': label, 'parent_pid': os.getpid(), 'root': str(root),
                                       'root_id': grant['root_id'], 'grant_sha256': grant_sha, 'allowed': allowed})
                    save(f'{label}.dispatch.json', dispatch)
                    argv = [*py, HERE / 'driver.py', '--worker', plan_path, grant_path, label, sha(dispatch)]
                path = root / f'{label}.jsonl'
                observed = guard.run_phase(label, list(map(str, argv)), cwd=root, env=env,
                                           stdout=root / f'{label}.stdout', stderr=root / f'{label}.stderr',
                                           receipt=path, allowed_exits=allowed)
                record = pin_control(path)
                phases[label] = record
                terminal = json.loads(check_control(record).splitlines()[-1])
                require(observed.status == terminal['status'] and observed.closure == terminal['closure']
                        and observed.raw_exit == terminal['raw_exit'] and observed.receipt_error is None,
                        f'Guard did not confirm persisted phase: {label}')
                return phase_status(terminal, allowed)
            finally:
                remaining[budget] -= time.monotonic() - started
                guard.deadline = whole_deadline

        with guard.shutdown_signals():
            note = grant['contention_note'].strip().encode() + b'\n'
            save('admission.txt', note)
            binding = encode({'plan': {'path': plan_path, 'sha256': grant['plan_sha256']},
                              'grant': {'path': grant_path, 'sha256': grant_sha},
                              'controls': grant['controls'], 'mode': mode, 'environment': env})
            save('admission-binding.json', binding)
            phase('inputs-open', 'identity')
            if mode == 'primary':
                # The frozen scripts import siblings by name. -E -s retains
                # their script directory and ignores Python environment/user site.
                outcome = run_primary(phase, root, [py[0], '-E', '-s', '-B'], BENCH, ORACLE, ASSETS,
                                      Path(require_functional(plan)['lookup']['path']))
                if not outcome['failed']:
                    phase('matrix-evidence', 'identity')
            else:
                outcome = run_warm(phase)
            science = None
            if outcome['failed']:
                code = 1
            elif mode == 'primary':
                science = json.loads(control_bytes(root / 'matrix-evidence.stdout'))
                raw = {'pass': 0, 'miss': 2, 'indeterminate': 3}[science['runner_target_status']]
                require(raw == outcome['measurement_exit'], 'Matrix raw exit and report disagree')
                rss = science['candidate_8GiB_status']
                code = 2 if raw == 2 or rss == 'miss' else 3 if raw == 3 or rss == 'indeterminate' else 0
            else:
                science = {size: json.loads(control_bytes(root / f'warm-verify-{size}.stdout'))
                           for size in WARM_ORDER}
                rss = {row['candidate_8GiB_status'] for row in science.values()}
                code = 2 if 'miss' in rss else 3 if 'indeterminate' in rss else 0
            publish(install, home, guard, {**outcome, 'mode': mode, 'science': science,
                    'functional': plan['functional'], 'outer_exit': code, 'phases': phases,
                    'admission_binding': {'path': str(root / 'admission-binding.json'), 'sha256': sha(binding)},
                    'outcome': 'incomplete' if code == 1 else 'scientifically_complete'})
            return code


if __name__ == '__main__':
    sys.modules['driver'] = sys.modules[__name__]
    os.umask(0o077)
    try:
        if len(sys.argv) == 6 and sys.argv[1] == '--worker':
            _, _, plan_arg, grant_arg, label_arg, digest_arg = sys.argv
            require(label_arg.isascii() and label_arg.replace('-', '').isalnum(), 'Invalid phase label')
            sys.exit(worker(label_arg, digest_arg, plan_arg, grant_arg))
        require(len(sys.argv) == 4 and sys.argv[1] == '--execute', 'Expected --execute PLAN GRANT')
        sys.exit(execute(*sys.argv[2:]))
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        sys.exit(f'MSS matrix stopped: {error}')
